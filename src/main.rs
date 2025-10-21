mod types;
mod config;
mod grid;
mod nanoparticle;
mod placement;
mod quasi_random;
mod background;
mod validation;
mod output;

// New custom dynamics modules
mod periodic;
mod potentials;
mod rigid_body;
mod surface_discretization;
mod forces;
mod minimization;

use clap::Parser;
use nalgebra::Vector3;
use rand::{rngs::StdRng, SeedableRng};

use config::Config;
use types::{AtomTypeMap, Nanoparticle};
use placement::place_nanoparticles_with_counts;
use background::generate_liquid_background;
use validation::analyze_nearest_neighbors;
use output::write_lammps_data;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Configuration file (TOML format)
    config: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_time = std::time::Instant::now();
    let args = Args::parse();

    // Load configuration
    let config_start = std::time::Instant::now();
    let config = Config::from_file(&args.config)?;
    let config_time = config_start.elapsed();
    println!("Loaded configuration from: {} ({:.2}s)", args.config, config_time.as_secs_f64());

    let box_size = Vector3::new(config.system.box_size, config.system.box_size, config.system.box_size);

    // Initialize RNG
    let mut rng = match config.system.random_seed {
        Some(seed) => {
            println!("Using random seed: {}", seed);
            StdRng::seed_from_u64(seed)
        }
        None => StdRng::from_entropy(),
    };

    // Initialize atom type mapping
    let mut type_map = AtomTypeMap::new();

    // First pass: collect all elements from nanoparticles and background
    let analysis_start = std::time::Instant::now();
    println!("Analyzing elements for consistent atom type assignment...");
    if config.nanoparticles.is_empty() {
        return Err("At least one nanoparticle must be specified in config".into());
    }

    let mut all_elements = Vec::new();
    for np_config in &config.nanoparticles {
        let elements = Nanoparticle::collect_elements_from_xyz(&np_config.file)?;
        all_elements.extend(elements);
    }

    // Add background element
    all_elements.push(config.background.element.clone());

    // Pre-populate atom types in mass order
    type_map.pre_populate_types(&all_elements);
    let analysis_time = analysis_start.elapsed();

    // Second pass: Load nanoparticle templates with pre-determined types
    let loading_start = std::time::Instant::now();
    println!("Loading nanoparticle templates...");
    let mut templates = Vec::new();
    for np_config in &config.nanoparticles {
        println!("Loading: {} (count: {})", np_config.file, np_config.count);
        let np = Nanoparticle::from_xyz_file(&np_config.file, &mut type_map)?;
        println!("  {} atoms, radius: {:.2} Å", np.atoms.len(), np.radius());
        templates.push((np, np_config.count));
    }
    let loading_time = loading_start.elapsed();

    // Place nanoparticles with individual counts
    let placement_start = std::time::Instant::now();
    let (placed_nanoparticles, occupied_grid) = place_nanoparticles_with_counts(
        &templates,
        box_size,
        config.background.separation,
        &mut rng,
        config.analysis.quasi_random,
        &type_map,
        config.background.density,
        &config.placement,
        &config.dynamics,
    );
    let placement_time = placement_start.elapsed();

    // Print voxel occupation summary
    occupied_grid.print_occupation_summary();

    // Generate background with proper element and type assignment
    let background_start = std::time::Instant::now();
    let background = match config.background.material {
        config::BackgroundMaterial::Liquid => {
            if config.background.density > 0.0 {
                let bg_type = type_map.get_or_create_type(&config.background.element);
                generate_liquid_background(
                    box_size,
                    config.background.density,
                    &occupied_grid,
                    &mut rng,
                    &config.background.element,
                    bg_type
                )
            } else {
                Vec::new()
            }
        }
        config::BackgroundMaterial::Fcc => {
            println!("FCC lattice generation not yet implemented, using liquid");
            if config.background.density > 0.0 {
                let bg_type = type_map.get_or_create_type(&config.background.element);
                generate_liquid_background(
                    box_size,
                    config.background.density,
                    &occupied_grid,
                    &mut rng,
                    &config.background.element,
                    bg_type
                )
            } else {
                Vec::new()
            }
        }
        config::BackgroundMaterial::None => Vec::new(),
    };
    let background_time = background_start.elapsed();

    // Write output
    let output_start = std::time::Instant::now();
    write_lammps_data(&config.system.output_file, &placed_nanoparticles, &background, box_size, &type_map)?;
    let output_time = output_start.elapsed();

    // Analyze nearest neighbor distances if requested
    let validation_start = std::time::Instant::now();
    if config.analysis.validate {
        if let Err(e) = analyze_nearest_neighbors(&placed_nanoparticles, &background) {
            println!("Warning: Could not analyze nearest neighbors: {}", e);
        }
    }
    let validation_time = validation_start.elapsed();

    // Summary
    let np_atoms: usize = placed_nanoparticles.iter().map(|np| np.atoms.len()).sum();
    let total_atoms = np_atoms + background.len();
    let volume = box_size.x * box_size.y * box_size.z;
    let actual_density = background.len() as f64 / volume;

    // Calculate mass-based density in g/cm³
    let mut total_mass_amu = 0.0;

    // Add nanoparticle masses
    for np in &placed_nanoparticles {
        for atom in &np.atoms {
            if let Some(&mass) = type_map.type_to_mass.get(&atom.atom_type) {
                total_mass_amu += mass;
            }
        }
    }

    // Add background masses
    for atom in &background {
        if let Some(&mass) = type_map.type_to_mass.get(&atom.atom_type) {
            total_mass_amu += mass;
        }
    }

    // Convert to g/cm³: 1 amu = 1.66054e-24 g, 1 Å³ = 1e-24 cm³
    let volume_cm3 = volume * 1e-24; // Å³ to cm³
    let mass_g = total_mass_amu * 1.66054e-24; // amu to g
    let density_g_cm3 = mass_g / volume_cm3;
    let total_time = start_time.elapsed();

    println!("\n=== Generation Complete ===");
    println!("Nanoparticles placed: {}", placed_nanoparticles.len());
    println!("Nanoparticle atoms: {}", np_atoms);
    println!("Background atoms: {}", background.len());
    println!("Total atoms: {}", total_atoms);
    println!("Box volume: {:.1} Å³", volume);
    println!("Actual liquid density: {:.6} atoms/Å³", actual_density);
    println!("Total mass density: {:.3} g/cm³", density_g_cm3);
    println!("Output file: {}", config.system.output_file);

    println!("\n=== Timing Breakdown ===");
    println!("Configuration loading: {:.2}s ({:.1}%)",
             config_time.as_secs_f64(),
             (config_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    println!("Element analysis: {:.2}s ({:.1}%)",
             analysis_time.as_secs_f64(),
             (analysis_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    println!("Template loading: {:.2}s ({:.1}%)",
             loading_time.as_secs_f64(),
             (loading_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    println!("Nanoparticle placement: {:.2}s ({:.1}%)",
             placement_time.as_secs_f64(),
             (placement_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    println!("Background generation: {:.2}s ({:.1}%)",
             background_time.as_secs_f64(),
             (background_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    println!("File output: {:.2}s ({:.1}%)",
             output_time.as_secs_f64(),
             (output_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    if config.analysis.validate {
        println!("Validation analysis: {:.2}s ({:.1}%)",
                 validation_time.as_secs_f64(),
                 (validation_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0);
    }
    println!("Total time: {:.2}s", total_time.as_secs_f64());

    Ok(())
}
