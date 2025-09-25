mod types;
mod grid;
mod nanoparticle;
mod placement;
mod quasi_random;
mod background;
mod validation;
mod output;

use clap::{Parser, ValueEnum};
use nalgebra::Vector3;
use rand::{rngs::StdRng, SeedableRng};

use types::{AtomTypeMap, Nanoparticle};
use placement::place_nanoparticles;
use background::generate_liquid_background;
use validation::analyze_nearest_neighbors;
use output::write_lammps_data;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input nanoparticle XYZ files (can specify multiple)
    #[arg(short, long)]
    nanoparticles: Vec<String>,

    /// Simulation box side length (cubic box)
    #[arg(short = 'b', long)]
    box_size: f64,

    /// Number of nanoparticles to place
    #[arg(short, long, default_value = "10")]
    count: usize,

    /// Background material
    #[arg(short = 'm', long, default_value = "liquid")]
    material: BackgroundMaterial,

    /// Output LAMMPS data file
    #[arg(short, long, default_value = "system.data")]
    output: String,

    /// Voxel size for space partitioning (Angstroms)
    #[arg(short, long, default_value = "1.0")]
    voxel_size: f64,

    /// Minimum separation distance between nanoparticles
    #[arg(short, long, default_value = "2.0")]
    separation: f64,

    /// Liquid density (atoms per cubic Angstrom)
    #[arg(short, long, default_value = "0.08")]
    density: f64,

    /// Random seed
    #[arg(short = 'r', long)]
    seed: Option<u64>,

    /// Validate system by analyzing nearest neighbor distances
    #[arg(long)]
    validate: bool,

    /// Use quasi-random (Korobov) positioning instead of pure random
    #[arg(long)]
    quasi_random: bool,
}

#[derive(Clone, ValueEnum)]
enum BackgroundMaterial {
    Liquid,
    Fcc,
    None,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let box_size = Vector3::new(args.box_size, args.box_size, args.box_size);

    // Initialize RNG
    let mut rng = match args.seed {
        Some(seed) => {
            println!("Using random seed: {}", seed);
            StdRng::seed_from_u64(seed)
        }
        None => StdRng::from_entropy(),
    };

    // Initialize atom type mapping
    let mut type_map = AtomTypeMap::new();

    // Load nanoparticle templates
    println!("Loading nanoparticle templates...");
    if args.nanoparticles.is_empty() {
        return Err("At least one nanoparticle file must be specified".into());
    }

    let mut templates = Vec::new();
    for filename in &args.nanoparticles {
        println!("Loading: {}", filename);
        let np = Nanoparticle::from_xyz_file(filename, &mut type_map)?;
        println!("  {} atoms, radius: {:.2} Å",
                np.atoms.len(), np.radius);
        templates.push(np);
    }

    // Place nanoparticles
    let (placed_nanoparticles, occupied_grid) = place_nanoparticles(
        &templates,
        box_size,
        args.count,
        args.separation,
        &mut rng,
        args.quasi_random,
        &type_map,
        args.density,
    );

    // Print voxel occupation summary
    occupied_grid.print_occupation_summary();

    // Generate background
    let mut background = match args.material {
        BackgroundMaterial::Liquid => {
            if args.density > 0.0 {
                generate_liquid_background(box_size, args.density, &occupied_grid, &mut rng)
            } else {
                Vec::new()
            }
        }
        BackgroundMaterial::Fcc => {
            println!("FCC lattice generation not yet implemented, using liquid");
            if args.density > 0.0 {
                generate_liquid_background(box_size, args.density, &occupied_grid, &mut rng)
            } else {
                Vec::new()
            }
        }
        BackgroundMaterial::None => Vec::new(),
    };

    // Assign proper atom types to background atoms
    if !background.is_empty() {
        let bg_type = type_map.get_or_create_type("Co");
        for atom in &mut background {
            atom.atom_type = bg_type;
        }
    }

    // Write output
    write_lammps_data(&args.output, &placed_nanoparticles, &background, box_size, &type_map)?;

    // Analyze nearest neighbor distances if requested
    if args.validate {
        if let Err(e) = analyze_nearest_neighbors(&placed_nanoparticles, &background) {
            println!("Warning: Could not analyze nearest neighbors: {}", e);
        }
    }

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

    println!("\n=== Generation Complete ===");
    println!("Nanoparticles placed: {}", placed_nanoparticles.len());
    println!("Nanoparticle atoms: {}", np_atoms);
    println!("Background atoms: {}", background.len());
    println!("Total atoms: {}", total_atoms);
    println!("Box volume: {:.1} Å³", volume);
    println!("Actual liquid density: {:.6} atoms/Å³", actual_density);
    println!("Total mass density: {:.3} g/cm³", density_g_cm3);
    println!("Output file: {}", args.output);

    Ok(())
}
