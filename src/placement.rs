use crate::grid::NanoparticleGrid;
use crate::types::Nanoparticle;
use crate::quasi_random::QuasiRandom;
use nalgebra::{Point3, Vector3};
use rand::{rngs::StdRng, Rng};

pub fn place_nanoparticles_with_counts(
    templates: &[(Nanoparticle, usize)], // (template, count) pairs
    box_size: Vector3<f64>,
    separation: f64,
    rng: &mut StdRng,
    use_quasi_random: bool,
    type_map: &crate::types::AtomTypeMap,
    target_density: f64,
    placement_mode: &crate::config::PlacementMode,
) -> (Vec<Nanoparticle>, NanoparticleGrid) {
    // Calculate total count and average radius
    let total_count: usize = templates.iter().map(|(_, count)| count).sum();
    let _avg_radius = templates.iter().map(|(np, _)| np.radius()).sum::<f64>() / templates.len() as f64;

    // Create expansion vector with proper counts (clone templates)
    let mut expanded_templates = Vec::new();
    for (template, count) in templates {
        for _ in 0..*count {
            expanded_templates.push(template.clone());
        }
    }

    println!("Placing {} nanoparticles ({} types) in {}x{}x{} Å³ box using {} positioning...",
             total_count, templates.len(), box_size.x, box_size.y, box_size.z,
             if use_quasi_random { "quasi-random" } else { "random" });

    // Print breakdown by type
    for (i, (template, count)) in templates.iter().enumerate() {
        println!("  Type {}: {} particles, radius {:.1} Å", i+1, count, template.radius());
    }

    // Use existing placement logic with expanded templates
    let (placed, grid) = place_nanoparticles(
        &expanded_templates,
        box_size,
        total_count,
        separation,
        rng,
        use_quasi_random,
        type_map,
        target_density,
        placement_mode,
    );

    (placed, grid)
}

pub fn place_nanoparticles(
    templates: &[Nanoparticle],
    box_size: Vector3<f64>,
    count: usize,
    separation: f64,
    rng: &mut StdRng,
    use_quasi_random: bool,
    type_map: &crate::types::AtomTypeMap,
    target_density: f64,
    placement_mode: &crate::config::PlacementMode,
) -> (Vec<Nanoparticle>, NanoparticleGrid) {
    // Calculate average radius for grid sizing
    let avg_radius = templates.iter().map(|np| np.radius()).sum::<f64>() / templates.len() as f64;
    let mut grid = NanoparticleGrid::new(avg_radius, box_size, type_map, target_density);
    let mut placed_particles = Vec::new();
    let mut all_placed_atoms = Vec::new(); // Track all placed atoms for overlap detection

    println!("Placing {} nanoparticles in {}x{}x{} Å³ box using {} positioning...",
             count, box_size.x, box_size.y, box_size.z,
             if use_quasi_random { "quasi-random" } else { "random" });

    // Analyze packing constraints
    let max_radius = templates.iter().map(|np| np.radius()).fold(0.0, f64::max);
    let min_box_dim = box_size.x.min(box_size.y).min(box_size.z);
    let max_diameter_with_sep = (max_radius + separation) * 2.0;
    let particles_per_dim = (min_box_dim / max_diameter_with_sep).floor() as usize;
    let theoretical_max = particles_per_dim.pow(3);

    println!("Packing analysis:");
    println!("  Largest nanoparticle radius: {:.1} Å", max_radius);
    println!("  Required spacing per particle: {:.1} Å", max_diameter_with_sep);
    println!("  Theoretical maximum particles: ~{} ({}³ grid)", theoretical_max, particles_per_dim);
    if count > theoretical_max {
        println!("  ⚠️  Warning: Requesting {} particles may exceed packing limit", count);
    }

    // Initialize quasi-random generator if needed (6D: 3 position + 3 orientation)
    let mut quasi_gen = if use_quasi_random {
        Some(QuasiRandom::new(6))
    } else {
        None
    };

    for i in 0..count {
        let template = &templates[i % templates.len()];
        let mut attempts = 0;
        let max_attempts = 50000;

        while attempts < max_attempts {
            let margin = template.radius() + separation;
            if margin * 2.0 > box_size.x.min(box_size.y).min(box_size.z) {
                println!("Warning: Nanoparticle too large for box with current separation");
                break;
            }

            // Get position and orientation
            let (test_center, rotation) = if let Some(ref mut qgen) = quasi_gen {
                let quasi_vals = qgen.next();
                let x = margin + quasi_vals[0] * (box_size.x - 2.0 * margin);
                let y = margin + quasi_vals[1] * (box_size.y - 2.0 * margin);
                let z = margin + quasi_vals[2] * (box_size.z - 2.0 * margin);

                // Use quasi-random values for orientation
                let u1 = quasi_vals[3];
                let u2 = quasi_vals[4];
                let u3 = quasi_vals[5];
                let rotation = Nanoparticle::quasi_random_orientation(u1, u2, u3);

                (Point3::new(x, y, z), rotation)
            } else {
                let x = rng.gen_range(margin..(box_size.x - margin));
                let y = rng.gen_range(margin..(box_size.y - margin));
                let z = rng.gen_range(margin..(box_size.z - margin));
                let rotation = Nanoparticle::random_orientation(rng);

                (Point3::new(x, y, z), rotation)
            };

            // Check for collisions with previously placed nanoparticles
            let should_place = match placement_mode {
                crate::config::PlacementMode::Collision => {
                    // Strict collision avoidance - no overlaps allowed
                    !grid.check_collision(template, &test_center, separation, type_map)
                },
                crate::config::PlacementMode::Overlap => {
                    // Allow overlaps - first placed wins, so always place if position is valid
                    true
                }
            };

            if should_place {
                // Create and place the nanoparticle
                let mut new_particle = template.clone();
                new_particle.rotate(&rotation);
                let offset = test_center - new_particle.center;
                new_particle.translate(&offset);

                grid.add_nanoparticle(&new_particle, test_center, type_map);

                // Handle atom overlap based on placement mode
                let final_atoms = grid.mark_atoms_occupied_with_overlap(&new_particle.atoms, placement_mode, &mut all_placed_atoms);

                // Update the nanoparticle with the final set of atoms (important for overlap mode)
                let mut final_particle = new_particle.clone();
                final_particle.atoms = final_atoms;

                placed_particles.push(final_particle);

                if (i + 1) % 5 == 0 || i == count - 1 {
                    println!("Placed {}/{} nanoparticles (attempt {})", i + 1, count, attempts + 1);
                }
                break;
            }

            attempts += 1;
        }

        if attempts == max_attempts {
            println!("Warning: Could only place {} out of {} nanoparticles",
                    placed_particles.len(), count);
            break;
        }
    }

    println!("Successfully placed {} nanoparticles", placed_particles.len());
    (placed_particles, grid)
}