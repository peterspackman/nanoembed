use crate::grid::NanoparticleGrid;
use crate::types::Nanoparticle;
use crate::quasi_random::QuasiRandom;
use nalgebra::{Point3, Vector3, Rotation3};
use rand::{rngs::StdRng, Rng};

// Dynamics imports
use crate::rigid_body::RigidBody;
use crate::surface_discretization;
use crate::potentials::{Potential, SoftSphereParams, HardSphereParams};
use crate::minimization::{minimize_energy, MinimizerParams};
use crate::config::{DynamicsConfig, PotentialType};

pub fn place_nanoparticles_with_counts(
    templates: &[(Nanoparticle, usize)], // (template, count) pairs
    box_size: Vector3<f64>,
    separation: f64,
    rng: &mut StdRng,
    use_quasi_random: bool,
    type_map: &crate::types::AtomTypeMap,
    target_density: f64,
    placement_config: &crate::config::PlacementConfig,
    dynamics_config: &DynamicsConfig,
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
        placement_config,
    );

    // Apply dynamics-based relaxation if enabled
    let relaxed = relax_with_dynamics(placed, box_size, dynamics_config);

    (relaxed, grid)
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
    placement_config: &crate::config::PlacementConfig,
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
            let should_place = match placement_config.mode {
                crate::config::PlacementMode::Collision => {
                    // Strict collision avoidance - no overlaps allowed
                    !grid.check_collision(template, &test_center, separation, type_map, placement_config.collision_buffer)
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

                grid.add_nanoparticle(&new_particle, test_center, type_map, placement_config.collision_buffer);

                // Handle atom overlap based on placement mode
                let final_atoms = grid.mark_atoms_occupied_with_overlap(&new_particle.atoms, &placement_config.mode, &mut all_placed_atoms, placement_config.min_atom_distance);

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

/// Convert a Nanoparticle to a RigidBody for dynamics simulation
fn nanoparticle_to_rigid_body(np: &Nanoparticle) -> RigidBody {
    // Extract hull points in body frame (centered at origin)
    let hull_points: Vec<Point3<f64>> = np.atoms
        .iter()
        .map(|atom| atom.position - np.center.coords)
        .collect();

    // Extract orientation from current nanoparticle state
    // For now, assume identity rotation (nanoparticle already rotated)
    let orientation = nalgebra::UnitQuaternion::identity();

    RigidBody::from_hull_points(hull_points, np.center, orientation)
}

/// Convert a RigidBody back to Nanoparticle (update positions)
fn rigid_body_to_nanoparticle(rb: &RigidBody, original_np: &Nanoparticle) -> Nanoparticle {
    let mut updated_np = original_np.clone();

    // Update center position
    updated_np.center = rb.position;

    // Update atom positions using rigid body transform
    for (i, atom) in updated_np.atoms.iter_mut().enumerate() {
        if i < rb.hull_points.len() {
            // Transform from body frame to world frame
            let world_pos = rb.position + rb.orientation * rb.hull_points[i].coords;
            atom.position = world_pos;
        }
    }

    updated_np
}

/// Apply dynamics-based relaxation to placed nanoparticles
/// This uses energy minimization with PBC to remove overlaps and optimize packing
pub fn relax_with_dynamics(
    placed_particles: Vec<Nanoparticle>,
    box_size: Vector3<f64>,
    dynamics_config: &DynamicsConfig,
) -> Vec<Nanoparticle> {
    if !dynamics_config.enabled {
        println!("Dynamics relaxation disabled, skipping...");
        return placed_particles;
    }

    println!("\n=== Starting Dynamics Relaxation ===");
    println!("Converting {} nanoparticles to rigid bodies...", placed_particles.len());

    // Convert nanoparticles to rigid bodies
    let mut rigid_bodies: Vec<RigidBody> = placed_particles
        .iter()
        .map(nanoparticle_to_rigid_body)
        .collect();

    // Discretize surfaces
    println!("Discretizing surfaces into collision spheres...");
    let mut total_spheres = 0;
    for body in &mut rigid_bodies {
        body.surface_spheres = surface_discretization::discretize_convex_hull(
            &body.hull_points,
            dynamics_config.surface.sphere_radius,
            dynamics_config.surface.target_spacing,
        );
        total_spheres += body.surface_spheres.len();
    }

    let avg_spheres = if !rigid_bodies.is_empty() {
        total_spheres / rigid_bodies.len()
    } else {
        0
    };
    println!("  Total surface spheres: {} (avg {:.1} per particle)",
             total_spheres, avg_spheres);

    // Create potential from config
    let potential = match &dynamics_config.potential {
        PotentialType::SoftSphere { epsilon, sigma } => {
            println!("Using soft sphere (WCA) potential:");
            println!("  epsilon = {:.3} kcal/mol", epsilon);
            println!("  sigma = {:.3} Å", sigma);
            Potential::SoftSphere(SoftSphereParams::new_wca(*epsilon, *sigma))
        }
        PotentialType::HardSphere { penalty_strength } => {
            println!("Using hard sphere potential:");
            println!("  penalty strength = {:.1} kcal/mol/Å²", penalty_strength);
            Potential::HardSphere(HardSphereParams::new(*penalty_strength))
        }
    };

    // Calculate initial energy
    let initial_energy = crate::forces::calculate_total_energy(&rigid_bodies, box_size, &potential);
    println!("Initial total energy: {:.6} kcal/mol", initial_energy);

    // Set up minimizer parameters
    let minimizer_params = MinimizerParams {
        max_iterations: dynamics_config.max_iterations,
        force_tolerance: dynamics_config.force_tolerance,
        energy_tolerance: dynamics_config.energy_tolerance,
        max_displacement: dynamics_config.max_displacement,
        max_rotation: dynamics_config.max_rotation,
        initial_step_size: 0.01,
        use_fire: dynamics_config.use_fire,
        print_interval: 100,  // Print every 100 iterations
    };

    println!("\nMinimization settings:");
    println!("  Algorithm: {}", if minimizer_params.use_fire { "FIRE" } else { "Steepest Descent" });
    println!("  Max iterations: {}", minimizer_params.max_iterations);
    println!("  Force tolerance: {:.6} kcal/mol/Å", minimizer_params.force_tolerance);
    println!("  Energy tolerance: {:.6} kcal/mol", minimizer_params.energy_tolerance);
    println!("  Max displacement: {:.3} Å", minimizer_params.max_displacement);
    println!("  Max rotation: {:.3} rad", minimizer_params.max_rotation);

    // Run minimization
    println!("\nRunning energy minimization...\n");
    let result = minimize_energy(&mut rigid_bodies, box_size, &potential, &minimizer_params);

    // Print results
    match result {
        crate::minimization::MinimizationResult::Converged {
            iterations,
            final_energy,
            final_max_force,
        } => {
            println!("\n✓ Minimization converged in {} iterations", iterations);
            println!("  Final energy: {:.6} kcal/mol (ΔE = {:.6})",
                     final_energy, final_energy - initial_energy);
            println!("  Final max force: {:.6} kcal/mol/Å", final_max_force);
        }
        crate::minimization::MinimizationResult::MaxIterations {
            iterations,
            final_energy,
            final_max_force,
        } => {
            println!("\n⚠ Reached maximum iterations ({})", iterations);
            println!("  Final energy: {:.6} kcal/mol (ΔE = {:.6})",
                     final_energy, final_energy - initial_energy);
            println!("  Final max force: {:.6} kcal/mol/Å", final_max_force);
            println!("  Note: System may not be fully relaxed");
        }
    }

    // Convert rigid bodies back to nanoparticles
    println!("\nConverting rigid bodies back to nanoparticles...");
    let relaxed_particles: Vec<Nanoparticle> = rigid_bodies
        .iter()
        .zip(placed_particles.iter())
        .map(|(rb, original_np)| rigid_body_to_nanoparticle(rb, original_np))
        .collect();

    println!("=== Dynamics Relaxation Complete ===\n");

    relaxed_particles
}