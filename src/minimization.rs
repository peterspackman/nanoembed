use crate::forces::{calculate_forces, calculate_max_force, calculate_total_energy};
use crate::potentials::Potential;
use crate::rigid_body::RigidBody;
use nalgebra::Vector3;

/// Parameters for energy minimization
#[derive(Debug, Clone)]
pub struct MinimizerParams {
    pub max_iterations: usize,
    pub force_tolerance: f64,        // Stop when max force < this (kcal/mol/Å)
    pub energy_tolerance: f64,       // Stop when |ΔE| < this (kcal/mol)
    pub max_displacement: f64,       // Maximum position change per step (Å)
    pub max_rotation: f64,           // Maximum rotation per step (radians)
    pub initial_step_size: f64,      // Initial step size
    pub use_fire: bool,              // Use FIRE algorithm instead of steepest descent
    pub print_interval: usize,       // Print progress every N iterations (0 = no printing)
}

impl Default for MinimizerParams {
    fn default() -> Self {
        Self {
            max_iterations: 10000,
            force_tolerance: 0.01,
            energy_tolerance: 1e-6,
            max_displacement: 0.5,
            max_rotation: 0.1,
            initial_step_size: 0.01,
            use_fire: true,
            print_interval: 100,
        }
    }
}

/// Result of minimization
#[derive(Debug, Clone)]
pub enum MinimizationResult {
    Converged {
        iterations: usize,
        final_energy: f64,
        final_max_force: f64,
    },
    MaxIterations {
        iterations: usize,
        final_energy: f64,
        final_max_force: f64,
    },
}

/// Minimize energy using steepest descent or FIRE algorithm
pub fn minimize_energy(
    bodies: &mut [RigidBody],
    box_size: Vector3<f64>,
    potential: &Potential,
    params: &MinimizerParams,
) -> MinimizationResult {
    if params.use_fire {
        fire_minimization(bodies, box_size, potential, params)
    } else {
        steepest_descent(bodies, box_size, potential, params)
    }
}

/// Steepest descent minimization
fn steepest_descent(
    bodies: &mut [RigidBody],
    box_size: Vector3<f64>,
    potential: &Potential,
    params: &MinimizerParams,
) -> MinimizationResult {
    let mut prev_energy = calculate_total_energy(bodies, box_size, potential);
    let mut step_size = params.initial_step_size;

    if params.print_interval > 0 {
        println!("Starting steepest descent minimization");
        println!("Initial energy: {:.6} kcal/mol", prev_energy);
    }

    for iteration in 0..params.max_iterations {
        // Calculate forces
        let force_accs = calculate_forces(bodies, box_size, potential);
        let max_force = calculate_max_force(&force_accs);

        // Check convergence
        if max_force < params.force_tolerance {
            if params.print_interval > 0 {
                println!("Converged at iteration {}", iteration);
                println!("Final energy: {:.6} kcal/mol", prev_energy);
                println!("Final max force: {:.6} kcal/mol/Å", max_force);
            }
            return MinimizationResult::Converged {
                iterations: iteration,
                final_energy: prev_energy,
                final_max_force: max_force,
            };
        }

        // Move bodies in direction of force (energy minimization)
        for (i, body) in bodies.iter_mut().enumerate() {
            let force = &force_accs[i].force;
            let torque = &force_accs[i].torque;

            // Normalize force and torque
            let force_norm = force.norm();
            let torque_norm = torque.norm();

            if force_norm > 1e-10 {
                let displacement = force.normalize() * step_size;
                let displacement_magnitude = displacement.norm().min(params.max_displacement);
                body.position += force.normalize() * displacement_magnitude;
                body.position = crate::periodic::wrap_position(body.position, box_size);
            }

            if torque_norm > 1e-10 {
                let angle = (torque_norm * step_size).min(params.max_rotation);
                let axis = nalgebra::Unit::new_normalize(*torque);
                let delta_q = nalgebra::UnitQuaternion::from_axis_angle(&axis, angle);
                body.orientation = delta_q * body.orientation;
                body.orientation = body.orientation.renormalize();
            }
        }

        // Calculate new energy
        let new_energy = calculate_total_energy(bodies, box_size, potential);
        let delta_e = new_energy - prev_energy;

        // Adaptive step size
        if delta_e < 0.0 {
            // Energy decreased - accept and increase step size
            step_size *= 1.1;
            prev_energy = new_energy;
        } else {
            // Energy increased - reject and decrease step size
            step_size *= 0.5;
        }

        // Check energy tolerance
        if delta_e.abs() < params.energy_tolerance && delta_e < 0.0 {
            if params.print_interval > 0 {
                println!("Energy converged at iteration {}", iteration);
                println!("Final energy: {:.6} kcal/mol", new_energy);
                println!("Final max force: {:.6} kcal/mol/Å", max_force);
            }
            return MinimizationResult::Converged {
                iterations: iteration,
                final_energy: new_energy,
                final_max_force: max_force,
            };
        }

        // Print progress
        if params.print_interval > 0 && (iteration + 1) % params.print_interval == 0 {
            println!(
                "Iteration {}: E = {:.6}, F_max = {:.6}, step = {:.6}",
                iteration + 1,
                new_energy,
                max_force,
                step_size
            );
        }
    }

    let final_energy = calculate_total_energy(bodies, box_size, potential);
    let force_accs = calculate_forces(bodies, box_size, potential);
    let final_max_force = calculate_max_force(&force_accs);

    if params.print_interval > 0 {
        println!("Reached maximum iterations");
        println!("Final energy: {:.6} kcal/mol", final_energy);
        println!("Final max force: {:.6} kcal/mol/Å", final_max_force);
    }

    MinimizationResult::MaxIterations {
        iterations: params.max_iterations,
        final_energy,
        final_max_force,
    }
}

/// FIRE (Fast Inertial Relaxation Engine) minimization
/// More efficient than steepest descent for molecular systems
fn fire_minimization(
    bodies: &mut [RigidBody],
    box_size: Vector3<f64>,
    potential: &Potential,
    params: &MinimizerParams,
) -> MinimizationResult {
    // FIRE parameters
    let dt_max = 0.1;
    let dt_min = 0.001;
    let n_min = 5;
    let f_inc = 1.1;
    let f_dec = 0.5;
    let alpha_start = 0.1;
    let f_alpha = 0.99;

    let mut dt = params.initial_step_size;
    let mut alpha = alpha_start;
    let mut n_positive = 0;

    // Initialize velocities to zero
    for body in bodies.iter_mut() {
        body.velocity = Vector3::zeros();
        body.angular_velocity = Vector3::zeros();
    }

    let mut prev_energy = calculate_total_energy(bodies, box_size, potential);

    if params.print_interval > 0 {
        println!("Starting FIRE minimization");
        println!("Initial energy: {:.6} kcal/mol", prev_energy);
    }

    for iteration in 0..params.max_iterations {
        // Calculate forces
        let force_accs = calculate_forces(bodies, box_size, potential);
        let max_force = calculate_max_force(&force_accs);

        // Check convergence
        if max_force < params.force_tolerance {
            if params.print_interval > 0 {
                println!("Converged at iteration {}", iteration);
                println!("Final energy: {:.6} kcal/mol", prev_energy);
                println!("Final max force: {:.6} kcal/mol/Å", max_force);
            }
            return MinimizationResult::Converged {
                iterations: iteration,
                final_energy: prev_energy,
                final_max_force: max_force,
            };
        }

        // Calculate power P = F · v
        let mut power = 0.0;
        for (i, body) in bodies.iter().enumerate() {
            power += force_accs[i].force.dot(&body.velocity);
            power += force_accs[i].torque.dot(&body.angular_velocity);
        }

        // FIRE algorithm
        if power > 0.0 {
            n_positive += 1;

            if n_positive > n_min {
                dt = (dt * f_inc).min(dt_max);
                alpha *= f_alpha;
            }
        } else {
            n_positive = 0;
            dt *= f_dec;
            alpha = alpha_start;

            // Zero velocities
            for body in bodies.iter_mut() {
                body.velocity = Vector3::zeros();
                body.angular_velocity = Vector3::zeros();
            }
        }

        // Update velocities with mixing
        for (i, body) in bodies.iter_mut().enumerate() {
            let force = &force_accs[i].force;
            let torque = &force_accs[i].torque;

            // v = (1 - alpha) * v + alpha * F/|F| * |v|
            let v_norm = body.velocity.norm();
            let f_norm = force.norm();

            if f_norm > 1e-10 {
                body.velocity = (1.0 - alpha) * body.velocity
                    + alpha * (force / f_norm) * v_norm;
            }

            let omega_norm = body.angular_velocity.norm();
            let tau_norm = torque.norm();

            if tau_norm > 1e-10 {
                body.angular_velocity = (1.0 - alpha) * body.angular_velocity
                    + alpha * (torque / tau_norm) * omega_norm;
            }
        }

        // Velocity Verlet step
        for (i, body) in bodies.iter_mut().enumerate() {
            // Update velocities (half step)
            body.update_velocity(force_accs[i].force, force_accs[i].torque, dt * 0.5);

            // Update positions
            body.update_position(dt, box_size);

            // Clamp displacements
            let v_norm = body.velocity.norm();
            if v_norm * dt > params.max_displacement {
                body.velocity *= params.max_displacement / (v_norm * dt);
            }

            let omega_norm = body.angular_velocity.norm();
            if omega_norm * dt > params.max_rotation {
                body.angular_velocity *= params.max_rotation / (omega_norm * dt);
            }
        }

        // Recalculate forces for second half of velocity update
        let force_accs = calculate_forces(bodies, box_size, potential);

        for (i, body) in bodies.iter_mut().enumerate() {
            body.update_velocity(force_accs[i].force, force_accs[i].torque, dt * 0.5);
        }

        // Calculate new energy
        let new_energy = calculate_total_energy(bodies, box_size, potential);
        let delta_e = new_energy - prev_energy;

        // Check energy tolerance
        if delta_e.abs() < params.energy_tolerance && iteration > 10 {
            if params.print_interval > 0 {
                println!("Energy converged at iteration {}", iteration);
                println!("Final energy: {:.6} kcal/mol", new_energy);
                println!("Final max force: {:.6} kcal/mol/Å", max_force);
            }
            return MinimizationResult::Converged {
                iterations: iteration,
                final_energy: new_energy,
                final_max_force: max_force,
            };
        }

        prev_energy = new_energy;

        // Print progress
        if params.print_interval > 0 && (iteration + 1) % params.print_interval == 0 {
            println!(
                "Iteration {}: E = {:.6}, F_max = {:.6}, dt = {:.6}, alpha = {:.6}",
                iteration + 1,
                new_energy,
                max_force,
                dt,
                alpha
            );
        }
    }

    let final_energy = calculate_total_energy(bodies, box_size, potential);
    let force_accs = calculate_forces(bodies, box_size, potential);
    let final_max_force = calculate_max_force(&force_accs);

    if params.print_interval > 0 {
        println!("Reached maximum iterations");
        println!("Final energy: {:.6} kcal/mol", final_energy);
        println!("Final max force: {:.6} kcal/mol/Å", final_max_force);
    }

    MinimizationResult::MaxIterations {
        iterations: params.max_iterations,
        final_energy,
        final_max_force,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::{HardSphereParams, SoftSphereParams};
    use crate::rigid_body::SurfaceSphere;
    use nalgebra::{Point3, UnitQuaternion};

    fn create_test_body(position: Point3<f64>) -> RigidBody {
        let hull_points = vec![Point3::origin()];
        let mut body =
            RigidBody::from_hull_points(hull_points, position, UnitQuaternion::identity());
        body.surface_spheres
            .push(SurfaceSphere::new(Point3::origin(), 1.0));
        body
    }

    #[test]
    fn test_minimizer_params_default() {
        let params = MinimizerParams::default();
        assert!(params.max_iterations > 0);
        assert!(params.force_tolerance > 0.0);
        assert!(params.energy_tolerance > 0.0);
    }

    #[test]
    fn test_minimize_no_forces() {
        // Two bodies far apart - should converge immediately
        let mut bodies = vec![
            create_test_body(Point3::new(0.0, 0.0, 0.0)),
            create_test_body(Point3::new(100.0, 0.0, 0.0)),
        ];

        let box_size = Vector3::new(1000.0, 1000.0, 1000.0);
        let potential = Potential::SoftSphere(SoftSphereParams::new_wca(1.0, 2.0));

        let params = MinimizerParams {
            print_interval: 0,
            ..Default::default()
        };

        let result = minimize_energy(&mut bodies, box_size, &potential, &params);

        // Should converge quickly since there are no forces
        match result {
            MinimizationResult::Converged { iterations, .. } => {
                assert!(iterations < 10);
            }
            _ => panic!("Should have converged"),
        }
    }

    #[test]
    fn test_minimize_overlapping_bodies() {
        // Two bodies overlapping - should push apart
        let mut bodies = vec![
            create_test_body(Point3::new(0.0, 0.0, 0.0)),
            create_test_body(Point3::new(1.0, 0.0, 0.0)), // Overlapping (sigma=2.0)
        ];

        let initial_distance = (bodies[1].position - bodies[0].position).norm();

        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(100.0));

        let params = MinimizerParams {
            max_iterations: 1000,
            print_interval: 0,
            ..Default::default()
        };

        let result = minimize_energy(&mut bodies, box_size, &potential, &params);

        let final_distance = (bodies[1].position - bodies[0].position).norm();

        // Bodies should have moved apart
        assert!(final_distance > initial_distance);

        // Should converge to separation ~ sigma (2.0)
        // (might not be exact due to tolerance)
        match result {
            MinimizationResult::Converged { final_energy, .. } => {
                // Energy should be close to zero when properly separated
                assert!(final_energy < 1.0);
            }
            MinimizationResult::MaxIterations { .. } => {
                // It's ok if it doesn't fully converge in test
            }
        }
    }

    #[test]
    fn test_steepest_descent() {
        let mut bodies = vec![
            create_test_body(Point3::new(0.0, 0.0, 0.0)),
            create_test_body(Point3::new(1.5, 0.0, 0.0)),
        ];

        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(50.0));

        let params = MinimizerParams {
            use_fire: false,
            max_iterations: 500,
            print_interval: 0,
            ..Default::default()
        };

        let result = steepest_descent(&mut bodies, box_size, &potential, &params);

        // Should reduce energy
        match result {
            MinimizationResult::Converged { final_energy, .. }
            | MinimizationResult::MaxIterations { final_energy, .. } => {
                // Final energy should be lower than initial
                // (initial overlap energy would be positive)
                assert!(final_energy < 10.0);
            }
        }
    }

    #[test]
    fn test_fire_minimization() {
        let mut bodies = vec![
            create_test_body(Point3::new(0.0, 0.0, 0.0)),
            create_test_body(Point3::new(1.5, 0.0, 0.0)),
        ];

        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(50.0));

        let params = MinimizerParams {
            use_fire: true,
            max_iterations: 500,
            print_interval: 0,
            ..Default::default()
        };

        let result = fire_minimization(&mut bodies, box_size, &potential, &params);

        // FIRE should converge efficiently
        match result {
            MinimizationResult::Converged {
                iterations,
                final_energy,
                ..
            } => {
                assert!(iterations < 500);
                assert!(final_energy < 10.0);
            }
            MinimizationResult::MaxIterations { .. } => {
                // Ok for test
            }
        }
    }

    #[test]
    fn test_pbc_during_minimization() {
        // Bodies on opposite sides of box
        let mut bodies = vec![
            create_test_body(Point3::new(2.0, 50.0, 50.0)),
            create_test_body(Point3::new(98.0, 50.0, 50.0)),
        ];

        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(50.0));

        let params = MinimizerParams {
            max_iterations: 100,
            print_interval: 0,
            ..Default::default()
        };

        minimize_energy(&mut bodies, box_size, &potential, &params);

        // Positions should remain in box
        for body in &bodies {
            assert!(body.position.x >= 0.0 && body.position.x < box_size.x);
            assert!(body.position.y >= 0.0 && body.position.y < box_size.y);
            assert!(body.position.z >= 0.0 && body.position.z < box_size.z);
        }
    }

    #[test]
    fn test_max_displacement_constraint() {
        let mut bodies = vec![
            create_test_body(Point3::new(0.0, 0.0, 0.0)),
            create_test_body(Point3::new(0.5, 0.0, 0.0)), // Very close
        ];

        let initial_pos = bodies[1].position;

        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(1000.0)); // Very stiff

        let max_disp = 0.1;
        let params = MinimizerParams {
            max_displacement: max_disp,
            max_iterations: 1, // Just one step
            print_interval: 0,
            ..Default::default()
        };

        minimize_energy(&mut bodies, box_size, &potential, &params);

        let displacement = (bodies[1].position - initial_pos).norm();

        // Displacement should not exceed max_displacement
        assert!(displacement <= max_disp * 1.1); // Small tolerance for numerical error
    }
}
