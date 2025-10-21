use crate::periodic::{minimum_image_distance, minimum_image_vector};
use crate::potentials::{evaluate_potential, Potential};
use crate::rigid_body::RigidBody;
use nalgebra::Vector3;

/// Accumulator for forces, torques, and energy
#[derive(Debug, Clone, Default)]
pub struct ForceAccumulator {
    pub force: Vector3<f64>,
    pub torque: Vector3<f64>,
    pub energy: f64,
}

impl ForceAccumulator {
    pub fn new() -> Self {
        Self {
            force: Vector3::zeros(),
            torque: Vector3::zeros(),
            energy: 0.0,
        }
    }

    pub fn add_force_at_point(
        &mut self,
        force: Vector3<f64>,
        point: nalgebra::Point3<f64>,
        body_position: nalgebra::Point3<f64>,
    ) {
        self.force += force;
        let r = point - body_position;
        self.torque += r.cross(&force);
    }

    pub fn add_energy(&mut self, energy: f64) {
        self.energy += energy;
    }
}

/// Calculate all forces between rigid bodies using PBC
pub fn calculate_forces(
    bodies: &[RigidBody],
    box_size: Vector3<f64>,
    potential: &Potential,
) -> Vec<ForceAccumulator> {
    let mut accumulators: Vec<ForceAccumulator> = vec![ForceAccumulator::new(); bodies.len()];

    // Double loop over all body pairs
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            calculate_pairwise_forces(
                &bodies[i],
                &bodies[j],
                i,
                j,
                box_size,
                potential,
                &mut accumulators,
            );
        }
    }

    accumulators
}

/// Calculate forces between two rigid bodies
fn calculate_pairwise_forces(
    body_i: &RigidBody,
    body_j: &RigidBody,
    idx_i: usize,
    idx_j: usize,
    box_size: Vector3<f64>,
    potential: &Potential,
    accumulators: &mut [ForceAccumulator],
) {
    // Early exit if bodies are far apart (broad phase)
    let (com_distance, _) = minimum_image_distance(body_i.position, body_j.position, box_size);

    // Conservative estimate of maximum interaction distance
    let max_radius_i = estimate_max_radius(&body_i.surface_spheres);
    let max_radius_j = estimate_max_radius(&body_j.surface_spheres);
    let cutoff = max_radius_i + max_radius_j + 10.0; // Add buffer

    if com_distance > cutoff {
        return; // Too far to interact
    }

    // Narrow phase: check all sphere pairs
    for sphere_i in &body_i.surface_spheres {
        let pos_i = body_i.get_sphere_global_position(sphere_i);

        for sphere_j in &body_j.surface_spheres {
            let pos_j = body_j.get_sphere_global_position(sphere_j);

            // Get minimum image distance with PBC
            let (distance, r_vec) = minimum_image_distance(pos_i, pos_j, box_size);

            // Sum of radii
            let sigma = sphere_i.radius + sphere_j.radius;

            // Evaluate potential
            let (energy, force_mag) = evaluate_potential(distance, sigma, potential);

            if force_mag.abs() > 1e-10 {
                // Force direction: from i to j
                let force_dir = if distance > 1e-10 {
                    r_vec.normalize()
                } else {
                    // If spheres exactly overlap, use arbitrary direction
                    Vector3::x()
                };

                let force_vec = force_dir * force_mag;

                // Apply to body i (force points away from j)
                accumulators[idx_i].add_force_at_point(force_vec, pos_i, body_i.position);
                accumulators[idx_i].add_energy(energy * 0.5);

                // Apply to body j (Newton's 3rd law)
                accumulators[idx_j].add_force_at_point(-force_vec, pos_j, body_j.position);
                accumulators[idx_j].add_energy(energy * 0.5);
            }
        }
    }
}

/// Estimate maximum radius of surface spheres from body center
fn estimate_max_radius(spheres: &[crate::rigid_body::SurfaceSphere]) -> f64 {
    spheres
        .iter()
        .map(|s| s.local_position.coords.norm() + s.radius)
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(0.0)
}

/// Get total energy of the system
pub fn calculate_total_energy(
    bodies: &[RigidBody],
    box_size: Vector3<f64>,
    potential: &Potential,
) -> f64 {
    let accumulators = calculate_forces(bodies, box_size, potential);
    accumulators.iter().map(|a| a.energy).sum()
}

/// Get maximum force magnitude on any body
pub fn calculate_max_force(accumulators: &[ForceAccumulator]) -> f64 {
    accumulators
        .iter()
        .map(|a| a.force.norm())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::{HardSphereParams, SoftSphereParams};
    use crate::rigid_body::SurfaceSphere;
    use approx::assert_abs_diff_eq;
    use nalgebra::{Point3, UnitQuaternion};

    fn create_simple_body(position: Point3<f64>) -> RigidBody {
        let hull_points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];

        let mut body =
            RigidBody::from_hull_points(hull_points, position, UnitQuaternion::identity());

        // Add a single surface sphere at origin (in body frame)
        body.surface_spheres
            .push(SurfaceSphere::new(Point3::origin(), 1.0));

        body
    }

    #[test]
    fn test_force_accumulator_default() {
        let acc = ForceAccumulator::new();

        assert_eq!(acc.force, Vector3::zeros());
        assert_eq!(acc.torque, Vector3::zeros());
        assert_eq!(acc.energy, 0.0);
    }

    #[test]
    fn test_force_accumulator_add_force() {
        let mut acc = ForceAccumulator::new();

        let force = Vector3::new(1.0, 2.0, 3.0);
        let point = Point3::new(0.0, 1.0, 0.0);
        let body_pos = Point3::origin();

        acc.add_force_at_point(force, point, body_pos);

        assert_eq!(acc.force, force);

        // Torque = r × F = (0,1,0) × (1,2,3) = (3,0,-1)
        assert_abs_diff_eq!(acc.torque.x, 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(acc.torque.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(acc.torque.z, -1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_force_accumulator_add_energy() {
        let mut acc = ForceAccumulator::new();
        acc.add_energy(10.0);
        assert_eq!(acc.energy, 10.0);

        acc.add_energy(5.0);
        assert_eq!(acc.energy, 15.0);
    }

    #[test]
    fn test_no_forces_when_far_apart() {
        let body1 = create_simple_body(Point3::new(0.0, 0.0, 0.0));
        let body2 = create_simple_body(Point3::new(100.0, 0.0, 0.0));

        let bodies = vec![body1, body2];
        let box_size = Vector3::new(1000.0, 1000.0, 1000.0);
        let potential = Potential::SoftSphere(SoftSphereParams::new_wca(1.0, 2.0));

        let forces = calculate_forces(&bodies, box_size, &potential);

        // Bodies are far apart (100 Å), should have no interaction
        assert_abs_diff_eq!(forces[0].force.norm(), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(forces[1].force.norm(), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(forces[0].energy, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(forces[1].energy, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_repulsive_force_when_overlapping() {
        // Two bodies with overlapping spheres
        let body1 = create_simple_body(Point3::new(0.0, 0.0, 0.0));
        let body2 = create_simple_body(Point3::new(1.5, 0.0, 0.0)); // Overlapping (sigma=2.0, r=1.5)

        let bodies = vec![body1, body2];
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(100.0));

        let forces = calculate_forces(&bodies, box_size, &potential);

        // Should have repulsive force
        assert!(forces[0].force.x > 0.0, "Body 1 should be pushed in +x");
        assert!(forces[1].force.x < 0.0, "Body 2 should be pushed in -x");

        // Forces should be equal and opposite
        assert_abs_diff_eq!(forces[0].force.x, -forces[1].force.x, epsilon = 1e-10);
        assert_abs_diff_eq!(forces[0].force.y, -forces[1].force.y, epsilon = 1e-10);
        assert_abs_diff_eq!(forces[0].force.z, -forces[1].force.z, epsilon = 1e-10);

        // Energy should be positive
        assert!(forces[0].energy > 0.0);
        assert!(forces[1].energy > 0.0);

        // Total energy should be conserved (sum = total)
        let total_energy = forces[0].energy + forces[1].energy;
        assert!(total_energy > 0.0);
    }

    #[test]
    fn test_forces_with_pbc() {
        // Two bodies on opposite sides of periodic box
        let body1 = create_simple_body(Point3::new(2.0, 50.0, 50.0));
        let body2 = create_simple_body(Point3::new(98.0, 50.0, 50.0));

        let bodies = vec![body1, body2];
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // With PBC, distance is 4.0 (through periodic boundary)
        // Spheres should be overlapping (sigma=2.0, r=4.0 → no overlap)
        let potential = Potential::HardSphere(HardSphereParams::new(100.0));

        let forces = calculate_forces(&bodies, box_size, &potential);

        // No overlap, so no force
        assert_abs_diff_eq!(forces[0].force.norm(), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(forces[1].force.norm(), 0.0, epsilon = 1e-10);

        // Now make them closer through PBC
        let body1 = create_simple_body(Point3::new(1.0, 50.0, 50.0));
        let body2 = create_simple_body(Point3::new(99.0, 50.0, 50.0));
        let bodies = vec![body1, body2];

        // Distance through PBC: 2.0, sigma=2.0 → exactly touching
        let forces = calculate_forces(&bodies, box_size, &potential);

        // Exactly touching, should have minimal or no force
        assert!(forces[0].force.norm() < 0.1);
        assert!(forces[1].force.norm() < 0.1);
    }

    #[test]
    fn test_multiple_surface_spheres() {
        let mut body1 = create_simple_body(Point3::new(0.0, 0.0, 0.0));
        let mut body2 = create_simple_body(Point3::new(3.0, 0.0, 0.0));

        // Add more surface spheres
        body1
            .surface_spheres
            .push(SurfaceSphere::new(Point3::new(1.0, 0.0, 0.0), 0.5));
        body2
            .surface_spheres
            .push(SurfaceSphere::new(Point3::new(-1.0, 0.0, 0.0), 0.5));

        let bodies = vec![body1, body2];
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::SoftSphere(SoftSphereParams::new_wca(1.0, 2.0));

        let forces = calculate_forces(&bodies, box_size, &potential);

        // Should have forces due to multiple sphere interactions
        // Can't predict exact values without detailed calculation,
        // but forces should exist
        assert!(forces[0].force.norm() > 0.0 || forces[1].force.norm() > 0.0);
    }

    #[test]
    fn test_calculate_total_energy() {
        let body1 = create_simple_body(Point3::new(0.0, 0.0, 0.0));
        let body2 = create_simple_body(Point3::new(1.5, 0.0, 0.0));

        let bodies = vec![body1, body2];
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(100.0));

        let total_energy = calculate_total_energy(&bodies, box_size, &potential);

        // Should have positive energy due to overlap
        assert!(total_energy > 0.0);
    }

    #[test]
    fn test_calculate_max_force() {
        let accumulators = vec![
            ForceAccumulator {
                force: Vector3::new(1.0, 0.0, 0.0),
                torque: Vector3::zeros(),
                energy: 0.0,
            },
            ForceAccumulator {
                force: Vector3::new(0.0, 5.0, 0.0),
                torque: Vector3::zeros(),
                energy: 0.0,
            },
            ForceAccumulator {
                force: Vector3::new(0.0, 0.0, 2.0),
                torque: Vector3::zeros(),
                energy: 0.0,
            },
        ];

        let max_force = calculate_max_force(&accumulators);

        // Maximum force should be 5.0 (second accumulator)
        assert_abs_diff_eq!(max_force, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_estimate_max_radius() {
        let spheres = vec![
            SurfaceSphere::new(Point3::new(1.0, 0.0, 0.0), 0.5),
            SurfaceSphere::new(Point3::new(0.0, 2.0, 0.0), 0.3),
            SurfaceSphere::new(Point3::new(0.0, 0.0, 1.5), 0.2),
        ];

        let max_radius = estimate_max_radius(&spheres);

        // Second sphere is at distance 2.0 with radius 0.3 = 2.3 total
        assert_abs_diff_eq!(max_radius, 2.3, epsilon = 1e-10);
    }

    #[test]
    fn test_no_self_interaction() {
        let body = create_simple_body(Point3::new(0.0, 0.0, 0.0));

        let bodies = vec![body];
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let potential = Potential::HardSphere(HardSphereParams::new(100.0));

        let forces = calculate_forces(&bodies, box_size, &potential);

        // Single body should not interact with itself
        assert_eq!(forces[0].force, Vector3::zeros());
        assert_eq!(forces[0].torque, Vector3::zeros());
        assert_eq!(forces[0].energy, 0.0);
    }
}
