use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

/// Surface sphere for collision detection and force calculation
/// Represents a sphere attached to the rigid body surface
#[derive(Debug, Clone)]
pub struct SurfaceSphere {
    /// Position in body frame (relative to center of mass)
    pub local_position: Point3<f64>,
    /// Radius of the sphere
    pub radius: f64,
}

impl SurfaceSphere {
    pub fn new(local_position: Point3<f64>, radius: f64) -> Self {
        Self {
            local_position,
            radius,
        }
    }
}

/// Rigid body representing a nanoparticle
#[derive(Debug, Clone)]
pub struct RigidBody {
    // Translational state
    pub position: Point3<f64>,
    pub velocity: Vector3<f64>,

    // Rotational state
    pub orientation: UnitQuaternion<f64>,
    pub angular_velocity: Vector3<f64>,

    // Physical properties
    pub mass: f64,
    pub inertia_tensor: Matrix3<f64>, // In body frame

    // Geometry
    pub surface_spheres: Vec<SurfaceSphere>,

    // Original geometry (in body frame, for reference)
    pub hull_points: Vec<Point3<f64>>,
}

impl RigidBody {
    /// Create a new rigid body from convex hull points
    pub fn from_hull_points(
        hull_points: Vec<Point3<f64>>,
        position: Point3<f64>,
        orientation: UnitQuaternion<f64>,
    ) -> Self {
        let mass = Self::compute_mass(&hull_points);
        let inertia_tensor = Self::compute_inertia_tensor(&hull_points, mass);

        Self {
            position,
            velocity: Vector3::zeros(),
            orientation,
            angular_velocity: Vector3::zeros(),
            mass,
            inertia_tensor,
            surface_spheres: Vec::new(),
            hull_points,
        }
    }

    /// Compute mass assuming uniform density
    /// For now, use simple approximation based on bounding sphere
    fn compute_mass(points: &[Point3<f64>]) -> f64 {
        if points.is_empty() {
            return 1.0;
        }

        // Compute bounding sphere radius
        let max_r2 = points
            .iter()
            .map(|p| p.coords.norm_squared())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(1.0);

        let radius = max_r2.sqrt();

        // Mass proportional to volume (assuming unit density)
        let volume = (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        volume
    }

    /// Compute inertia tensor for a collection of points
    /// Assumes points represent a uniform density object
    fn compute_inertia_tensor(points: &[Point3<f64>], total_mass: f64) -> Matrix3<f64> {
        if points.is_empty() {
            return Matrix3::identity();
        }

        let point_mass = total_mass / points.len() as f64;

        let mut ixx = 0.0;
        let mut iyy = 0.0;
        let mut izz = 0.0;
        let mut ixy = 0.0;
        let mut ixz = 0.0;
        let mut iyz = 0.0;

        for p in points {
            let x = p.x;
            let y = p.y;
            let z = p.z;

            ixx += point_mass * (y * y + z * z);
            iyy += point_mass * (x * x + z * z);
            izz += point_mass * (x * x + y * y);

            ixy -= point_mass * x * y;
            ixz -= point_mass * x * z;
            iyz -= point_mass * y * z;
        }

        Matrix3::new(
            ixx, ixy, ixz,
            ixy, iyy, iyz,
            ixz, iyz, izz,
        )
    }

    /// Get global position of a surface sphere
    pub fn get_sphere_global_position(&self, sphere: &SurfaceSphere) -> Point3<f64> {
        self.position + self.orientation * sphere.local_position.coords
    }

    /// Get all surface sphere positions in global frame
    pub fn get_all_sphere_positions(&self) -> Vec<Point3<f64>> {
        self.surface_spheres
            .iter()
            .map(|s| self.get_sphere_global_position(s))
            .collect()
    }

    /// Apply a force at a point (in global frame)
    /// Returns (force, torque) to be applied to the body
    pub fn apply_force_at_point(
        &self,
        force: Vector3<f64>,
        point: Point3<f64>,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let r = point - self.position;
        let torque = r.cross(&force);
        (force, torque)
    }

    /// Update position and orientation (for integration)
    pub fn update_position(&mut self, dt: f64, box_size: Vector3<f64>) {
        // Update position
        self.position += self.velocity * dt;

        // Wrap position with PBC
        self.position = crate::periodic::wrap_position(self.position, box_size);

        // Update orientation using angular velocity
        let angle = self.angular_velocity.norm() * dt;
        if angle > 1e-10 {
            let axis = nalgebra::Unit::new_normalize(self.angular_velocity);
            let delta_q = UnitQuaternion::from_axis_angle(&axis, angle);
            self.orientation = delta_q * self.orientation;
            self.orientation = self.orientation.renormalize(); // Prevent drift
        }
    }

    /// Update velocities given forces and torques
    pub fn update_velocity(&mut self, force: Vector3<f64>, torque: Vector3<f64>, dt: f64) {
        // Linear acceleration
        let acceleration = force / self.mass;
        self.velocity += acceleration * dt;

        // Angular acceleration
        // α = I^(-1) * τ (in world frame, need to rotate to body frame)
        // For simplicity, assume diagonal inertia tensor aligned with principal axes
        let inertia_inv = self.inertia_tensor.try_inverse()
            .unwrap_or_else(|| Matrix3::identity());
        let angular_acceleration = inertia_inv * torque;
        self.angular_velocity += angular_acceleration * dt;
    }

    /// Apply damping to velocities
    pub fn apply_damping(&mut self, damping: f64) {
        self.velocity *= damping;
        self.angular_velocity *= damping;
    }

    /// Get kinetic energy
    pub fn kinetic_energy(&self) -> f64 {
        let translational = 0.5 * self.mass * self.velocity.norm_squared();
        let rotational = 0.5 * self.angular_velocity.dot(&(self.inertia_tensor * self.angular_velocity));
        translational + rotational
    }

    /// Get total momentum
    pub fn linear_momentum(&self) -> Vector3<f64> {
        self.velocity * self.mass
    }

    /// Get angular momentum
    pub fn angular_momentum(&self) -> Vector3<f64> {
        self.inertia_tensor * self.angular_velocity
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn create_test_cube_points() -> Vec<Point3<f64>> {
        // Unit cube centered at origin
        vec![
            Point3::new(-0.5, -0.5, -0.5),
            Point3::new(0.5, -0.5, -0.5),
            Point3::new(-0.5, 0.5, -0.5),
            Point3::new(0.5, 0.5, -0.5),
            Point3::new(-0.5, -0.5, 0.5),
            Point3::new(0.5, -0.5, 0.5),
            Point3::new(-0.5, 0.5, 0.5),
            Point3::new(0.5, 0.5, 0.5),
        ]
    }

    #[test]
    fn test_rigid_body_creation() {
        let points = create_test_cube_points();
        let pos = Point3::new(10.0, 20.0, 30.0);
        let orientation = UnitQuaternion::identity();

        let body = RigidBody::from_hull_points(points.clone(), pos, orientation);

        assert_eq!(body.position, pos);
        assert_eq!(body.orientation, UnitQuaternion::identity());
        assert_eq!(body.velocity, Vector3::zeros());
        assert_eq!(body.angular_velocity, Vector3::zeros());
        assert_eq!(body.hull_points.len(), points.len());
        assert!(body.mass > 0.0);
    }

    #[test]
    fn test_mass_computation() {
        let points = create_test_cube_points();
        let mass = RigidBody::compute_mass(&points);

        // Mass should be positive and reasonable
        assert!(mass > 0.0);
        assert!(mass < 10.0); // Sanity check for unit cube
    }

    #[test]
    fn test_inertia_tensor_symmetry() {
        let points = create_test_cube_points();
        let mass = RigidBody::compute_mass(&points);
        let inertia = RigidBody::compute_inertia_tensor(&points, mass);

        // Inertia tensor should be symmetric
        assert_abs_diff_eq!(inertia[(0, 1)], inertia[(1, 0)], epsilon = 1e-10);
        assert_abs_diff_eq!(inertia[(0, 2)], inertia[(2, 0)], epsilon = 1e-10);
        assert_abs_diff_eq!(inertia[(1, 2)], inertia[(2, 1)], epsilon = 1e-10);
    }

    #[test]
    fn test_inertia_tensor_cube() {
        let points = create_test_cube_points();
        let mass = RigidBody::compute_mass(&points);
        let inertia = RigidBody::compute_inertia_tensor(&points, mass);

        // For a symmetric cube, diagonal elements should be equal
        // and off-diagonal elements should be zero
        assert_abs_diff_eq!(inertia[(0, 0)], inertia[(1, 1)], epsilon = 1e-10);
        assert_abs_diff_eq!(inertia[(1, 1)], inertia[(2, 2)], epsilon = 1e-10);

        assert_abs_diff_eq!(inertia[(0, 1)], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(inertia[(0, 2)], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(inertia[(1, 2)], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_surface_sphere_global_position() {
        let points = create_test_cube_points();
        let pos = Point3::new(10.0, 20.0, 30.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);

        // Add a surface sphere at local position (1, 0, 0)
        let sphere = SurfaceSphere::new(Point3::new(1.0, 0.0, 0.0), 0.5);
        body.surface_spheres.push(sphere.clone());

        let global_pos = body.get_sphere_global_position(&sphere);

        // With identity rotation, global position should be body position + local position
        assert_abs_diff_eq!(global_pos.x, 11.0, epsilon = 1e-10);
        assert_abs_diff_eq!(global_pos.y, 20.0, epsilon = 1e-10);
        assert_abs_diff_eq!(global_pos.z, 30.0, epsilon = 1e-10);
    }

    #[test]
    fn test_surface_sphere_global_position_rotated() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);

        // 90 degree rotation around z-axis
        let axis = Vector3::z_axis();
        let angle = std::f64::consts::PI / 2.0;
        let orientation = UnitQuaternion::from_axis_angle(&axis, angle);

        let mut body = RigidBody::from_hull_points(points, pos, orientation);

        // Add a surface sphere at local position (1, 0, 0)
        let sphere = SurfaceSphere::new(Point3::new(1.0, 0.0, 0.0), 0.5);
        body.surface_spheres.push(sphere.clone());

        let global_pos = body.get_sphere_global_position(&sphere);

        // After 90° rotation around z, (1,0,0) -> (0,1,0)
        assert_abs_diff_eq!(global_pos.x, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(global_pos.y, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(global_pos.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_apply_force_at_point() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let body = RigidBody::from_hull_points(points, pos, orientation);

        // Apply force at point offset from center
        let force = Vector3::new(1.0, 0.0, 0.0);
        let point = Point3::new(0.0, 1.0, 0.0); // 1 unit in +y from center

        let (returned_force, torque) = body.apply_force_at_point(force, point);

        assert_eq!(returned_force, force);

        // r × F = (0,1,0) × (1,0,0) = (0,0,-1)
        assert_abs_diff_eq!(torque.x, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(torque.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(torque.z, -1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_update_position_no_pbc() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);
        body.velocity = Vector3::new(1.0, 2.0, 3.0);

        let dt = 0.1;
        let box_size = Vector3::new(1000.0, 1000.0, 1000.0);

        body.update_position(dt, box_size);

        assert_abs_diff_eq!(body.position.x, 0.1, epsilon = 1e-10);
        assert_abs_diff_eq!(body.position.y, 0.2, epsilon = 1e-10);
        assert_abs_diff_eq!(body.position.z, 0.3, epsilon = 1e-10);
    }

    #[test]
    fn test_update_position_with_pbc() {
        let points = create_test_cube_points();
        let pos = Point3::new(99.0, 50.0, 50.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);
        body.velocity = Vector3::new(20.0, 0.0, 0.0); // Move 2.0 Å in x

        let dt = 0.1;
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        body.update_position(dt, box_size);

        // Should wrap from 101.0 to 1.0
        assert_abs_diff_eq!(body.position.x, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(body.position.y, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(body.position.z, 50.0, epsilon = 1e-10);
    }

    #[test]
    fn test_update_velocity() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);

        let force = Vector3::new(10.0, 0.0, 0.0);
        let torque = Vector3::new(0.0, 0.0, 1.0);
        let dt = 0.1;

        let initial_velocity = body.velocity;
        let initial_angular_velocity = body.angular_velocity;

        body.update_velocity(force, torque, dt);

        // Linear velocity should increase
        assert!(body.velocity.x > initial_velocity.x);

        // Angular velocity should increase
        assert!(body.angular_velocity.norm() > initial_angular_velocity.norm());
    }

    #[test]
    fn test_damping() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);
        body.velocity = Vector3::new(10.0, 10.0, 10.0);
        body.angular_velocity = Vector3::new(1.0, 1.0, 1.0);

        let damping = 0.9;
        body.apply_damping(damping);

        assert_abs_diff_eq!(body.velocity.x, 9.0, epsilon = 1e-10);
        assert_abs_diff_eq!(body.velocity.y, 9.0, epsilon = 1e-10);
        assert_abs_diff_eq!(body.velocity.z, 9.0, epsilon = 1e-10);

        assert_abs_diff_eq!(body.angular_velocity.x, 0.9, epsilon = 1e-10);
        assert_abs_diff_eq!(body.angular_velocity.y, 0.9, epsilon = 1e-10);
        assert_abs_diff_eq!(body.angular_velocity.z, 0.9, epsilon = 1e-10);
    }

    #[test]
    fn test_kinetic_energy() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);

        // Zero velocity = zero energy
        assert_abs_diff_eq!(body.kinetic_energy(), 0.0, epsilon = 1e-10);

        // Add linear velocity
        body.velocity = Vector3::new(1.0, 0.0, 0.0);
        let ke = body.kinetic_energy();
        assert!(ke > 0.0);

        // KE = 0.5 * m * v^2
        let expected_ke = 0.5 * body.mass * 1.0;
        assert_abs_diff_eq!(ke, expected_ke, epsilon = 1e-10);
    }

    #[test]
    fn test_linear_momentum() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);
        body.velocity = Vector3::new(2.0, 3.0, 4.0);

        let momentum = body.linear_momentum();

        assert_abs_diff_eq!(momentum.x, body.mass * 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(momentum.y, body.mass * 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(momentum.z, body.mass * 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_quaternion_renormalization() {
        let points = create_test_cube_points();
        let pos = Point3::new(0.0, 0.0, 0.0);
        let orientation = UnitQuaternion::identity();

        let mut body = RigidBody::from_hull_points(points, pos, orientation);
        body.angular_velocity = Vector3::new(0.1, 0.2, 0.3);

        let dt = 0.1;
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // Update several times
        for _ in 0..100 {
            body.update_position(dt, box_size);
        }

        // Quaternion should still be normalized
        let quat_norm = body.orientation.norm();
        assert_abs_diff_eq!(quat_norm, 1.0, epsilon = 1e-6);
    }
}
