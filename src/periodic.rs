use nalgebra::{Point3, Vector3};

/// Wrap a position into the primary periodic box [0, box_size)
///
/// # Arguments
/// * `pos` - Position that may be outside the box
/// * `box_size` - Dimensions of the periodic box
///
/// # Returns
/// Position wrapped into [0, box_size) in each dimension
pub fn wrap_position(pos: Point3<f64>, box_size: Vector3<f64>) -> Point3<f64> {
    Point3::new(
        pos.x - (pos.x / box_size.x).floor() * box_size.x,
        pos.y - (pos.y / box_size.y).floor() * box_size.y,
        pos.z - (pos.z / box_size.z).floor() * box_size.z,
    )
}

/// Apply minimum image convention to a displacement vector
/// Returns the shortest vector between two points accounting for periodic boundaries
///
/// # Arguments
/// * `r` - Displacement vector (may cross periodic boundary)
/// * `box_size` - Dimensions of the periodic box
///
/// # Returns
/// Shortest displacement vector accounting for PBC
pub fn minimum_image_vector(r: Vector3<f64>, box_size: Vector3<f64>) -> Vector3<f64> {
    Vector3::new(
        r.x - (r.x / box_size.x).round() * box_size.x,
        r.y - (r.y / box_size.y).round() * box_size.y,
        r.z - (r.z / box_size.z).round() * box_size.z,
    )
}

/// Calculate minimum image distance between two points
///
/// # Arguments
/// * `pos1` - First position
/// * `pos2` - Second position
/// * `box_size` - Dimensions of the periodic box
///
/// # Returns
/// Tuple of (distance, minimum image vector from pos1 to pos2)
pub fn minimum_image_distance(
    pos1: Point3<f64>,
    pos2: Point3<f64>,
    box_size: Vector3<f64>,
) -> (f64, Vector3<f64>) {
    let r = pos2 - pos1;
    let r_min = minimum_image_vector(r, box_size);
    let distance = r_min.norm();
    (distance, r_min)
}

/// Check if a distance is within a cutoff using minimum image convention
pub fn within_cutoff(
    pos1: Point3<f64>,
    pos2: Point3<f64>,
    box_size: Vector3<f64>,
    cutoff: f64,
) -> bool {
    let (distance, _) = minimum_image_distance(pos1, pos2, box_size);
    distance < cutoff
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_wrap_position_inside_box() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let pos = Point3::new(50.0, 50.0, 50.0);
        let wrapped = wrap_position(pos, box_size);

        assert_abs_diff_eq!(wrapped.x, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.y, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.z, 50.0, epsilon = 1e-10);
    }

    #[test]
    fn test_wrap_position_outside_box() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // Test positive overflow
        let pos = Point3::new(150.0, 250.0, 99.0);
        let wrapped = wrap_position(pos, box_size);
        assert_abs_diff_eq!(wrapped.x, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.y, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.z, 99.0, epsilon = 1e-10);

        // Test negative overflow
        let pos = Point3::new(-10.0, -150.0, 1.0);
        let wrapped = wrap_position(pos, box_size);
        assert_abs_diff_eq!(wrapped.x, 90.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.y, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.z, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_wrap_position_exact_boundary() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // At zero boundary
        let pos = Point3::new(0.0, 0.0, 0.0);
        let wrapped = wrap_position(pos, box_size);
        assert_abs_diff_eq!(wrapped.x, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.z, 0.0, epsilon = 1e-10);

        // At upper boundary (should wrap to 0)
        let pos = Point3::new(100.0, 100.0, 100.0);
        let wrapped = wrap_position(pos, box_size);
        assert_abs_diff_eq!(wrapped.x, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_minimum_image_vector_no_wrapping() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let r = Vector3::new(10.0, 20.0, 30.0);
        let r_min = minimum_image_vector(r, box_size);

        assert_abs_diff_eq!(r_min.x, 10.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_min.y, 20.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_min.z, 30.0, epsilon = 1e-10);
    }

    #[test]
    fn test_minimum_image_vector_with_wrapping() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // Vector going the "long way" - should wrap to shorter path
        let r = Vector3::new(80.0, 0.0, 0.0); // 80 Å to the right
        let r_min = minimum_image_vector(r, box_size);
        assert_abs_diff_eq!(r_min.x, -20.0, epsilon = 1e-10); // Should be 20 Å to the left
        assert_abs_diff_eq!(r_min.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_min.z, 0.0, epsilon = 1e-10);

        // Negative direction
        let r = Vector3::new(-80.0, 0.0, 0.0); // 80 Å to the left
        let r_min = minimum_image_vector(r, box_size);
        assert_abs_diff_eq!(r_min.x, 20.0, epsilon = 1e-10); // Should be 20 Å to the right
        assert_abs_diff_eq!(r_min.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_min.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_minimum_image_vector_exactly_half() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // Exactly at half box length - should stay same (round chooses closer to zero)
        let r = Vector3::new(50.0, 0.0, 0.0);
        let r_min = minimum_image_vector(r, box_size);
        assert_abs_diff_eq!(r_min.x, 50.0, epsilon = 1e-10);

        // Slightly over half - should wrap
        let r = Vector3::new(50.1, 0.0, 0.0);
        let r_min = minimum_image_vector(r, box_size);
        assert_abs_diff_eq!(r_min.x, -49.9, epsilon = 1e-10);
    }

    #[test]
    fn test_minimum_image_distance() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // Two points close together
        let pos1 = Point3::new(10.0, 10.0, 10.0);
        let pos2 = Point3::new(15.0, 10.0, 10.0);
        let (dist, r_vec) = minimum_image_distance(pos1, pos2, box_size);

        assert_abs_diff_eq!(dist, 5.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.x, 5.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_minimum_image_distance_across_boundary() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // Two points on opposite sides of boundary
        let pos1 = Point3::new(5.0, 50.0, 50.0);
        let pos2 = Point3::new(95.0, 50.0, 50.0);
        let (dist, r_vec) = minimum_image_distance(pos1, pos2, box_size);

        // Direct distance would be 90, but through PBC it's 10
        assert_abs_diff_eq!(dist, 10.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.x, -10.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.y, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_minimum_image_distance_3d() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);

        // 3D case with wrapping in all dimensions
        let pos1 = Point3::new(5.0, 5.0, 5.0);
        let pos2 = Point3::new(95.0, 95.0, 95.0);
        let (dist, r_vec) = minimum_image_distance(pos1, pos2, box_size);

        // Each component should wrap to -10
        assert_abs_diff_eq!(r_vec.x, -10.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.y, -10.0, epsilon = 1e-10);
        assert_abs_diff_eq!(r_vec.z, -10.0, epsilon = 1e-10);

        // Distance should be sqrt(10^2 + 10^2 + 10^2) = sqrt(300)
        assert_abs_diff_eq!(dist, (300.0_f64).sqrt(), epsilon = 1e-10);
    }

    #[test]
    fn test_within_cutoff() {
        let box_size = Vector3::new(100.0, 100.0, 100.0);
        let cutoff = 15.0;

        // Within cutoff
        let pos1 = Point3::new(10.0, 10.0, 10.0);
        let pos2 = Point3::new(20.0, 10.0, 10.0);
        assert!(within_cutoff(pos1, pos2, box_size, cutoff));

        // Outside cutoff
        let pos2 = Point3::new(30.0, 10.0, 10.0);
        assert!(!within_cutoff(pos1, pos2, box_size, cutoff));

        // Within cutoff through PBC
        let pos1 = Point3::new(5.0, 50.0, 50.0);
        let pos2 = Point3::new(95.0, 50.0, 50.0);
        assert!(within_cutoff(pos1, pos2, box_size, cutoff)); // Distance is 10 through PBC
    }

    #[test]
    fn test_non_cubic_box() {
        let box_size = Vector3::new(100.0, 50.0, 200.0);

        // Test wrapping in non-cubic box
        let pos = Point3::new(150.0, 75.0, 350.0);
        let wrapped = wrap_position(pos, box_size);

        assert_abs_diff_eq!(wrapped.x, 50.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.y, 25.0, epsilon = 1e-10);
        assert_abs_diff_eq!(wrapped.z, 150.0, epsilon = 1e-10);

        // Test minimum image in non-cubic box
        let pos1 = Point3::new(5.0, 5.0, 10.0);
        let pos2 = Point3::new(95.0, 45.0, 190.0);
        let (dist, r_vec) = minimum_image_distance(pos1, pos2, box_size);

        assert_abs_diff_eq!(r_vec.x, -10.0, epsilon = 1e-10); // Wraps in x (100 box)
        assert_abs_diff_eq!(r_vec.y, -10.0, epsilon = 1e-10); // Wraps in y (50 box)
        assert_abs_diff_eq!(r_vec.z, -20.0, epsilon = 1e-10); // Wraps in z (200 box)
    }
}
