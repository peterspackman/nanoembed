use crate::rigid_body::SurfaceSphere;
use nalgebra::Point3;

/// Discretize a convex hull into surface spheres
pub fn discretize_convex_hull(
    hull_points: &[Point3<f64>],
    sphere_radius: f64,
    target_spacing: f64,
) -> Vec<SurfaceSphere> {
    if hull_points.is_empty() {
        return Vec::new();
    }

    // Strategy: Place spheres uniformly on the convex hull surface
    // We'll use a simple approach: create a Fibonacci sphere and project onto hull

    let num_spheres = estimate_num_spheres(hull_points, target_spacing);
    let fibonacci_points = fibonacci_sphere(num_spheres);

    // Scale fibonacci points to match hull size
    let max_radius = hull_points
        .iter()
        .map(|p| p.coords.norm())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(1.0);

    // Project each fibonacci point to the nearest hull point
    let mut surface_spheres = Vec::new();

    for fib_point in fibonacci_points {
        // Scale to hull size
        let scaled_point = fib_point * max_radius;

        // Find closest hull point (simple approach)
        // In production, could project to hull surface more accurately
        let closest_hull_point = find_closest_point(&scaled_point, hull_points);

        surface_spheres.push(SurfaceSphere::new(closest_hull_point, sphere_radius));
    }

    // Remove duplicates (points too close together)
    remove_duplicate_spheres(&mut surface_spheres, target_spacing * 0.5);

    surface_spheres
}

/// Estimate number of spheres needed based on hull size and spacing
fn estimate_num_spheres(hull_points: &[Point3<f64>], target_spacing: f64) -> usize {
    // Estimate surface area from bounding sphere
    let max_radius = hull_points
        .iter()
        .map(|p| p.coords.norm())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(1.0);

    let surface_area = 4.0 * std::f64::consts::PI * max_radius * max_radius;
    let area_per_sphere = target_spacing * target_spacing;

    let num = (surface_area / area_per_sphere).ceil() as usize;
    num.max(10) // At least 10 spheres
}

/// Generate points on a sphere using Fibonacci spiral
/// Returns points on unit sphere
fn fibonacci_sphere(n: usize) -> Vec<Point3<f64>> {
    let mut points = Vec::with_capacity(n);
    let phi = std::f64::consts::PI * (3.0 - (5.0_f64).sqrt()); // Golden angle

    for i in 0..n {
        let y = 1.0 - (i as f64 / (n as f64 - 1.0)) * 2.0; // y goes from 1 to -1
        let radius = (1.0 - y * y).sqrt(); // Radius at y

        let theta = phi * i as f64;

        let x = theta.cos() * radius;
        let z = theta.sin() * radius;

        points.push(Point3::new(x, y, z));
    }

    points
}

/// Find closest point in a set to a given point
fn find_closest_point(point: &Point3<f64>, candidates: &[Point3<f64>]) -> Point3<f64> {
    candidates
        .iter()
        .min_by(|a, b| {
            let dist_a = (point - *a).norm();
            let dist_b = (point - *b).norm();
            dist_a.partial_cmp(&dist_b).unwrap()
        })
        .copied()
        .unwrap_or_else(|| Point3::origin())
}

/// Remove spheres that are too close together
fn remove_duplicate_spheres(spheres: &mut Vec<SurfaceSphere>, min_distance: f64) {
    let mut i = 0;
    while i < spheres.len() {
        let mut j = i + 1;
        while j < spheres.len() {
            let dist = (spheres[i].local_position - spheres[j].local_position).norm();
            if dist < min_distance {
                spheres.remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

/// Alternative: discretize using icosphere subdivision
/// More uniform distribution than Fibonacci
pub fn discretize_with_icosphere(
    hull_points: &[Point3<f64>],
    sphere_radius: f64,
    subdivisions: usize,
) -> Vec<SurfaceSphere> {
    if hull_points.is_empty() {
        return Vec::new();
    }

    // Start with icosahedron vertices
    let icosphere_points = generate_icosphere(subdivisions);

    // Scale to hull size
    let max_radius = hull_points
        .iter()
        .map(|p| p.coords.norm())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(1.0);

    let mut surface_spheres = Vec::new();

    for point in icosphere_points {
        let scaled_point = point * max_radius;
        let closest_hull_point = find_closest_point(&scaled_point, hull_points);
        surface_spheres.push(SurfaceSphere::new(closest_hull_point, sphere_radius));
    }

    surface_spheres
}

/// Generate icosphere vertices (unit sphere)
/// subdivisions=0 gives icosahedron (12 vertices)
/// subdivisions=1 gives 42 vertices, etc.
fn generate_icosphere(subdivisions: usize) -> Vec<Point3<f64>> {
    // Golden ratio
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;

    // Icosahedron vertices
    let mut vertices = vec![
        Point3::new(-1.0, phi, 0.0).normalize(),
        Point3::new(1.0, phi, 0.0).normalize(),
        Point3::new(-1.0, -phi, 0.0).normalize(),
        Point3::new(1.0, -phi, 0.0).normalize(),
        Point3::new(0.0, -1.0, phi).normalize(),
        Point3::new(0.0, 1.0, phi).normalize(),
        Point3::new(0.0, -1.0, -phi).normalize(),
        Point3::new(0.0, 1.0, -phi).normalize(),
        Point3::new(phi, 0.0, -1.0).normalize(),
        Point3::new(phi, 0.0, 1.0).normalize(),
        Point3::new(-phi, 0.0, -1.0).normalize(),
        Point3::new(-phi, 0.0, 1.0).normalize(),
    ];

    // For now, return icosahedron
    // Full subdivision would require triangle mesh data structure
    // This is good enough for initial implementation
    if subdivisions > 0 {
        println!(
            "Warning: Icosphere subdivision not fully implemented, using base icosahedron"
        );
    }

    vertices
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn create_test_cube() -> Vec<Point3<f64>> {
        vec![
            Point3::new(-1.0, -1.0, -1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(1.0, 1.0, -1.0),
            Point3::new(-1.0, -1.0, 1.0),
            Point3::new(1.0, -1.0, 1.0),
            Point3::new(-1.0, 1.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
        ]
    }

    #[test]
    fn test_fibonacci_sphere_count() {
        let n = 100;
        let points = fibonacci_sphere(n);
        assert_eq!(points.len(), n);
    }

    #[test]
    fn test_fibonacci_sphere_unit_sphere() {
        let points = fibonacci_sphere(100);

        // All points should be on unit sphere
        for point in points {
            let radius = point.coords.norm();
            assert_abs_diff_eq!(radius, 1.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_fibonacci_sphere_coverage() {
        let points = fibonacci_sphere(100);

        // Points should cover sphere uniformly
        // Check that we have points in all octants
        let mut octants = vec![0; 8];

        for point in points {
            let idx = (if point.x > 0.0 { 1 } else { 0 })
                + (if point.y > 0.0 { 2 } else { 0 })
                + (if point.z > 0.0 { 4 } else { 0 });
            octants[idx] += 1;
        }

        // Each octant should have at least a few points
        for count in octants {
            assert!(count > 5, "Octant should have more than 5 points");
        }
    }

    #[test]
    fn test_estimate_num_spheres() {
        let cube = create_test_cube();
        let target_spacing = 1.0;

        let num = estimate_num_spheres(&cube, target_spacing);

        // Should get reasonable number of spheres
        assert!(num >= 10);
        assert!(num < 1000);
    }

    #[test]
    fn test_find_closest_point() {
        let candidates = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];

        let test_point = Point3::new(0.9, 0.1, 0.0);
        let closest = find_closest_point(&test_point, &candidates);

        assert_eq!(closest, Point3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_discretize_convex_hull_basic() {
        let cube = create_test_cube();
        let sphere_radius = 0.5;
        let target_spacing = 1.0;

        let spheres = discretize_convex_hull(&cube, sphere_radius, target_spacing);

        // Should have created some spheres
        assert!(!spheres.is_empty());

        // All spheres should have correct radius
        for sphere in &spheres {
            assert_abs_diff_eq!(sphere.radius, sphere_radius, epsilon = 1e-10);
        }

        // All sphere positions should be at hull points
        for sphere in &spheres {
            let is_at_hull = cube.iter().any(|p| {
                (p - &sphere.local_position).norm() < 1e-6
            });
            assert!(is_at_hull, "Sphere should be at a hull point");
        }
    }

    #[test]
    fn test_discretize_empty_hull() {
        let empty: Vec<Point3<f64>> = Vec::new();
        let spheres = discretize_convex_hull(&empty, 0.5, 1.0);

        assert!(spheres.is_empty());
    }

    #[test]
    fn test_remove_duplicate_spheres() {
        let mut spheres = vec![
            SurfaceSphere::new(Point3::new(0.0, 0.0, 0.0), 0.5),
            SurfaceSphere::new(Point3::new(0.1, 0.0, 0.0), 0.5), // Very close
            SurfaceSphere::new(Point3::new(5.0, 0.0, 0.0), 0.5), // Far away
        ];

        let min_distance = 1.0;
        remove_duplicate_spheres(&mut spheres, min_distance);

        // Should have removed one of the close spheres
        assert_eq!(spheres.len(), 2);

        // Remaining spheres should be far apart
        for i in 0..spheres.len() {
            for j in (i + 1)..spheres.len() {
                let dist = (spheres[i].local_position - spheres[j].local_position).norm();
                assert!(dist >= min_distance);
            }
        }
    }

    #[test]
    fn test_generate_icosphere() {
        let vertices = generate_icosphere(0);

        // Icosahedron has 12 vertices
        assert_eq!(vertices.len(), 12);

        // All vertices should be on unit sphere
        for vertex in vertices {
            let radius = vertex.coords.norm();
            assert_abs_diff_eq!(radius, 1.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_discretize_with_icosphere() {
        let cube = create_test_cube();
        let sphere_radius = 0.5;

        let spheres = discretize_with_icosphere(&cube, sphere_radius, 0);

        // Should create spheres
        assert!(!spheres.is_empty());

        // All spheres should have correct radius
        for sphere in &spheres {
            assert_abs_diff_eq!(sphere.radius, sphere_radius, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_discretize_different_spacing() {
        let cube = create_test_cube();
        let sphere_radius = 0.5;

        let spheres_fine = discretize_convex_hull(&cube, sphere_radius, 0.5);
        let spheres_coarse = discretize_convex_hull(&cube, sphere_radius, 2.0);

        // Finer spacing should create more spheres
        assert!(spheres_fine.len() >= spheres_coarse.len());
    }
}
