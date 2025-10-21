/// Potential energy functions for particle interactions

#[derive(Debug, Clone)]
pub enum Potential {
    SoftSphere(SoftSphereParams),
    HardSphere(HardSphereParams),
}

/// Parameters for soft sphere (WCA) potential
/// WCA = Weeks-Chandler-Andersen potential (purely repulsive Lennard-Jones)
#[derive(Debug, Clone)]
pub struct SoftSphereParams {
    /// Energy scale (kcal/mol or kJ/mol)
    pub epsilon: f64,
    /// Length scale (Angstroms) - equilibrium distance
    pub sigma: f64,
    /// Interaction cutoff (typically 2^(1/6) * sigma for WCA)
    pub cutoff: f64,
}

impl SoftSphereParams {
    /// Create WCA parameters with standard cutoff
    pub fn new_wca(epsilon: f64, sigma: f64) -> Self {
        let cutoff = sigma * 2.0_f64.powf(1.0 / 6.0); // ~1.122 * sigma
        Self {
            epsilon,
            sigma,
            cutoff,
        }
    }

    /// Create soft sphere with custom cutoff
    pub fn new_custom(epsilon: f64, sigma: f64, cutoff: f64) -> Self {
        Self {
            epsilon,
            sigma,
            cutoff,
        }
    }
}

/// Parameters for hard sphere potential
/// Uses harmonic repulsion when spheres overlap
#[derive(Debug, Clone)]
pub struct HardSphereParams {
    /// Spring constant for overlap penalty (kcal/mol/Å²)
    pub penalty_strength: f64,
}

impl HardSphereParams {
    pub fn new(penalty_strength: f64) -> Self {
        Self { penalty_strength }
    }
}

/// Compute WCA (Weeks-Chandler-Andersen) potential energy and force
/// This is a purely repulsive potential based on Lennard-Jones
///
/// For r < cutoff:
///   U(r) = 4ε[(σ/r)^12 - (σ/r)^6] + ε
///   F(r) = 24ε/r * [2(σ/r)^12 - (σ/r)^6]
///
/// # Arguments
/// * `r` - Distance between particles
/// * `params` - WCA parameters
///
/// # Returns
/// Tuple of (energy, force_magnitude)
/// Force points from particle 1 to particle 2, so multiply by unit vector
pub fn wca_potential(r: f64, params: &SoftSphereParams) -> (f64, f64) {
    if r < params.cutoff && r > 0.0 {
        let sigma_over_r = params.sigma / r;
        let sr6 = sigma_over_r.powi(6);
        let sr12 = sr6 * sr6;

        let energy = 4.0 * params.epsilon * (sr12 - sr6) + params.epsilon;
        let force_mag = 24.0 * params.epsilon / r * (2.0 * sr12 - sr6);

        (energy, force_mag)
    } else {
        (0.0, 0.0)
    }
}

/// Compute hard sphere potential energy and force using harmonic repulsion
/// Only acts when spheres overlap (r < sigma)
///
/// For r < sigma:
///   overlap = sigma - r
///   U(r) = 0.5 * k * overlap^2
///   F(r) = k * overlap
///
/// # Arguments
/// * `r` - Distance between sphere centers
/// * `sigma` - Sum of sphere radii (contact distance)
/// * `params` - Hard sphere parameters
///
/// # Returns
/// Tuple of (energy, force_magnitude)
pub fn hard_sphere_potential(r: f64, sigma: f64, params: &HardSphereParams) -> (f64, f64) {
    let overlap = sigma - r;
    if overlap > 0.0 {
        let energy = 0.5 * params.penalty_strength * overlap * overlap;
        let force_mag = params.penalty_strength * overlap;
        (energy, force_mag)
    } else {
        (0.0, 0.0)
    }
}

/// Evaluate the configured potential
pub fn evaluate_potential(r: f64, sigma: f64, potential: &Potential) -> (f64, f64) {
    match potential {
        Potential::SoftSphere(params) => wca_potential(r, params),
        Potential::HardSphere(params) => hard_sphere_potential(r, sigma, params),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_wca_at_cutoff() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // At cutoff, energy should be ~zero and force should be ~zero
        let (energy, force) = wca_potential(params.cutoff, &params);

        assert_abs_diff_eq!(energy, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(force, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_wca_beyond_cutoff() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // Beyond cutoff, no interaction
        let (energy, force) = wca_potential(params.cutoff + 1.0, &params);

        assert_eq!(energy, 0.0);
        assert_eq!(force, 0.0);
    }

    #[test]
    fn test_wca_at_minimum() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // At r = sigma, the LJ potential has its minimum
        // For WCA (shifted), at sigma we have specific energy
        let (energy, force) = wca_potential(params.sigma, &params);

        // At r = sigma, force should be zero (it's the minimum before shifting)
        // Energy should be -epsilon + epsilon = 0 due to WCA shift
        assert_abs_diff_eq!(force, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(energy, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_wca_repulsive_at_close_distance() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // Very close distance should give strong repulsion
        let (energy1, force1) = wca_potential(0.5, &params);
        let (energy2, force2) = wca_potential(1.0, &params);

        // Energy and force should be positive (repulsive)
        assert!(energy1 > 0.0);
        assert!(force1 > 0.0);
        assert!(energy2 > 0.0);
        assert!(force2 > 0.0);

        // Closer distance should have stronger interaction
        assert!(energy1 > energy2);
        assert!(force1 > force2);
    }

    #[test]
    fn test_wca_force_direction() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // Within cutoff, force should be positive (repulsive)
        let (_, force) = wca_potential(2.0, &params);
        assert!(force > 0.0, "WCA force should be repulsive (positive)");
    }

    #[test]
    fn test_hard_sphere_no_overlap() {
        let params = HardSphereParams::new(100.0);
        let sigma = 5.0;

        // No overlap - no interaction
        let (energy, force) = hard_sphere_potential(6.0, sigma, &params);

        assert_eq!(energy, 0.0);
        assert_eq!(force, 0.0);
    }

    #[test]
    fn test_hard_sphere_with_overlap() {
        let params = HardSphereParams::new(100.0);
        let sigma = 5.0;

        // Overlap of 1.0 Å
        let r = 4.0;
        let overlap = sigma - r; // 1.0
        let (energy, force) = hard_sphere_potential(r, sigma, &params);

        // Energy = 0.5 * k * overlap^2
        let expected_energy = 0.5 * 100.0 * 1.0 * 1.0;
        assert_abs_diff_eq!(energy, expected_energy, epsilon = 1e-10);

        // Force = k * overlap
        let expected_force = 100.0 * 1.0;
        assert_abs_diff_eq!(force, expected_force, epsilon = 1e-10);
    }

    #[test]
    fn test_hard_sphere_larger_overlap() {
        let params = HardSphereParams::new(100.0);
        let sigma = 5.0;

        // Larger overlap
        let (energy1, force1) = hard_sphere_potential(4.0, sigma, &params); // 1.0 overlap
        let (energy2, force2) = hard_sphere_potential(3.0, sigma, &params); // 2.0 overlap

        // Force should scale linearly with overlap
        assert_abs_diff_eq!(force2 / force1, 2.0, epsilon = 1e-10);

        // Energy should scale quadratically with overlap
        assert_abs_diff_eq!(energy2 / energy1, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hard_sphere_exactly_touching() {
        let params = HardSphereParams::new(100.0);
        let sigma = 5.0;

        // Exactly touching (r = sigma)
        let (energy, force) = hard_sphere_potential(sigma, sigma, &params);

        assert_abs_diff_eq!(energy, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(force, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_evaluate_potential_soft_sphere() {
        let potential = Potential::SoftSphere(SoftSphereParams::new_wca(1.0, 3.0));

        let (energy, force) = evaluate_potential(2.0, 5.0, &potential);

        // Should call WCA potential
        let (expected_energy, expected_force) = wca_potential(2.0, &SoftSphereParams::new_wca(1.0, 3.0));

        assert_abs_diff_eq!(energy, expected_energy, epsilon = 1e-10);
        assert_abs_diff_eq!(force, expected_force, epsilon = 1e-10);
    }

    #[test]
    fn test_evaluate_potential_hard_sphere() {
        let potential = Potential::HardSphere(HardSphereParams::new(100.0));
        let sigma = 5.0;

        let (energy, force) = evaluate_potential(4.0, sigma, &potential);

        // Should call hard sphere potential
        let (expected_energy, expected_force) = hard_sphere_potential(4.0, sigma, &HardSphereParams::new(100.0));

        assert_abs_diff_eq!(energy, expected_energy, epsilon = 1e-10);
        assert_abs_diff_eq!(force, expected_force, epsilon = 1e-10);
    }

    #[test]
    fn test_wca_cutoff_value() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // Cutoff should be 2^(1/6) * sigma ≈ 1.122 * 3.0
        let expected_cutoff = 3.0 * 2.0_f64.powf(1.0 / 6.0);
        assert_abs_diff_eq!(params.cutoff, expected_cutoff, epsilon = 1e-10);
    }

    #[test]
    fn test_zero_distance_handling() {
        let params = SoftSphereParams::new_wca(1.0, 3.0);

        // At r=0, force should be well-defined (not infinity due to cutoff check)
        let (energy, force) = wca_potential(0.0, &params);

        // Should return (0, 0) due to r > 0.0 check
        assert_eq!(energy, 0.0);
        assert_eq!(force, 0.0);
    }

    #[test]
    fn test_penalty_strength_scaling() {
        let params1 = HardSphereParams::new(100.0);
        let params2 = HardSphereParams::new(200.0);
        let sigma = 5.0;
        let r = 4.0;

        let (energy1, force1) = hard_sphere_potential(r, sigma, &params1);
        let (energy2, force2) = hard_sphere_potential(r, sigma, &params2);

        // Doubling penalty strength should double force and energy
        assert_abs_diff_eq!(force2 / force1, 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(energy2 / energy1, 2.0, epsilon = 1e-10);
    }
}
