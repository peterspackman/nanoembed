use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

#[derive(Debug, Deserialize, Serialize)]
pub struct Config {
    pub system: SystemConfig,
    pub background: BackgroundConfig,
    pub nanoparticles: Vec<NanoparticleConfig>,
    #[serde(default)]
    pub analysis: AnalysisConfig,
    #[serde(default)]
    pub placement: PlacementConfig,
    #[serde(default)]
    pub dynamics: DynamicsConfig,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SystemConfig {
    pub box_size: f64,
    #[serde(default = "default_output_file")]
    pub output_file: String,
    #[serde(default = "default_voxel_size")]
    pub voxel_size: f64,
    pub random_seed: Option<u64>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct BackgroundConfig {
    #[serde(default = "default_background_material")]
    pub material: BackgroundMaterial,
    #[serde(default = "default_background_element")]
    pub element: String,
    #[serde(default = "default_density")]
    pub density: f64,
    #[serde(default = "default_separation")]
    pub separation: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct NanoparticleConfig {
    pub file: String,
    pub count: usize,
}

#[derive(Debug, Deserialize, Serialize, Default)]
pub struct AnalysisConfig {
    #[serde(default)]
    pub validate: bool,
    #[serde(default)]
    pub quasi_random: bool,
}

#[derive(Debug, Deserialize, Serialize, Default)]
pub struct PlacementConfig {
    #[serde(default = "default_placement_mode")]
    pub mode: PlacementMode,
    /// Minimum distance between atom centers in overlap mode (Å)
    #[serde(default = "default_min_atom_distance")]
    pub min_atom_distance: f64,
    /// Buffer distance added to collision detection shapes (Å)
    #[serde(default = "default_collision_buffer")]
    pub collision_buffer: f64,
}

#[derive(Debug, Deserialize, Serialize, Clone, Default)]
#[serde(rename_all = "lowercase")]
pub enum PlacementMode {
    #[default]
    Collision,  // No overlaps allowed - strict collision avoidance
    Overlap,    // Allow overlaps - first placed wins in overlap regions
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "lowercase")]
pub enum BackgroundMaterial {
    Liquid,
    Fcc,
    None,
}

/// Configuration for energy minimization and dynamics
#[derive(Debug, Deserialize, Serialize, Default)]
pub struct DynamicsConfig {
    /// Enable dynamics-based relaxation
    #[serde(default)]
    pub enabled: bool,

    /// Maximum number of minimization iterations
    #[serde(default = "default_max_iterations")]
    pub max_iterations: usize,

    /// Force tolerance for convergence (kcal/mol/Å)
    #[serde(default = "default_force_tolerance")]
    pub force_tolerance: f64,

    /// Energy tolerance for convergence (kcal/mol)
    #[serde(default = "default_energy_tolerance")]
    pub energy_tolerance: f64,

    /// Maximum displacement per step (Å)
    #[serde(default = "default_max_displacement")]
    pub max_displacement: f64,

    /// Maximum rotation per step (radians)
    #[serde(default = "default_max_rotation")]
    pub max_rotation: f64,

    /// Use FIRE algorithm (faster) vs steepest descent
    #[serde(default = "default_use_fire")]
    pub use_fire: bool,

    /// Potential type
    #[serde(default = "default_potential_type")]
    pub potential: PotentialType,

    /// Surface sphere discretization
    #[serde(default)]
    pub surface: SurfaceConfig,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "lowercase")]
pub enum PotentialType {
    SoftSphere {
        epsilon: f64,  // Energy scale (kcal/mol)
        sigma: f64,    // Length scale (Å)
    },
    HardSphere {
        penalty_strength: f64,  // Spring constant (kcal/mol/Å²)
    },
}

impl Default for PotentialType {
    fn default() -> Self {
        PotentialType::SoftSphere {
            epsilon: 1.0,
            sigma: 3.0,
        }
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SurfaceConfig {
    /// Radius of surface spheres (Å)
    #[serde(default = "default_sphere_radius")]
    pub sphere_radius: f64,

    /// Target spacing between surface spheres (Å)
    #[serde(default = "default_target_spacing")]
    pub target_spacing: f64,
}

impl Default for SurfaceConfig {
    fn default() -> Self {
        Self {
            sphere_radius: default_sphere_radius(),
            target_spacing: default_target_spacing(),
        }
    }
}

impl Config {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(path)?;
        let config: Config = toml::from_str(&content)?;
        Ok(config)
    }

    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let content = toml::to_string_pretty(self)?;
        fs::write(path, content)?;
        Ok(())
    }

    pub fn total_nanoparticle_count(&self) -> usize {
        self.nanoparticles.iter().map(|np| np.count).sum()
    }
}

// Default value functions
fn default_output_file() -> String {
    "system.data".to_string()
}

fn default_voxel_size() -> f64 {
    1.0
}

fn default_background_material() -> BackgroundMaterial {
    BackgroundMaterial::Liquid
}

fn default_background_element() -> String {
    "Co".to_string()
}

fn default_density() -> f64 {
    0.08
}

fn default_separation() -> f64 {
    2.0  // Separation between nanoparticle surfaces (atomic buffer added automatically)
}

fn default_placement_mode() -> PlacementMode {
    PlacementMode::Collision
}

fn default_min_atom_distance() -> f64 {
    2.0  // Minimum distance between atom centers in overlap mode
}

fn default_collision_buffer() -> f64 {
    3.0  // Buffer distance added to collision detection shapes
}

// Dynamics config defaults
fn default_max_iterations() -> usize {
    10000
}

fn default_force_tolerance() -> f64 {
    0.01  // kcal/mol/Å
}

fn default_energy_tolerance() -> f64 {
    1e-6  // kcal/mol
}

fn default_max_displacement() -> f64 {
    0.5  // Å
}

fn default_max_rotation() -> f64 {
    0.1  // radians (~5.7 degrees)
}

fn default_use_fire() -> bool {
    true  // FIRE is generally faster than steepest descent
}

fn default_potential_type() -> PotentialType {
    PotentialType::SoftSphere {
        epsilon: 1.0,
        sigma: 3.0,
    }
}

fn default_sphere_radius() -> f64 {
    2.0  // Å
}

fn default_target_spacing() -> f64 {
    3.0  // Å
}