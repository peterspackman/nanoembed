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