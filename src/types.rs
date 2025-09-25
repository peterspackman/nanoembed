use nalgebra::Point3;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Atom {
    pub element: String,
    pub position: Point3<f64>,
    pub atom_type: u32,
}

#[derive(Debug, Clone)]
pub struct Nanoparticle {
    pub atoms: Vec<Atom>,
    pub center: Point3<f64>,
    pub radius: f64, // Bounding sphere radius
}

#[derive(Debug)]
pub struct AtomTypeMap {
    pub element_to_type: HashMap<String, u32>,
    pub type_to_element: HashMap<u32, String>,
    pub type_to_mass: HashMap<u32, f64>,
    pub type_to_covalent_radius: HashMap<u32, f64>,
    pub next_type: u32,
}

impl AtomTypeMap {
    pub fn new() -> Self {
        Self {
            element_to_type: HashMap::new(),
            type_to_element: HashMap::new(),
            type_to_mass: HashMap::new(),
            type_to_covalent_radius: HashMap::new(),
            next_type: 1,
        }
    }

    pub fn get_or_create_type(&mut self, element: &str) -> u32 {
        if let Some(&atom_type) = self.element_to_type.get(element) {
            atom_type
        } else {
            let atom_type = self.next_type;
            self.next_type += 1;

            // Approximate atomic masses (amu) and covalent radii (Å)
            let (mass, cov_radius) = match element {
                "H" => (1.008, 0.31),
                "C" => (12.011, 0.76),
                "N" => (14.007, 0.71),
                "O" => (15.999, 0.66),
                "F" => (18.998, 0.57),
                "Si" => (28.085, 1.11),
                "P" => (30.974, 1.07),
                "S" => (32.065, 1.05),
                "Cl" => (35.453, 0.99),
                "Ar" => (39.948, 1.06),
                "K" => (39.098, 2.03),
                "Ca" => (40.078, 1.76),
                "Ti" => (47.867, 1.60),
                "Cr" => (51.996, 1.39),
                "Mn" => (54.938, 1.39),
                "Fe" => (55.845, 1.32),
                "Co" => (58.933, 1.26),
                "Ni" => (58.693, 1.24),
                "Cu" => (63.546, 1.32),
                "Zn" => (65.409, 1.22),
                "Ga" => (69.723, 1.22),
                "As" => (74.922, 1.19),
                "Se" => (78.971, 1.20),
                "Br" => (79.904, 1.20),
                "Zr" => (91.224, 1.75),
                "Mo" => (95.96, 1.54),
                "Tc" => (98.0, 1.47),
                "Ru" => (101.07, 1.46),
                "Pd" => (106.42, 1.39),
                "Ag" => (107.868, 1.45),
                "Cd" => (112.411, 1.44),
                "In" => (114.818, 1.42),
                "Sn" => (118.710, 1.39),
                "Sb" => (121.760, 1.39),
                "Te" => (127.60, 1.38),
                "I" => (126.904, 1.39),
                "W" => (183.84, 1.62),
                "Re" => (186.207, 1.51),
                "Os" => (190.23, 1.44),
                "Ir" => (192.217, 1.41),
                "Pt" => (195.084, 1.36),
                "Au" => (196.966, 1.36),
                "Hg" => (200.592, 1.32),
                "Tl" => (204.383, 1.45),
                "Pb" => (207.2, 1.46),
                "Bi" => (208.980, 1.48),
                _ => (12.011, 0.76), // Default to carbon
            };

            self.element_to_type.insert(element.to_string(), atom_type);
            self.type_to_element.insert(atom_type, element.to_string());
            self.type_to_mass.insert(atom_type, mass);
            self.type_to_covalent_radius.insert(atom_type, cov_radius);

            atom_type
        }
    }

    pub fn get_minimum_voxel_size(&self) -> f64 {
        // Find the largest covalent radius among all atom types
        let max_cov_radius = self.type_to_covalent_radius.values()
            .fold(0.76_f64, |max, &radius| max.max(radius)); // Default to carbon's radius

        // Use 1.5x the largest covalent radius to prevent unphysical overlaps
        let min_separation = max_cov_radius * 1.5;

        println!("Maximum covalent radius: {:.2} Å, minimum voxel size: {:.2} Å",
                 max_cov_radius, min_separation);

        min_separation
    }
}