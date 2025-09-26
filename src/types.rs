use nalgebra::Point3;
use std::collections::HashMap;
use std::str::FromStr;
use rapier3d::geometry::SharedShape;
use rapier3d::math::{Point, Isometry};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Element {
    H, C, N, O, F, Si, P, S, Cl, Ar, K, Ca, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn,
    Ga, As, Se, Br, Zr, Mo, Tc, Ru, Pd, Ag, Cd, In, Sn, Sb, Te, I, W, Re,
    Os, Ir, Pt, Au, Hg, Tl, Pb, Bi,
}

impl Element {
    pub fn symbol(&self) -> &'static str {
        match self {
            Element::H => "H", Element::C => "C", Element::N => "N", Element::O => "O",
            Element::F => "F", Element::Si => "Si", Element::P => "P", Element::S => "S",
            Element::Cl => "Cl", Element::Ar => "Ar", Element::K => "K", Element::Ca => "Ca",
            Element::Ti => "Ti", Element::Cr => "Cr", Element::Mn => "Mn", Element::Fe => "Fe",
            Element::Co => "Co", Element::Ni => "Ni", Element::Cu => "Cu", Element::Zn => "Zn",
            Element::Ga => "Ga", Element::As => "As", Element::Se => "Se", Element::Br => "Br",
            Element::Zr => "Zr", Element::Mo => "Mo", Element::Tc => "Tc", Element::Ru => "Ru",
            Element::Pd => "Pd", Element::Ag => "Ag", Element::Cd => "Cd", Element::In => "In",
            Element::Sn => "Sn", Element::Sb => "Sb", Element::Te => "Te", Element::I => "I",
            Element::W => "W", Element::Re => "Re", Element::Os => "Os", Element::Ir => "Ir",
            Element::Pt => "Pt", Element::Au => "Au", Element::Hg => "Hg", Element::Tl => "Tl",
            Element::Pb => "Pb", Element::Bi => "Bi",
        }
    }

    pub fn atomic_mass(&self) -> f64 {
        match self {
            Element::H => 1.008, Element::C => 12.011, Element::N => 14.007, Element::O => 15.999,
            Element::F => 18.998, Element::Si => 28.085, Element::P => 30.974, Element::S => 32.065,
            Element::Cl => 35.453, Element::Ar => 39.948, Element::K => 39.098, Element::Ca => 40.078,
            Element::Ti => 47.867, Element::Cr => 51.996, Element::Mn => 54.938, Element::Fe => 55.845,
            Element::Co => 58.933, Element::Ni => 58.693, Element::Cu => 63.546, Element::Zn => 65.409,
            Element::Ga => 69.723, Element::As => 74.922, Element::Se => 78.971, Element::Br => 79.904,
            Element::Zr => 91.224, Element::Mo => 95.96, Element::Tc => 98.0, Element::Ru => 101.07,
            Element::Pd => 106.42, Element::Ag => 107.868, Element::Cd => 112.411, Element::In => 114.818,
            Element::Sn => 118.710, Element::Sb => 121.760, Element::Te => 127.60, Element::I => 126.904,
            Element::W => 183.84, Element::Re => 186.207, Element::Os => 190.23, Element::Ir => 192.217,
            Element::Pt => 195.084, Element::Au => 196.966, Element::Hg => 200.592, Element::Tl => 204.383,
            Element::Pb => 207.2, Element::Bi => 208.980,
        }
    }

    pub fn covalent_radius(&self) -> f64 {
        match self {
            Element::H => 0.31, Element::C => 0.76, Element::N => 0.71, Element::O => 0.66,
            Element::F => 0.57, Element::Si => 1.11, Element::P => 1.07, Element::S => 1.05,
            Element::Cl => 0.99, Element::Ar => 1.06, Element::K => 2.03, Element::Ca => 1.76,
            Element::Ti => 1.60, Element::Cr => 1.39, Element::Mn => 1.39, Element::Fe => 1.32,
            Element::Co => 1.26, Element::Ni => 1.24, Element::Cu => 1.32, Element::Zn => 1.22,
            Element::Ga => 1.22, Element::As => 1.19, Element::Se => 1.20, Element::Br => 1.20,
            Element::Zr => 1.75, Element::Mo => 1.54, Element::Tc => 1.47, Element::Ru => 1.46,
            Element::Pd => 1.39, Element::Ag => 1.45, Element::Cd => 1.44, Element::In => 1.42,
            Element::Sn => 1.39, Element::Sb => 1.39, Element::Te => 1.38, Element::I => 1.39,
            Element::W => 1.62, Element::Re => 1.51, Element::Os => 1.44, Element::Ir => 1.41,
            Element::Pt => 1.36, Element::Au => 1.36, Element::Hg => 1.32, Element::Tl => 1.45,
            Element::Pb => 1.46, Element::Bi => 1.48,
        }
    }
}

impl FromStr for Element {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Element::H), "C" => Ok(Element::C), "N" => Ok(Element::N), "O" => Ok(Element::O),
            "F" => Ok(Element::F), "Si" => Ok(Element::Si), "P" => Ok(Element::P), "S" => Ok(Element::S),
            "Cl" => Ok(Element::Cl), "Ar" => Ok(Element::Ar), "K" => Ok(Element::K), "Ca" => Ok(Element::Ca),
            "Ti" => Ok(Element::Ti), "Cr" => Ok(Element::Cr), "Mn" => Ok(Element::Mn), "Fe" => Ok(Element::Fe),
            "Co" => Ok(Element::Co), "Ni" => Ok(Element::Ni), "Cu" => Ok(Element::Cu), "Zn" => Ok(Element::Zn),
            "Ga" => Ok(Element::Ga), "As" => Ok(Element::As), "Se" => Ok(Element::Se), "Br" => Ok(Element::Br),
            "Zr" => Ok(Element::Zr), "Mo" => Ok(Element::Mo), "Tc" => Ok(Element::Tc), "Ru" => Ok(Element::Ru),
            "Pd" => Ok(Element::Pd), "Ag" => Ok(Element::Ag), "Cd" => Ok(Element::Cd), "In" => Ok(Element::In),
            "Sn" => Ok(Element::Sn), "Sb" => Ok(Element::Sb), "Te" => Ok(Element::Te), "I" => Ok(Element::I),
            "W" => Ok(Element::W), "Re" => Ok(Element::Re), "Os" => Ok(Element::Os), "Ir" => Ok(Element::Ir),
            "Pt" => Ok(Element::Pt), "Au" => Ok(Element::Au), "Hg" => Ok(Element::Hg), "Tl" => Ok(Element::Tl),
            "Pb" => Ok(Element::Pb), "Bi" => Ok(Element::Bi),
            _ => Err(format!("Unknown element: {}", s)),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub element: String,
    pub position: Point3<f64>,
    pub atom_type: u32,
}

#[derive(Debug, Clone)]
pub struct ConvexHull {
    pub shape: SharedShape,
}

impl ConvexHull {
    pub fn from_points(points: &[Point3<f64>]) -> Self {
        // Convert nalgebra points to rapier format
        let rapier_points: Vec<Point<f32>> = points
            .iter()
            .map(|p| Point::new(p.x as f32, p.y as f32, p.z as f32))
            .collect();

        // Create convex hull using Rapier - this is the actual precise boundary
        let shape = SharedShape::convex_hull(&rapier_points).unwrap_or_else(|| {
            // Fallback to bounding sphere only if convex hull computation fails
            println!("Warning: Convex hull computation failed, falling back to sphere");
            let center = rapier_points.iter().fold(Point::origin(), |acc, p| Point::new(
                acc.x + p.x / rapier_points.len() as f32,
                acc.y + p.y / rapier_points.len() as f32,
                acc.z + p.z / rapier_points.len() as f32,
            ));
            let radius = rapier_points.iter().map(|p| (*p - center).norm()).fold(0.0, f32::max);
            SharedShape::ball(radius)
        });

        ConvexHull { shape }
    }

    pub fn bounding_radius(&self) -> f64 {
        // Get the AABB (axis-aligned bounding box) and use it to estimate radius
        let aabb = self.shape.compute_aabb(&Isometry::identity());
        let half_extents = aabb.half_extents();
        (half_extents.x.max(half_extents.y).max(half_extents.z)) as f64
    }
}

#[derive(Debug, Clone)]
pub struct Nanoparticle {
    pub atoms: Vec<Atom>,
    pub center: Point3<f64>,
    pub hull: ConvexHull, // Precise convex hull boundary
}

impl Nanoparticle {
    pub fn radius(&self) -> f64 {
        self.hull.bounding_radius()
    }
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

    pub fn pre_populate_types(&mut self, elements: &[String]) {
        // Convert strings to Element enum and sort by atomic mass
        let mut unique_elements: Vec<Element> = elements
            .iter()
            .filter_map(|s| Element::from_str(s).ok())
            .collect();

        unique_elements.sort_by(|a, b| a.atomic_mass().partial_cmp(&b.atomic_mass()).unwrap());
        unique_elements.dedup();

        // Pre-assign types in mass order (lightest to heaviest)
        for element in unique_elements {
            self.get_or_create_type(element.symbol());
        }
    }

    pub fn get_or_create_type(&mut self, element: &str) -> u32 {
        if let Some(&atom_type) = self.element_to_type.get(element) {
            atom_type
        } else {
            let atom_type = self.next_type;
            self.next_type += 1;

            // Get atomic properties from Element enum, fallback to carbon for unknown elements
            let (mass, cov_radius) = match Element::from_str(element) {
                Ok(elem) => (elem.atomic_mass(), elem.covalent_radius()),
                Err(_) => (Element::C.atomic_mass(), Element::C.covalent_radius()),
            };

            self.element_to_type.insert(element.to_string(), atom_type);
            self.type_to_element.insert(atom_type, element.to_string());
            self.type_to_mass.insert(atom_type, mass);
            self.type_to_covalent_radius.insert(atom_type, cov_radius);

            atom_type
        }
    }

    pub fn get_adaptive_voxel_size(&self, target_density: f64) -> f64 {
        // Find the largest covalent radius among all atom types
        let max_cov_radius = self.type_to_covalent_radius.values()
            .fold(0.76_f64, |max, &radius| max.max(radius)); // Default to carbon's radius

        // Physical minimum: 1.5x largest covalent radius to prevent unphysical overlaps
        let min_separation = max_cov_radius * 1.5;

        // Density-based optimal size: cube root of inverse density gives natural spacing
        let volume_per_atom = 1.0 / target_density; // Å³ per atom
        let density_spacing = volume_per_atom.powf(1.0 / 3.0); // Å

        // Use the larger of physical minimum or density-based spacing
        let optimal_size = density_spacing.max(min_separation);

        // Clamp to reasonable range
        let voxel_size = optimal_size.max(2.0).min(20.0);

        println!("Density-adaptive voxel sizing:");
        println!("  Target density: {:.3} atoms/Å³", target_density);
        println!("  Volume per atom: {:.1} Å³", volume_per_atom);
        println!("  Density spacing: {:.2} Å", density_spacing);
        println!("  Physical minimum: {:.2} Å (1.5 × {:.2} Å cov radius)",
                 min_separation, max_cov_radius);
        println!("  Selected voxel size: {:.2} Å", voxel_size);

        voxel_size
    }
}