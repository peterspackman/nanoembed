use crate::types::{Atom, AtomTypeMap, Nanoparticle, ConvexHull};
use nalgebra::{Point3, Vector3, Rotation3, Unit};
use rand::{rngs::StdRng, Rng};
use std::fs::File;
use std::io::{BufRead, BufReader};

impl Nanoparticle {
    pub fn collect_elements_from_xyz(filename: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Parse header
        let num_atoms: usize = lines.next().unwrap()?.parse()?;
        lines.next(); // Skip comment line

        let mut elements = Vec::with_capacity(num_atoms);

        for line in lines {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                elements.push(parts[0].to_string());
            }
        }

        Ok(elements)
    }

    pub fn from_xyz_file(filename: &str, type_map: &mut AtomTypeMap) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Parse header
        let num_atoms: usize = lines.next().unwrap()?.parse()?;
        lines.next(); // Skip comment line

        let mut atoms = Vec::with_capacity(num_atoms);
        let mut positions = Vec::with_capacity(num_atoms);

        for line in lines {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                let element = parts[0].to_string();
                let x: f64 = parts[1].parse()?;
                let y: f64 = parts[2].parse()?;
                let z: f64 = parts[3].parse()?;

                let position = Point3::new(x, y, z);
                positions.push(position);

                let atom_type = type_map.get_or_create_type(&element);
                atoms.push(Atom {
                    element,
                    position,
                    atom_type,
                });
            }
        }

        // Calculate center
        let center = positions.iter().fold(Point3::origin(), |acc, p| acc + p.coords) / positions.len() as f64;

        // Center the nanoparticle at origin
        for (atom, pos) in atoms.iter_mut().zip(positions.iter_mut()) {
            atom.position -= center.coords;
            *pos -= center.coords;
        }

        // Compute convex hull from atom positions - this defines the precise boundary
        let atom_positions: Vec<Point3<f64>> = atoms.iter().map(|a| a.position).collect();
        let hull = ConvexHull::from_points(&atom_positions);

        Ok(Nanoparticle {
            atoms,
            center: Point3::origin(), // Now centered at origin
            hull,
        })
    }

    pub fn translate(&mut self, offset: &Vector3<f64>) {
        self.center += offset;
        for atom in &mut self.atoms {
            atom.position += offset;
        }
    }

    pub fn rotate(&mut self, rotation: &Rotation3<f64>) {
        // Rotate atoms around current center
        for atom in &mut self.atoms {
            let relative_pos = atom.position - self.center;
            let rotated_pos = rotation * relative_pos;
            atom.position = self.center + rotated_pos;
        }
    }

    pub fn random_orientation(rng: &mut StdRng) -> Rotation3<f64> {
        // Generate random rotation using random quaternion
        // Sample uniformly from unit quaternion space
        let u1 = rng.gen::<f64>();
        let u2 = rng.gen::<f64>();
        let u3 = rng.gen::<f64>();
        Self::quasi_random_orientation(u1, u2, u3)
    }

    pub fn quasi_random_orientation(u1: f64, u2: f64, u3: f64) -> Rotation3<f64> {
        // Generate rotation using quaternion from uniform [0,1] values
        // Based on Shepperd's method for uniform quaternion sampling
        let sqrt_1_u1 = (1.0 - u1).sqrt();
        let sqrt_u1 = u1.sqrt();

        let q_w = sqrt_1_u1 * (2.0 * std::f64::consts::PI * u2).sin();
        let q_x = sqrt_1_u1 * (2.0 * std::f64::consts::PI * u2).cos();
        let q_y = sqrt_u1 * (2.0 * std::f64::consts::PI * u3).sin();
        let q_z = sqrt_u1 * (2.0 * std::f64::consts::PI * u3).cos();

        // Create rotation from quaternion
        let axis = Vector3::new(q_x, q_y, q_z);
        if let Some(unit_axis) = Unit::try_new(axis, 1e-10) {
            Rotation3::from_axis_angle(&unit_axis, 2.0 * q_w.acos())
        } else {
            Rotation3::identity()
        }
    }
}