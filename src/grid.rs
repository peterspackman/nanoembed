use nalgebra::{Point3, Vector3};

#[derive(Debug, Clone)]
pub struct VoxelGrid {
    pub voxel_size: f64,
    pub box_size: Vector3<f64>,
    pub grid_dims: (usize, usize, usize),
    pub occupied: Vec<bool>,
}

impl VoxelGrid {
    pub fn new(voxel_size: f64, box_size: Vector3<f64>) -> Self {
        let nx = (box_size.x / voxel_size).ceil() as usize;
        let ny = (box_size.y / voxel_size).ceil() as usize;
        let nz = (box_size.z / voxel_size).ceil() as usize;

        Self {
            voxel_size,
            box_size,
            grid_dims: (nx, ny, nz),
            occupied: vec![false; nx * ny * nz],
        }
    }

    pub fn get_index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.grid_dims.0 * self.grid_dims.1 + j * self.grid_dims.0 + i
    }

    pub fn get_voxel_coords(&self, pos: &Point3<f64>) -> Option<(usize, usize, usize)> {
        if pos.x < 0.0 || pos.x >= self.box_size.x ||
           pos.y < 0.0 || pos.y >= self.box_size.y ||
           pos.z < 0.0 || pos.z >= self.box_size.z {
            return None;
        }

        let i = (pos.x / self.voxel_size).floor() as usize;
        let j = (pos.y / self.voxel_size).floor() as usize;
        let k = (pos.z / self.voxel_size).floor() as usize;

        Some((
            i.min(self.grid_dims.0 - 1),
            j.min(self.grid_dims.1 - 1),
            k.min(self.grid_dims.2 - 1),
        ))
    }

    pub fn get_voxel_center(&self, i: usize, j: usize, k: usize) -> Point3<f64> {
        Point3::new(
            (i as f64 + 0.5) * self.voxel_size,
            (j as f64 + 0.5) * self.voxel_size,
            (k as f64 + 0.5) * self.voxel_size,
        )
    }

    pub fn is_occupied(&self, i: usize, j: usize, k: usize) -> bool {
        if i >= self.grid_dims.0 || j >= self.grid_dims.1 || k >= self.grid_dims.2 {
            return true; // Out of bounds = occupied
        }
        let idx = self.get_index(i, j, k);
        self.occupied[idx]
    }
}

#[derive(Debug, Clone)]
pub struct NanoparticleGrid {
    pub voxel_size: f64,
    pub box_size: Vector3<f64>,
    pub grid_dims: (usize, usize, usize),
    pub occupied: Vec<bool>,
    pub nanoparticle_centers: Vec<Point3<f64>>,
    pub nanoparticle_radii: Vec<f64>,
}

impl NanoparticleGrid {
    pub fn new(_average_np_radius: f64, box_size: Vector3<f64>, type_map: &crate::types::AtomTypeMap, target_density: f64) -> Self {
        // Use density-adaptive voxel sizing
        let voxel_size = type_map.get_adaptive_voxel_size(target_density);

        let nx = (box_size.x / voxel_size).ceil() as usize;
        let ny = (box_size.y / voxel_size).ceil() as usize;
        let nz = (box_size.z / voxel_size).ceil() as usize;

        println!("Creating nanoparticle grid: {}x{}x{} voxels ({:.1} Å each)",
                 nx, ny, nz, voxel_size);

        Self {
            voxel_size,
            box_size,
            grid_dims: (nx, ny, nz),
            occupied: vec![false; nx * ny * nz],
            nanoparticle_centers: Vec::new(),
            nanoparticle_radii: Vec::new(),
        }
    }

    pub fn get_index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.grid_dims.0 * self.grid_dims.1 + j * self.grid_dims.0 + i
    }

    pub fn get_voxel_coords(&self, pos: &Point3<f64>) -> Option<(usize, usize, usize)> {
        if pos.x < 0.0 || pos.x >= self.box_size.x ||
           pos.y < 0.0 || pos.y >= self.box_size.y ||
           pos.z < 0.0 || pos.z >= self.box_size.z {
            return None;
        }

        let i = (pos.x / self.voxel_size).floor() as usize;
        let j = (pos.y / self.voxel_size).floor() as usize;
        let k = (pos.z / self.voxel_size).floor() as usize;

        Some((
            i.min(self.grid_dims.0 - 1),
            j.min(self.grid_dims.1 - 1),
            k.min(self.grid_dims.2 - 1),
        ))
    }

    pub fn check_collision(&self, center: &Point3<f64>, radius: f64, separation: f64) -> bool {
        let search_radius = radius + separation;

        // Quick sphere-sphere check against all placed nanoparticles
        for (i, &np_center) in self.nanoparticle_centers.iter().enumerate() {
            let distance = (center - np_center).norm();
            let min_distance = search_radius + self.nanoparticle_radii[i];

            if distance < min_distance {
                return true; // Collision detected
            }
        }

        false
    }

    pub fn add_nanoparticle(&mut self, center: Point3<f64>, radius: f64) {
        // Store nanoparticle info for distance checks
        self.nanoparticle_centers.push(center);
        self.nanoparticle_radii.push(radius);
    }

    pub fn mark_atoms_occupied(&mut self, atoms: &[crate::types::Atom]) {
        // Mark voxels that contain actual nanoparticle atoms
        let mut marked_voxels = std::collections::HashSet::new();
        for atom in atoms {
            if let Some((i, j, k)) = self.get_voxel_coords(&atom.position) {
                let idx = self.get_index(i, j, k);
                if !self.occupied[idx] {
                    self.occupied[idx] = true;
                    marked_voxels.insert((i, j, k));
                }
            }
        }
        println!("  Marked {} voxels as occupied for {} atoms", marked_voxels.len(), atoms.len());
    }

    pub fn is_voxel_occupied(&self, i: usize, j: usize, k: usize) -> bool {
        if i >= self.grid_dims.0 || j >= self.grid_dims.1 || k >= self.grid_dims.2 {
            return true; // Out of bounds = occupied
        }
        let idx = self.get_index(i, j, k);
        self.occupied[idx]
    }

    pub fn print_occupation_summary(&self) {
        let total_voxels = self.occupied.len();
        let occupied_voxels = self.occupied.iter().filter(|&&x| x).count();
        let free_voxels = total_voxels - occupied_voxels;
        let free_percentage = (free_voxels as f64 / total_voxels as f64) * 100.0;

        println!("Voxel occupation: {}/{} occupied, {} free ({:.1}% free)",
                 occupied_voxels, total_voxels, free_voxels, free_percentage);
    }
}