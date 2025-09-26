use nalgebra::{Point3, Vector3};
use rapier3d::geometry::{ColliderBuilder, ColliderSet, SharedShape};
use rapier3d::math::Isometry;
use rapier3d::pipeline::{QueryPipeline, QueryFilter};
use rapier3d::dynamics::{RigidBodyBuilder, RigidBodySet};

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

pub struct NanoparticleGrid {
    pub voxel_size: f64,
    pub box_size: Vector3<f64>,
    pub grid_dims: (usize, usize, usize),
    pub occupied: Vec<bool>,
    pub rigid_body_set: RigidBodySet,
    pub collider_set: ColliderSet,
    pub query_pipeline: QueryPipeline,
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
            rigid_body_set: RigidBodySet::new(),
            collider_set: ColliderSet::new(),
            query_pipeline: QueryPipeline::new(),
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

    pub fn check_collision(&self, nanoparticle: &crate::types::Nanoparticle, center: &Point3<f64>, separation: f64, type_map: &crate::types::AtomTypeMap) -> bool {
        // Create the candidate shape at the proposed position
        let candidate_pos = Isometry::translation(center.x as f32, center.y as f32, center.z as f32);

        // Use 3 Å expansion - this was tested and works to prevent close contacts
        let max_cov_radius = 3.0;

        // Since both shapes are expanded by max_cov_radius, we need additional separation
        // The distance query returns distance between expanded surfaces, so we need the requested separation
        let total_buffer = separation;

        // Create an expanded convex hull using Minkowski sum with a sphere
        let test_shape = if let Some(polyhedron) = nanoparticle.hull.shape.as_convex_polyhedron() {
            rapier3d::geometry::SharedShape::round_convex_hull(
                polyhedron.points(),
                max_cov_radius as f32
            ).unwrap_or_else(|| {
                let radius = nanoparticle.radius() as f32 + max_cov_radius as f32;
                rapier3d::geometry::SharedShape::ball(radius)
            })
        } else {
            let radius = nanoparticle.radius() as f32 + max_cov_radius as f32;
            rapier3d::geometry::SharedShape::ball(radius)
        };

        // Check if any existing nanoparticle is too close
        for (_, collider) in self.collider_set.iter() {
            // Calculate distance between the candidate shape and this existing collider
            let distance_result = rapier3d::parry::query::distance(
                &candidate_pos,
                test_shape.as_ref(),
                collider.position(),
                collider.shape()
            );

            if let Ok(distance) = distance_result {
                if distance < total_buffer as f32 {
                    return true; // Too close!
                }
            }
        }

        false
    }

    pub fn add_nanoparticle(&mut self, nanoparticle: &crate::types::Nanoparticle, center: Point3<f64>, type_map: &crate::types::AtomTypeMap) {
        // Create a rigid body for the nanoparticle (static/fixed)
        let rigid_body = RigidBodyBuilder::fixed()
            .translation(nalgebra::Vector3::new(center.x as f32, center.y as f32, center.z as f32))
            .build();
        let rb_handle = self.rigid_body_set.insert(rigid_body);

        // Create an expanded collider that accounts for atomic radii - use same expansion as in check_collision
        let max_cov_radius = 3.0;
        let expanded_shape = rapier3d::geometry::SharedShape::round_convex_hull(
            nanoparticle.hull.shape.as_convex_polyhedron().unwrap().points(),
            max_cov_radius as f32
        ).unwrap_or_else(|| {
            // Fallback to expanded sphere if round_convex_hull fails
            let radius = nanoparticle.radius() as f32 + max_cov_radius as f32;
            rapier3d::geometry::SharedShape::ball(radius)
        });
        let collider = ColliderBuilder::new(expanded_shape).build();
        self.collider_set.insert_with_parent(collider, rb_handle, &mut self.rigid_body_set);

        // Update the query pipeline to include the new collider
        self.query_pipeline.update(&self.rigid_body_set, &self.collider_set);
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

    pub fn mark_atoms_occupied_with_overlap(&mut self, atoms: &[crate::types::Atom], placement_mode: &crate::config::PlacementMode, placed_atoms: &mut Vec<crate::types::Atom>) -> Vec<crate::types::Atom> {
        match placement_mode {
            crate::config::PlacementMode::Collision => {
                // Same as before - mark all atoms and add to placed_atoms list
                self.mark_atoms_occupied(atoms);
                placed_atoms.extend_from_slice(atoms);
                atoms.to_vec()
            },
            crate::config::PlacementMode::Overlap => {
                // Check each new atom against all previously placed atoms
                let mut kept_atoms = Vec::new();
                let min_distance = 2.0; // Minimum distance between atom centers (2x covalent radius)

                for new_atom in atoms {
                    let mut overlaps = false;

                    // Check against all previously placed atoms
                    for existing_atom in placed_atoms.iter() {
                        let distance = (new_atom.position - existing_atom.position).norm();
                        if distance < min_distance {
                            overlaps = true;
                            break; // First placed wins - skip this atom
                        }
                    }

                    if !overlaps {
                        kept_atoms.push(new_atom.clone());
                        placed_atoms.push(new_atom.clone()); // Add to global placed atoms list
                    }
                }

                // Still mark voxels for background generation
                self.mark_atoms_occupied(&kept_atoms);

                println!("  Kept {} atoms out of {} (overlap mode: {:.1}% rejected)",
                         kept_atoms.len(), atoms.len(),
                         100.0 * (atoms.len() - kept_atoms.len()) as f64 / atoms.len() as f64);
                kept_atoms
            }
        }
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