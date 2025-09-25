use crate::grid::NanoparticleGrid;
use crate::types::Atom;
use nalgebra::{Point3, Vector3};
use rand::{rngs::StdRng, Rng};
use rayon::prelude::*;

pub fn generate_liquid_background(
    box_size: Vector3<f64>,
    density: f64,
    grid: &NanoparticleGrid,
    _rng: &mut StdRng,
) -> Vec<Atom> {
    println!("Generating liquid background with density {} atoms/Å³...", density);

    let mut atoms = Vec::new();
    let voxel_volume = grid.voxel_size.powi(3);
    let atoms_per_voxel = density * voxel_volume;

    println!("Using {}x{}x{} background voxels ({:.1} Å each)",
             grid.grid_dims.0, grid.grid_dims.1, grid.grid_dims.2, grid.voxel_size);

    // Process voxels in parallel for better performance
    let voxel_results: Vec<Vec<Atom>> = (0..grid.grid_dims.2).into_par_iter().map(|k| {
        let mut local_rng = rand::thread_rng();
        let mut local_atoms = Vec::new();

        for j in 0..grid.grid_dims.1 {
            for i in 0..grid.grid_dims.0 {
                // Skip voxels that contain nanoparticle atoms
                if !grid.is_voxel_occupied(i, j, k) {
                    let voxel_center = Point3::new(
                        (i as f64 + 0.5) * grid.voxel_size,
                        (j as f64 + 0.5) * grid.voxel_size,
                        (k as f64 + 0.5) * grid.voxel_size,
                    );

                    // Decide whether to place an atom in this voxel
                    let place_atom = if atoms_per_voxel < 1.0 {
                        local_rng.gen::<f64>() < atoms_per_voxel
                    } else {
                        true // If density is high enough, always try to place one atom per voxel
                    };

                    if place_atom {
                        // Place atom at voxel center with small random offset to avoid k-d tree collisions
                        use rand::Rng;
                        let offset = 0.1; // Small offset to break ties while maintaining ~1 Å spacing
                        let test_pos = Point3::new(
                            voxel_center.x + local_rng.gen_range(-offset..offset),
                            voxel_center.y + local_rng.gen_range(-offset..offset),
                            voxel_center.z + local_rng.gen_range(-offset..offset),
                        );

                        // Ensure within box bounds
                        if test_pos.x >= 0.0 && test_pos.x < box_size.x &&
                           test_pos.y >= 0.0 && test_pos.y < box_size.y &&
                           test_pos.z >= 0.0 && test_pos.z < box_size.z {
                            local_atoms.push(Atom {
                                element: "Co".to_string(),
                                position: test_pos,
                                atom_type: 999, // Will be updated later with proper type
                            });
                        }
                    }
                }
            }
        }
        local_atoms
    }).collect();

    // Collect all atoms
    for mut voxel_atoms in voxel_results {
        atoms.append(&mut voxel_atoms);
    }

    println!("Generated {} liquid atoms", atoms.len());
    atoms
}