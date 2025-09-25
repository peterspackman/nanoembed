use crate::grid::NanoparticleGrid;
use crate::types::Atom;
use nalgebra::{Point3, Vector3};
use rand::{rngs::StdRng, Rng};
use rayon::prelude::*;

pub fn generate_liquid_background(
    box_size: Vector3<f64>,
    density: f64,
    grid: &NanoparticleGrid,
    rng: &mut StdRng,
) -> Vec<Atom> {
    println!("Generating liquid background with density {} atoms/Å³...", density);

    let mut all_atoms = Vec::new();
    let voxel_volume = grid.voxel_size.powi(3);

    // Calculate free volume fraction
    let total_voxels = grid.occupied.len();
    let occupied_voxels = grid.occupied.iter().filter(|&&x| x).count();
    let free_voxels = total_voxels - occupied_voxels;
    let free_volume_fraction = free_voxels as f64 / total_voxels as f64;

    // Adjust density to account for free volume
    let adjusted_density = density / free_volume_fraction;
    let atoms_per_voxel = adjusted_density * voxel_volume;

    println!("Using {}x{}x{} background voxels ({:.1} Å each)",
             grid.grid_dims.0, grid.grid_dims.1, grid.grid_dims.2, grid.voxel_size);
    println!("Free volume: {:.1}% ({}/{} voxels)",
             free_volume_fraction * 100.0, free_voxels, total_voxels);
    println!("Adjusted density: {:.3} → {:.3} atoms/Å³ for free volume",
             density, adjusted_density);
    println!("Processing in waves along Z dimension for better locality...");

    // Process in waves along Z dimension
    let mut previous_wave_atoms: Vec<Vec<Vec<Option<Point3<f64>>>>> = vec![
        vec![vec![None; grid.grid_dims.0]; grid.grid_dims.1];
        2 // Keep last 2 Z layers for collision checking
    ];

    for k in 0..grid.grid_dims.2 {
        let wave_start = std::time::Instant::now();

        // Process current Z slice in parallel (by Y strips)
        let y_strip_results: Vec<Vec<Atom>> = (0..grid.grid_dims.1).into_par_iter().map(|j| {
            let mut local_rng = rand::thread_rng();
            let mut local_atoms = Vec::new();

            for i in 0..grid.grid_dims.0 {
                // Skip voxels that contain nanoparticle atoms
                if !grid.is_voxel_occupied(i, j, k) {
                    // Decide whether to place an atom in this voxel
                    let place_atom = if atoms_per_voxel < 1.0 {
                        local_rng.gen::<f64>() < atoms_per_voxel
                    } else {
                        true
                    };

                    if place_atom {
                        let voxel_center = Point3::new(
                            (i as f64 + 0.5) * grid.voxel_size,
                            (j as f64 + 0.5) * grid.voxel_size,
                            (k as f64 + 0.5) * grid.voxel_size,
                        );

                        // Small random offset to avoid k-d tree collisions
                        let offset = 0.1;
                        let candidate_pos = Point3::new(
                            voxel_center.x + local_rng.gen_range(-offset..offset),
                            voxel_center.y + local_rng.gen_range(-offset..offset),
                            voxel_center.z + local_rng.gen_range(-offset..offset),
                        );

                        // Check collision with previous waves and nanoparticles
                        let mut collision = false;

                        // Check against previous waves
                        if k > 0 {
                            collision = check_collision_with_previous_waves(
                                i, j, &previous_wave_atoms, grid.voxel_size, &candidate_pos
                            );
                        }

                        // Check against nanoparticles using voxel occupancy
                        if !collision {
                            collision = check_collision_with_nanoparticle_voxels(
                                &candidate_pos, grid
                            );
                        }

                        if !collision &&
                           candidate_pos.x >= 0.0 && candidate_pos.x < box_size.x &&
                           candidate_pos.y >= 0.0 && candidate_pos.y < box_size.y &&
                           candidate_pos.z >= 0.0 && candidate_pos.z < box_size.z {

                            local_atoms.push(Atom {
                                element: "Co".to_string(),
                                position: candidate_pos,
                                atom_type: 999,
                            });
                        }
                    }
                }
            }
            local_atoms
        }).collect();

        // Update previous wave tracking
        let current_wave_idx = k % 2;
        let mut current_wave: Vec<Vec<Option<Point3<f64>>>> = vec![vec![None; grid.grid_dims.0]; grid.grid_dims.1];

        // Collect atoms from this wave and update tracking
        for (j, y_strip_atoms) in y_strip_results.into_iter().enumerate() {
            for atom in y_strip_atoms.iter() {
                // Convert position back to voxel coordinates
                let i = ((atom.position.x / grid.voxel_size).floor() as usize).min(grid.grid_dims.0 - 1);
                current_wave[j][i] = Some(atom.position);
            }
            all_atoms.extend(y_strip_atoms);
        }

        previous_wave_atoms[current_wave_idx] = current_wave;

        // Progress reporting
        if (k + 1) % 20 == 0 || k + 1 == grid.grid_dims.2 {
            let elapsed = wave_start.elapsed();
            println!("Processed wave {}/{} in {:.2}s ({} atoms so far)",
                     k + 1, grid.grid_dims.2, elapsed.as_secs_f64(), all_atoms.len());
        }
    }

    println!("Generated {} liquid atoms", all_atoms.len());
    all_atoms
}

// Check collision with atoms in previous 1-2 Z layers
fn check_collision_with_previous_waves(
    i: usize,
    j: usize,
    previous_waves: &[Vec<Vec<Option<Point3<f64>>>>],
    voxel_size: f64,
    candidate_pos: &Point3<f64>,
) -> bool {
    let min_distance = voxel_size * 0.8; // Minimum separation

    // Check 3x3 neighborhood in previous 1-2 layers
    for wave in previous_waves {
        for dj in -1..=1_i32 {
            for di in -1..=1_i32 {
                let ni = (i as i32 + di).max(0) as usize;
                let nj = (j as i32 + dj).max(0) as usize;

                if ni < wave.len() && nj < wave[0].len() {
                    if let Some(other_pos) = wave[nj][ni] {
                        let distance = (candidate_pos - other_pos).norm();
                        if distance < min_distance {
                            return true;
                        }
                    }
                }
            }
        }
    }

    false
}

// Check collision with nanoparticle voxels (much more accurate than sphere)
fn check_collision_with_nanoparticle_voxels(
    candidate_pos: &Point3<f64>,
    grid: &NanoparticleGrid,
) -> bool {
    // Convert position to voxel coordinates
    let voxel_i = (candidate_pos.x / grid.voxel_size).floor() as usize;
    let voxel_j = (candidate_pos.y / grid.voxel_size).floor() as usize;
    let voxel_k = (candidate_pos.z / grid.voxel_size).floor() as usize;

    // Check if this voxel or nearby voxels contain nanoparticle atoms
    for dk in -1..=1_i32 {
        for dj in -1..=1_i32 {
            for di in -1..=1_i32 {
                let ni = (voxel_i as i32 + di).max(0) as usize;
                let nj = (voxel_j as i32 + dj).max(0) as usize;
                let nk = (voxel_k as i32 + dk).max(0) as usize;

                if ni < grid.grid_dims.0 && nj < grid.grid_dims.1 && nk < grid.grid_dims.2 {
                    if grid.is_voxel_occupied(ni, nj, nk) {
                        return true;
                    }
                }
            }
        }
    }

    false
}