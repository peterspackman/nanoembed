use crate::types::{Atom, Nanoparticle};
use kiddo::{KdTree, SquaredEuclidean};

pub fn analyze_nearest_neighbors(
    nanoparticles: &[Nanoparticle],
    background: &[Atom],
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Analyzing nearest neighbor distances...");

    // Build k-d tree with all atom positions
    let mut kdtree: KdTree<f64, 3> = KdTree::new();
    let mut all_positions = Vec::new();
    let mut atom_types = Vec::new();

    let start_time = std::time::Instant::now();

    // Add nanoparticle atoms
    for np in nanoparticles {
        for atom in &np.atoms {
            let pos = [atom.position.x, atom.position.y, atom.position.z];
            kdtree.add(&pos, all_positions.len() as u64);
            all_positions.push(pos);
            atom_types.push(atom.atom_type);
        }
    }

    // Add background atoms
    for atom in background {
        let pos = [atom.position.x, atom.position.y, atom.position.z];
        kdtree.add(&pos, all_positions.len() as u64);
        all_positions.push(pos);
        atom_types.push(atom.atom_type);
    }

    println!("Built k-d tree with {} atoms in {:.2}s",
             all_positions.len(), start_time.elapsed().as_secs_f64());

    // Analyze nearest neighbors
    let analysis_start = std::time::Instant::now();
    let mut min_distances = Vec::new();
    let mut np_min_distances = Vec::new();
    let mut bg_min_distances = Vec::new();
    let mut close_pairs = Vec::new(); // Store details of very close atoms

    let total_np_atoms: usize = nanoparticles.iter().map(|np| np.atoms.len()).sum();

    for (i, pos) in all_positions.iter().enumerate() {
        let nearest = kdtree.nearest_n::<SquaredEuclidean>(pos, 2); // Get 2 nearest (self + closest)

        if nearest.len() > 1 {
            let distance = nearest[1].distance.sqrt(); // Convert from squared distance
            let nearest_idx = nearest[1].item as usize;

            min_distances.push(distance);

            // Store details of very close atoms
            if distance < 0.7 {
                let atom_type_1 = atom_types[i];
                let atom_type_2 = atom_types[nearest_idx];
                let is_np_1 = i < total_np_atoms;
                let is_np_2 = nearest_idx < total_np_atoms;

                close_pairs.push((
                    i, nearest_idx, distance,
                    atom_type_1, atom_type_2,
                    is_np_1, is_np_2,
                    *pos, all_positions[nearest_idx]
                ));
            }

            // Categorize by atom location
            if i < total_np_atoms {
                np_min_distances.push(distance);
            } else {
                bg_min_distances.push(distance);
            }
        }
    }

    println!("Analyzed {} nearest neighbors in {:.2}s",
             min_distances.len(), analysis_start.elapsed().as_secs_f64());

    // Calculate statistics
    if !min_distances.is_empty() {
        min_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let min_dist = min_distances[0];
        let max_dist = min_distances[min_distances.len() - 1];
        let avg_dist = min_distances.iter().sum::<f64>() / min_distances.len() as f64;
        let median_dist = min_distances[min_distances.len() / 2];
        let p95_dist = min_distances[(min_distances.len() as f64 * 0.95) as usize];

        println!("\n=== Nearest Neighbor Analysis ===");
        println!("Overall statistics:");
        println!("  Minimum distance: {:.3} Å", min_dist);
        println!("  Maximum distance: {:.3} Å", max_dist);
        println!("  Average distance: {:.3} Å", avg_dist);
        println!("  Median distance:  {:.3} Å", median_dist);
        println!("  95th percentile:  {:.3} Å", p95_dist);

        if !np_min_distances.is_empty() {
            np_min_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let np_avg = np_min_distances.iter().sum::<f64>() / np_min_distances.len() as f64;
            let np_min = np_min_distances[0];
            println!("Nanoparticle atoms:");
            println!("  Minimum distance: {:.3} Å", np_min);
            println!("  Average distance: {:.3} Å", np_avg);
        }

        if !bg_min_distances.is_empty() {
            bg_min_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let bg_avg = bg_min_distances.iter().sum::<f64>() / bg_min_distances.len() as f64;
            let bg_min = bg_min_distances[0];
            println!("Background atoms:");
            println!("  Minimum distance: {:.3} Å", bg_min);
            println!("  Average distance: {:.3} Å", bg_avg);
        }

        // Quality warnings and detailed analysis
        if min_dist < 0.5 {
            println!("\n⚠️  WARNING: Very close atoms detected (< 0.5 Å)");
            println!("   Consider increasing separation or reducing density");
        } else if min_dist < 0.7 {
            println!("\n⚠️  WARNING: Close atoms detected (< 0.7 Å)");
            println!("   System may need equilibration before MD");
        } else {
            println!("\n✓ System appears to have reasonable atomic spacing");
        }

        // Show details of close atom pairs
        if !close_pairs.is_empty() {
            println!("\n=== Close Atom Pairs (< 0.7 Å) ===");
            println!("Found {} close pairs:", close_pairs.len());

            // Sort by distance and show worst cases
            close_pairs.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
            let num_to_show = close_pairs.len().min(10);

            for (i, &(idx1, idx2, dist, type1, type2, is_np1, is_np2, pos1, pos2)) in
                close_pairs.iter().take(num_to_show).enumerate() {

                let loc1 = if is_np1 { "NP" } else { "BG" };
                let loc2 = if is_np2 { "NP" } else { "BG" };

                println!("  {}: Atoms {} ({}) and {} ({}) at {:.3} Å",
                        i + 1, idx1, loc1, idx2, loc2, dist);
                println!("      Types: {} - {}", type1, type2);
                println!("      Pos1: ({:.3}, {:.3}, {:.3})", pos1[0], pos1[1], pos1[2]);
                println!("      Pos2: ({:.3}, {:.3}, {:.3})", pos2[0], pos2[1], pos2[2]);
                println!();
            }

            if close_pairs.len() > num_to_show {
                println!("  ... and {} more close pairs", close_pairs.len() - num_to_show);
            }
        }
    }

    Ok(())
}
