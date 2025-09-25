use crate::types::{AtomTypeMap, Nanoparticle, Atom};
use nalgebra::Vector3;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_lammps_data(
    filename: &str,
    nanoparticles: &[Nanoparticle],
    background: &[Atom],
    box_size: Vector3<f64>,
    type_map: &AtomTypeMap,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Writing LAMMPS data file: {}", filename);

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    let total_atoms: usize = nanoparticles.iter().map(|np| np.atoms.len()).sum::<usize>()
                           + background.len();

    let num_types = type_map.next_type - 1;

    // LAMMPS header
    writeln!(writer, "LAMMPS data file - nanoparticle composite system")?;
    writeln!(writer)?;
    writeln!(writer, "{} atoms", total_atoms)?;
    writeln!(writer, "{} atom types", num_types)?;
    writeln!(writer)?;
    writeln!(writer, "0.0000000000 {:.10} xlo xhi", box_size.x)?;
    writeln!(writer, "0.0000000000 {:.10} ylo yhi", box_size.y)?;
    writeln!(writer, "0.0000000000 {:.10} zlo zhi", box_size.z)?;
    writeln!(writer)?;

    // Masses (ordered by ascending mass)
    writeln!(writer, "Masses")?;
    writeln!(writer)?;

    let mut type_mass_pairs: Vec<(u32, f64, String)> = Vec::new();
    for atom_type in 1..type_map.next_type {
        if let (Some(element), Some(mass)) = (
            type_map.type_to_element.get(&atom_type),
            type_map.type_to_mass.get(&atom_type)
        ) {
            type_mass_pairs.push((atom_type, *mass, element.clone()));
        }
    }

    // Sort by mass
    type_mass_pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    for (atom_type, mass, element) in type_mass_pairs {
        writeln!(writer, "{} {:.6}   # {}", atom_type, mass, element)?;
    }
    writeln!(writer)?;

    // Atoms
    writeln!(writer, "Atoms")?;
    writeln!(writer)?;

    let mut atom_id = 1;

    // Write nanoparticle atoms (atomic style)
    for np in nanoparticles {
        for atom in &np.atoms {
            writeln!(writer, "{} {} {:.10} {:.10} {:.10}",
                    atom_id, atom.atom_type,
                    atom.position.x, atom.position.y, atom.position.z)?;
            atom_id += 1;
        }
    }

    // Write background atoms (atomic style)
    for atom in background {
        writeln!(writer, "{} {} {:.10} {:.10} {:.10}",
                atom_id, atom.atom_type,
                atom.position.x, atom.position.y, atom.position.z)?;
        atom_id += 1;
    }

    writer.flush()?;
    println!("Successfully wrote {} atoms to LAMMPS file", total_atoms);
    Ok(())
}