# NanoEmbed - Nanoparticle Composite Generator

A Rust program for generating LAMMPS input files with nanoparticles embedded in liquid matrices. Uses precise convex hull collision detection and supports both collision-free and overlapping placement modes.

## Features

- **Convex hull collision detection**: Precise nanoparticle boundaries using Rapier3D physics engine
- **Flexible placement modes**: Collision (no overlaps) or overlap (first placed wins)
- **TOML configuration**: Simple configuration files instead of complex command-line arguments
- **Multiple nanoparticle types**: Support for different particle types with individual counts
- **Quasi-random positioning**: Even distribution using Halton sequences
- **Atomic radius expansion**: Proper separation based on covalent radii
- **Mass-based atom types**: Consistent atom type assignment ordered by atomic mass

## Quick Start

```bash
# Build the project
cargo build --release

# Run with a configuration file
./target/release/nanoembed config.toml
```

## Configuration File Format

Create a TOML configuration file (e.g., `config.toml`):

```toml
[system]
box_size = 250.0          # Cubic box size in Angstroms
output_file = "system.data"
random_seed = 12345       # Optional: for reproducible results

[background]
material = "liquid"       # liquid, fcc, or none
element = "Co"           # Background element
density = 0.04           # Atoms per Å³
separation = 2.0         # Minimum separation between nanoparticles (Å)

[[nanoparticles]]
file = "np/np_1.xyz"
count = 5

[[nanoparticles]]
file = "np/np_2.xyz"
count = 10

[analysis]
validate = true          # Run nearest-neighbor analysis
quasi_random = true      # Use quasi-random positioning for even distribution

[placement]
mode = "collision"       # "collision" (no overlaps) or "overlap" (first wins)
```

## Placement Modes

### Collision Mode (Default)
- **No overlaps allowed**: Nanoparticles cannot overlap
- **Strict separation**: Maintains minimum distance between surfaces
- **Higher packing challenge**: May fail to place all particles in dense systems

### Overlap Mode
- **Overlaps allowed**: Nanoparticles can overlap
- **First placed wins**: Atoms from earlier-placed nanoparticles take precedence
- **Higher packing density**: Can achieve denser packings
- **Atom-level resolution**: Only overlapping atoms are removed, not entire particles

## Input Format

Nanoparticle files should be in standard XYZ format:
```
8240
Generated nanoparticle
C    20.205    27.316    26.547
W    21.664    26.474    25.127
...
```

## Key Technical Features

- **Convex hull collision detection**: Uses Rapier3D physics engine for precise geometric calculations
- **Rounded convex hulls**: Expands boundaries by atomic radii to prevent unphysical overlaps
- **Density-adaptive voxel sizing**: Automatically adjusts spatial resolution based on target density
- **Mass-ordered atom types**: Assigns atom types consistently ordered by atomic mass
- **Timing analysis**: Detailed performance breakdown for optimization

## Output

The program generates:
1. **LAMMPS data file**: Standard format ready for molecular dynamics simulation
2. **Timing breakdown**: Performance analysis of each generation stage
3. **Validation analysis**: Nearest-neighbor distance statistics (optional)
4. **Placement statistics**: Success rates and packing information