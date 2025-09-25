# NanoEmbed - Efficient Nanoparticle Composite Generator

A high-performance Rust program for generating LAMMPS input files with randomly distributed nanoparticles in liquid or solid matrices. Designed to scale efficiently to millions of atoms.

## Features

- **Voxel-based collision detection**: O(1) collision checking using 3D spatial hashing
- **Convex hull geometry**: Accurate exclusion volumes for complex nanoparticle shapes
- **Multiple atom types**: Preserves original atom types from nanoparticle files
- **Parallel processing**: Multi-threaded liquid generation for better performance
- **Flexible backgrounds**: Liquid (quasi-random) or FCC lattice filling
- **Rotation support**: Random orientation of nanoparticles for realistic packing
- **LAMMPS output**: Standard LAMMPS data file format

## Quick Start

```bash
# Build the project
cargo build --release

# Basic usage - place 5 nanoparticles in a 200x200x200 Å box
./target/release/nanoembed \
    -n np/np_1.xyz \
    -b 200 200 200 \
    -c 5 \
    -o system.data

# Advanced usage with custom density and voxel size
./target/release/nanoembed \
    -n np/np_1.xyz np/np_2.xyz \
    -b 300 300 400 \
    -c 20 \
    -d 0.05 \
    -v 0.5 \
    -s 3.0 \
    -m liquid \
    -r 12345
```

## Architecture

### Voxel-Based Spatial Partitioning
- Simulation box divided into uniform voxels (default 1.0 Å)
- Each voxel marked as occupied/empty based on convex hull tests
- Enables O(1) collision detection for efficient placement

### Convex Hull Collision Detection
- Nanoparticles represented by convex hulls computed from atom positions
- Point-in-hull tests using hyperplane inequalities
- Accurate exclusion volumes for complex geometries

### Efficient Background Generation
- Parallel voxel processing for liquid generation
- Poisson sampling within empty voxels at specified density
- Scales to millions of background atoms

## Parameters

- `--nanoparticles`: XYZ files containing nanoparticle geometries
- `--box-size`: Orthorhombic simulation box dimensions [Lx Ly Lz] in Angstroms
- `--count`: Number of nanoparticles to randomly place
- `--density`: Liquid density in atoms/Å³ (typical: 0.05-0.1)
- `--voxel-size`: Spatial resolution for collision detection (smaller = more accurate)
- `--separation`: Minimum distance between nanoparticle surfaces
- `--material`: Background type (liquid, fcc, or none)

## Input Format

Nanoparticle files should be in standard XYZ format:
```
62607

C    20.205    27.316    26.547
W    21.664    26.474    25.127
...
```

## Performance

- Tested with systems containing millions of atoms
- Parallel processing scales with CPU cores
- Memory usage scales linearly with system size
- Voxel-based approach much faster than naive O(n²) collision detection

## Output

LAMMPS data file with:
- Proper atom IDs and molecule groupings
- Multiple atom types (nanoparticle + background)
- Periodic boundary conditions
- Ready for MD simulation

## Example Systems

Generate a 500³ Å system with tungsten carbide nanoparticles in cobalt matrix:
```bash
./target/release/nanoembed \
    -n nanoparticles/wc_np1.xyz nanoparticles/wc_np2.xyz \
    -b 500 500 500 \
    -c 50 \
    -d 0.08 \
    -s 2.0 \
    -m liquid \
    -o wc_co_composite.data
```