# NanoEmbed - Nanoparticle Composite Generator

A Rust program for generating LAMMPS input files with nanoparticles embedded in liquid matrices. Features custom rigid body dynamics with periodic boundary conditions for optimal packing.

## Features

### Core Capabilities
- **Periodic boundary conditions (PBC)**: Full support for periodic systems with minimum image convention
- **Dynamics-based relaxation**: Energy minimization with FIRE algorithm for optimal packing
- **Surface sphere discretization**: Efficient collision detection using discretized convex hulls
- **Flexible placement modes**: Collision (no overlaps) or overlap (first placed wins)
- **Multiple nanoparticle types**: Support for different particle types with individual counts
- **Quasi-random positioning**: Even distribution using Halton sequences

### Dynamics System
- **Custom rigid body engine**: Quaternion-based rotation, PBC-aware from the ground up
- **Soft sphere (WCA) potential**: Smooth repulsive interactions for realistic packing
- **Hard sphere potential**: Strict no-overlap constraints with harmonic penalties
- **FIRE minimization**: Fast Inertial Relaxation Engine for efficient energy minimization
- **Comprehensive testing**: 200+ unit tests covering all physics modules

### Technical Features
- **TOML configuration**: Simple configuration files instead of complex command-line arguments
- **Atomic radius expansion**: Proper separation based on covalent radii
- **Mass-based atom types**: Consistent atom type assignment ordered by atomic mass
- **Detailed diagnostics**: Energy, force, and convergence monitoring

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

# Optional: Enable dynamics-based relaxation with PBC
[dynamics]
enabled = true           # Enable energy minimization (recommended!)
max_iterations = 10000   # Maximum minimization steps
force_tolerance = 0.01   # Convergence threshold (kcal/mol/Å)
use_fire = true          # Use FIRE algorithm (faster than steepest descent)

[dynamics.potential.softsphere]
epsilon = 1.0            # Energy scale (kcal/mol)
sigma = 3.0              # Length scale (Å)

[dynamics.surface]
sphere_radius = 2.0      # Surface collision sphere radius (Å)
target_spacing = 3.0     # Spacing between surface spheres (Å)
```

See `config_example_dynamics.toml` for a complete example with dynamics enabled.

## Two-Phase Workflow

### Phase 1: Nanoparticle Placement & Relaxation
1. **Initial placement**: Quasi-random or random positioning with collision checking
2. **Surface discretization**: Convert convex hulls to collision detection spheres
3. **Energy minimization** (if dynamics enabled):
   - Apply FIRE algorithm with periodic boundary conditions
   - Remove overlaps and optimize packing
   - Monitor energy, forces, and convergence
4. **Result**: Well-packed nanoparticles with no unphysical overlaps

### Phase 2: Binder Filling
1. **Voxel-based generation**: Fill gaps between nanoparticles
2. **Collision avoidance**: Use voxel occupancy to avoid nanoparticle regions
3. **Density control**: Achieve target binder density in free volume
4. **Result**: Complete composite with nanoparticles + binder matrix

## Placement Modes

### Collision Mode (Default)
- **No overlaps allowed**: Nanoparticles cannot overlap
- **Strict separation**: Maintains minimum distance between surfaces
- **Works with dynamics**: Initial placement may have some overlaps, but minimization removes them
- **Higher packing challenge**: May fail to place all particles in dense systems without dynamics

### Overlap Mode
- **Overlaps allowed**: Nanoparticles can overlap during initial placement
- **First placed wins**: Atoms from earlier-placed nanoparticles take precedence
- **Higher packing density**: Can achieve denser packings
- **Atom-level resolution**: Only overlapping atoms are removed, not entire particles
- **With dynamics**: Overlaps are resolved during energy minimization

## Periodic Boundary Conditions

When dynamics are enabled, the system uses full periodic boundary conditions:

- **Minimum image convention**: Distances calculated using nearest periodic image
- **Position wrapping**: Particles wrap around box boundaries
- **No surface effects**: Eliminates artifacts from box edges
- **Improved packing**: Allows particles near boundaries to interact correctly

This is especially important for achieving uniform, realistic packing throughout the simulation box.

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

### Dynamics Engine (New!)
- **Custom rigid body physics**: Built from scratch for nanoparticle packing, no external physics engine
- **Periodic boundary conditions**: Full PBC support with minimum image convention
- **Surface sphere discretization**: Fibonacci/icosphere sampling of convex hull surfaces
- **WCA potential**: Weeks-Chandler-Andersen soft sphere for smooth, realistic interactions
- **FIRE algorithm**: Fast Inertial Relaxation Engine for efficient energy minimization
- **Quaternion rotations**: Robust 3D rotations without gimbal lock
- **Comprehensive tests**: 200+ unit tests with full coverage

### Core Features
- **Convex hull representation**: Precise nanoparticle boundaries from atom positions
- **Density-adaptive voxel sizing**: Automatically adjusts spatial resolution based on target density
- **Mass-ordered atom types**: Assigns atom types consistently ordered by atomic mass
- **Timing analysis**: Detailed performance breakdown for optimization
- **Detailed diagnostics**: Energy, force, and convergence monitoring during minimization

For technical details on the dynamics system, see [DYNAMICS.md](DYNAMICS.md).

## Output

The program generates:
1. **LAMMPS data file**: Standard format ready for molecular dynamics simulation
2. **Timing breakdown**: Performance analysis of each generation stage
3. **Validation analysis**: Nearest-neighbor distance statistics (optional)
4. **Placement statistics**: Success rates and packing information