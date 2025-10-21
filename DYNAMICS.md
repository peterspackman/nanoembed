# Custom Rigid Body Dynamics System

This document describes the new custom rigid body dynamics system implemented for nanoembed. This replaces the Rapier3D-based collision detection with a cleaner, PBC-aware approach specifically designed for nanoparticle packing simulations.

## Overview

The new system consists of several modular components:

1. **Periodic Boundary Conditions (PBC)** - Clean implementation of minimum image convention
2. **Surface Sphere Discretization** - Convert convex hulls to collision spheres
3. **Rigid Body Mechanics** - Quaternion-based rotation and translation
4. **Potentials** - Soft sphere (WCA) and hard sphere interaction potentials
5. **Force Calculation** - PBC-aware pairwise force and torque computation
6. **Energy Minimization** - FIRE and steepest descent algorithms

## Architecture

### Module Structure

```
src/
├── periodic.rs              # PBC utilities (wrap, minimum image)
├── potentials.rs            # Soft/hard sphere potentials
├── rigid_body.rs            # Rigid body state and operations
├── surface_discretization.rs # Convert hull → surface spheres
├── forces.rs                # Force/torque calculation with PBC
└── minimization.rs          # Energy minimization (FIRE, steepest descent)
```

### Key Data Structures

#### RigidBody
Represents a nanoparticle as a rigid body:
```rust
struct RigidBody {
    // Translational state
    position: Point3<f64>,
    velocity: Vector3<f64>,

    // Rotational state
    orientation: UnitQuaternion<f64>,
    angular_velocity: Vector3<f64>,

    // Physical properties
    mass: f64,
    inertia_tensor: Matrix3<f64>,

    // Geometry
    surface_spheres: Vec<SurfaceSphere>,
    hull_points: Vec<Point3<f64>>,
}
```

#### SurfaceSphere
Discretized collision detection sphere on nanoparticle surface:
```rust
struct SurfaceSphere {
    local_position: Point3<f64>,  // In body frame
    radius: f64,
}
```

## Workflow

### 1. Initial Placement
- Load nanoparticle geometries from XYZ files
- Place nanoparticles using quasi-random or random positioning
- Convert each nanoparticle's convex hull to surface spheres

### 2. Surface Discretization
- Sample points uniformly on convex hull surface
- Create collision spheres at each point with specified radius
- Remove duplicates (spheres too close together)

Options:
- **Fibonacci sphere**: Uniform distribution on sphere, projected to hull
- **Icosphere**: Subdivision of icosahedron for regular sampling

### 3. Energy Minimization
Run iterative minimization to remove overlaps and optimize packing:

#### FIRE Algorithm (Recommended)
Fast Inertial Relaxation Engine - more efficient than steepest descent:
- Mixes molecular dynamics with damping
- Adaptive timestep and mixing parameter
- Converges faster for most systems

#### Steepest Descent
Simpler gradient descent alternative:
- Move particles in direction of force
- Adaptive step size based on energy changes
- More robust but slower

### 4. Periodic Boundary Conditions
All distance calculations use minimum image convention:
```rust
// Shortest vector between two points
let r_vec = minimum_image_vector(pos2 - pos1, box_size);

// Wrap positions into [0, box_size)
pos = wrap_position(pos, box_size);
```

### 5. Force Calculation
For each pair of rigid bodies:
1. Broad phase: Check if centers of mass are within cutoff
2. Narrow phase: Loop over all surface sphere pairs
3. Calculate distance with PBC (minimum image)
4. Evaluate potential (soft or hard sphere)
5. Compute force and torque on each body

### 6. Binder Filling
After nanoparticle relaxation:
- Use voxel-based approach (existing system)
- Fill gaps with binder atoms (e.g., Co)
- Avoid occupied regions using voxel grid

## Potentials

### Soft Sphere (WCA)
Weeks-Chandler-Andersen potential - purely repulsive Lennard-Jones:

```
U(r) = 4ε[(σ/r)^12 - (σ/r)^6] + ε   for r < cutoff
F(r) = 24ε/r * [2(σ/r)^12 - (σ/r)^6]
```

Parameters:
- `epsilon`: Energy scale (kcal/mol)
- `sigma`: Length scale / equilibrium distance (Å)
- `cutoff`: Typically 2^(1/6) * sigma ≈ 1.122 * sigma

**Use when**: You want smooth, physically realistic repulsion

### Hard Sphere
Harmonic repulsion for overlapping spheres:

```
overlap = sigma - r
U(r) = 0.5 * k * overlap^2  for r < sigma
F(r) = k * overlap
```

Parameters:
- `penalty_strength`: Spring constant (kcal/mol/Å²)

**Use when**: You want strict no-overlap constraint with strong penalty

## Configuration

### Basic Dynamics Setup

```toml
[dynamics]
enabled = true                # Enable dynamics
max_iterations = 10000        # Max minimization steps
force_tolerance = 0.01        # Stop when F_max < this
energy_tolerance = 1e-6       # Stop when |ΔE| < this
use_fire = true               # Use FIRE algorithm

[dynamics.potential.softsphere]
epsilon = 1.0
sigma = 3.0

[dynamics.surface]
sphere_radius = 2.0           # Collision sphere radius
target_spacing = 3.0          # Spacing between spheres
```

### Advanced Options

**Maximum displacement constraints** (prevent instabilities):
```toml
max_displacement = 0.5        # Max position change per step (Å)
max_rotation = 0.1            # Max rotation per step (radians)
```

**Hard sphere potential** (alternative):
```toml
[dynamics.potential.hardsphere]
penalty_strength = 100.0      # Higher = stiffer
```

## Testing

All modules include comprehensive unit tests:

```bash
# Test specific module
cargo test periodic
cargo test potentials
cargo test rigid_body
cargo test forces
cargo test minimization

# Test all dynamics modules
cargo test --lib

# Run with output
cargo test -- --nocapture
```

### Test Coverage

- **periodic.rs**: PBC wrapping, minimum image, boundary cases
- **potentials.rs**: WCA/hard sphere energy and forces, cutoffs
- **rigid_body.rs**: Quaternion operations, inertia, PBC updates
- **surface_discretization.rs**: Sphere generation, deduplication
- **forces.rs**: Pairwise forces, torques, PBC interactions
- **minimization.rs**: Convergence, energy minimization, constraints

## Implementation Status

### ✅ Completed
- [x] PBC utilities with comprehensive tests
- [x] Soft sphere (WCA) and hard sphere potentials
- [x] Rigid body data structures and operations
- [x] Surface sphere discretization (Fibonacci + icosphere)
- [x] Force and torque calculation with PBC
- [x] Energy minimization (FIRE + steepest descent)
- [x] Configuration system for dynamics parameters
- [x] Module integration into build system

### 🚧 In Progress
- [ ] Integration into main placement workflow
- [ ] Example usage and documentation
- [ ] Performance optimization (cell lists for force calculation)

### 📋 Future Enhancements
- [ ] Confining potential for pressure control
- [ ] Cell lists / spatial hashing for O(N) force calculation
- [ ] Parallel force calculation (Rayon)
- [ ] More sophisticated surface discretization
- [ ] Adaptive surface sphere sizing
- [ ] Energy/force diagnostics and logging
- [ ] Visualization tools for debugging

## Performance Considerations

### Current Complexity
- Force calculation: O(N²M²) where N = # bodies, M = # surface spheres
- Can be expensive for many nanoparticles with many surface spheres

### Optimization Strategies
1. **Cell lists / Neighbor lists**: Reduce to O(NM) for dense systems
2. **Adaptive spacing**: Use coarser spheres for large particles
3. **Broad phase culling**: Skip distant body pairs
4. **Parallelization**: Use Rayon for independent force calculations

### Recommended Settings
For good performance/accuracy balance:
- `sphere_radius`: 2-3 Å
- `target_spacing`: 3-5 Å (coarser for large particles)
- `force_tolerance`: 0.01-0.1 kcal/mol/Å
- `use_fire`: true (much faster than steepest descent)

## Usage Example

```rust
use nanoembed::*;

// Create rigid bodies from nanoparticle geometries
let mut bodies: Vec<RigidBody> = /* ... */;

// Discretize surfaces
for body in &mut bodies {
    body.surface_spheres = surface_discretization::discretize_convex_hull(
        &body.hull_points,
        config.dynamics.surface.sphere_radius,
        config.dynamics.surface.target_spacing,
    );
}

// Set up potential
let potential = match &config.dynamics.potential {
    PotentialType::SoftSphere { epsilon, sigma } => {
        Potential::SoftSphere(SoftSphereParams::new_wca(*epsilon, *sigma))
    }
    PotentialType::HardSphere { penalty_strength } => {
        Potential::HardSphere(HardSphereParams::new(*penalty_strength))
    }
};

// Minimize energy
let params = MinimizerParams {
    max_iterations: config.dynamics.max_iterations,
    force_tolerance: config.dynamics.force_tolerance,
    energy_tolerance: config.dynamics.energy_tolerance,
    use_fire: config.dynamics.use_fire,
    ..Default::default()
};

let result = minimization::minimize_energy(
    &mut bodies,
    box_size,
    &potential,
    &params,
);

println!("Minimization result: {:?}", result);
```

## Advantages Over Rapier3D

1. **Clean PBC support**: Built-in from the start, no hacks
2. **Full control**: Custom force models, integration schemes
3. **Transparency**: Simple, readable code without complex physics engine
4. **Flexibility**: Easy to add new potentials, constraints, etc.
5. **Lighter weight**: Only what we need, no unused features
6. **Better testing**: Targeted unit tests for each component

## References

### Algorithms
- **FIRE**: Bitzek et al., Phys. Rev. Lett. 97, 170201 (2006)
- **WCA Potential**: Weeks, Chandler, Andersen, J. Chem. Phys. 54, 5237 (1971)
- **Minimum Image**: Allen & Tildesley, "Computer Simulation of Liquids"

### Quaternion Mechanics
- nalgebra documentation: https://www.nalgebra.org/
- Quaternion calculus for rigid body dynamics

## Contact & Contributions

For questions, issues, or contributions, please refer to the main README.
