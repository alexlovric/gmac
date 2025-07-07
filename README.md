# GMAC
[![Build & Test](https://github.com/alexlovric/gmac/actions/workflows/build&test.yml/badge.svg?branch=main)](https://github.com/alexlovric/gmac/actions/workflows/build&test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A fast geometry manipulation and creation library made in rust, with a convenient python interface, and very few dependencies. Primary features include:
- Deform geometries using RBF and FFD
- Transform geometries (or selection just a selection of nodes)
- Large range of selection and transformation tools
- Import/export stl, obj and vtk-type files
- Convenient python interface (gmac_py)
- Create primitives
- Great performance

Here's a demonstration of a plane tail deformed using the Free Form deformer (FFD):

| Variation 1 | Variation 2 |
|-------------|-------------|
| <img src="https://github.com/alexlovric/gmac/blob/main/assets/plane_tail_variationa1.png?raw=true" /> | <img src="https://github.com/alexlovric/gmac/blob/main/assets/plane_tail_variationa2.png?raw=true" /> |


## Add to your rust project
Add the following to your `Cargo.toml`:

```toml
[dependencies]
gmac = "^0.2.0" # includes rayon by default
```

For those who want a lightweight dependency free version use `default-features = false`

For **Rbf** specifically, if you still need more performance try adding features `openblas` or `intel-mkl`.
Make sure you have the required dependencies installed for the features you choose. For openblas openssl is required.

## Examples in Rust
### Deforming with Free Form Deformer (FFD)
Heres a demonstration of deformation using the `FreeFormDeformer`:

```rust
use gmac::{
    core::{
        primitives::generate_box,
        transformation::{build_transformation_matrix, transform_node},
    },
    io::stl::StlFormat,
    morph::{ffd::FreeFormDeformer, design_block::DesignBlock},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a simple primitive or import something else
    let mut geometry = generate_box(
        [2.0, 1.0, 1.0],  // Lengths (x, y, z)
        [0.0, 0.0, 0.0],  // Center coordinates
        [0.0, 0.0, 0.0],  // Rotation angles (degrees)
        [12, 12, 12],     // Elements in each direction
    )?;
    
    // Alternative: geometry = Mesh::from_stl("path_to_stl")? or from_obj

    // Create a design block (control lattice) for FFD
    let design_block = DesignBlock::new(
        [1.8, 1.2, 1.2],  // Lengths (x, y, z)
        [0.2, 0.0, 0.0],  // Center offset
        [0.0, 0.0, 0.0],  // Rotation angles (degrees)
        [3, 2, 2],        // Control points each direction
    )?;

    // Select which control points will be free to move during 
    // deformation, the second parameter specifies the number of 
    // fixed layers of control points at the boundaries
    let free_design_ids = design_block
        .select_free_design_nodes(&geometry, Some(2))?;

    // Create a transformation matrix
    let transform_matrix = build_transformation_matrix(
        [0.25, 0.0, 0.0],   // Translation vector (x, y, z)
        [125.0, 0.0, 0.0],  // Rotation angles (degrees)
        [1.0, 1.25, 1.25],  // Scaling factors (x, y, z)
    );

    // Create a copy of the original control points to modify
    let mut deformed_design_nodes = design_block.nodes.clone();

    // Apply the transformation to each free control point
    free_design_ids.iter().for_each(|&id| {
        transform_node(
            &mut deformed_design_nodes[id], // Control points
            &transform_matrix,              // Transformation
            &[0.2, 0., 0.],                 // Origin or pivot point
        )
    });

    // Create a Free-Form Deformer with original design block
    let ffd = FreeFormDeformer::new(design_block);

    // Apply the deformation to the original geometry
    geometry.nodes = ffd.deform(&geometry.nodes, &deformed_design_nodes)?;

    // Save the final deformed geometry as an STL/OBJ/VTK file
    geometry.write_stl(Some("deformed.stl"), Some(StlFormat::Binary))?;

    Ok(())
}
```
| Original mesh & control points | Deformed mesh & control points |
|-------------------------|-------------------------|
| <img src="https://github.com/alexlovric/gmac/blob/main/assets/example_1_original.png?raw=true" /> | <img src="https://github.com/alexlovric/gmac/blob/main/assets/example_1_deformed.png?raw=true" width="99%" /> |


### Deforming with Radial Basis Functions (RBF)
Using the `RbfDeformer` is very similar to using the FFD, but instead of using a design block (control lattice), you use a set of control points:

```rust
use gmac::core::{
    clusters::generate_block_cluster, transformation::transform_node,
    primitives::generate_box, transformation::build_transformation_matrix,
    selection::select_nodes_in_plane_direction,
};
use gmac::io::{stl::StlFormat, vtk::write_vtp};
use gmac::morph::rbf::RbfDeformer;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a simple geometry or import one
    let mut geometry = generate_torus(1.0, 0.3, [0.0, 0.0, 0.0], 64, 64)?;

    // Create a set of control points
    let original_control_points = generate_block_cluster(
        [2.8, 1.0, 2.8], // Slightly larger than the original box
        [0.0, 0.0, 0.0], // Centered at the origin
        [0.0, 0.0, 0.0], // No rotation
        [3, 3, 3],       // Number of control points in each direction
    )?;

    // Select control points that lie in a plane defined 
    let target_control_point_ids = select_nodes_in_plane_direction(
        &original_control_points,
        [0.0, 0.0, 0.0], // A point in the plane
        [1.0, 0.0, 0.0], // Normal vector (x-axis in this case)
    );

    // Apply transformation to the selected control points
    let transform_matrix = build_transformation_matrix(
        [0.0, 0.0, 0.0],   
        [45.0, 20.0, 0.0], 
        [1.0, 1.0, 1.0],   
    );

    // Apply the transformation to each selected control point
    let mut deformed_control_points = original_control_points.clone();
    for &id in &target_control_point_ids {
        let mut point = deformed_control_points[id];
        transform_node(
            &mut point,
            &transform_matrix,
            &[0.0, 0.0, 0.0], // Origin point
        );
        deformed_control_points[id] = point;
    }

    // Create an RBF deformer with the control point configurations
    let rbf = RbfDeformer::new(
        original_control_points,
        deformed_control_points,
        Some("thin_plate_spline"), // Type of RBF kernel
        Some(1.0),                 // Shape parameter
    )?;

    // Apply the RBF deformation to the original geometry
    geometry.nodes = rbf.deform(&geometry.nodes)?;

    // Save the final deformed geometry as an STL/OBJ/VTK file
    geometry.write_stl(Some("deformed.stl"), None)?;

    Ok(())
}
```

Here you can see the original control points and mesh, as well as the deformed control points and mesh:

| Original mesh & control points | Deformed mesh & control points |
|-------------------------|-------------------------|
| <img src="https://github.com/alexlovric/gmac/blob/main/assets/example_2_original.png?raw=true" /> | <img src="https://github.com/alexlovric/gmac/blob/main/assets/example_2_deformed.png?raw=true" width="99%"/> |


## Examples in Python
### Using GMAC to deform a box

Heres a simple demonstration of using the Free Form Deformer (FFD) to deform a generated box (or an imported STL file).

```python
import gmac
import gmac.morph as morph
import gmac.io as io
import numpy as np

# Create a simple box geometry
geometry = gmac.generate_box(
    [1.0, 1.0, 1.0],  # Dimensions (length, width, height)
    [0.0, 0.0, 0.0],  # Center coordinates
    [0.0, 0.0, 0.0],  # Rotation angles (degrees)
    [5, 5, 5]         # Number of divisions in each direction
)

# Or import one from stl
# geometry = gmac.Mesh.from_stl("path_to_stl")

# Create a design block (control lattice) for FFD
design_block = morph.DesignBlock(
    [0.8, 1.2, 1.2],  # Dimensions of the control lattice
    [0.2, 0.0, 0.0],  # Center offset
    [0.0, 0.0, 0.0],  # Rotation angles (degrees)
    [2, 2, 2]         # Number of control points in each direction
)

# Select which control points will be free to move
free_design_ids = design_block.select_free_design_nodes(geometry, 2)

# Create a transformation matrix
transformation_matrix = gmac.build_transformation_matrix(
    [0.25, 0.0, 0.0],   # Translation vector (x, y, z)
    [45.0, 0.0, 0.0],   # Rotation angles (degrees)
    [1.0, 1.5, 1.5]     # Scaling factors (x, y, z)
)

# Transform only the free control points
deformed_design_nodes = np.array(design_block.nodes)
deformed_design_nodes[free_design_ids] = gmac.transform_nodes(
    deformed_design_nodes[free_design_ids],
    transformation_matrix,
    [0.2, 0., 0.],
)

# Save the deformed control points for visualization
io.write_vtp(deformed_design_nodes, "deformed_design_nodes.vtp")

# Create a Free-Form Deformer with the original design block
ffd = morph.FreeFormDeformer(design_block)

# Apply the deformation to the original geometry
geometry.nodes = ffd.deform(geometry.nodes, deformed_design_nodes)

# Save the final deformed geometry as an STL file (default binary)
geometry.write_stl("deformed_geometry.stl") # binary default
```

For Radial Basis Function (RBF) deformation, see the RbfDeformer example in examples.

## Build python from source
These instructions assume that Python3 and Cargo are installed on your system. To set up this project, follow these steps:
1. Clone the repository:
    ```bash
    git clone https://github.com/alexlovric/gmac.git
    cd gmac/gmac_py
    ```
2. Create a virtual environment and install build system:
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate # In windows /Scripts/activate
    python3 -m pip install -r requirements.txt
    ```
3. Build the release binary:
    ```bash
    maturin develop --release
    ```
4. Build the python wheel:
    ```bash
    maturin build --release
    ```
5. Running examples:
    ```bash
    python3 -m pip install <path to wheel (target/wheels/*.whl)>
    cd examples
    python3 -m pip install -r requirements.txt
    python3 *simple_ffd_deformation*.py
    ```

## References

The gmac_morph is heavily influenced by PyGEM (https://github.com/mathLab/PyGeM), and the following

Sieger, Menzel, Botsch. *On Shape Deformation Techniques for Simulation-based Design Optimization.* SEMA SIMAI Springer Series, 2015.

Lombardi, Parolini, Quarteroni, Rozza. *Numerical Simulation of Sailing Boats: Dynamics, FSI, and Shape Optimization.* Springer Optimization and Its Applications, 2012.

## License

MIT License - See [LICENSE](LICENSE) for details.

## Support
If you'd like to support the project consider:
- Identifying the features you'd like to see implemented or bugs you'd like to fix and open an issue.
- Contributing to the code by resolving existing issues, I'm happy to have you.
- Donating to help me continue development, [Buy Me a Coffee](https://coff.ee/alexlovric)