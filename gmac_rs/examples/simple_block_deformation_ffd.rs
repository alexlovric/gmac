//! Simple Block Deformation using Free-Form Deformation (FFD)
//!
//! This example demonstrates how to use the GMAC library to perform Free-Form Deformation (FFD)
//! on a simple box geometry. FFD is a technique used to deform solid geometry in a smooth way
//! by manipulating a grid of control points that enclose the geometry.

use gmac::{
    core::{
        primitives::generate_box,
        transformation::{build_transformation_matrix, transform_node},
    },
    io::{stl::StlFormat, vtk::write_vtu},
    morph::{ffd::FreeFormDeformer, design_block::DesignBlock},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a simple box geometry with specified dimensions, center, orientation, and resolution
    let mut geometry = generate_box(
        [1.0, 1.0, 1.0], // Dimensions (length, width, height)
        [0.0, 0.0, 0.0], // Center coordinates
        [0.0, 0.0, 0.0], // Rotation angles (degrees)
        [5, 5, 5],       // Number of divisions in each direction
    );

    // Alternative: Load geometry from an STL file
    // let mut geometry = Mesh::from_stl("original.stl")?;

    // Save the original geometry as an STL file
    geometry.write_stl(Some("target/original.stl"), Some(StlFormat::Binary))?;

    // Create a design block (control lattice) for FFD
    // The design block defines the control points that will be used to deform the geometry
    let design_block = DesignBlock::new(
        [0.8, 1.2, 1.2], // Dimensions of the control lattice
        [0.2, 0.0, 0.0], // Center offset
        [0.0, 0.0, 0.0], // Rotation angles (degrees)
        [2, 2, 2],       // Number of control points in each direction
    );

    // Select which control points will be free to move during deformation
    // The second parameter (Some(2)) specifies the number of fixed layers of control points at the boundaries
    let free_design_ids = design_block.select_free_design_nodes(&geometry, Some(2))?;

    // Create a transformation matrix that combines translation, rotation, and scaling
    let transform_matrix = build_transformation_matrix(
        [0.25, 0.0, 0.0], // Translation vector (x, y, z)
        [45.0, 0.0, 0.0], // Rotation angles (degrees) around x, y, z axes
        [1.0, 1.5, 1.5],  // Scaling factors in x, y, z directions
    );

    // Create a copy of the original control points to modify
    let mut deformed_design_nodes = design_block.nodes.clone();

    // Apply the transformation to each free control point
    free_design_ids.iter().for_each(|&id| {
        transform_node(
            &mut deformed_design_nodes[id], // The control point to transform
            &transform_matrix,              // The transformation to apply
            &[0.2, 0., 0.],                 // Pivot point for transformations
        )
    });

    // Create a Free-Form Deformer with the original design block
    let ffd = FreeFormDeformer::new(design_block);

    // Apply the deformation to the original geometry using the deformed control points
    geometry.nodes = ffd.deform(&geometry.nodes, &deformed_design_nodes)?;

    // Save the deformed geometry as a VTK file for visualization
    write_vtu(
        &geometry.nodes,
        &geometry.cells,
        Some("target/deformed.vtu"),
    )?;

    // Save the final deformed geometry as an STL file

    geometry.write_stl(Some("target/deformed.stl"), Some(StlFormat::Ascii))?;
    // or
    // write_stl(
    //     &geometry.nodes,
    //     &geometry.cells,
    //     Some("target/deformed.stl"),
    //     None,
    // )?;

    Ok(())
}
