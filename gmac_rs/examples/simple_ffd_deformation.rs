//! Simple Block Deformation using Free-Form Deformation (FFD)
//!
//! This example demonstrates how to use the GMAC library to perform Free-Form Deformation (FFD)
//! on a simple box geometry. FFD is a technique used to deform solid geometry in a smooth way
//! by manipulating a grid of control points that enclose the geometry.

use std::time::Instant;

use gmac::{
    core::{
        primitives::generate_box,
        transformation::{build_transformation_matrix, transform_selected_nodes},
    },
    io::{stl::StlFormat, vtk::write_vtp},
    morph::{design_block::DesignBlock, ffd::FreeFormDeformer},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_timer = Instant::now();

    // Create a simple box geometry
    let mut mesh = generate_box(
        [2.0, 1.0, 1.0], // Dimensions (length, width, height)
        [0.0, 0.0, 0.0], // Center coordinates
        [0.0, 0.0, 0.0], // Rotation angles (degrees)
        [12, 12, 12],    // Number of divisions in each direction
    )?;

    // Alternative: Mesh::from_stl("path_to_stl")? or Mesh::from_obj("path_to_obj")?

    // Save the original geometry as an STL file
    mesh.write_stl(Some("target/original.stl"), Some(StlFormat::Binary))?;

    // Create a design block (control lattice) for FFD
    // The design block defines the control points that will be used to deform the geometry
    let design_block = DesignBlock::new(
        [1.8, 1.2, 1.2], // Dimensions of the control lattice
        [0.2, 0.0, 0.0], // Center offset
        [0.0, 0.0, 0.0], // Rotation angles (degrees)
        [3, 2, 2],       // Number of control points in each direction
    )?;

    write_vtp(&design_block.nodes, Some("target/design_block.vtp"))?;

    // Select which control points will be free to move during deformation
    // The second parameter (Some(2)) specifies the number of fixed layers of control points at the boundaries
    let free_design_ids = design_block.select_free_design_nodes(&mesh, Some(2))?;

    // Create a transformation matrix that combines translation, rotation, and scaling
    let transform_matrix = build_transformation_matrix(
        [0.25, 0.0, 0.0],
        [125.0, 0.0, 0.0],
        [1.0, 1.25, 1.25],
    );

    // Apply the transformation to each free control point
    let deformed_design_nodes = transform_selected_nodes(
        &design_block.nodes,
        &free_design_ids,
        &transform_matrix,
        &[0.2, 0.0, 0.0], // Pivot point
    );

    // Create a Free-Form Deformer with the original design block
    let deformation_timer = Instant::now();

    let ffd = FreeFormDeformer::new(design_block)?;

    // Apply the deformation to the original geometry using the deformed control points
    mesh.nodes = ffd.deform(&mesh.nodes, &deformed_design_nodes)?;

    let elapsed = deformation_timer.elapsed();
    println!("Deformation took: {}ms", elapsed.as_millis());

    write_vtp(
        &deformed_design_nodes,
        Some("target/deformed_design_block.vtp"),
    )?;

    // Save the final deformed geometry as an STL file
    mesh.write_stl(Some("target/deformed.stl"), Some(StlFormat::Binary))?;

    let elapsed = start_timer.elapsed();
    println!("Total took: {}ms", elapsed.as_millis());

    Ok(())
}
