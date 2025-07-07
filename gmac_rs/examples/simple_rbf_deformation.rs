//! Simple Block Deformation using Radial Basis Functions (RBF)
//!
//! This example demonstrates how to use the GMAC library to perform mesh deformation
//! using Radial Basis Functions (RBF). RBF is a powerful technique for creating smooth
//! deformations by interpolating the movement of control points throughout the mesh.

use std::time::Instant;

use gmac::core::{
    clusters::generate_block_cluster,
    primitives::generate_torus,
    selection::select_nodes_in_plane_direction,
    transformation::{build_transformation_matrix, transform_selected_nodes},
};
use gmac::io::{stl::StlFormat, vtk::write_vtp};
use gmac::morph::rbf::RbfDeformer;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_timer = Instant::now();

    // Create a simple geometry or import one
    let mut mesh = generate_torus(1.0, 0.3, [0.0, 0.0, 0.0], 64, 64)?;

    // Alternative: Mesh::from_stl("path_to_stl")? or Mesh::from_obj("path_to_obj")?

    // Save the original geometry as an STL file
    mesh.write_stl(Some("target/original_torus.stl"), None)?;

    // Create a set of control points that will be used to drive the deformation
    // These points form a lattice around the geometry
    let original_control_points = generate_block_cluster(
        [2.8, 1.0, 2.8], // Slightly larger than the original box
        [0.0, 0.0, 0.0], // Centered at the origin
        [0.0, 0.0, 0.0], // No rotation
        [3, 3, 3],       // Number of control points in each direction
    )?;

    // Save the original control points for visualization
    write_vtp(
        &original_control_points,
        Some("target/original_control_points.vtp"),
    )?;

    // Select control points that lie in a plane defined by an origin point and normal vector
    let target_control_point_ids = select_nodes_in_plane_direction(
        &original_control_points,
        [0.0, 0.0, 0.0], // A point in the plane
        [1.0, 0.0, 0.0], // Normal vector (x-axis in this case)
    );

    // Apply transformation to the selected control points
    let transform_matrix = build_transformation_matrix(
        [0.0, 0.0, 0.0],   // Translation
        [45.0, 20.0, 0.0], // Rotation
        [1.0, 1.0, 1.0],   // Scale
    );

    let deformed_control_points = transform_selected_nodes(
        &original_control_points,
        &target_control_point_ids,
        &transform_matrix,
        &[0.0, 0.0, 0.0],
    );

    // Save the deformed control points for visualization
    write_vtp(
        &deformed_control_points,
        Some("target/deformed_control_points.vtp"),
    )?;

    // Create an RBF deformer with the control point configurations
    let deformation_timer = Instant::now();

    let rbf = RbfDeformer::new(
        original_control_points,
        deformed_control_points,
        Some("thin_plate_spline"), // Type of RBF kernel
        Some(1.0),                 // Shape parameter for the RBF kernel
    )?;

    // Apply the RBF deformation to the original geometry
    mesh.nodes = rbf.deform(&mesh.nodes)?;

    let elapsed = deformation_timer.elapsed();
    println!("Deformation took: {}ms", elapsed.as_millis());

    // Save the final deformed geometry as an STL file
    mesh.write_stl(Some("target/deformed_torus.stl"), Some(StlFormat::Binary))?;

    let elapsed = start_timer.elapsed();
    println!("Total took: {}ms", elapsed.as_millis());

    Ok(())
}
