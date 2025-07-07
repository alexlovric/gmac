//! Simple Block Deformation using Radial Basis Functions (RBF)
//!
//! This example demonstrates how to use the GMAC library to perform mesh deformation
//! using Radial Basis Functions (RBF). RBF is a powerful technique for creating smooth
//! deformations by interpolating the movement of control points throughout the mesh.

use gmac::core::primitives::generate_torus;
use gmac::core::{
    clusters::generate_block_cluster, transformation::transform_node,
    transformation::build_transformation_matrix,
    selection::select_nodes_in_plane_direction,
};
use gmac::io::{stl::StlFormat, vtk::write_vtp};
use gmac::morph::rbf::RbfDeformer;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a simple geometry or import one
    let mut geometry = generate_torus(1.0, 0.3, [0.0, 0.0, 0.0], 64, 64)?;

    // Alternative: Mesh::from_stl("path_to_stl")? or Mesh::from_obj("path_to_obj")?

    // Save the original geometry as an STL file
    geometry.write_stl(Some("target/original_torus.stl"), None)?;

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

    // Create a copy of the original control points to modify
    let mut deformed_control_points = original_control_points.clone();

    // Apply transformation to the selected control points
    let transform_matrix = build_transformation_matrix(
        [0.0, 0.0, 0.0],   // Move 1 unit in x-direction
        [45.0, 20.0, 0.0],  // Rotate 45 degrees around x-axis
        [1.0, 1.0, 1.0], // Scale y and z dimensions to 75%
    );

    // Apply the transformation to each selected control point
    for &id in &target_control_point_ids {
        let mut point = deformed_control_points[id];
        transform_node(
            &mut point,
            &transform_matrix,
            &[0.0, 0.0, 0.0], // Origin point for the transformation
        );
        deformed_control_points[id] = point;
    }

    // Save the deformed control points for visualization
    write_vtp(
        &deformed_control_points,
        Some("target/deformed_control_points.vtp"),
    )?;

    // Create an RBF deformer with the control point configurations
    let rbf = RbfDeformer::new(
        original_control_points,
        deformed_control_points,
        Some("thin_plate_spline"), // Type of RBF kernel
        Some(1.0),        // Shape parameter for the RBF kernel
    )?;

    // Apply the RBF deformation to the original geometry
    geometry.nodes = rbf.deform(&geometry.nodes)?;

    // Save the final deformed geometry as an STL file
    geometry.write_stl(Some("target/deformed_torus.stl"), Some(StlFormat::Binary))?;

    Ok(())
}
