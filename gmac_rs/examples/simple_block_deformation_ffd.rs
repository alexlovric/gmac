use gmac::core::{primitives::generate_box};
use gmac::core::transformation::{build_transformation_matrix, transform_node};
use gmac::io::stl::write_stl;
use gmac::morph::{ffd::FreeFormDeformer, design_block::DesignBlock};
use gmac::io::vtk::write_vtu;

fn main() {
    // 1. Create a base geometry - a 3D box centered at origin
    //    Dimensions: 1.0 x 1.0 x 1.0
    //    Discretization: 5x5x5 points
    let mut geometry =
        generate_box([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [5, 5, 5]);

    // 2. Create a control lattice (design block) that will be used to deform the geometry
    //    The design block is slightly larger than the geometry to ensure full coverage
    //    Dimensions: 0.8 x 1.2 x 1.2
    //    Offset from origin: [0.2, 0.0, 0.0]
    //    Control points: 2x2x2 grid
    let design_block =
        DesignBlock::new([0.8, 1.2, 1.2], [0.2, 0.0, 0.0], [0.0, 0.0, 0.0], [2, 2, 2]);

    // 3. Select which control points will be free to move during deformation
    //    The second parameter (Some(2)) specifies the number of layers of control points
    //    that will be kept fixed at the intersecting boundaries
    let free_design_ids = design_block
        .select_free_design_nodes(&geometry, Some(2))
        .unwrap();

    // 4. Create a transformation matrix that will be applied to the control points
    //    - Translation: [0.25, 0.0, 0.0]
    //    - Rotation: 45 degrees around X-axis
    //    - Scale: [1.0, 1.5, 1.5]
    let transform_matrix =
        build_transformation_matrix([0.25, 0.0, 0.0], [45.0, 0.0, 0.0], [1.0, 1.5, 1.5]);

    // 5. Create a copy of the original control points to modify
    let mut deformed_design_nodes = design_block.nodes.clone();

    // 6. Apply the transformation to each free control point
    free_design_ids.iter().for_each(|&id| {
        transform_node(
            &mut deformed_design_nodes[id],
            &transform_matrix,
            &[0.2, 0., 0.],
        )
    });

    // 7. Create a deformed version of the design block with the transformed control points
    let ffd = FreeFormDeformer::new(design_block);

    // 8. Apply the FFD deformation to the original geometry using the modified control lattice
    geometry.nodes = ffd.deform(&geometry.nodes, &deformed_design_nodes).unwrap();

    // 9. Save the deformed geometry to an STL file or points
    write_vtu(
        &geometry.nodes,
        &geometry.cells,
        Some("target/deformed.vtu"),
    )
    .unwrap();

    write_stl(
        &geometry.nodes,
        &geometry.cells,
        Some("target/deformed.stl"),
    )
    .unwrap();
}
