use gmac_core::{primitives::generate_box};
use gmac_core::transformation::{build_transformation_matrix, transform_node};
use gmac_io::stl::write_stl;
use gmac_morph::{ffd::FreeFormDeformer, design_block::DesignBlock};
use gmac_io::vtk::write_vtu;

fn main() {
    let mut geometry =
        generate_box([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [5, 5, 5]);

    let design_block =
        DesignBlock::new([0.8, 1.2, 1.2], [0.2, 0.0, 0.0], [0.0, 0.0, 0.0], [2, 2, 2]);

    let mut deformed_design_nodes = design_block.nodes.clone();

    let free_design_ids = design_block
        .select_free_design_nodes(&geometry, Some(2))
        .unwrap();

    let transform_matrix =
        build_transformation_matrix([0.25, 0.0, 0.0], [45.0, 0.0, 0.0], [1.0, 1.5, 1.5]);

    free_design_ids.iter().for_each(|&id| {
        transform_node(
            &mut deformed_design_nodes[id],
            &transform_matrix,
            &[0.2, 0., 0.],
        )
    });

    let ffd = FreeFormDeformer::new(design_block);

    geometry.nodes = ffd.deform(&geometry.nodes, &deformed_design_nodes).unwrap();

    write_vtu(
        &geometry.nodes,
        &geometry.cells,
        Some("target/deformed.vtu"),
    )
    .unwrap();

    write_stl(&geometry.nodes, &geometry.cells, Some("target/deformed.stl")).unwrap();
}
