use approx::assert_abs_diff_eq;
use gmac::{
    core::{
        clusters::generate_block_cluster,
        primitives::{generate_box, generate_torus},
        selection::select_nodes_in_plane_direction,
        transformation::{build_transformation_matrix, transform_selected_nodes},
    },
    error::Result,
    morph::{design_block::DesignBlock, ffd::FreeFormDeformer, rbf::RbfDeformer},
};

#[test]
fn test_full_ffd_workflow() -> Result<()> {
    let mesh =
        generate_box([2.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [6, 6, 6])?;

    let design_block =
        DesignBlock::new([1.8, 1.2, 1.2], [0.2, 0.0, 0.0], [0.0, 0.0, 0.0], [3, 2, 2])?;

    let free_design_ids = design_block.select_free_design_nodes(&mesh, Some(2))?;

    let transform_matrix = build_transformation_matrix(
        [0.25, 0.0, 0.0],
        [125.0, 0.0, 0.0],
        [1.0, 1.25, 1.25],
    );

    let deformed_design_nodes = transform_selected_nodes(
        &design_block.nodes,
        &free_design_ids,
        &transform_matrix,
        &[0.2, 0.0, 0.0],
    );

    let ffd = FreeFormDeformer::new(design_block)?;
    let deformed_nodes = ffd.deform(&mesh.nodes, &deformed_design_nodes)?;

    // Check that the workflow ran without errors and produced a result.
    assert_eq!(
        mesh.nodes.len(),
        deformed_nodes.len(),
        "Deformed mesh should have the same number of nodes as the original"
    );

    // Check that the deformation actually changed the nodes.
    // This is a crucial sanity check.
    assert_ne!(
        mesh.nodes, deformed_nodes,
        "Deformed nodes should not be identical to original nodes"
    );

    // Check that the number of deformed control points is correct.
    assert_eq!(
        deformed_design_nodes.len(),
        ffd.original_design_block.nodes.len()
    );

    // Check the results of the deformation
    let expected = [
        [1.247770919067215, -0.455213671485261, -0.0647521143485552],
        [1.247770919067215, -0.624348717209654, -0.1816956943850295],
        [1.247770919067215, 0.273517977100231, 0.689100831558209],
        [1.247770919067215, 0.104382931375838, 0.572157251521735],
        [1.247770919067215, -0.064752114348555, 0.455213671485260],
        [1.247770919067215, -0.233887160072948, 0.338270091448786],
        [1.247770919067215, -0.403022205797342, 0.221326511412312],
        [1.247770919067215, -0.572157251521735, 0.104382931375838],
        [1.247770919067215, -0.741292297246128, -0.0125606486606363],
        [1.247770919067215, 0.156574397063757, 0.858235877282603],
        [1.247770919067215, -0.0125606486606361, 0.741292297246128],
        [1.247770919067215, -0.181695694385029, 0.624348717209654],
        [1.247770919067215, -0.350830740109423, 0.507405137173179],
        [1.247770919067215, -0.519965785833816, 0.390461557136706],
        [1.247770919067215, -0.689100831558209, 0.273517977100231],
        [1.247770919067215, -0.858235877282603, 0.156574397063757],
    ];

    // Create a slice that views only the last 16 nodes.
    let last_nodes = &deformed_nodes[deformed_nodes.len() - 16..];

    // Loop through the slice and the expected values.
    for (actual_node, expected_node) in last_nodes.iter().zip(expected.iter()) {
        assert_abs_diff_eq!(actual_node[0], expected_node[0], epsilon = 1e-6);
        assert_abs_diff_eq!(actual_node[1], expected_node[1], epsilon = 1e-6);
        assert_abs_diff_eq!(actual_node[2], expected_node[2], epsilon = 1e-6);
    }

    Ok(())
}

#[test]
fn test_full_rbf_workflow() -> Result<()> {
    let mesh = generate_torus(1.0, 0.3, [0.0, 0.0, 0.0], 16, 16)?;

    let original_control_points = generate_block_cluster(
        [2.8, 1.0, 2.8],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [3, 3, 3],
    )?;

    let target_control_point_ids = select_nodes_in_plane_direction(
        &original_control_points,
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    );

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

    let rbf = RbfDeformer::new(
        original_control_points.clone(),
        deformed_control_points,
        Some("thin_plate_spline"),
        None,
    )?;

    let deformed_nodes = rbf.deform(&mesh.nodes)?;

    assert_eq!(
        mesh.nodes.len(),
        deformed_nodes.len(),
        "Deformed mesh should have the same number of nodes as the original"
    );

    // A crucial sanity check: the deformation should have actually moved the points.
    assert_ne!(
        mesh.nodes, deformed_nodes,
        "Deformed nodes should not be identical to original nodes"
    );

    // Check the results of the deformation
    let expected = [
        [1.201026741426032, 0.380021673595413, -0.214829056838749],
        [1.141532148939275, 0.425340988838948, -0.115588654691320],
        [1.051384534466215, 0.438513431817737, -0.027920917559857],
        [0.943929453262135, 0.421054044627878, 0.0273690851087761],
        [0.835688218383731, 0.377937044106749, 0.0368490174994807],
        [0.743627819277862, 0.315298215076788, -0.000120481024923],
        [0.682103703154963, 0.239098872069635, -0.070472682956937],
        [0.660485561666371, 0.157818421694127, -0.157818421694128],
        [0.682103703154958, 0.085877360537827, -0.254503549650527],
        [0.743627819277855, 0.037258033085639, -0.352435767137506],
        [0.835688218383720, 0.018489563067904, -0.433275624674130],
        [0.943929453262126, 0.030858091509537, -0.479281221246205],
        [1.051384534466198, 0.072850002397361, -0.483442516655241],
        [1.141532148939266, 0.139597019481615, -0.449349353629244],
        [1.201026741426018, 0.222344786447333, -0.387537403204021],
        [1.221727043405326, 0.307899641641036, -0.307899641641045],
    ];

    // Create a slice that views only the last 16 nodes.
    let last_nodes = &deformed_nodes[deformed_nodes.len() - 16..];

    // Loop through the slice and the expected values.
    for (actual_node, expected_node) in last_nodes.iter().zip(expected.iter()) {
        assert_abs_diff_eq!(actual_node[0], expected_node[0], epsilon = 1e-6);
        assert_abs_diff_eq!(actual_node[1], expected_node[1], epsilon = 1e-6);
        assert_abs_diff_eq!(actual_node[2], expected_node[2], epsilon = 1e-6);
    }

    Ok(())
}
