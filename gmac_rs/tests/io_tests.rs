use gmac::{
    core::primitives::{generate_cone, generate_torus},
    io::{
        obj::{read_obj, write_obj},
        stl::{read_stl, write_stl, StlFormat},
    },
};
use std::fs::remove_file;

#[test]
fn test_obj_write_read_roundtrip() {
    // 1. Setup: Generate a known piece of geometry.
    let original_mesh = generate_cone(1.0, 2.0, [0.0, 0.0, 0.0], 16).unwrap();
    let filename = "test_cone.obj";

    // 2. Action: Write the mesh to an OBJ file.
    let write_result =
        write_obj(&original_mesh.nodes, &original_mesh.cells, Some(filename));
    assert!(write_result.is_ok(), "Failed to write OBJ file");

    // 3. Action: Read the mesh back from the file.
    let (read_nodes, read_cells) =
        read_obj(filename).expect("Failed to read OBJ file back");

    // 4. Verification: Check if the data survived the roundtrip.
    assert_eq!(
        original_mesh.nodes.len(),
        read_nodes.len(),
        "Node count should match after read/write cycle"
    );
    assert_eq!(
        original_mesh.cells.len(),
        read_cells.len(),
        "Cell count should match after read/write cycle"
    );

    // Clean up the test file.
    remove_file(filename).unwrap();
}

#[test]
fn test_stl_write_read_roundtrip() {
    // 1. Setup: Generate a known piece of geometry.
    let original_mesh = generate_torus(1.0, 0.5, [0.0, 0.0, 0.0], 16, 16).unwrap();
    let filename = "test_cone.stl";

    // 2. Action: Write the mesh to an OBJ file.
    let write_result = write_stl(
        &original_mesh.nodes,
        &original_mesh.cells,
        Some(filename),
        Some(StlFormat::Binary),
    );
    assert!(write_result.is_ok(), "Failed to write STL file");

    // 3. Action: Read the mesh back from the file.
    let (read_nodes, read_cells) =
        read_stl(filename).expect("Failed to read STL file back");

    // 4. Verification: Check if the data survived the roundtrip.
    assert_eq!(
        original_mesh.nodes.len(),
        read_nodes.len(),
        "Node count should match after read/write cycle"
    );
    assert_eq!(
        original_mesh.cells.len(),
        read_cells.len(),
        "Cell count should match after read/write cycle"
    );

    // Clean up the test file.
    remove_file(filename).unwrap();
}
