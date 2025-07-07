#![allow(unused_imports)]

use std::time::Instant;

use gmac::{
    core::{
        clusters::generate_sphere_cluster,
        mesh::Mesh,
        primitives::{
            generate_box, generate_capsule, generate_cone, generate_cylinder,
            generate_torus, generate_uvsphere,
        },
    },
    io::{
        obj::{read_obj, write_obj},
        stl::{read_stl, StlFormat},
        vtk::{read_vtu, write_vtp, write_vtu},
    },
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_timer = Instant::now();

    // Create a simple box geometry with specified dimensions, center, orientation, and resolution
    // let _geometry = generate_cylinder(0.5, 1.0, [0.0, 0.0, 0.0], 16, 16)?;
    // let _geometry = generate_torus(1.0, 0.5, [0.0, 0.0, 0.0], 16, 16)?;
    // let _geometry = generate_box(
    //     [1.0, 1.0, 1.0],
    //     [1.0, 1.0, 1.0],
    //     [0.0, 0.0, 0.0],
    //     [100, 100, 100],
    // )?;

    // let obj = read_stl("assets/temp/car/Audi_R8_2017.stl")?;
    // let obj = read_obj("assets/temp/airplane/11803_Airplane_v1_l1.obj")?;
    // let obj = read_obj("assets/temp/hand/16834_hand_v1_NEW.obj")?;

    let mesh = Mesh::from_stl("assets/temp/car/Audi_R8_2017.stl")?;

    // mesh.write_stl(Some("target/hand.stl"), Some(StlFormat::Ascii))?;

    // write_obj(&mesh.nodes, &mesh.cells, Some("target/hand.obj"))?;

    // let _sphere = generate_sphere_cluster(1.0, [0.0, 0.0, 0.0], 16)?;

    // _geometry.write_stl(Some("target/box.stl"), Some(StlFormat::Ascii))?;

    write_vtu(&mesh.nodes, &mesh.cells, Some("target/car.vtu"))?;

    let vtu = read_vtu("target/car.vtu")?;

    let mesh2 = Mesh::new(vtu.0, vtu.1);

    mesh2.write_stl(Some("target/car2.stl"), Some(StlFormat::Ascii))?;

    write_vtp(&mesh.nodes, Some("target/car.vtp"))?;
    let elapsed = start_timer.elapsed();
    println!("Total took: {}ms", elapsed.as_millis());

    Ok(())
}
