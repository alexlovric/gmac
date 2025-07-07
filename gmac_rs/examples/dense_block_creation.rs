#![allow(unused_imports)]

use std::time::Instant;

use gmac::{
    core::{clusters::generate_sphere_cluster, primitives::{generate_box, generate_capsule, generate_cone, generate_cylinder, generate_torus, generate_uvsphere}},
    io::{stl::StlFormat, vtk::write_vtp},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_timer = Instant::now();

    // Create a simple box geometry with specified dimensions, center, orientation, and resolution
    let _geometry = generate_cylinder(0.5, 1.0, [0.0, 0.0, 0.0], 16, 16)?;

    let _sphere = generate_sphere_cluster(1.0, [0.0, 0.0, 0.0], 16)?;

    let elapsed = start_timer.elapsed();
    println!("Total took: {}ms", elapsed.as_millis());

    _geometry.write_stl(Some("target/box.stl"), Some(StlFormat::Binary))?;

    write_vtp(&_sphere, Some("target/sphere.vtp"))?;

    Ok(())
}
