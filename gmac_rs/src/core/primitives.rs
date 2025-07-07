use std::collections::HashMap;

use crate::{
    error::{Error, Result},
    core::mesh::Mesh,
    core::transformation::{build_rotation_matrix, rotate_node},
};

use std::f64::consts::PI;

/// Generates a high-quality sphere mesh using icosahedron subdivision.
///
/// # Arguments
/// * `radius`: The desired radius of the sphere.
/// * `centre`: The desired center of the sphere.
/// * `subdivisions`: The number of times to subdivide the mesh. 0 = base icosahedron.
///                  A value of 4 or 5 provides a nicely detailed sphere.
///
/// # Returns
/// Returns a `Result<Mesh>` instance containing the generated sphere.
pub fn generate_icosphere(
    radius: f64,
    centre: [f64; 3],
    subdivisions: u32,
) -> Result<Mesh> {
    if radius <= 0.0 {
        return Err(Error::MeshGeneration(
            "Radius must be positive.".to_string(),
        ));
    }

    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut nodes = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];

    let mut cells = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];

    let mut midpoint_cache = HashMap::new();

    let get_midpoint =
        |i: usize,
         j: usize,
         nodes: &mut Vec<[f64; 3]>,
         cache: &mut HashMap<(usize, usize), usize>| {
            let key = if i < j { (i, j) } else { (j, i) };
            *cache.entry(key).or_insert_with(|| {
                let a = nodes[i];
                let b = nodes[j];
                let new = [
                    (a[0] + b[0]) * 0.5,
                    (a[1] + b[1]) * 0.5,
                    (a[2] + b[2]) * 0.5,
                ];
                nodes.push(new);
                nodes.len() - 1
            })
        };

    for _ in 0..subdivisions {
        let mut new_cells = Vec::with_capacity(cells.len() * 4);
        for &[a, b, c] in &cells {
            let ab = get_midpoint(a, b, &mut nodes, &mut midpoint_cache);
            let bc = get_midpoint(b, c, &mut nodes, &mut midpoint_cache);
            let ca = get_midpoint(c, a, &mut nodes, &mut midpoint_cache);
            new_cells.push([a, ab, ca]);
            new_cells.push([b, bc, ab]);
            new_cells.push([c, ca, bc]);
            new_cells.push([ab, bc, ca]);
        }
        cells = new_cells;
    }

    // Normalize and scale each vertex, then offset by center
    for p in &mut nodes {
        let norm = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
            .sqrt()
            .max(f64::EPSILON);
        p[0] = (p[0] / norm) * radius + centre[0];
        p[1] = (p[1] / norm) * radius + centre[1];
        p[2] = (p[2] / norm) * radius + centre[2];
    }

    Ok(Mesh::new(nodes, cells))
}

/// Generates a sphere mesh using the UV mapping method.
///
/// # Arguments
/// * `radius`: The desired radius of the sphere.
/// * `centre`: The desired center of the sphere.
/// * `sectors`: The number of vertical divisions (like longitude lines). Must be 3 or more.
/// * `stacks`: The number of horizontal divisions (like latitude lines). Must be 2 or more.
///
/// # Returns
/// Returns a `Result<Mesh>` instance containing the generated sphere.
pub fn generate_uvsphere(
    radius: f64,
    centre: [f64; 3],
    sectors: u32,
    stacks: u32,
) -> Result<Mesh> {
    if sectors < 3 || stacks < 2 {
        return Err(Error::MeshGeneration(
            "Sectors must be 3 or more, and stacks must be 2 or more.".to_string(),
        ));
    }

    let mut nodes = Vec::new();
    let mut cells = Vec::new();

    let offset = |x: f64, y: f64, z: f64| -> [f64; 3] {
        [x + centre[0], y + centre[1], z + centre[2]]
    };

    // North pole
    nodes.push(offset(0.0, radius, 0.0));

    // Latitude rings
    for i in 1..stacks {
        let stack_angle = PI / 2.0 - (i as f64 / stacks as f64) * PI; // From +π/2 to -π/2
        let xy = radius * stack_angle.cos();
        let y = radius * stack_angle.sin();

        for j in 0..sectors {
            let sector_angle = (j as f64 / sectors as f64) * 2.0 * PI; // 0 to 2π
            let x = xy * sector_angle.cos();
            let z = xy * sector_angle.sin();
            nodes.push(offset(x, y, z));
        }
    }

    // South pole
    nodes.push(offset(0.0, -radius, 0.0));

    let n_pole_idx = 0;
    let s_pole_idx = nodes.len() - 1;

    // Connect top fan
    create_triangle_fan(&mut cells, n_pole_idx, 1, sectors, true);

    // Connect middle quads
    for i in 0..stacks - 2 {
        let ring1_idx = 1 + (i * sectors) as usize;
        let ring2_idx = 1 + ((i + 1) * sectors) as usize;
        stitch_vertex_rings(&mut cells, ring1_idx, ring2_idx, sectors, false);
    }

    // Connect bottom fan
    let bottom_ring_idx = (1 + (stacks - 2) * sectors) as usize;
    create_triangle_fan(&mut cells, s_pole_idx, bottom_ring_idx, sectors, false);

    Ok(Mesh::new(nodes, cells))
}

/// Creates a standard 6 sided box mesh.
///
/// # Arguments
/// * `length`: An array `[x, y, z]` that specifies the dimensions of the box.
/// * `centre`: An array `[x, y, z]` that specifies the origin of the box.
/// * `theta`: An array `[x, y, z]` that specifies the rotation of the box.
/// * `resolution`: An array `[x, y, z]` that specifies the number of nodes
///                 along each axis.
///
/// # Returns
/// Returns a `Result<Mesh>` containing the `Mesh` instance or an error message.
pub fn generate_box(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> Result<Mesh> {
    // Validation checks
    if length.iter().any(|&ll| ll <= 0.0) {
        return Err(Error::MeshGeneration(
            "Invalid length: dimensions must be positive.".to_string(),
        ));
    }
    if resolution.iter().any(|&rr| rr == 0) {
        return Err(Error::MeshGeneration(
            "Invalid resolution: dimensions must be non-zero.".to_string(),
        ));
    }

    let [rx, ry, rz] = resolution;

    let step = [
        length[0] / (rx as f64),
        length[1] / (ry as f64),
        length[2] / (rz as f64),
    ];

    let start = [
        centre[0] - length[0] * 0.5,
        centre[1] - length[1] * 0.5,
        centre[2] - length[2] * 0.5,
    ];

    let rotation_matrix = build_rotation_matrix(&theta);

    let n_nodes = (rx + 1) * (ry + 1) * (rz + 1) - (rx - 1) * (ry - 1) * (rz - 1);
    let mut nodes = vec![[0.0; 3]; n_nodes];

    let mut sides = [
        Vec::with_capacity((ry + 1) * (rz + 1)), // Min x
        Vec::with_capacity((ry + 1) * (rz + 1)), // Max x
        Vec::with_capacity((rx + 1) * (rz + 1)), // Min y
        Vec::with_capacity((rx + 1) * (rz + 1)), // Max y
        Vec::with_capacity((rx + 1) * (ry + 1)), // Min z
        Vec::with_capacity((rx + 1) * (ry + 1)), // Max z
    ];

    let mut node_count = 0;

    for ii in 0..=rx {
        for jj in 0..=ry {
            for kk in 0..=rz {
                // Add nodes on the boundary
                if ii == 0 || ii == rx || jj == 0 || jj == ry || kk == 0 || kk == rz {
                    nodes[node_count] = [
                        ii as f64 * step[0] + start[0],
                        jj as f64 * step[1] + start[1],
                        kk as f64 * step[2] + start[2],
                    ];

                    if theta.iter().any(|angle| *angle != 0.0) {
                        rotate_node(&mut nodes[node_count], &rotation_matrix, &centre);
                    }

                    if ii == 0 {
                        sides[0].push(node_count);
                    }
                    if ii == rx {
                        sides[1].push(node_count);
                    }
                    if jj == 0 {
                        sides[2].push(node_count);
                    } else if jj == ry {
                        sides[3].push(node_count);
                    }
                    if kk == 0 {
                        sides[4].push(node_count);
                    } else if kk == rz {
                        sides[5].push(node_count);
                    }

                    node_count += 1;
                }
            }
        }
    }

    // Cells
    let n_cells = 4 * (rx * ry + ry * rz + rx * rz);
    let mut cells: Vec<[usize; 3]> = Vec::with_capacity(n_cells);

    // X - min/max
    for jj in 0..resolution[1] {
        for kk in 0..resolution[2] {
            let p1 = sides[0][jj * (rz + 1) + kk];
            let p2 = sides[0][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[0][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[0][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p2, p3]);
            cells.push([p3, p2, p4]);

            let p1 = sides[1][jj * (rz + 1) + kk];
            let p2 = sides[1][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[1][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[1][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);
        }
    }

    // Y - min/max
    for jj in 0..resolution[0] {
        for kk in 0..resolution[2] {
            let p1 = sides[2][jj * (rz + 1) + kk];
            let p2 = sides[2][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[2][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[2][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);

            let p1 = sides[3][jj * (rz + 1) + kk];
            let p2 = sides[3][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[3][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[3][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p3, p2]);
            cells.push([p3, p4, p2]);
        }
    }

    // Z - min/max
    for jj in 0..resolution[0] {
        for kk in 0..resolution[1] {
            let p1 = sides[4][jj * (ry + 1) + kk];
            let p2 = sides[4][jj * (ry + 1) + (kk + 1)];
            let p3 = sides[4][(jj + 1) * (ry + 1) + kk];
            let p4 = sides[4][(jj + 1) * (ry + 1) + (kk + 1)];

            cells.push([p1, p2, p3]);
            cells.push([p3, p2, p4]);

            let p1 = sides[5][jj * (ry + 1) + kk];
            let p2 = sides[5][jj * (ry + 1) + (kk + 1)];
            let p3 = sides[5][(jj + 1) * (ry + 1) + kk];
            let p4 = sides[5][(jj + 1) * (ry + 1) + (kk + 1)];

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);
        }
    }

    Ok(Mesh::new(nodes, cells))
}

/// Generates a cylinder mesh with closed caps.
///
/// # Arguments
/// * `radius`: The radius of the cylinder.
/// * `height`: The height of the cylinder.
/// * `centre`: The center of the cylinder.
/// * `sectors`: The number of divisions around the circumference. Must be 3 or more.
/// * `stacks`: The number of divisions along the height. Must be 1 or more.
///
/// # Returns
/// A `Result<Mesh>` instance containing the generated cylinder.
pub fn generate_cylinder(
    radius: f64,
    height: f64,
    centre: [f64; 3],
    sectors: u32,
    stacks: u32,
) -> Result<Mesh> {
    if sectors < 3 || stacks < 1 {
        return Err(Error::MeshGeneration(
            "Sectors must be 3 or more, and stacks must be 1 or more.".to_string(),
        ));
    }

    let mut nodes = Vec::new();
    let mut cells = Vec::new();
    let half_height = height / 2.0;

    let offset = |x: f64, y: f64, z: f64| -> [f64; 3] {
        [x + centre[0], y + centre[1], z + centre[2]]
    };

    // Generate Nodes
    let bottom_center_idx = nodes.len();
    nodes.push(offset(0.0, -half_height, 0.0));
    let top_center_idx = nodes.len();
    nodes.push(offset(0.0, half_height, 0.0));

    for i in 0..=stacks {
        let y = -half_height + (i as f64 / stacks as f64) * height;
        for j in 0..sectors {
            let angle = (j as f64 / sectors as f64) * 2.0 * PI;
            let x = radius * angle.cos();
            let z = radius * angle.sin();
            nodes.push(offset(x, y, z));
        }
    }

    // Generate Cells
    let wall_node_idx = 2;

    create_triangle_fan(&mut cells, bottom_center_idx, wall_node_idx, sectors, true);

    let top_ring_idx = wall_node_idx + (stacks * sectors) as usize;
    create_triangle_fan(&mut cells, top_center_idx, top_ring_idx, sectors, false);

    for i in 0..stacks {
        let ring1_start = wall_node_idx + (i * sectors) as usize;
        let ring2_start = wall_node_idx + ((i + 1) * sectors) as usize;
        stitch_vertex_rings(&mut cells, ring1_start, ring2_start, sectors, false);
    }

    Ok(Mesh::new(nodes, cells))
}

/// Generates a torus (donut) mesh.
///
/// # Arguments
/// * `major_radius`: The distance from the center of the torus to the center of the tube.
/// * `minor_radius`: The radius of the tube itself.
/// * `centre`: The center of the torus.
/// * `major_segments`: The number of divisions around the main ring. Must be 3 or more.
/// * `minor_segments`: The number of divisions around the tube. Must be 3 or more.
///
/// # Returns
/// Returns a `Result<Mesh>` instance containing the generated torus.
pub fn generate_torus(
    major_radius: f64,
    minor_radius: f64,
    centre: [f64; 3],
    major_segments: u32,
    minor_segments: u32,
) -> Result<Mesh> {
    if major_segments < 3 || minor_segments < 3 {
        return Err(Error::MeshGeneration(
            "Major and minor segments must be 3 or more.".to_string(),
        ));
    }
    if major_radius <= 0.0 || minor_radius <= 0.0 {
        return Err(Error::MeshGeneration("Radii must be positive.".to_string()));
    }
    if minor_radius > major_radius {
        return Err(Error::MeshGeneration(
            "Minor radius cannot be greater than major radius.".to_string(),
        ));
    }

    let mut nodes = Vec::new();
    let mut cells = Vec::new();

    // Generate Nodes around the origin (0,0,0)
    for i in 0..=major_segments {
        let major_angle = (i as f64 / major_segments as f64) * 2.0 * PI;
        let maj_cos = major_angle.cos();
        let maj_sin = major_angle.sin();

        for j in 0..=minor_segments {
            let minor_angle = (j as f64 / minor_segments as f64) * 2.0 * PI;
            let min_cos = minor_angle.cos();
            let min_sin = minor_angle.sin();

            let x = (major_radius + minor_radius * min_cos) * maj_cos;
            let y = minor_radius * min_sin;
            let z = (major_radius + minor_radius * min_cos) * maj_sin;

            nodes.push([x, y, z]);
        }
    }

    // Generate Cells (Triangles)
    for i in 0..major_segments {
        for j in 0..minor_segments {
            let p1 = (i * (minor_segments + 1) + j) as usize;
            let p2 = (p1 + 1) as usize;
            let p3 = ((i + 1) * (minor_segments + 1) + j) as usize;
            let p4 = (p3 + 1) as usize;

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);
        }
    }

    for node in &mut nodes {
        node[0] += centre[0];
        node[1] += centre[1];
        node[2] += centre[2];
    }

    Ok(Mesh::new(nodes, cells))
}

/// Generates a cone mesh with a closed circular base.
///
/// # Arguments
/// * `radius`: The radius of the base of the cone.
/// * `height`: The height of the cone from the base to the apex.
/// * `centre`: The center of the cone.
/// * `sectors`: The number of divisions around the circumference. Must be 3 or more.
///
/// # Returns
/// A `Result<Mesh>` instance containing the generated cone.
pub fn generate_cone(
    radius: f64,
    height: f64,
    centre: [f64; 3],
    sectors: u32,
) -> Result<Mesh> {
    if sectors < 3 {
        return Err(Error::MeshGeneration(
            "Sectors must be 3 or more.".to_string(),
        ));
    }
    if radius <= 0.0 || height <= 0.0 {
        return Err(Error::MeshGeneration(
            "Radius and height must be positive.".to_string(),
        ));
    }

    let mut nodes = Vec::new();
    let mut cells = Vec::new();
    let half_height = height / 2.0;

    // Generate Nodes around the origin (0,0,0)
    nodes.push([0.0, half_height, 0.0]); // Apex, index 0
    let apex_idx = 0;

    nodes.push([0.0, -half_height, 0.0]); // Base center, index 1
    let base_center_idx = 1;

    let base_ring_idx = nodes.len();
    for i in 0..sectors {
        let angle = (i as f64 / sectors as f64) * 2.0 * PI;
        let x = radius * angle.cos();
        let z = radius * angle.sin();
        nodes.push([x, -half_height, z]);
    }

    // Generate Cells
    create_triangle_fan(&mut cells, base_center_idx, base_ring_idx, sectors, true);
    create_triangle_fan(&mut cells, apex_idx, base_ring_idx, sectors, false);

    for node in &mut nodes {
        node[0] += centre[0];
        node[1] += centre[1];
        node[2] += centre[2];
    }

    Ok(Mesh::new(nodes, cells))
}

/// Generates a capsule mesh, which is a cylinder with hemispherical ends.
///
/// # Arguments
/// * `radius`: The radius of the cylinder and hemispheres.
/// * `cylinder_height`: The height of the central cylindrical part.
/// * `centre`: The center of the capsule.
/// * `sectors`: The number of divisions around the circumference. Must be 3 or more.
/// * `stacks`: The number of horizontal divisions for EACH hemisphere. Must be 1 or more.
///
/// # Returns
/// A `Result<Mesh>` instance containing the generated capsule.
pub fn generate_capsule(
    radius: f64,
    cylinder_height: f64,
    centre: [f64; 3],
    sectors: u32,
    stacks: u32,
) -> Result<Mesh> {
    if sectors < 3 || stacks < 1 {
        return Err(Error::MeshGeneration(
            "Sectors must be 3 or more, and stacks must be 1 or more.".to_string(),
        ));
    }
    if radius <= 0.0 || cylinder_height < 0.0 {
        return Err(Error::MeshGeneration(
            "Radius must be positive and cylinder height must be non-negative."
                .to_string(),
        ));
    }

    let total_stacks = stacks * 2;
    let mut nodes = Vec::new();
    let mut cells = Vec::new();
    let half_cyl_height = cylinder_height / 2.0;

    // Generate Nodes (Stretched Sphere) around the origin
    for i in 0..=total_stacks {
        let stack_angle = PI / 2.0 - (i as f64 / total_stacks as f64) * PI;
        let r = radius * stack_angle.cos();
        let mut y = radius * stack_angle.sin();

        if y > 0.0 {
            y += half_cyl_height;
        } else {
            y -= half_cyl_height;
        }

        for j in 0..sectors {
            let sector_angle = (j as f64 / sectors as f64) * 2.0 * PI;
            let x = r * sector_angle.cos();
            let z = r * sector_angle.sin();
            nodes.push([x, y, z]);
        }
    }

    // Generate Cells
    for i in 0..total_stacks {
        for j in 0..sectors {
            let p1 = (i * sectors + j) as usize;
            let p2 = (i * sectors + (j + 1) % sectors) as usize;
            let p3 = ((i + 1) * sectors + j) as usize;
            let p4 = ((i + 1) * sectors + (j + 1) % sectors) as usize;

            cells.push([p1, p2, p4]);
            cells.push([p4, p3, p1]);
        }
    }

    for node in &mut nodes {
        node[0] += centre[0];
        node[1] += centre[1];
        node[2] += centre[2];
    }

    Ok(Mesh::new(nodes, cells))
}

/// Creates a new NACA wing mesh.
///
/// # Arguments
/// * `maximum_camber`: The maximum camber of the airfoil, usually expressed as a
///                     percentage of the chord length.
/// * `camber_distance`: The distance from the leading edge to the location of
///                      maximum camber, usually as a percentage of chord length.
/// * `maximum_thickness`: The maximum thickness of the airfoil, usually as a
///                        percentage of chord length.
/// * `n_points`: The number of points to generate for each half of the airfoil
///               (upper and lower).
/// * `wing_span`: A tuple representing the z-coordinates for the beginning and
///                end of the wing span.
///
/// # Returns
/// A `Result<Mesh>` struct containing the generated nodes and cells for the NACA wing mesh.
pub fn generate_naca_wing(
    maximum_camber: f64,
    camber_distance: f64,
    maximum_thickness: f64,
    n_points: usize,
    wing_span: (f64, f64),
) -> Result<Mesh> {
    let nodes_2d = generate_naca_profile(
        maximum_camber,
        camber_distance,
        maximum_thickness,
        n_points,
        false,
    );

    let nodes = nodes_2d
        .iter()
        .flat_map(|&[x, y, _]| vec![[x, y, wing_span.0], [x, y, wing_span.1]])
        .collect::<Vec<[f64; 3]>>();

    let mut cells = Vec::with_capacity(4 * (n_points - 1));

    for i in 0..(n_points) {
        let i0 = 2 * i;
        let i1 = i0 + 1;
        let i2 = i0 + 2;
        let i3 = i0 + 3;

        // Upper surface
        cells.push([i0, i2, i1]); // First triangle
        cells.push([i3, i1, i2]); // Second triangle
    }

    for i in (n_points + 1)..2 * (n_points) {
        let i0 = 2 * i;
        let i1 = i0 + 1;
        let i2 = i0 + 2;
        let i3 = i0 + 3;

        // Upper surface
        cells.push([i0, i2, i1]); // First triangle
        cells.push([i3, i1, i2]); // Second triangle
    }

    // Trailing edge
    cells.push([2 * (n_points + 1), 0, 2 * (n_points + 1) + 1]);
    cells.push([1, 2 * (n_points + 1) + 1, 0]);

    Ok(Mesh::new(nodes, cells))
}

fn generate_naca_profile(
    maximum_camber: f64,
    camber_distance: f64,
    maximum_thickness: f64,
    n_points: usize,
    finite_trailing_edge: bool,
) -> Vec<[f64; 3]> {
    let mm = 0.01 * maximum_camber;
    let pp = 0.10 * camber_distance;
    let tt = 0.01 * maximum_thickness;

    let a0 = 0.2969;
    let a1 = -0.1260;
    let a2 = -0.3516;
    let a3 = 0.2843;
    let a4 = if finite_trailing_edge {
        -0.1015
    } else {
        -0.1036
    };

    let mut xu = Vec::with_capacity(n_points + 1);
    let mut yu = Vec::with_capacity(n_points + 1);
    let mut xl = Vec::new();
    let mut yl = Vec::new();

    for i in (0..=n_points).rev() {
        let x = i as f64 / n_points as f64;
        let yt = 5.0
            * tt
            * (a0 * x.sqrt() + a1 * x + a2 * x.powi(2) + a3 * x.powi(3) + a4 * x.powi(4));

        let (yc, dyc_dx) = if pp == 0.0 {
            (0.0, 0.0) // Replace with whatever is appropriate when pp is zero
        } else if x <= pp {
            (
                mm / pp.powi(2) * x * (2.0 * pp - x),
                mm / pp.powi(2) * (2.0 * pp - 2.0 * x),
            )
        } else {
            (
                mm / (1.0 - pp).powi(2) * (1.0 - 2.0 * pp + x) * (1.0 - x),
                mm / (1.0 - pp).powi(2) * (2.0 * pp - 2.0 * x),
            )
        };

        let theta = dyc_dx.atan();

        xu.push(x - yt * theta.sin());
        yu.push(yc + yt * theta.cos());

        xl.push(x + yt * theta.sin());
        yl.push(yc - yt * theta.cos());
    }

    let combined: Vec<[f64; 3]> = xu
        .iter()
        .chain(xl.iter().skip(1))
        .zip(yu.iter().chain(yl.iter().skip(1)))
        .map(|(&x, &y)| [x - 0.25, y, 0.0])
        .collect();

    combined
}

/// Connects two rings of vertices into a quad mesh wall (e.g., for a cylinder or torus).
fn stitch_vertex_rings(
    cells: &mut Vec<[usize; 3]>,
    start_index_ring1: usize,
    start_index_ring2: usize,
    segments: u32,
    reverse_winding: bool,
) {
    for i in 0..segments {
        let p1 = start_index_ring1 + i as usize;
        let p2 = start_index_ring1 + ((i + 1) % segments) as usize;
        let p3 = start_index_ring2 + i as usize;
        let p4 = start_index_ring2 + ((i + 1) % segments) as usize;

        if reverse_winding {
            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);
        } else {
            cells.push([p1, p2, p4]);
            cells.push([p4, p3, p1]);
        }
    }
}

/// Connects a ring of vertices to a central point to form a fan of triangles (e.g., for a cone).
fn create_triangle_fan(
    cells: &mut Vec<[usize; 3]>,
    center_vertex_index: usize,
    ring_start_index: usize,
    segments: u32,
    reverse_winding: bool,
) {
    for i in 0..segments {
        let p1 = center_vertex_index;
        let p2 = ring_start_index + i as usize;
        let p3 = ring_start_index + ((i + 1) % segments) as usize;

        if reverse_winding {
            cells.push([p1, p3, p2]);
        } else {
            cells.push([p1, p2, p3]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_generate_icosphere() {
        // Test base icosahedron (0 subdivisions)
        let mesh = generate_icosphere(1.0, [0.0, 0.0, 0.0], 0).unwrap();
        assert_eq!(
            mesh.nodes.len(),
            12,
            "Base icosahedron should have 12 vertices"
        );
        assert_eq!(
            mesh.cells.len(),
            20,
            "Base icosahedron should have 20 faces"
        );

        // Test first subdivision
        let mesh = generate_icosphere(1.0, [0.0, 0.0, 0.0], 1).unwrap();
        let expected_nodes = 12 + 30; // 12 original + 30 edge midpoints
        let expected_cells = 20 * 4; // Each face becomes 4 smaller faces
        assert_eq!(mesh.nodes.len(), expected_nodes);
        assert_eq!(mesh.cells.len(), expected_cells);

        // Verify all nodes are on the sphere's surface
        for node in mesh.nodes {
            let dist = (node[0].powi(2) + node[1].powi(2) + node[2].powi(2)).sqrt();
            assert!(
                (dist - 1.0).abs() < EPSILON,
                "Node is not on the sphere surface"
            );
        }
    }

    #[test]
    fn test_generate_uvsphere() {
        let sectors = 16;
        let stacks = 8;
        let mesh = generate_uvsphere(1.0, [0.0, 0.0, 0.0], sectors, stacks).unwrap();

        // Expected counts: 2 poles + (stacks - 1) rings * sectors per ring
        let expected_nodes = 2 + (stacks - 1) * sectors;
        // Expected cells: 2 caps * sectors + middle quads * 2
        let expected_cells = 2 * sectors + (stacks - 2) * sectors * 2;

        assert_eq!(mesh.nodes.len() as u32, expected_nodes);
        assert_eq!(mesh.cells.len() as u32, expected_cells);

        // Verify all nodes are on the sphere's surface
        for node in mesh.nodes {
            let dist = (node[0].powi(2) + node[1].powi(2) + node[2].powi(2)).sqrt();
            assert!(
                (dist - 1.0).abs() < EPSILON,
                "Node is not on the sphere surface"
            );
        }
    }

    #[test]
    fn test_generate_cylinder() {
        let radius = 1.0;
        let height = 2.0;
        let centre = [0.0, 0.0, 0.0];
        let sectors = 16;
        let stacks = 4;
        let mesh = generate_cylinder(radius, height, centre, sectors, stacks).unwrap();

        // Expected counts: 2 center points + (stacks + 1) rings * sectors per ring
        let expected_nodes = 2 + (stacks + 1) * sectors;
        // Expected cells: 2 caps * sectors + wall quads * 2
        let expected_cells = 2 * sectors + stacks * sectors * 2;

        assert_eq!(mesh.nodes.len() as u32, expected_nodes);
        assert_eq!(mesh.cells.len() as u32, expected_cells);

        // Check that the top and bottom center nodes are at the correct height
        assert!((mesh.nodes[0][1] - (-height / 2.0)).abs() < EPSILON);
        assert!((mesh.nodes[1][1] - (height / 2.0)).abs() < EPSILON);
    }

    #[test]
    fn test_generate_torus() {
        let major_radius = 2.0;
        let minor_radius = 0.5;
        let centre = [0.0, 0.0, 0.0];
        let major_segments = 32;
        let minor_segments = 16;
        let mesh = generate_torus(
            major_radius,
            minor_radius,
            centre,
            major_segments,
            minor_segments,
        )
        .unwrap();

        // Expected counts for a torus with a seam
        let expected_nodes = (major_segments + 1) * (minor_segments + 1);
        let expected_cells = major_segments * minor_segments * 2;

        assert_eq!(mesh.nodes.len() as u32, expected_nodes);
        assert_eq!(mesh.cells.len() as u32, expected_cells);
    }

    #[test]
    fn test_generate_cone() {
        let radius = 1.0;
        let height = 2.0;
        let centre = [0.0, 0.0, 0.0];
        let sectors = 16;
        let mesh = generate_cone(radius, height, centre, sectors).unwrap();

        // Expected counts: 1 apex + 1 base center + `sectors` base ring vertices
        let expected_nodes = 2 + sectors;
        // Expected cells: `sectors` for the base + `sectors` for the walls
        let expected_cells = 2 * sectors;

        assert_eq!(mesh.nodes.len() as u32, expected_nodes);
        assert_eq!(mesh.cells.len() as u32, expected_cells);
    }

    #[test]
    fn test_generate_capsule() {
        let radius = 1.0;
        let cylinder_height = 2.0;
        let centre = [0.0, 0.0, 0.0];
        let sectors = 16;
        let stacks = 8; // Stacks per hemisphere
        let mesh =
            generate_capsule(radius, cylinder_height, centre, sectors, stacks).unwrap();

        let total_stacks = stacks * 2;

        // Expected counts based on the "stretched sphere" algorithm
        let expected_nodes = (total_stacks + 1) * sectors;
        let expected_cells = total_stacks * sectors * 2;

        assert_eq!(mesh.nodes.len() as u32, expected_nodes);
        assert_eq!(mesh.cells.len() as u32, expected_cells);
    }

    #[test]
    fn test_generate_box() {
        let uniform_box =
            generate_box([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1, 1, 1])
                .unwrap();

        assert_eq!(uniform_box.nodes.len(), 8);

        for node in &uniform_box.nodes {
            assert_eq!(node.len(), 3);
        }

        assert_eq!(uniform_box.cells.len(), 12);

        for cell in &uniform_box.cells {
            assert_eq!(cell.len(), 3);
        }
    }
}
