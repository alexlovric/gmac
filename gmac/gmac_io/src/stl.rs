use std::fs::File;
use std::io::{Write, Result};

pub fn write_stl(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    filename: Option<&str>,
) -> Result<()> {
    let mut file = File::create(filename.unwrap_or("mesh.stl"))?;

    writeln!(file, "solid exported_grid")?;

    for cell in cells.iter() {
        let p1 = nodes[cell[0]];
        let p2 = nodes[cell[1]];
        let p3 = nodes[cell[2]];

        // Calculate normal
        let u = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
        let v = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];

        let n = [
            (u[1] * v[2]) - (u[2] * v[1]),
            (u[2] * v[0]) - (u[0] * v[2]),
            (u[0] * v[1]) - (u[1] * v[0]),
        ];

        let norm = (n[0].powi(2) + n[1].powi(2) + n[2].powi(2)).sqrt();

        writeln!(
            file,
            "  facet normal {} {} {}",
            n[0] / norm,
            n[1] / norm,
            n[2] / norm
        )?;
        writeln!(file, "    outer loop")?;
        writeln!(file, "      vertex {} {} {}", p1[0], p1[1], p1[2])?;
        writeln!(file, "      vertex {} {} {}", p2[0], p2[1], p2[2])?;
        writeln!(file, "      vertex {} {} {}", p3[0], p3[1], p3[2])?;
        writeln!(file, "    endloop")?;
        writeln!(file, "  endfacet")?;
    }

    writeln!(file, "endsolid exported_grid")?;

    Ok(())
}
