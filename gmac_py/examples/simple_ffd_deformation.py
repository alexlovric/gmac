"""
Simple Block Deformation using Free-Form Deformation (FFD)

This example demonstrates how to use the GMAC library to perform Free-Form Deformation (FFD)
on a simple box geometry. FFD is a technique used to deform solid geometry in a smooth way
by manipulating a grid of control points.
"""

import gmac
import gmac.morph
import gmac.io
import numpy as np

# Create a simple geometry or import one
geometry = gmac.generate_box(
    [2.0, 1.0, 1.0],  # Dimensions (length, width, height)
    [0.0, 0.0, 0.0],  # Center coordinates
    [0.0, 0.0, 0.0],  # Rotation angles (degrees)
    [12, 12, 12]      # Number of divisions in each direction
)

# Alternative: geometry = Mesh::from_stl("path_to_stl")? or from_obj

# geometry.write_stl("original.stl", "binary")

# Create a design block (control lattice) for FFD
# The design block defines the control points that will be used to deform the geometry
design_block = gmac.morph.DesignBlock(
    [1.8, 1.2, 1.2],  # Dimensions of the control lattice
    [0.2, 0.0, 0.0],   # Center offset
    [0.0, 0.0, 0.0],   # Rotation angles (degrees)
    [3, 2, 2]          # Number of control points in each direction
)

# Save the initial design block (control points) for visualization
gmac.io.write_vtp(design_block.nodes, "design_nodes.vtp")

# Select which control points will be free to move during deformation
# The second parameter (2) specifies the number of fixed layers of control points at the boundaries
free_design_ids = design_block.select_free_design_nodes(geometry, 2)

# Create a transformation matrix that combines translation, rotation, and scaling
transformation_matrix = gmac.build_transformation_matrix(
    [0.25, 0.0, 0.0],   # Translation vector (x, y, z)
    [125.0, 0.0, 0.0],   # Rotation angles (degrees) around x, y, z axes
    [1.0, 1.25, 1.25]     # Scaling factors in x, y, z directions
)

# Transform only the free control points using the transformation matrix
deformed_design_nodes = np.array(design_block.nodes)  # Create a copy of the original nodes
deformed_design_nodes[free_design_ids] = gmac.transform_nodes(
    deformed_design_nodes[free_design_ids],  # Points to transform
    transformation_matrix,                   # Combined transformation matrix
    [0.2, 0., 0.],                          # Pivot point for transformations
)

# Save the deformed control points for visualisation
gmac.io.write_vtp(deformed_design_nodes, "deformed_design_nodes.vtp")

# Create a Free-Form Deformer with the original design block
ffd = gmac.morph.FreeFormDeformer(design_block)

# Apply the deformation to the original geometry using the deformed control points
geometry.nodes = ffd.deform(geometry.nodes, deformed_design_nodes)

# Save the final deformed geometry as an STL/OBJ/VTK file
geometry.write_stl("deformed_geometry.stl", "ascii")