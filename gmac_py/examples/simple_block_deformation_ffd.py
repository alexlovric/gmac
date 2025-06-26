import gmac
import gmac.morph as morph
import gmac.io as io
import numpy as np

geometry = gmac.generate_box([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [5, 5, 5])

design_block = morph.DesignBlock([0.8, 1.2, 1.2], [0.2, 0.0, 0.0], [0.0, 0.0, 0.0], [2, 2, 2])

io.write_vtp(design_block.nodes, "design_block.vtp")

free_design_ids = design_block.select_free_design_nodes(geometry, 2)

transformation_matrix = gmac.build_transformation_matrix([0.25, 0.0, 0.0], [45.0, 0.0, 0.0], [1.0, 1.5, 1.5])

deformed_design_nodes = np.array(design_block.nodes)
deformed_design_nodes[free_design_ids] = gmac.transform_nodes(
    deformed_design_nodes[free_design_ids],
    transformation_matrix,
    [0.2, 0., 0.],
)

io.write_vtp(deformed_design_nodes, "deformed_design_nodes.vtp")

ffd = morph.FreeFormDeformer(design_block)

geometry.nodes = ffd.deform(geometry.nodes, deformed_design_nodes)

io.write_stl(geometry.nodes, geometry.cells, "deformed_geometry.stl")