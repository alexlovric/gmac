import gmac
import gmac.morph as morph
import gmac.io as io
import numpy as np

box = gmac.generate_box(
    length=[1.0, 1.0, 1.0],
    centre=[0.0, 0.0, 0.0],
    theta=[0.0, 0.0, 0.0],
    resolution=[5, 5, 5],
)

io.write_stl(nodes=box.nodes, cells=box.cells, filename="original_box.stl")

original_control_points = np.array(
    gmac.generate_block_cluster(
        length=[1.2, 1.2, 1.2],
        centre=[0.0, 0.0, 0.0],
        theta=[0.0, 0.0, 0.0],
        resolution=[2, 2, 2],
    )
)

io.write_vtp(original_control_points, "original_control_points.vtp")

target_control_point_ids = gmac.select_nodes_in_plane_direction(
    nodes=original_control_points, origin=[0.3, 0.0, 0.0], normal=[1.0, 0.0, 0.0]
)

io.write_vtp(
    nodes=original_control_points[target_control_point_ids],
    filename="target_points.vtp",
)

deformed_control_points = original_control_points.copy()
deformed_control_points[target_control_point_ids] = gmac.transform_nodes(
    nodes=deformed_control_points[target_control_point_ids],
    transformation_matrix=gmac.build_transformation_matrix(
        translation=[1.0, 0.0, 0.0],
        rotation=[45.0, 0.0, 0.0],
        scaling=[1.0, 0.75, 0.75],
    ),
    origin=[0.0, 0.0, 0.0],
)

io.write_vtp(nodes=deformed_control_points, filename="deformed_control_points.vtp")

rbf = morph.RbfDeformer(
    original_control_points=original_control_points,
    deformed_control_points=deformed_control_points,
    kernel="thin_plate_spline",
    epsilon=1.0,
)

box.nodes = rbf.deform(points=box.nodes)

io.write_stl(nodes=box.nodes, cells=box.cells, filename="deformed_box.stl")