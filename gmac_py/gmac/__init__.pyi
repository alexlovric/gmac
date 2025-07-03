class Mesh:
    """
    Mesh dataclass.
    """

    def __init__(self, nodes: list[list[float]], cells: list[list[int]]) -> None:
        """
        Instantiate mesh from nodes and cells.

        Parameters
        ----------
        nodes : list[list[float]]
            Node list.
        cells : list[list[int]]
            Cell list.
        """
        self.nodes = nodes
        self.cells = cells
    @classmethod
    def from_stl(cls, filename: str) -> Mesh:
        """
        Instantiate mesh from an stl file.

        Parameters
        ----------
        filename : str
            Path to the stl file.

        Returns
        -------
        Mesh
            The mesh.
        """
    def triangles(self) -> list[list[tuple[float, float, float]]]:
        """
        Get triangles that make up the mesh cells.

        Returns
        -------
        list[list[tuple[float, float, float]]]
            Nodes making up the triangles.
        """
    def cell_normals(self) -> list[float]:
        """
        Get cell normals for the mesh.
        """
    def slice_from_plane(
        self, origin: list[float], normal: list[float]
    ) -> list[list[float]]:
        """
        Find interpolated points for the intersection of the mesh and slicing plane.

        Parameters
        ----------
        origin : list[float]
            Origin of slicing plane.
        normal : list[float]
            Normal of slicing plane.

        Returns
        -------
        list[list[float]]
            New nodes that interpolate onto slice plane.
        """
    def clip_from_plane(self, origin: list[float], normal: list[float]) -> Mesh:
        """
        Get mesh in the direction of a clipping plane.

        Parameters
        ----------
        origin : list[float]
            Origin of clipping plane.
        normal : list[float]
            Normal of clipping plane.

        Returns
        -------
        Mesh
            New clipped mesh.
        """
    def write_stl(self, filename: str | None, format: str = "binary") -> None:
        """
        Write an stl file from the nodes and cells.

        Parameters
        ----------
        filename : str | None
            Optional filename. Defaults to mesh.stl
        format : str
            Format of the stl file. Defaults to "binary".
        """
    def write_vtp(self, filename: str | None) -> None:
        """
        Write an VTP file from the nodes.

        Parameters
        ----------
        filename : str | None
            Optional filename. Defaults to mesh.stl
        """

def select_nodes_closest_to_point(
    nodes: list[list[float]],
    point: list[float],
) -> list[int]:
    """
    Select the indices of nodes closest to a point.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    point : list[float]
        Point target.

    Returns
    -------
    list[int]
        The indices of the nodes that match the condition.
    """

def select_nodes_on_line(
    nodes: list[list[float]],
    point_a: list[float],
    point_b: list[float],
) -> list[int]:
    """
    Select the indices of nodes that pass through a line.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    point_a : list[float]
        Starting point of the line.
    point_b : list[float]
        Ending point of the line.

    Returns
    -------
    list[int]
        The indices of the nodes that match the condition.
    """

def select_nodes_in_sphere(
    nodes: list[list[float]],
    radius: float,
    centre: list[float],
) -> list[int]:
    """
    Select the indices of nodes in a sphere.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    radius : float
        The radius of the sphere.
    centre : list[float]
        The centre of the sphere.

    Returns
    -------
    list[int]
        The indices of the nodes that match the condition.
    """

def select_nodes_in_box(
    nodes: list[list[float]],
    length: list[float],
    centre: list[float],
    theta: list[float],
) -> list[int]:
    """
    Select the indices of nodes in a box.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    length : list[float]
        The dimensions of the box [length, width, height].
    centre : list[float]
        The centre of the box.
    theta : list[float]
        The rotation angles of the box.
        
    Returns
    -------
    list[int]
        The indices of the nodes that match the condition.
    """

def select_nodes_closest_to_plane(
    nodes: list[list[float]],
    origin: list[float],
    normal: list[float],
) -> list[int]:
    """
    Select the indices of nodes closest to a plane.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    origin : list[float]
        An array of x, y, z coordinates representing a point on the plane.
    normal : list[float]
        An array representing the normal vector to the plane.

    Returns
    -------
    list[int]
        The indices of the nodes that match the condition.
    """

def select_nodes_in_plane_direction(
    nodes: list[list[float]],
    origin: list[float],
    normal: list[float],
) -> list[int]:
    """
    Select the indices of nodes in the direction of a plane's normal.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    origin : list[float]
        An array of x, y, z coordinates representing a point on the plane.
    normal : list[float]
        An array representing the normal vector to the plane.

    Returns
    -------
    list[int]
        The indices of the nodes that match the condition.
    """

def translate_nodes(
    nodes: list[list[float]],
    translation: list[float],
) -> list[list[float]]:
    """
    Translates nodes by a given translation.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    translation : list[float]
        A list of floats (x, y, z) that represents the translation.

    Returns
    -------
    list[list[float]]
        The translated nodes.
    """

def rotate_nodes(
    nodes: list[list[float]],
    rotation: list[float],
) -> list[list[float]]:
    """
    Rotates nodes by given rotation angles.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    rotation : list[float]
        A list of floats (x, y, z) that represent the angles of rotation.
    origin : list[float]
        A list of floats (x, y, z) that represent the origin of the rotation.

    Returns
    -------
    list[list[float]]
        The rotated nodes.
    """

def scale_nodes(
    nodes: list[list[float]],
    scaling: list[float],
) -> list[list[float]]:
    """
    Scales nodes by a given scaling.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    scaling : list[float]
        A list of floats (x, y, z) that represent the scaling.
    origin : list[float]
        A list of floats (x, y, z) that represent the origin of the scaling.

    Returns
    -------
    list[list[float]]
        The scaled nodes.
    """

def transform_nodes(
    nodes: list[list[float]],
    transformation_matrix: list[float],
    origin: list[float] | None,
) -> list[list[float]]:
    """
    Transforms nodes using a transformation matrix.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to select from.
    transformation_matrix : list[list[float]]
        A 4x4 transformation matrix.
    origin : list[float]
        A list of floats (x, y, z) that represent the origin of the transformation.

    Returns
    -------
    list[list[float]]
        The transformed nodes.
    """

def build_transformation_matrix(
    translation: list[float],
    rotation: list[float],
    scaling: list[float],
) -> list[list[float]]:
    """
    Builds 4x4 transformation matrix.

    Parameters
    ----------
    translation : list[float]
        A list of floats (x, y, z) that represents the translation.
    rotation : list[float]
        A list of floats (x, y, z) that represent the angles of rotation.
    scaling : list[float]
        A list of floats (x, y, z) that represent the scaling.

    Returns
    -------
    list[list[float]]
        4x4 transformation matrix.
    """

def generate_box(
    length: list[float],
    centre: list[float],
    theta: list[float],
    resolution: list[int],
) -> Mesh:
    """
    Generate 3D box mesh shell.

    Parameters
    ----------
    length : list[float]
        Lengths in x, y and z directions.
    centre : list[float]
        Centre of the block.
    theta : list[float]
        Rotation angles in x, y and z directions.
    resolution : list[int]
        Number of cells in x, y and z directions.

    Returns
    -------
    Mesh
        The mesh containing the nodes and cells.
    """

def generate_naca_wing(
    maximum_camber: float,
    camber_distance: float,
    maximum_thickness: float,
    n_points: int,
    wing_span: tuple[float, float],
) -> Mesh:
    """
    Generate a 3D Naca wing mesh shell.

    Parameters
    ----------
    maximum_camber : float
        The maximum camber of the airfoil, usually expressed as a percentage of the chord length.
    camber_distance : float
        The distance from the leading edge to the location of maximum camber,
        usually as a percentage of chord length.
    maximum_thickness : float
        The maximum thickness of the airfoil, usually as a percentage of chord length.
    n_points : int
        The number of points to generate for each half of the airfoil (upper and lower).
    wing_span : tuple[float, float]
        A tuple representing the z-coordinates for the beginning and end of the wing span.

    Returns
    -------
    Mesh
        The mesh containing the nodes and cells.
    """

def generate_block_cluster(
    length: list[float],
    centre: list[float],
    theta: list[float],
    resolution: list[int],
) -> list[list[float]]:
    """
    Generate 3D node cluster block.

    Parameters
    ----------
    length : list[float]
        Lengths in x, y and z directions.
    centre : list[float]
        Centre of the block.
    theta : list[float]
        Rotation angles in x, y and z directions.
    resolution : list[int]
        Number of cells in x, y and z directions.

    Returns
    -------
    list[list[float]]
        The nodes making up the block.
    """
