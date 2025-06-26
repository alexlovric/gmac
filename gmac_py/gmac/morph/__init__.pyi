class RbfDeformer:
    """
    Radial basis function mesh deformer.
    """

    def __init__(
        self,
        original_control_points: list[list[float]],
        deformed_control_points: list[list[float]],
        kernel: str | None,
        epsilon: float | None,
    ) -> None:
        """Instantiate a radial basis function deformer.

        Parameters
        ----------
        original_control_points : list[list[float]]
            A n*3 array containing the original node positions.
        deformed_control_points : list[list[float]]
            A n*3 array containing the deformed node positions.
        kernel : str | None
            An optional kernel function name.
        epsilon : float | None
            An optional bandwidth parameter for the kernel.
            Defaults to 1. if `None` given.
        """

    def deform(self, points: list[list[float]]) -> list[list[float]]:
        """Predict the output node positions given input nodes `points`.

        Parameters
        ----------
        points : list[list[float]]
            New input nodes to predict the deformed positions for.

        Returns
        -------
        list[list[float]]
            Deformed node positions.
        """

class FreeFormDeformer:
    """
    Free form deformation mesh deformer.
    """

    def __init__(
        self,
        design_block: DesignBlock,
    ) -> None:
        """Instantiate a free form deformer.

        Parameters
        ----------
        design_block : DesignBlock
            The design block to use for the deformation."""

    def deform(
        self,
        points: list[list[float]],
        deformed_design_nodes: list[list[float]],
    ) -> list[list[float]]:
        """Predict the output node positions given input nodes `points`.

        Parameters
        ----------
        points : list[list[float]]
            New input nodes to predict the deformed positions for.
        deformed_design_nodes : list[list[float]]
            Deformed design nodes.

        Returns
        -------
        list[list[float]]
            Deformed node positions.
        """

class DesignBlock:
    """
    Design block in 3D space for free form deformation.
    """

    def __init__(
        self,
        length: list[float],
        centre: list[float],
        theta: list[float],
        resolution: list[int],
    ) -> None:
        """Instantiate a design block.

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
        """

    def nodes(self) -> list[list[float]]:
        """Get the nodes of the design block."""

    def length(self) -> list[float]:
        """Get the length of the design block."""

    def centre(self) -> list[float]:
        """Get the centre of the design block."""

    def theta(self) -> list[float]:
        """Get the rotation angles of the design block."""

    def resolution(self) -> list[int]:
        """Get the resolution of the design block."""

    def select_free_design_nodes(
        self,
        target_mesh: Mesh,
        fixed_layers: int | None,
    ) -> list[int]:
        """Select the free design nodes.

        Parameters
        ----------
        target_mesh : Mesh
            Mesh to check intersections with.
        fixed_layers : int | None
            How many layers to exclude from the intersection, for instance quadratic=2, linear=1.

        Returns
        -------
        list[int]
            The indices of the free design nodes.
        """
