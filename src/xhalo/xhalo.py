class Halo:
    """
    Class representing the x-ray Halo created by dust scattering.

    Attributes
    """

    def __init__(self, source, dust_model, scatter_model=None):
        """Initialize an x-ray dust halo instance."""
        self.source = source
        self.dusts = dust_model
        self.scatter_model = scatter_model

        self.thetas = None
        self.intensities = None

    def __call__(self, theta):
        """
        Calculate light intensity at a given angular distance from source.
        """
        return self.I(theta)
    
    def I(self, theta):
        """
        Calculate light intensity at a given angular distance from source.
        """