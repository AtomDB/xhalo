
class Source:
    """
    Class representing a source of x-ray emission.
    """ # AUTUMN: increase specificity

    def __init__(self, N_H=1e22, E=None, f=1, R=None, ): # AUTUMN add defaults
        """
        Initialize an x-ray Source instance.
        """

        self.N_H = N_H
        self.R = R
        self.E = E
        self.f = f

        self.halo = None
    
    def get_halo(self, dust_model=None, scatter_model=None): #AUTUMN add defaults
        """
        Create an x-ray dust halo for this source.
        """

        self.halo = Halo(self, dust_model, scatter_model) # AUTUMN make multiple halos possible
        return self.halo


class Halo:
    """
    Class representing the x-ray Halo caused by dust between observer and source.
    """ # AUTUMN: increase specificity

    def __init__(self, source, dust_model=None, scatter_model=None):
        """Initialize an x-ray dust halo instance."""
        self.source = source
        self.dust_model = dust_model
        self.scatter_model = scatter_model

        self.thetas = None
        self.intensities = None

    def __call__(self, thetas):
        """
        Calculate light intensity at a given distance from source.
        """

class DustModel:
    """
    Class representing an interstellar dust model.
    """ # AUTUMN: increase specificity

    def __init__(self, name=None, compositions = None):
        """
        
        """
        self.name = name
        self.compositions = compositions

class Composition:
    """
    Class representing a specific dust composition.
    """

    def __init__(self, material=None, density=None, n=None):
        """Initialize a dust composition instance."""
        self.material = material
        self.density = density
        self.n = n

# create children of DustModel for popular dust models

# create children of Composition for popular compositions