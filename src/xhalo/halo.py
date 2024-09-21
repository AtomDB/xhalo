import scipy.integrate as integrate
import numpy as np
import math
class Halo:
    """
    Class representing the halo created by dust scattering xray light. 
    ------------------------------------------
    Attributes:
    E (float or Quantity):
        The energy of the xray light.
        Assumed units: [keV]
    N_H (float or Quantity):
        The column density of hydrogen atoms. 
        Assumed units: [cm^-2]
    dust_model (list or numpy.ndarray):
        The list of dust components causing scattering. 
        Each item within the list is of type Dust.
    scatter_model (str):
        Name of the approximation used for modeling scattering.
        Choices:
            - "GaussRG"
            - "InfCylRG"
            - "ExactRG"
            - "Mie"
    x (float):
        Location of dust sheet(s) as a value between 0 and 1.
        Ex: x=.3 signifies all dust is in a dust sheet 30%
            of the way to the source.
    """
    def __init__(self, N_H, E, dust_model, scatter_model = "GaussRG", x = None):
        self.E = E #AUTUMN: allow for multiple E segments
        self.N_H = N_H
        self.dust_model = dust_model
        self.scatter_model = scatter_model
        
        if x is not None:
            print("WARNING: The implementation of nonuniform dust \
                  distributions is not yet complete.")

        if self.scatter_model == "GaussRG" and x is None:
            self.dsigma_dOmega = self._gaussRG_dsigma_dOmega
        elif self.scatter_model == "GaussRG" and x is not None:
            self.dsigma_dOmega = self.gaussRGx_dsigma_dOmega
        elif self.scatter_model == "ExactRG":
            print("The exact Rayleigh Gans Approximation has not yet been introduced to this program.")
        elif self.scatter_model == "":
            pass
        
    # #AUTUMN : allow other scatter models
    # def gaussRG_dsigma_dOmega(self, a, dust, theta):
    #     c = 1.1*8.4616e-8
    #     exponent = math.exp(-.4575 * self.E**2 * a**2 * (theta/60)**2)
    #     F = dust.henke_F(self.E)
    #     return (c * math.pi**.5 * F * (dust.rho/3)**2 * a**6 * exponent)

    def _gaussRG_dsigma_dOmega(self, a, dust, theta):
        # modeled from shalo.sl
        import math
        c = 1.1*8.4616e-8
        F = dust.henke_F(self.E)

        beta = self.E * a * (theta/60) * 0.4575**.5

        return c * a**6 * math.pi**.5 * F * (dust.rho/3)**2 * (math.erfc(beta)/(2*beta))

        
    def I(self, theta):
        I = 0
        for dust in self.dust_model:
            I += self.dust_I(dust, theta)
        return I
    
    def dust_I(self, dust, theta):
        import scipy.integrate as integrate
        def integrand(a):
            return dust.n(a) * self.dsigma_dOmega(a, dust, theta)
        integration = integrate.quad(integrand, dust.a_min, dust.a_max)
        return self.N_H * integration[0]
    
    def plot_I(self):
        from matplotlib import pyplot as plt
        thetas = np.linspace(5,30,100)
        I_vals = np.empty(thetas.size)
        for index,theta in enumerate(thetas):
            I_vals[index] = self.I(theta)
        plt.plot(thetas, I_vals)