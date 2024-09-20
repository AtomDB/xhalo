import math
import scipy.integrate as integrate
import scipy.constants as constants
import numpy as np
import pandas as pd

class Dust:
    """
    Class representing a specific dust composition.

    Attributes
    ----------

    """

    def __init__(self, name, a_min, a_max, size_func, rho, ppm, atom_weight):
        self.name = name
        self.a_min = a_min #um
        self.a_max = a_max #um
        self.rho = rho # [g/cm^3]
        self.ppm = ppm * 1e-6 # grains per H atom [grains/H_atoms]
        self.atom_weight = atom_weight # [g/mol]
        self.henke_E = None
        self.henke_F = None

        self.size_func = size_func # function of a
 
        self.norm = (self.ppm * self.atom_weight) / (constants.N_A * self._total_mass()) # unitless  
    
    def __str__(self):
            summary = f"""
            {self.name}
            a_min = {self.a_min} um
            a_max = {self.a_max} um
            norm = {self.norm}
            rho = {self.rho} g/cm^3
            """
            return summary
    
    def n(self,a):
        return self.norm * self.size_func(a)

    def plot_size_func(self):
        from matplotlib import pyplot as plt
        import numpy as np
        a_vals = np.linspace(self.a_min, self.a_max, 100)
        plt.yscale("log")
        plt.plot(a_vals, self.size_func(a_vals))

    def _total_mass(self):
        """
        Return total mass of dust per H_atom [g/H_atom]
        """
        def mass_dist(a):
            return self.size_func(a) * (4/3) * math.pi * (a/1e4)**3 * self.rho
        
        return integrate.quad(mass_dist, self.a_min, self.a_max)[0]
    
class Silicate(Dust):
    def __init__(self):
        Dust.__init__(
            self,
            name = "Silicate",
            a_min = .005, #um
            a_max = .250, #um
            size_func = self.size_func,
            rho = 3.3, #g/cm^3
            ppm = 33,
            atom_weight = 172)

    def size_func(self, a):
        # return grains/H_atom/um for grain size a in um
        return a**-3.5
    
    def F(self,E):
        if self.henke_E is None or self.henke_F is None:
            data = pd.read_csv("dust_data/silicate_f.csv", header=None, names=["E","F"], sep='\s+')
            self.henke_E = np.array(data["E"])
            self.henke_F = np.array(data["F"])
        return np.interp(E, self.henke_E, self.henke_F)
    
class Graphite(Dust):
    def __init__(self):
        Dust.__init__(
            self,
            name = "Graphite",
            a_min = .005, #um
            a_max = .250, #um
            size_func = self.size_func,
            rho = 2.2, #g/cm^3
            ppm = 270,
            atom_weight = 12)

    def size_func(self, a):
        return a**-3.5
    
    def F(self,E):
        if self.henke_E is None or self.henke_F is None:
            data = pd.read_csv("dust_data/graphite_f.csv", header=None, names=["E","F"], sep='\s+')
            self.henke_E = np.array(data["E"])
            self.henke_F = np.array(data["F"])
        return np.interp(E, self.henke_E, self.henke_F)