import math
import scipy.integrate as integrate
import scipy.constants as constants
from scipy import interpolate 
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

class Dust:
    """
    Class representing a dust composition.
    ------------------------------------------
    Attributes:
        name (str):
            The name of this dust composition.
        a_min (float or astropy Quantity):
            The minimum radius of dust grains. 
            Assumed units: [um]
        a_max (float or Quantity):
            The maximum radius of dust grains.
            Assumed units: [um]
        rho (float of Quantity):
            The density of the material within a dust grain.
            Assumed units: [g/cm^2]
        ppm (float):
            Number of dust grains per million hydrogen atoms.
            Assumed Units: ppm = [grains/(10^6 H_atoms)]
            Ex: ppm=33 means there are 33 dust grains per 
                million hydrogen atoms.
        atom_weight (float or quantity):
            The weight of one mole of this material.
            Assumed units: [g/mol]
        norm (float or Quantity):
            The factor used for normalizing the size function.
            Assumed units: Unitless
        file (str)
    """

    def __init__(self, name, size_func, a_min, a_max, \
                 rho, ppm, atom_weight, henke_F, mie_file = None):
        self.name = name
        self.a_min = a_min #um
        self.a_max = a_max #um
        self.rho = rho # [g/cm^3]
        self.ppm = ppm # grains per million H atoms [grains/H_atoms]
        self.atom_weight = atom_weight # [g/mol]
        
        self.size_func = size_func # function of a
        self.henke_F = henke_F # function of E
 
        self.norm = (self.ppm * 1e-6 * self.atom_weight) / (constants.N_A * self._total_mass()) # unitless 
        
        self._mie_file = mie_file
        self._mie_loaded = False
        self._mie_table = None
    
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
        fig,ax = plt.subplots(1, 1, figsize=(8,4))
        ax.set_title("Size Function")
        ax.set_ylabel("Proportion (Unnormalized)")
        ax.set_xlabel("Dust Size [$\mu m$]")
        ax.set_yscale("log")


        a_vals = np.linspace(self.a_min, self.a_max, 100)
        plt.plot(a_vals, self.size_func(a_vals))

    def get_mie_dsigma_dOmega(self, a, theta, E):
        if not self._mie_loaded:
            self._load_mie_table()
        
        # check bounds
        var_dict = {"Energy":E, "Size":a, "ThetaObs":theta}
        for var in var_dict:
            min = self._mie_table[var].min()
            max = self._mie_table[var].max()
            if var_dict[var] < min or var_dict[var] > max:
                print("Mie data only available for {var} values between {min} and {max}.")
                raise(LookupError)

        #get closest vals and create mask
        mask = self._mie_table["Energy"] >= 0
        for var in var_dict:
            closest_index = abs(self._mie_table[var] - var_dict[var]).argmin()
            closest_val = self._mie_table.iloc[closest_index][var]
            mask = mask & (self._mie_table[var] == closest_val)

        intensity = np.array(self._mie_table[mask]["Intensity"])[0] 
        table_norm = (self.ppm * 1e-6 * self.atom_weight)/(constants.N_A * self._mass_grain(a))
        return intensity / table_norm #/ self._mass_grain(a)

    def _load_mie_table(self):
        if self._mie_file is None:
            raise(FileNotFoundError) #AUTUMN: FIX ERRORS
        
        hdu = fits.open(self._mie_file)
        self._mie_table = Table(hdu[1].data).to_pandas()

    def _mass_grain(self, a):
        return (4/3) * math.pi * (a/1e4)**3 * self.rho
    
    def _total_mass(self):
        """
        Return total mass of dust per H_atom [g/H_atom]
        """
        def integrand(a):
            return self.size_func(a) * self._mass_grain(a)
        return integrate.quad(integrand, self.a_min, self.a_max)[0]

class Silicate(Dust):
    def __init__(self):
        data = pd.read_csv("dust_data/silicate_f.csv", header=None, names=["E","F"], sep='\s+')
        henke_E = np.array(data["E"])
        henke_F = np.array(data["F"])

        Dust.__init__(
            self,
            name = "Silicate",
            a_min = .005, #um
            a_max = .250, #um
            size_func = self.size_func,
            rho = 3.3, #g/cm^3
            ppm = 33,
            atom_weight = 172,
            henke_F = interpolate.interp1d(henke_E, henke_F),
            mie_file = "dust_data/DsdO_3.30.fits")

    def size_func(self, a):
        # return grains/H_atom/um for grain size a in um
        return a**-3.5
    
class Graphite(Dust):
    def __init__(self):
        data = pd.read_csv("dust_data/graphite_f.csv", header=None, names=["E","F"], sep='\s+')
        henke_E = np.array(data["E"])
        henke_F = np.array(data["F"])

        Dust.__init__(
            self,
            name = "Graphite",
            a_min = .005, #um
            a_max = .250, #um
            size_func = self.size_func,
            rho = 2.2, #g/cm^3
            ppm = 270,
            atom_weight = 12,
            henke_F = interpolate.interp1d(henke_E, henke_F),
            mie_file = "dust_data/DsdO_2.20.fits")

    def size_func(self, a):
        return a**-3.5