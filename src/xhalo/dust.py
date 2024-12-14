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
    
class WD(Dust):
    """
    Class representing a dust component following the models from 
    Weingartner & Draine, ApJ, 2001, 548, 296
    """

    WD_dict = {
            "alphag" :  [-2.25,-2.17,-2.04,-1.91,-1.84,-1.72,-1.54,-2.26,-2.16,-2.01,
                        -1.83,-1.64,-2.35,-2.12,-1.94,-1.61,-2.62,-2.52,-2.36,-2.09,
                        -1.96,-2.80,-2.67,-2.45,-1.90],
            "betag" :   [-0.0648,-0.0382,-0.111,-0.125,-0.132,-0.322,-0.165,-0.199,-0.0862,
                        -0.0973,-0.175,-0.247,-0.668,-0.67,-0.853,-0.722,-0.0144,-0.0541,
                        -0.0957,-0.193,-0.813,0.0356,0.0129,-0.00132,-0.0517],
            "atg" :     [0.00745,0.00373,0.00828,0.00837,0.00898,0.0254,.0107,0.0241,
                        0.00867,0.00811,0.0117,0.0152,0.148,0.0686,0.0786,0.0418,0.0187,
                        0.0366,0.0305,0.0199,0.0693,0.0203,0.0134,0.0275,0.012],
            "acg" :     [0.606,0.586,0.543,0.499,0.489,0.438,0.428,0.861,0.803,0.696,
                        0.604,0.536,1.96,1.35,0.921,0.72,5.74,6.65,6.44,4.6,3.48,3.43,
                        3.44,5.14,7.28],
            "cg" :      [9.94e-11,3.79e-10,5.57e-11,4.15e-11,2.90e-11,3.20e-12,9.99e-12,
                        5.47e-12,4.58e-11,3.96e-11,1.42e-11,5.83e-12,4.82e-14,3.65e-13,
                        2.57e-13,7.58e-13,6.46e-12,1.08e-12,1.62e-12,4.21e-12,2.95e-13,
                        2.74e-12,7.25e-12,8.79e-13,2.86e-12],
            "alphas" :  [-1.48,-1.46,-1.43,-1.41,-2.1,-2.1,-2.21,-2.03,-2.05,-2.06,-2.08,
                        -2.09,-1.57,-1.57,-1.55,-1.59,-2.01,-2.11,-2.05,-2.1,-2.11,-1.09,
                        -1.14,-1.08,-1.13],
            "betas" :   [-9.34,-10.3,-11.7,-11.5,-0.114,-0.0407,0.3,0.668,0.832,0.995,
                        1.29,1.58,1.1,1.25,1.33,2.12,0.894,1.58,1.19,1.64,2.1,-0.37,
                        -0.195,-0.336,-0.109],
            "ats" :     [0.172,0.174,0.173,0.171,0.169,0.166,0.164,0.189,0.188,0.185,
                        0.184,0.183,0.198,0.197,0.195,0.193,0.198,0.197,0.197,0.198,0.198,
                        0.218,0.216,0.216,0.211],
            "cs" :      [1.02e-12,1.09e-12,1.27e-12,1.33e-12,1.26e-13,1.27e-13,1.e-13,
                        5.2e-14,4.81e-14,4.7e-14,4.26e-14,3.94e-14,4.24e-14,4.e-14,
                        4.05e-14,3.2e-14,4.95e-14,3.69e-14,4.37e-14,3.63e-14,3.13e-14,
                        1.17e-13,1.05e-13,1.17e-13,1.04e-13],
            "bc5" :     [0.,1.,2.,3.,4.,5.,6.,0.,1.,2.,3.,4.,0.,1.,2.,3.,0.,1.,2.,3.,4.,
                        0.,1.,2.,3.]
    }

    def wd_graindist(self, index, dust_type, a):
        index -= 1

        WD_dict = {
            "alphag" :  [-2.25,-2.17,-2.04,-1.91,-1.84,-1.72,-1.54,-2.26,-2.16,-2.01,
                        -1.83,-1.64,-2.35,-2.12,-1.94,-1.61,-2.62,-2.52,-2.36,-2.09,
                        -1.96,-2.80,-2.67,-2.45,-1.90],
            "betag" :   [-0.0648,-0.0382,-0.111,-0.125,-0.132,-0.322,-0.165,-0.199,-0.0862,
                        -0.0973,-0.175,-0.247,-0.668,-0.67,-0.853,-0.722,-0.0144,-0.0541,
                        -0.0957,-0.193,-0.813,0.0356,0.0129,-0.00132,-0.0517],
            "atg" :     [0.00745,0.00373,0.00828,0.00837,0.00898,0.0254,.0107,0.0241,
                        0.00867,0.00811,0.0117,0.0152,0.148,0.0686,0.0786,0.0418,0.0187,
                        0.0366,0.0305,0.0199,0.0693,0.0203,0.0134,0.0275,0.012],
            "acg" :     [0.606,0.586,0.543,0.499,0.489,0.438,0.428,0.861,0.803,0.696,
                        0.604,0.536,1.96,1.35,0.921,0.72,5.74,6.65,6.44,4.6,3.48,3.43,
                        3.44,5.14,7.28],
            "cg" :      [9.94e-11,3.79e-10,5.57e-11,4.15e-11,2.90e-11,3.20e-12,9.99e-12,
                        5.47e-12,4.58e-11,3.96e-11,1.42e-11,5.83e-12,4.82e-14,3.65e-13,
                        2.57e-13,7.58e-13,6.46e-12,1.08e-12,1.62e-12,4.21e-12,2.95e-13,
                        2.74e-12,7.25e-12,8.79e-13,2.86e-12],
            "alphas" :  [-1.48,-1.46,-1.43,-1.41,-2.1,-2.1,-2.21,-2.03,-2.05,-2.06,-2.08,
                        -2.09,-1.57,-1.57,-1.55,-1.59,-2.01,-2.11,-2.05,-2.1,-2.11,-1.09,
                        -1.14,-1.08,-1.13],
            "betas" :   [-9.34,-10.3,-11.7,-11.5,-0.114,-0.0407,0.3,0.668,0.832,0.995,
                        1.29,1.58,1.1,1.25,1.33,2.12,0.894,1.58,1.19,1.64,2.1,-0.37,
                        -0.195,-0.336,-0.109],
            "ats" :     [0.172,0.174,0.173,0.171,0.169,0.166,0.164,0.189,0.188,0.185,
                        0.184,0.183,0.198,0.197,0.195,0.193,0.198,0.197,0.197,0.198,0.198,
                        0.218,0.216,0.216,0.211],
            "cs" :      [1.02e-12,1.09e-12,1.27e-12,1.33e-12,1.26e-13,1.27e-13,1.e-13,
                        5.2e-14,4.81e-14,4.7e-14,4.26e-14,3.94e-14,4.24e-14,4.e-14,
                        4.05e-14,3.2e-14,4.95e-14,3.69e-14,4.37e-14,3.63e-14,3.13e-14,
                        1.17e-13,1.05e-13,1.17e-13,1.04e-13],
            "bc5" :     [0.,1.,2.,3.,4.,5.,6.,0.,1.,2.,3.,4.,0.,1.,2.,3.,0.,1.,2.,3.,4.,
                        0.,1.,2.,3.]
        }

        alphag = WD_dict["alphag"][index]
        betag = WD_dict["betag"][index]
        atg = WD_dict["atg"][index]*1.e-4
        acg = WD_dict["acg"][index]*1.e-4
        cg = WD_dict["cg"][index]
        alphas = WD_dict["alphas"][index]
        betas = WD_dict["betas"][index]
        ats = WD_dict["ats"][index]*1.e-4
        acs = 1e-5
        cs = WD_dict["cs"][index]
        bc5 = WD_dict["bc5"][index]

        if dust_type == "Silicate":
            dnda = (cs/a) * (a/ats)**alphas

            if betas >= 0:
	            dnda = dnda * (1 + betas*(a/ats)) 
            else:
                dnda = dnda / (1 - betas*(a/ats))

            if a > ats:
                dnda = dnda * math.exp(((ats - a)/acs)**3)
        
        elif dust_type == "Graphite":
            dnda = (cg/a) * (a/atg)**alphag

            if betag >= 0:
                dnda = dnda * (1 + betag*(a/atg))
            else:
                dnda = dnda / (1 - betag*(a/atg))

            if a > atg:
                dnda = dnda * math.exp(((atg - a)/acg)**3)

            a01 = 3.5e-8
            a02 = 3e-7
            sig = 0.4
            b1 = 2.0496e-7
            b2 = 9.6005e-11
            dndavsg = (b1/a) * math.exp(-0.5 * (math.log(a/a01)/sig)**2) \
                        + (b2/a) * math.exp(-0.5 * (math.log(a/a02)/sig)**2)
            
            if dndavsg >= 0.0001*dnda:
                dnda = dnda + (bc5 * dndavsg)

        return dnda

class Silicate_WD(WD):
    def __init__(self):
        data = pd.read_csv("dust_data/silicate_f.csv", header=None, names=["E","F"], sep='\s+')
        henke_E = np.array(data["E"])
        henke_F = np.array(data["F"])

        Dust.__init__(
            self,
            name = "WD Silicate",
            a_min = .0003, #um
            a_max = .5, #um
            size_func = self.size_func,
            rho = 3.5, #g/cm^3
            ppm = 36.3,
            atom_weight = 172,
            henke_F = interpolate.interp1d(henke_E, henke_F),
            mie_file = "dust_data/DsdO_3.30.fits")
        
    def size_func(self, a):
        return 1e-4 * self.wd_graindist(7, "Silicate", a*1e-4)
    
class Graphite_WD(WD):
    def __init__(self):
        data = pd.read_csv("dust_data/graphite_f.csv", header=None, names=["E","F"], sep='\s+')
        henke_E = np.array(data["E"])
        henke_F = np.array(data["F"])

        Dust.__init__(
            self,
            name = "WD Graphite",
            a_min = .0003, #um
            a_max = 1.25, #um
            size_func = self.size_func,
            rho = 2.24, #g/cm^3
            ppm = 330 * .7,
            atom_weight = 12,
            henke_F = interpolate.interp1d(henke_E, henke_F),
            mie_file = "dust_data/DsdO_2.20.fits")
        
    def size_func(self, a):
        return 1e-4 * self.wd_graindist(7, "Graphite", a*1e-4)