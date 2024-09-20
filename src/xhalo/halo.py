class Halo:
    import scipy.integrate as integrate
    import numpy as np
    import math

    def __init__(self, N_H, E, dusts = None, scatter_model = "GaussRG"):
        self.E = E #AUTUMN: allow for multiple E segments
        self.N_H = N_H
        self.dusts = dusts
        self.scatter_model = scatter_model

        if self.scatter_model == "GaussRG":
            self.dsigma_dOmega = self.gaussRG_dsigma_dOmega
        elif self.scatter_model == "ExactRG":
            print("The exact Rayleigh Gans Approximation has not yet been introduced to this program.")
        elif self.scatter_model == "":
        
    # #AUTUMN : allow other scatter models
    # def dsigma_dOmega(self, a, dust, theta):
    #     c = 1.1*8.4616e-8
    #     exponent = math.exp(-.4575 * self.E**2 * a**2 * (theta/60)**2)
    #     F = dust.F(self.E)
    #     return (c * math.pi**.5 * F * (dust.rho/3)**2 * a**6 * exponent)

    def gaussRG_dsigma_dOmega(self, a, dust, theta):
        # modeled from shalo.sl
        import math
        c = 1.1*8.4616e-8
        F = dust.F(self.E)

        beta = self.E * a * (theta/60) * 0.4575**.5

        return c * a**6 * math.pi**.5 * F * (dust.rho/3)**2 * (math.erfc(beta)/(2*beta))

        
    def I(self, theta):
        I = 0
        for dust in self.dusts:
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
        thetas = np.linspace(0,30,100)
        I_vals = np.empty(thetas.size)
        for index,theta in enumerate(thetas):
            I_vals[index] = self.I(theta)
        plt.plot(thetas, I_vals)