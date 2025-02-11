import numpy as np
from neqsim.thermo import TPflash, fluid

### START TEST VISCOSITY ###
fluid1 = fluid("srk")
fluid1.addComponent("helium", 1.0)
fluid1.setTemperature(293.15, "K")
fluid1.setPressure(30.175, "MPa")
TPflash(fluid1)
fluid1.getPhase(0).getPhysicalProperties().setViscosityModel("LBC")
fluid1.initProperties()
density = fluid1.getPhase(0).getDensity_GERG2008() / (fluid1.getPhase(0).getMolarMass() *1000)
viscosity = fluid1.getPhase(0).getViscosity("Pas")*10**(6)
print("viscosity LBC: ", viscosity, " [µPa*s]")
### END TEST VISCOSITY ###

def visc_1972(fluid):
    T = fluid.getTemperature("K")
    dens = fluid1.getPhase(0).getDensity_GERG2008() * 0.001     #[g/cm^3]
    x = np.log(T)

    B = -47.5295259/x + 87.6799309 - 42.0741589*x + 8.33128289*x**2 - 0.589252385*x**3
    C = 547.309267/x - 904.870586 + 431.404928*x - 81.4504854*x**2 + 5.37008433*x**3
    D = -1684.39324/x + 3331.08630 - 1632.19172*x + 308.804413*x**2 - 20.2936367*x**3


    eta0_m = -0.135311743/x + 1.00347841 + 1.20654649*x - 0.149564551*x**2 + 0.0125208416*x**3
    etaE_m = dens*B + dens**2*C + dens**3*D

    eta0 = 196*T**(0.71938) * np.exp(12.451/T - 295.67/T**2 - 4.1249)
    etae = np.exp(eta0_m + etaE_m) - np.exp(eta0_m + 0)

    if T > 3.5 and T <= 100:
        eta = (np.exp(eta0_m + etaE_m)) * 0.1       #[µPa*s]
    if T > 100:
        eta = (eta0 + etae) * 0.1       #[µPa*s]
    '''
    if T > 300:
        T = 300
        x = np.log(T)

        B = -47.5295259/x + 87.6799309 - 42.0741589*x + 8.33128289*x**2 - 0.589252385*x**3
        C = 547.309267/x - 904.870586 + 431.404928*x - 81.4504854*x**2 + 5.37008433*x**2
        D = -1684.39324/x + 3331.08630 - 1632.19172*x + 308.804413*x**2 - 20.2936367*x**3


        eta0_m = -0.135311743/x + 1.00347841 + 1.20654649*x - 0.149564551*x**2 + 0.0125208416*x**3
        etaE_m = dens*B + dens**2*C + dens**3*D

        eta0 = 196*T**(0.71938) * np.exp(12.451/T - 295.67/T**2 - 4.1249)
        etae = np.exp(eta0_m + etaE_m) - np.exp(eta0_m + 0)

        eta = (eta0 + etae) * 0.1       #[µPa*s]
    '''
    return eta

#from https://www.sciencedirect.com/science/article/pii/S0149197024004670#bib28:
def visc_liu(fluid):
    T = fluid.getTemperature("K")
    eta = (3.817*10**(-7) * T**(0.6938)) * 10**6    #[µPa*s]
    return eta

print("viscosity 1972 model:", visc_1972(fluid1), " [µPa*s]")
print("viscosity Liu model:", visc_liu(fluid1), " [µPa*s]")