""" POP II STELLAR FEEDBACK MODULE

    This module contain everything needed to calculate stellar
    feedback from Pop II stars. """

# Values of Psi/c (in cm s^-2):

psi_ion_cII = 3.27e-8               # Direct stellar emission in ionizing band
psi_FUV_cII = 2.40e-8               # Direct stellar emission in FUV band
psi_NUV_cII = 1.31e-8               # Direct stellar emission in NUV band
psi_opt_cII = 6.73e-9               # Direct stellar emission in optical band
psi_Lya_cII = 9.13e-9               # From indirect Ly-alpha emission

# Band opacities (in cm^2 g^-1):

k_ionII     = 1.5e6                 # Opacity in ionizing band (from HI) 
k_FUVII     = lambda Z: 2000.0*Z    # Opacity in FUV band (from dust)
k_NUVII     = lambda Z: 1800.0*Z    # Opacity in NUV band (from dust)
k_optII     = lambda Z: 180.0*Z     # Opacity in optical/near-IR band (from dust)
k_IRII      = lambda Z: 5.0*Z       # Rosseland mean opacity in the infrared (from dust)



