""" STARBURST MODULE

    This module calculate the stellar mass and half-mass radius
    following a single starburst in a halo of mass M at redshift
    z, for a given spin parameter. """

import numpy as np

###############################################################
#######                                                 #######
#######          D I S K   P R O P E R T I E S          #######
#######                                                 #######
###############################################################

def Mdisk(M, z, t):

    """ Mass of disk (Msun) a time t (yrs) after
        collapse. Utilizes the self-similar solution
        for the gas accretion rate. """

    # Calculation of disk accretion rate, using results
    # of self-similar model:

    M_div = 2.92e6*(((1+z)/10)**(3/2))

    if M > M_div:

        # Temperature and sound speed of gas in the halo:

        T_h  = min( T_vir(M,z), 8000.0 )
        c_sh = np.sqrt(kB*T_h/(mu_h*mH))

        # Calculation of g:

        beta = (1-fB)*((v_vir(M,z)/c_sh)**2) - 2.0
        y    = np.log10(beta+1.693)**0.52

        g = 0.347*(1.0 + 2.01*np.exp(-6*y))*np.sqrt(0.53 + beta/(beta+2))

    else:

        # Calculation of zeta (requires virial temperature and IGM temperature):

        a1    = 1/119.0
        T_IGM = T_CMB0*a1*((1+z)**2)
        zeta  = (T_vir(M,z)/T_IGM)**(3/2)

        # Calculation of g:

        x = np.log10(zeta/32)**3

        g = 2.00e-5*np.exp(0.00281*x)*zeta*np.sqrt(0.563 + 3.14e-5*zeta)

    # Resulting gas accretion rate (g s^-1):

    Mdot  = g*(v_vir(M,z)**3)/G

    # Disk mass (in Msun):

    Mdisk = Mdot*t*yr/Msun 

    return Mdisk

def f_disk(spin):

    """ The spin parameter-dependent f_disk factor. """

    A = (fB**2)/(16*(spin**2))
    B = 32*(spin**2)/(fB**2)

    f_disk = A*(np.sqrt(1 + B) - 1)

    return f_disk 

def Rdisk(M, z, spin, t):

    """ Disk radius (cm) for a halo of mass M (Msun) at
        redshift z with spin parameter spin,
        a time t (yrs) following collapse. """

    # Fraction of baryons in the disk:

    F = Mdisk(M,z,t)/(fB*M)

    # Disk radius:

    Rdisk = 8*F*f_disk(spin)*(spin**2)*R_vir(M,z)/fB

    return Rdisk

def sigma0(M, z, spin, t):

    """ The disk surface density (g cm^-2)
        at the outer edge of the disk
        (i.e. R = R_disk) at time t (yrs). """

    sigma0 = Mdisk(M,z,t)*Msun/(2*np.pi*(Rdisk(M,z,spin,t)**2))

    return sigma0

def T(z,Z):

    """ Temperature of the disk gas (K) as a function of
        redshift z and gas metallicity Z. """

    # Gas temperature:

    if Z < 1e-4 and H2 == False:

        # In this case Ly-alpha cooling yields T = 8000 K:

        T = 8000.0

    elif Z < 1e-4 and H2 == True:

        # In very metal-poor gas, H2 cools the gas to 200 K:

        T = 200.0

    else:

        # If a sufficient amount of metals is present, gas
        # cools to either the CMB temperature or 10 K:

        T = max( 10, T_CMB(z) )

    return T

def c_s(z,Z):

    """ Sound speed in the disk (cm s^-1)
        as a function of redshift z. """

    c_s = np.sqrt(kB*T(z)/(mu*mH))

    return c_s

def Mach(M, z, spin, Z):

    """ Mach number of disk gas as a function
        of halo mass (M), redshift z, and gas metallicty Z. """

    # Note that the Mach number is independent of time since
    # the time-dependence of sigma0 ~ 1/t and t cancel. For
    # numerical reasons (avoiding problem at t = 0) we just set
    # t = 1.0 (i.e. 1 year):

    x    = 2*e_t*(v_vir(M,z)**2)/(np.pi*G*sigma0(M, z, spin, 1.0)*c_s(z,Z)*yr)

    # Resulting Mach number:

    Mach = ( x**(2/3) + (x/3)**2 )**(1/2)

    return Mach

##############################################################
#######                                                 #######
#######        S T E L L A R   F E E D B A C K          #######
#######                                                 #######
###############################################################


def Transmission(k_band, sigma0):

    """ The factor exp(-tau_band) for a given band
        as a function of band opacity (k_band) and disk
        surface density (sigma0). """

    Transmission = np.exp( - k_band*sigma0 - 0.883*((k_band*sigma0)**0.48))

    return Transmission

def Pdot_m_rad_noion(sigma0, Z):

    """ Momentum injection rate per Solar mass formed
        (in cm s^-2) from direct non-ionizing radiation
        from stars. """

    if Z < Z_crit:
        FUV    =  0.0
        NUV    =  0.0
        OPT    =  0.0
    else:
        FUV    =  psi_FUV_cII*(1 - Transmission(k_FUVII(Z), sigma0))
        UV     =  psi_NUV_cII*(1 - Transmission(k_NUVII(Z), sigma0))
        OPT    =  psi_opt_cII*(1 - Transmission(k_optII(Z), sigma0))

    # PUT IT ALL TOGETHER:

    Pdot_m_rad_noion = FUV + NUV + OPT

    return Pdot_m_rad_noion

def Pdot_m_rest(sigma0, z, Z):

    """ Momentum injection rate per Solar mass formed
        (in cm s^-2) from ionizing radiation, Ly-alpha
        scattering, and stellar winds. """

    # CONTRIBUTION FROM IONIZING BAND:

    if Z < Z_crit:
        ION      = psi_ion_cIII*(1 - Transmission(k_ionIII, sigma0))
    else:
        ION = psi_ion_cII*(1 - Transmission(k_ionII, sigma0))

    # CONTRIBUTION FROM STELLAR WIND/MASS LOSS:

    if Z > 1e-6:

        WINDS = 1.6e-8*np.exp(8e-4*(np.log10(Z/1e-6)**4.3))*(0.01 + Z)

    else:

        # For lower metallicities the fit above yields complex values,
        # and in any case is negligible. So in this case we ignore winds:

        WINDS = 0.0

    # CONTRIBUTION FROM LYMAN-ALPHA SCATTERING:

    sigma0_Msun = sigma0*(pc**2)/Msun    # Surface density in Msun/(pc^2)
    T4          = T(z,Z)/1e4             # Normalized gas temperature
    a_V         = 4.7e-4*(T4**(-1/2))    # Voigt parameter of Lyman-alpha line
    f_weight    = 3/2                    # Mass-weighting factor for the disk

    # Force multiplier without dust:
    
    M_F_nodust   = lambda tau0: 3.51*((a_V*tau0)**(1/3))
    tau0         = 5.33e6*(T4**(-1/2))*sigma0_Msun

    # Force multiplier in general case:

    if Z > 0.0:

        # Maximum force multiplier in the presence of dust:

        tau0_star =  3.58e6*(T4**(-1/4))*(Z**(-3/4))
        
        M_F       =  min(f_weight*M_F_nodust(tau0), M_F_nodust(tau0_star))

    else:

        M_F       =  f_weight*M_F_nodust(tau0)

    # Resulting contribution from Lyman-alpha scattering:

    if Z < Z_crit:
        LY_ALPHA = M_F*psi_Lya_cIII*(1 - Transmission(k_ionIII, sigma0))
    else:
        LY_ALPHA = M_F*psi_Lya_cII*(1 - Transmission(k_ionII, sigma0))

    # PUT IT ALL TOGETHER:

    Pdot_m_rest = ION + WINDS + LY_ALPHA

    return Pdot_m_rest

Pdot_m_rest = np.vectorize( Pdot_m_rest )
