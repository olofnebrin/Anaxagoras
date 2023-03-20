""" COSMOLOGY MODULE

    This module contain all the relevant cosmological parameters used,
    as well as derived halo properties (e.g. virial radius). """

import numpy as np

###############################################################
#######                                                 #######
#######   C O S M O L O G I C A L  P A R A M E T E R S  #######
#######                                                 #######
###############################################################

h       = 0.674            # Normalized Hubble constant
omega_B = 0.0493           # Baryon density parameter
omega_M = 0.315            # Matter density parameter
fB      = omega_B/omega_M  # Universal baryon fraction
T_CMB0  = 2.726            # Present CMB temperature (K)
Y       = 0.247            # Helium mass fraction
X       = 1 - Y            # Hydrogen mass fraction

def T_CMB(z):

    """ The CMB temperature (K) at redshift z. """

    return T_CMB0*(1.0 + z)

###############################################################
#######                                                 #######
#######          H A L O  P R O P E R T I E S           #######
#######                                                 #######
###############################################################

def R_vir(M,z):

    """ Virial radius (in cm) of a halo of mass M (in Msun)
        at redshift z. """

    # Hubble constant and mean matter density in the Universe:

    H0    = h*3.2404e-18
    rho_m = omega_M*((3*(H0**2)/(8*np.pi*G)))*((1+z)**3)

    # Virial overdensity:

    delta_vir = 18*(np.pi**2)
    rho_vir   = delta_vir*rho_m

    # Resulting virial radius:

    R_vir = (3*M*Msun/(4*np.pi*rho_vir))**(1/3)

    return R_vir

def v_vir(M,z):

    """ Virial velocity (in cm/s) of a halo of
        mass M (in Msun) at redshift z. """

    v_vir = np.sqrt(G*M*Msun/R_vir(M,z))

    return v_vir

def T_vir(M, z):

    """ Virial temperature (K) of a halo of mass
        M (in Msun) at redshift z. """

    # Virial temperature:

    T_vir = 0.75*mu_h*mH*(v_vir(M,z)**2)/(2*kB)

    return T_vir


def c(M,z):

    """ NFW halo concentration c as a function
        of halo mass M (in Msun) and redshift z.

        The expression is valid for all halo masses M
        and redshifts z > 4, and is taken from:

        C. A. Correa et al. (2015). The accretion history of dark matter
        halos - III. A physical model for the concentration-mass relation,
        MNRAS, vol. 452, 2, pp. 1217-1232. """

    a = 1.3081 - 0.1078*(1+z) + 0.00398*((1+z)**2)
    b = 0.0223 - 0.0944*((1+z)**(-0.3907))

    # Halo concentration:

    c = (10**a)*(M**b)

    return c


def M_NFWencl(R, M, z):

    """ Mass eneclosed (Msun) within a radius R (pc)
        of a DM halo of mass M at redshift z,
        with an NFW density profile. """

    # Useful definitions:

    NFW_conc = c(M,z)  # Halo concentration (call just once)

    f = lambda x: np.log(1+x) - x/(1+x)

    M_encl = M*f( NFW_conc*R*pc/R_vir(M,z) )/f( NFW_conc )

    return M_encl
