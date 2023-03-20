""" COOLING THRESHOLD MODULE

    This module contains the calculated cooling threshold. """

import numpy as np

def J21_bg(z):

    """ Background LW background. """

    # My thesis:

    J21_bg_thesis = 10*np.exp(-(z-6)/3.3)

    # Result from Incatasciato+ 2023:

    A = 2.119; B = - 1.117e-1; C = -2.782e-3

    logJ21 = A + B*(1+z) + C*((1+z)**2)

    J21_bg_Incatasciato = 10**logJ21
    
    return J21_bg_Incatasciato

J21_bg = np.vectorize( J21_bg )

###############################################################
#######                                                 #######
#######       C O O L I N G  T H R E S H O L D          #######
#######                                                 #######
###############################################################

def M(vel,z):

    """ Halo mass (Msun) corresponding to a given
        virial velocity vel (km/s) at redshift z. """

    Z = (1+z)/10.0

    M = 1e6*((vel/3.65)**3)*(Z**(-3/2))

    return M

def vel(M,z):

    """ Virial velocity (in km/s) associated with
        a halo of mass M (Msun) at redshift z. """

    Z = (1+z)/10.0; M6 = M/1e6

    vel = 3.65*(M6**(1/3))*(Z**0.5)

    return vel

def M_H2(z,J21):

    """ Critical halo mass for H2-cooling at redshift z in
        a Lyman-Werner background
        J21 = J_LW/10^(-21) erg/s/cm^2/Hz/st. """

    Z = (1+z)/10.0

    J21 = np.fmax(1e-10, J21)

    # Threshold when J21 = 0:

    eta  = 6.0

    top  = np.exp(0.016*(eta**(-0.18))*(Z**0.271))
    bott = np.log(1 + 0.27*(eta**0.859)*(Z**(-1.29)))**0.2703

    M_H2_noLW = (1.45e6/(eta**0.2703))*(top/bott)*(Z**(-1.094))

    # When J21 > 0:

    f_LW = 1 + 3.5*(J21**0.31)*np.exp( - (0.06/J21)**(3/4) )
    ztr  = 6.5*(J21**0.23)
    F    = ((z/ztr)**7)/(1 + (z/ztr)**7)

    # Atomic-cooling threshold:

    M_Lya = 5.11e7*(Z**(-1.5))

    # Put it all together:

    M_H2 = min(M_H2_noLW*f_LW*F + M_Lya*(1 - F), M_Lya)

    return M_H2

M_H2 = np.vectorize( M_H2 )

def M_reion(z, Fion):

    """ Critical halo mass above which reionization feedback
        is ineffective. """

    Z = (1+z)/10.0

    # Convenient unit of flux (from Visbal+ 2017):

    F0 = 6.7e6

    # Threshold due to R-type ionization fronts:

    M_R = 5.5e9*((Fion/F0)**3)*(Z**(-15))

    # Threshold due to D-type ionization fronts:

    stuff = (1.0 + 120.0*((Fion/F0)**0.975)*(Z**(-1.95)))**0.923

    M_D = (1.37e8/stuff)*((Fion/F0)**(3/2))*(Z**(-6))

    # Threshold above which heated gas remains bound:

    M_bound = 1.9e8*(Z**(-3/2))

    # Put it all together:

    M_reion = min( M_bound, max(M_R, M_D) )

    return M_reion

M_reion = np.vectorize( M_reion )


def M_crit(z, sigma_str):

    """ Critical halo mass for cooling, taking baryon
        streaming velocities into account. """

    # Effect of streaming velocities:

    vstr   = 30.0*sigma_str*((1+z)/1e3)
    
    Mcool  = max( M_H2(z, J21_bg(z)), M_reion(z, Fion_bg(z)) )
    vcool  = vel(Mcool, z) 

    alpha  = 6.0
    veff   = np.sqrt( vcool**2 + (alpha*vstr)**2 )

    M_crit = M(veff,z)

    return M_crit

M_crit = np.vectorize( M_crit )



