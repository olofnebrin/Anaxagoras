""" THE COOLING THRESHOLD 

    This module contain everything related to cosmological parameters
    and derived quantities (e.g. age of Universe), halo properties (e.g. 
    virial radius), and the the cooling threshold Mcrit. """

module CosmologyAndCooling
    
###############################################################
#######                                                 #######
#######        C O N S T A N T S  &  U N I T            #######
#######        C O N V E R S I O N                      #######
#######                                                 #######
###############################################################
    
G       = 6.67e-8          # Gravitational constant (cgs units)
kB      = 1.38e-16         # Boltzmann constant (erg/K)
mH      = 1.67e-24         # Mass of a hydrogen atom (g)
Msun    = 1.989e33         # Solar mass (g)
pc      = 3.08e18          # Parsec (cm)
mpc     = 1e6*pc           # Megaparsec (cm)
kms     = 1e5              # 1 km/s (cm/s)
yr      = 3.1556926e7      # 1 year (s)
    
###############################################################
#######                                                 #######
#######   C O S M O L O G I C A L  P A R A M E T E R S  #######
#######                                                 #######
###############################################################
    
h       = 0.674            # Normalized Hubble constant
H0      = h*100*kms/mpc    # Hubble constant (s^-1)
ΩB      = 0.0493           # Baryon density parameter
ΩM      = 0.315            # Matter density parameter
ΩΛ      = 1.0-ΩM           # Dark energy density parameter (we assume flat Universe)
fB      = ΩB/ΩM            # Universal baryon fraction
T_CMB0  = 2.726            # Present CMB temperature (K)
Y       = 0.247            # Helium mass fraction
X       = 1 - Y            # Hydrogen mass fraction
Δvir    = 18*(pi^2)        # Virial overdensity of a halo   
    
# Derived quantities: The CMB tempoerature (K), the mean matter density (g/cm^3),
# and the age of the Universe (yr):
    
T_CMB(z) = T_CMB0*(1+z)    
ρM(z)    = 3*ΩM*(H0^2)*((1+z)^3)/(8*pi*G)
    
function tuni(z)
    
    """ Function that calculate the age of the Universe (yrs)
        at redshift z. See e.g. Eq. (2.177) in Baumann (2022).
        Only valid after matter-radiation equality. """
    
    a    = 1/(1+z)
    tuni = (2/(3*sqrt(ΩΛ)*H0))*asinh( (a*((ΩΛ/ΩM)^(1/3)))^(3/2) )
    
     return tuni/yr

end

export T_CMB
export ρM

###############################################################
#######                                                 #######
#######          H A L O  P R O P E R T I E S           #######
#######                                                 #######
###############################################################

function Rvir(M,z)
    """ The virial radius (cm) of a halo of mass M (Msun)
        at redshift z.  """

    ρvir = Δvir*ρM(z)
    Rvir = (3*M*Msun/(4*pi*ρvir))^(1/3)

    return Rvir
end

export Rvir 

function vvir(M,z)
    """ The virial velocity (km/s) of a halo of mass M (Msun)
        at redshift z.  """

    vvir = sqrt(G*M*Msun/Rvir(M,z))

    return vvir/kms
end

function Mvel(vvir,z)
    """ The halo mass (Msun) as a function of the
        virial velocity (km/s) and redshift. Needed to 
        compute Mcrit with streaming velocities below. """

    
    ρvir = Δvir*ρM(z)
    M = sqrt(3/(4*pi))*((vvir*kms)^3)/((G^(3/2))*sqrt(ρvir))

    return M/Msun
end

function c_NFW(M,z)

    """ Navarro-Frenk-White (NFW) halo concentration as a function
        of halo mass (Msun) and redshift z. The expression is valid 
        for all halo masses M and redshifts z > 4, and is taken from:
    
        C. A. Correa et al. (2015). The accretion history of dark matter
        halos - III. A physical model for the concentration-mass relation,
        MNRAS, vol. 452, 2, pp. 1217-1232."""

    a = 1.3081 - 0.1078*(1+z) + 0.00398*((1+z)^2)
    b = 0.0223 - 0.0944*((1+z)^(-0.3907))
    
    # Halo concentration:
    
    c_NFW = (10^a)*(M^b)
    return c_NFW
end


###############################################################
#######                                                 #######
#######       C O O L I N G  T H R E S H O L D          #######
#######                                                 #######
###############################################################

function Mcrit_H2LW(z, JLW21)
    """ This function calculate the cooling threshold (Msun) only
        considering H2 and Lyα-cooling and Lyman-Werner feedback.
        JLW21 is the average normalized Lyman-Werner intensity. """

    Z       = (1+z)/10
    J       = max(1e-10, JLW21)
    η       = 6.0

    # Threshold in absence of LW feedback:

    Top     = exp(0.016*(η^(-0.18))*(Z^0.271))
    Bottom  = log(1 + 0.27*(η^0.859)*(Z^(-1.29)))^0.2703

    McritH2 = (1.45e6/(η^0.2703))*(Top/Bottom)*(Z^(-1.094))

    # Threshold with LW feedback:

    f_LW    = 1.0 + 3.5*(J^0.31)*exp( - (0.06/J)^(3/4))
    ztr     = 6.5*(J^0.23)
    F       = ((z/ztr)^7)/(1 + (z/ztr)^7)
    M1e4    = 5.11e7*(Z^(-3/2))

    Mcrit_H2LW = min(M1e4, McritH2*f_LW*F + M1e4*(1-F))

    return Mcrit_H2LW
end

function M_reion(z, Fion)
    """ Critical halo mass (Msun) above which reionization feedback
        is ineffective. Fion is the ionizing photon flux (cm^-2 s^-1). """

    Z     = (1+z)/10
    F     = Fion/6.7e6

    # Threshold due to R-type IF:

    M_R   = 5.5e9*(F^3)*(Z^(-15))

    # Threshold due to D-type IF:

    stuff = (1 + 120*(F^0.975)*(Z^(-1.95)))^0.923
    M_D   = (1.37e8/stuff)*(F^(3/2))*(Z^(-6))

    # Threshold above which heated gas remains bound:

    M_bound = 1.9e8*(Z^(-3/2))

    # Put it all together:

    M_reion = min( M_bound, max(M_R, M_D) )

    return M_reion
end

function Mcrit(z, JLW21, Fion, vstr)
    """ The full expression for the cooling threshold (Msun),
        taking into account both radiative feedback and streaming
        velocities. The streaming velocity (vstr) here is in units
        of the RMS streaming velocity. """

    # The streaming velocity after recombination/decoupling (km/s):

    VStr = vstr*30.0*(1+z)/1000.0

    # Threshold in absence of streaming:

    Mcrit_nostr = max( Mcrit_H2LW(z, JLW21), M_reion(z, Fion) )
    vcrit_nostr = vvir(Mcrit_nostr,z)

    # Threshold with streaming:

    α           = 6.0
    veff        = sqrt(vcrit_nostr^2 + (α*VStr)^2)
    Mcrit       = Mvel(veff,z)

    return Mcrit
end

export Mcrit 

end # End module