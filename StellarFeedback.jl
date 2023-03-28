""" STELLAR FEEDBACK MODULE

    This module contain everything related to stellar feedback
    from Pop II and Pop III stars (mass-to-light ratios, SN rates,
    metal yields etc). """

module StellarFeedback

# The assumed mass of Pop III stars (Msun) & the critical 
# metallicity for Pop II star formation:

mPopIII = 25.0; export mPopIII
Zcrit   = 1e-5; export Zcrit 

###############################################################
#######                                                 #######
#######   N E E D E D  M O D U L E S / P A C K A G E S  #######
#######                                                 #######
###############################################################

include("./CosmologyAndCooling.jl")
using .CosmologyAndCooling

using Dierckx            # For making interpolations of data 
using DelimitedFiles     # For saving & reading data in text files

###############################################################
#######                                                 #######
#######     L I G H T -  T O - M A S S  R A T I O S     #######
#######                                                 #######
###############################################################

# Below the light-to-mass ratios in different bands are estimated
# for both Pop II and Pop III stars. The Pop II band parameters 
# come directly from the FIRE-3 fits (Hopkins et al. 2023), which
# are based on STARBURST99 for a Kroupa IMF. For Pop III stars we
# assume a Dirac delta IMF at some specified mass, and the 
# light-to-mass ratio (in the ionizing band) is calculated 
# accordingly. 

################ POP III LIGHT-TO-MASS RATIO ##################

# To calculate the light-to-mass ratio in the ionizing band
# we need to know the ionizing photon emission rate (N_LyCdot),
# and the average energy of LyC photons (E_LyCIII).

# Data from Mas-Ribas (2016):

mPopIII_data     = [9.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 
80.0, 100.0, 120.0, 200.0, 300.0, 400.0, 500.0, 1000.0]
                
E_LyCIII_data    = [1.66, 1.85, 1.95, 2.04, 2.09, 2.11, 2.15, 2.23, 2.30, 2.35, 2.44, 
2.48, 2.53, 2.67, 2.77, 2.91, 2.95, 2.99]*13.6*eV

# Make interpolation from the above data:

E_LyCIII     = Spline1D(mPopIII_data, E_LyCIII_data)

# Fit to N_LyCdotIII from Schaerer (2002), averaged over the Pop III lifetime:

N_LyCdotIII(m) = 10^( 43.61 + 4.9*log10(m) - 0.83*(log10(m)^2) )

# Resulting light-to-mass ratio (divided by c) in the ionizing band:

Ψ_ion_cIII = N_LyCdotIII(mPopIII)*E_LyCIII(mPopIII)/(mPopIII*Msun*c)

################ POP II LIGHT-TO-MASS RATIOS ##################

function Ψbol(t,zFe)

end

###############################################################
#######                                                 #######
#######              S U P E R N O V A E                #######
#######                                                 #######
###############################################################

# Rates of SNe, metal yields, and SN energies are calculated for
# Pop III stars below using data from Nomoto et al. (2013). For
# Pop III stars we assume that stars with 11-25 Msun explode
# as normal SNe, 25-40 Msun as hypernovae, and 140-300 Msun as 
# pair-instability SNe. For Pop II stars we take the same yields
# and expressions from FIRE-3 (Hopkins et al. 2023) that assume
# a Kroupa IMF.

##################### POP III SUPERNOVAE ######################

# Read data of yields and energies for CCSNe with m = 11-40 Msun
# and create an interpolation for the SN energy:

CCSN_data     = readdlm("Nomoto2013_SNyieldsAndEnergies_CCSN_11to40Msun.txt")
m_CCSN_data   = CCSN_data[2,2:9]
E_CCSN_data   = CCSN_data[3,2:9]; E_SNIII_CCSN = Spline1D(m_CCSN_data, E_CCSN_data)

# Read metal yields from CCSNe with m = 11-40 Msun and create interpolations:

C_CCSN_data   = CCSN_data[14,3:10] + CCSN_data[15,3:10]
C_SNIII_CCSN  = Spline1D(m_CCSN_data, C_CCSN_data)

N_CCSN_data   = CCSN_data[16,3:10] + CCSN_data[17,3:10]
N_SNIII_CCSN  = Spline1D(m_CCSN_data, N_CCSN_data)

O_CCSN_data   = CCSN_data[18,3:10] + CCSN_data[19,3:10] + CCSN_data[20,3:10]
O_SNIII_CCSN  = Spline1D(m_CCSN_data, O_CCSN_data)

Ne_CCSN_data  = CCSN_data[22,3:10] + CCSN_data[23,3:10] + CCSN_data[24,3:10]
Ne_SNIII_CCSN = Spline1D(m_CCSN_data, Ne_CCSN_data)

Mg_CCSN_data  = CCSN_data[26,3:10] + CCSN_data[27,3:10] + CCSN_data[28,3:10]
Mg_SNIII_CCSN = Spline1D(m_CCSN_data, Mg_CCSN_data)

Si_CCSN_data  = CCSN_data[30,3:10] + CCSN_data[31,3:10] + CCSN_data[32,3:10]
Si_SNIII_CCSN = Spline1D(m_CCSN_data, Si_CCSN_data)

S_CCSN_data   = (CCSN_data[34,3:10] + CCSN_data[35,3:10] + CCSN_data[36,3:10]
              + CCSN_data[37,3:10])
S_SNIII_CCSN  = Spline1D(m_CCSN_data, S_CCSN_data)

Ca_CCSN_data  = (CCSN_data[46,3:10] + CCSN_data[47,3:10] + CCSN_data[48,3:10]
              + CCSN_data[49,3:10] + CCSN_data[50,3:10] + CCSN_data[51,3:10])
Ca_SNIII_CCSN = Spline1D(m_CCSN_data, Ca_CCSN_data)

Fe_CCSN_data  = (CCSN_data[65,3:10] + CCSN_data[66,3:10] + CCSN_data[67,3:10]
              + CCSN_data[68,3:10])
Fe_SNIII_CCSN = Spline1D(m_CCSN_data, Fe_CCSN_data)

# Read data of yields and energies for PISNe with m = 140-300 Msun
# and create an interpolation for the SN energy:

PISN_data     = readdlm("Nomoto2013_SNyieldsAndEnergies_PISN_140to300Msun.txt")
m_PISN_data   = PISN_data[2,2:7]
E_PISN_data   = PISN_data[3,2:7]; E_SNIII_PISN = Spline1D(m_PISN_data, E_PISN_data)

# Read metal yields from PISNe with m = 140-300 Msun and create interpolations:

C_PISN_data   = PISN_data[14,3:8] + PISN_data[15,3:8]
C_SNIII_PISN  = Spline1D(m_PISN_data, C_PISN_data)

N_PISN_data   = PISN_data[16,3:8] + PISN_data[17,3:8]
N_SNIII_PISN  = Spline1D(m_PISN_data, N_PISN_data)

O_PISN_data   = PISN_data[18,3:8] + PISN_data[19,3:8] + PISN_data[20,3:8]
O_SNIII_PISN  = Spline1D(m_PISN_data, O_PISN_data)

Ne_PISN_data  = PISN_data[22,3:8] + PISN_data[23,3:8] + PISN_data[24,3:8]
Ne_SNIII_PISN = Spline1D(m_PISN_data, Ne_PISN_data)

Mg_PISN_data  = PISN_data[26,3:8] + PISN_data[27,3:8] + PISN_data[28,3:8]
Mg_SNIII_PISN = Spline1D(m_PISN_data, Mg_PISN_data)

Si_PISN_data  = PISN_data[30,3:8] + PISN_data[31,3:8] + PISN_data[32,3:8]
Si_SNIII_PISN = Spline1D(m_PISN_data, Si_PISN_data)

S_PISN_data   = (PISN_data[34,3:8] + PISN_data[35,3:8] + PISN_data[36,3:8]
              + PISN_data[37,3:8])
S_SNIII_PISN  = Spline1D(m_PISN_data, S_PISN_data)

Ca_PISN_data  = (PISN_data[46,3:8] + PISN_data[47,3:8] + PISN_data[48,3:8]
              + PISN_data[49,3:8] + PISN_data[50,3:8] + PISN_data[51,3:8])
Ca_SNIII_PISN = Spline1D(m_PISN_data, Ca_PISN_data)

Fe_PISN_data  = (PISN_data[65,3:8] + PISN_data[66,3:8] + PISN_data[67,3:8]
              + PISN_data[68,3:8])
Fe_SNIII_PISN = Spline1D(m_PISN_data, Fe_PISN_data)

# Calculate the SN energy and yields for the assumed masses of Pop III stars:

if mPopIII >= minimum(m_CCSN_data) && mPopIII <= maximum(m_CCSN_data)
    E_SNIII    = E_SNIII_CCSN(mPopIII)
    metals_III = [C_SNIII_CCSN(mPopIII); N_SNIII_CCSN(mPopIII); O_SNIII_CCSN(mPopIII); 
                  Ne_SNIII_CCSN(mPopIII); Mg_SNIII_CCSN(mPopIII); Si_SNIII_CCSN(mPopIII); 
                  S_SNIII_CCSN(mPopIII); Ca_SNIII_CCSN(mPopIII); Fe_SNIII_CCSN(mPopIII)]
elseif mPopIII >= minimum(m_PISN_data) && mPopIII <= maximum(m_PISN_data)
    E_SNIII    = E_SNIII_PISN(mPopIII)  
    metals_III = [C_SNIII_PISN(mPopIII); N_SNIII_PISN(mPopIII); O_SNIII_PISN(mPopIII); 
                  Ne_SNIII_PISN(mPopIII); Mg_SNIII_PISN(mPopIII); Si_SNIII_PISN(mPopIII); 
                  S_SNIII_PISN(mPopIII); Ca_SNIII_PISN(mPopIII); Fe_SNIII_PISN(mPopIII)]
else
    E_SNIII    = 0.0
    metals_III = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
end

export E_SNIII
export metals_III

##################### POP II SUPERNOVAE #######################

function R_CCSN_II(t)
    """ The Pop II core-collapse SN rate per Solar mass formed
        per year (see Hopkins et al. 2023). """

    # Core-collapse SNe:

    aS1 = 0.39/1e9; aS2 = 0.51/1e9; aS3 = 0.18/1e9
    tS1 = 3.7e6; tS2 = 7.0e6; tS3 = 44e6
    ψS1 = log(aS2/aS1)/log(tS2/tS1)
    ψS2 = log(aS3/aS2)/log(tS3/tS2)

    if t < tS1 || t > tS3
        R_CC = 0.0
    elseif t >= tS1 && t <= tS2
        R_CC = aS1*((t/tS1)^ψS1)
    elseif t > tS2 && t <= tS3
        R_CC = aS2*((t/tS2)^ψS2)
    end

    return R_CC
end

export R_CCSN_II

function Mej_CCSN_II(t)
    """ The ejecta mass (Msun) from core-collapse Pop II SN at time
        t (yrs) assuming a Kroupa IMF (see Hopkins et al. 2023). """

    if t > 0.0 && t <= 6.5e6
        Mej_CCSN_II = 10.0*((t/6.5e6)^(-2.22))
    elseif t > 6.5e6
        Mej_CCSN_II = 10.0*((t/6.5e6)^(-0.267))
    end

    return Mej_CCSN_II
end

function metals_CCSNII(t)
    """ The metal yields (Msun) at time t (yrs) from
        Pop II core-collapse SNe (see Hopkins et al. 2023). """

    # Time division:

    t_cc1 = 3.7e6; t_cc2 = 8.0e6; t_cc3 = 18e6
    t_cc4 = 30e6; t_cc5 = 44e6

    # a-coefficients for each metal species (C,N,O,Ne,Mg,Si,S,Ca,Fe):

    a_cc1 = [0.237; 1.07e-2; 9.53e-2; 2.60e-2; 2.89e-2; 4.12e-4;
              3.63e-4; 4.28e-5; 5.46e-4]
    a_cc2 = [8.57e-3; 3.48e-3; 0.102; 2.20e-2; 1.25e-2; 7.69e-3;
              5.61e-3; 3.21e-4; 2.18e-3]
    a_cc3 = [1.69e-2; 3.44e-4; 9.85e-2; 1.93e-2; 5.77e-3; 8.73e-3;
              5.49e-3; 6.00e-4; 1.08e-2]
    a_cc4 = [9.33e-3; 3.72e-3; 1.73e-2; 2.70e-3; 1.03e-3; 2.23e-3;
              1.26e-3; 1.84e-4; 4.57e-3]
    a_cc5 = [4.47e-3; 3.50e-3; 8.20e-3; 2.75e-3; 1.03e-3; 1.18e-3;
              5.75e-4; 9.64e-5; 1.83e-3]

    # Fraction of ejecta mass in each species:

    if t <= t_cc1
        y_CCSNII = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] 
    elseif t > t_cc1 && t <= t_cc2
        ψ_cc1    = log.(a_cc2./a_cc1)/log(t_cc2/t_cc1)
        y_CCSNII = a_cc1.*((t/t_cc1).^ψ_cc1)
    elseif t > t_cc2 && t <= t_cc3
        ψ_cc2    = log.(a_cc3./a_cc2)/log(t_cc3/t_cc2)
        y_CCSNII = a_cc2.*((t/t_cc2).^ψ_cc2)
    elseif t > t_cc3 && t <= t_cc4
        ψ_cc3    = log.(a_cc4./a_cc3)/log(t_cc4/t_cc3)
        y_CCSNII = a_cc3.*((t/t_cc3).^ψ_cc3)
    elseif t > t_cc4 && t < t_cc5
        ψ_cc4    = log.(a_cc5./a_cc4)/log(t_cc5/t_cc4)
        y_CCSNII = a_cc4.*((t/t_cc4).^ψ_cc4)
    else
        y_CCSNII = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
    end

    return y_CCSNII*Mej_CCSN_II(t)
end

export metals_CCSNII

end # end module