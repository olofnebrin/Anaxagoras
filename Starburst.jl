""" STARBURST MODULE

    This module contain everything related to calculating the stellar mass
    formed in a single starburst. It also contains the needed calculations
    of the half-mass radius and more. """

include("./CosmologyAndCooling.jl")
using .CosmologyAndCooling

###############################################################
#######                                                 #######
#######          D I S K   P R O P E R T I E S          #######
#######                                                 #######
###############################################################

μh = 1.23
μ  = 1.23

function Mdisk(M, z, t)
    """ Mass of the disk (Msun) at a time t (yrs) after
        gas collapse. Utalizes the self-similar solution
        for the gas accretion rate. """

    Mdiv = 2.92e6*(((1+z)/10)^(-3/2))

    if M > Mdiv
        # Temperature and sound speed of gas in the halo:

        Th  = min( Tvir(M,z),  8000.0 )
        csh = sqrt(kB*Th/(μh*mH))

        # Calculation of g
        
        β    = (1-fB)*((vvir(M,z)*kms/csh)^2) - 2.0
        y    = log10(β+1.693)^0.52

        g    = 0.347*(1.0 + 2.01*exp(-6*y))*sqrt(0.53 + β/(β+2))

    else
        # Calculation of ζ (requires virial & IGM temperature):

        a1    = 1/119.0
        T_IGM = T_CMB0*a1*((1+z)^2)
        ζ     = (Tvir(M,z)/T_IGM)^(3/2)

        # Calculation of g:

        x = log10(ζ/32)^3
        g = 2.00e-5*exp(0.00281*x)*ζ*sqrt(0.563 + 3.14e-5*ζ)

    end

    # Resulting gas accretion rate (g s^-1):

    Mdot  = g*((vvir(M,z)*kms)^3)/G

    # Disk mass (in Msun):

    Mdisk = Mdot*t*yr/Msun 

    return Mdisk
end

function Rdisk(M,z,λ,t)
    """ Disk radius (cm) for a halo of mass M (Msun) at
        redshift z with spin parameter spin,
        a time t (yrs) following collapse. """

    # The λ-dependent factor f_disk:

    A = (fB^2)/(16*(λ^2))
    B = 32*(λ^2)/(fB^2)

    f_disk = A*(sqrt(1+B) - 1)

    # Fraction of baryons in the disk:

    F = Mdisk(M,z,t)/(fB*M)

    # Disk radius:

    Rdisk = 8*F*f_disk*(λ^2)*Rvir(M,z)/fB

    return Rdisk
end

function vdisk(M,z,λ,t)
    """ The disk velocity (cm/s). Note that it is independent
        of time because of cancellation of time-dependent factors.  """

    # To ensure that Rdisk > 0 in the division:

    t = max( 1.0*yr, t )

    # The disk velocity:

    vdisk = sqrt(G*Mdisk(M,z,t)*Msun/Rdisk(M,z,λ,t))

    return vdisk
end

function Σdisk0(M, z, λ, t)
    """ The disk surface density (g/cm^2) at the outer
        edge of the disk (R = Rdisk) at time t (yrs). """

    Σdisk0 = Mdisk(M,z,t)*Msun/(2*pi*(Rdisk(M,z,λ,t)^2))

    return Σdisk0
end


function T(z,Z)
    """ The temperature of dense gas in the disk (K) 
        as a function of redshift and gas metallicity Z. """

    T = 8000.0

    return T
end

function cs(z,Z)
    """ The sound speed in the disk (cm/s) as
        a function of redshift z and gas metallicity
        Z (relative to Solar). """

    return sqrt(kB*T(z,Z)/(μ*mH))
end

function Mach(M,z,λ,Z,t)
    """ Mach number of disk gas as a function
        of halo mass M (Msun), redshift z, 
        gas metallicty (relative to Solar),
        and time t (yrs). Note that the time-dependence
        cancels."""

    e_t = 0.1

    t = max(1.0*yr, t)

    x = 2*e_t*((vvir(M,z)*kms)^2)/(pi*G*Σdisk0(M,z,λ,t)*cs(z,Z)*t*yr)

    Mach = sqrt( x^(2/3) + (x/3)^3 )

    return Mach
end

function Q_Toomre(M,z,λ,Z,t)
    """ The Toomre parameter for the disk. """

    # Epicycle frequency:

    κ    = sqrt(2)*vdisk(M,z,λ,t)/Rdisk(M,z,λ,t)

    # Gas velocity dispersion:

    σgas = cs(z,Z)*sqrt(1 + (Mach(M,z,λ,Z,t)^2)/3)

    # Toomre parameter:

    Q_Toomre = κ*σgas/(pi*G*Σdisk0(M,z,λ,t))

    return Q_Toomre
end

###############################################################
#######                                                 #######
#######        S T E L L A R   F E E D B A C K          #######
#######                                                 #######
###############################################################

Z_crit = 1e-5

###### POP II STELLAR FEEDBACK ######

# For Pop II stars we assume the Kroupa IMF (for 0.1-100 Msun)
# and use the band parameters from FIRE-2:

Ψ_ion_cII  = 3.27e-8     # Direct stellar emission in ionizing band
k_ionII    = 1.5e6       # Hydrogen opacity in ionizing band (for Pop II SED)

Ψ_FUV_cII  = 2.40e-8     # Direct stellar emission in FUV band
k_FUV(Z)   = 2.0e3*Z     # Dust opacity in FUV band

Ψ_NUV_cII  = 1.31e-8     # Direct stellar emission in NUV band
k_NUV(Z)   = 1.8e3*Z     # Dust opacity in NUV band

Ψ_opt_cII  = 6.73e-9     # Direct stellar emission in optical band
k_opt(Z)   = 180.0*Z     # Dust opacity in optical band

k_IR(Z)    = 5.0*Z       # Rosseland mean opacity in the IR band

# Ionizing photon output and Lyman-α band:

Qdot_ionII = 5e46/Msun
E_Lyα      = 10.2*eV
Ψ_Lyα_cII  = (2/3)*Qdot_ionII*E_Lyα/c

###### POP III STELLAR FEEDBACK ######

# For Pop III stars the dust abundance is so low that we 
# can neglect radiation pressure in non-ionizing bands.

Ψ_ion_cIII = 2.10e-7     # Direct stellar emission in ionizing band
k_ionIII   = 7.8e5       # Hydrogen opacity in ionizing band (for Pop III SED)

# Ionizing photon output and Lyman-α band:

Qdot_ionIII = 5e46/Msun
Ψ_Lyα_cIII  = 5.02e-8
 
function Transmission(k_band, Σ0)
    """ The factor exp(-τ_band) for a given band
        as a function of band opacity (k_band) and disk
        surface density (Σ0). """

    Transmission = exp( - k_band*Σ0 - 0.883*((k_band*Σ0)^0.48))

    return Transmission
end

function Pdot_m_rad_noion(Σ0, Z)
    """ Momentum injection rate per Solar mass formed
        (in cm s^-2) from direct non-ionizing radiation
        from stars. """

    if Z < Z_crit
        Pdot_m_rad_noion = 0.0
    else
        FUV    =  Ψ_FUV_cII*(1 - Transmission(k_FUV(Z), Σ0))
        NUV    =  Ψ_NUV_cII*(1 - Transmission(k_NUV(Z), Σ0))
        OPT    =  Ψ_opt_cII*(1 - Transmission(k_opt(Z), Σ0))

        Pdot_m_rad_noion = FUV + NUV + OPT
    end

    return Pdot_m_rad_noion
end

function Pdot_m_rest(Σ0, z, Z)
    """ Momentum injection rate per Solar mass formed
        (in cm s^-2) from ionizing radiation, Ly-alpha
        scattering, and stellar winds. """

    # CONTRIBUTION FROM IONIZING BAND:

    if Z < Z_crit
        ION = Ψ_ion_cIII*(1 - Transmission(k_ionIII, Σ0))
    else
        ION = Ψ_ion_cII*(1 - Transmission(k_ionII, Σ0))
    end

    # CONTRIBUTION FROM STELLAR WIND/MASS LOSS:

    if Z > max(1e-6, Z_crit)
        WINDS = 1.6e-8*exp(8e-4*(log10(Z/1e-6)^4.3))*(0.01 + Z)
    else
        # For lower metallicities the fit above yields complex values,
        # and in any case is negligible. So in this case we ignore winds:

        WINDS = 0.0
    end

    # CONTRIBUTION FROM LYMAN-ALPHA SCATTERING:

    Σ0_Msun     = Σ0*(pc^2)/Msun         # Surface density in Msun/(pc^2)
    T4          = T(z,Z)/1e4             # Normalized gas temperature
    a_V         = 4.7e-4*(T4^(-1/2))     # Voigt parameter of Lyman-alpha line
    f_weight    = 3/2                    # Mass-weighting factor for the disk

    # Force multiplier without dust:
    
    M_F_nodust(τ0) = 3.51*((a_V*τ0)^(1/3))
    τ0             = 5.33e6*(T4^(-1/2))*Σ0_Msun

    # Force multiplier in general case:

    if Z > 0.0
        # Maximum force multiplier in the presence of dust:

        τ0_star   =  3.58e6*(T4^(-1/4))*(Z^(-3/4))      
        M_F       =  min(f_weight*M_F_nodust(τ0), M_F_nodust(τ0_star))
    else
        M_F       =  f_weight*M_F_nodust(τ0)
    end

    # Resulting contribution from Lyman-α scattering:

    if Z < Z_crit
        LY_ALPHA = M_F*Ψ_Lyα_cIII*(1 - Transmission(k_ionIII, Σ0))
    else
        LY_ALPHA = M_F*Ψ_Lyα_cII*(1 - Transmission(k_ionII, Σ0))
    end

    # PUT IT ALL TOGETHER:

    Pdot_m_rest = ION + WINDS + LY_ALPHA

    return Pdot_m_rest
end

function τ_IR(Σ0,Z)
    """ The factor that controls the feedback from multiple
        scattering of IR radiation. """

    η_IR = 0.9

    # Calculated value of tau_IR:

    tau_IR = η_IR*k_IR(Z)*Pdot_m_rad_noion(Σ0, Z)/(2*pi*G)

    return tau_IR
end



###############################################################
#######                                                 #######
#######           S T E L L A R   M A S S               #######
#######                                                 #######
###############################################################

ϵff = 1.0

function Mstar(M,z,λ,Z,t)
    """ The stellar mass (Msun) formed by time t (yrs)
        in a halo of mass M (Msun) at redshift z,
        with spin parameter λ, and gas metallicity
        Z (relative to Solar). """

    # Surface density at the disk outer edge (for
    # convenience):

    Σ0 = Σdisk0(M,z,λ,t)

    # Critical surface density

    Σcrit = (Pdot_m_rad_noion(Σ0,Z) +  Pdot_m_rest(Σ0,z,Z))/(2*pi*G)

    # tbar (independent of time):

    t_ff  = sqrt(3)*cs(z,Z)*sqrt(1 + (Mach(M,z,λ,Z,t)^2)/3)/(4*G*Σ0)
    tbar  = ϵff*t*yr/t_ff

    # Convenient definitions:

    τIR   = τ_IR(Σ0,Z)
    Σbar  = Σ0*(1 + τIR)/Σcrit
    A     = (Σbar+1)*tbar*log(1 + 1/Σbar)/((Σbar+1)*tbar - Σbar)
    B     = tbar*log(tbar)/((tbar - 1)*((Σbar + 1)*tbar - Σbar))

    # Only if the disk is Toomre unstable will stars form:

    if Q_Toomre(M,z,λ,Z,t) < 1
        f_Q = 1.0
    else
        f_Q = 0.0
    end

    # Put it all together:

    Mstar = f_Q*Σbar*Mdisk(M,z,t)*(A-B)/(1 + τIR)

    return Mstar
end

function dMstardt(M,z,λ,Z,t,dt)
    """ The star formation rate (Msun/yr) during
        the starburst at time t. """

    if t-dt > 0
        dMstardt = (Mstar(M,z,λ,Z,t+dt) - Mstar(M,z,λ,Z,t-dt))/(2*dt)
    else
        dMstardt = (Mstar(M,z,λ,Z,t+dt) - Mstar(M,z,λ,Z,t))/dt

    end

    return dMstardt

end

function Starburst(M,z,λ,Z)
    """ The starburst history and properties of the
        object(s) produced is calculated here. """

    # Time iteration to find the total stellar 
    # mass produced and the starburst time-scale:

    t = 0.0; dt = 1e3 # Starting time & time-step (yr)

    while Efb < Efb_crit & Pdotfb < Fgrav
        
        StellarMass += dMstardt(M,z,λ,Z,t,dt)*dt

    end

end
