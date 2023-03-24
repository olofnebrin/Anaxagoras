""" STARBURST MODULE

    This module contain everything related to calculating the stellar mass
    formed in a single starburst. It also contains the needed calculations
    of the half-mass radius and more. """


###############################################################
#######                                                 #######
#######      S T E L L A R  F E E D B A C K             #######
#######                                                 #######
###############################################################

###### POP II STELLAR FEEDBACK ######

Ψ_ion_cII = 3.27e-8     # Direct stellar emission in ionizing band
Ψ_FUV_cII = 2.40e-8     # Direct stellar emission in FUV band
Ψ_NUV_cII = 1.31e-8     # Direct stellar emission in NUV band
Ψ_opt_cII = 6.73e-9     # Direct stellar emission in optical band

# Lyman-α band:

E_Lyα      = 10.2*eV
Qdot_ionII = 5e46/MSun
Ψ_Lyα_cII  = (2/3)*Qdot_ionII*E_Lyα/c

###### POP III STELLAR FEEDBACK ######

Ψ_ion_cIII = 2.10e-7     # Direct stellar emission in ionizing band
Ψ_FUV_cIII = 0.0         # Direct stellar emission in FUV band
Ψ_NUV_cIII = 0.0         # Direct stellar emission in NUV band
Ψ_opt_cIII = 0.0         # Direct stellar emission in optical band

# Lyman-α band:

Qdot_ionIII = 5e46/MSun
Ψ_Lyα_cIII  = 5.02e-8

###############################################################
#######                                                 #######
#######          D I S K   P R O P E R T I E S          #######
#######                                                 #######
###############################################################

function Mdisk(M, z, t)
    """ Mass of the disk (Msun) at a time t (yrs) after
        gas collapse. Utalizes the self-simimlar solution
        for the gas accretion rate. """

    Mdiv = 2.92e6*(((1+z)/10)^(-3/2))

    if M > Mdiv
        # Temperature and sound speed of gas in the halo:

        Th  = min( Tvir(M,z),  T_Lyα )
        csh = sqrt(kB*Th/(μh*mH))

        # Calculation of g
        
        beta = (1-fB)*((vvir(M,z)/csh)^2) - 2.0
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

    Mdot  = g*(vvir(M,z)^3)/G

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

function Σdisk0(M, z, λ, t)
    """ The disk surface density (g/cm^2) at the outer
        edge of the disk (R = Rdisk) at time t (yrs). """

    Σdisk0 = Mdisk(M,z,t)*Msun/(2*pi*(Rdisk(M,z,λ,t)^2))

    return Σdisk0
end

function T(z,Z)
    """ The temperature of dense gas in the disk (K) 
        as a function of redshift and gas metallicity Z. """

    T = 200.0

    return T
end

function cs(z,Z)
    """ The sound speed in the disk (cm/s) as
        a function of redshift z and gas metallicity
        Z (relative to Solar). """

    return sqrt(kB*T(z,Z)/(μ*mH))
end

###############################################################
#######                                                 #######
#######        S T E L L A R   F E E D B A C K          #######
#######                                                 #######
###############################################################

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
        FUV    =  0.0
        NUV    =  0.0
        OPT    =  0.0
    else
        FUV    =  Ψ_FUV_cII*(1 - Transmission(k_FUVII(Z), Σ0))
        NUV    =  Ψ_NUV_cII*(1 - Transmission(k_NUVII(Z), Σ0))
        OPT    =  Ψ_opt_cII*(1 - Transmission(k_optII(Z), Σ0))
    end

    # PUT IT ALL TOGETHER:

    Pdot_m_rad_noion = FUV + NUV + OPT

    return Pdot_m_rad_noion
end

function Pdot_m_rest(Σ0, z, Z)

    """ Momentum injection rate per Solar mass formed
        (in cm s^-2) from ionizing radiation, Ly-alpha
        scattering, and stellar winds. """

    # CONTRIBUTION FROM IONIZING BAND:

    if Z < Z_crit
        ION = Ψ_ion_cIII*(1 - Transmission(k_ionIII(Z), Σ0))
    else
        ION = Ψ_ion_cII*(1 - Transmission(k_ionII(Z), Σ0))
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

    # Calculated value of tau_IR:

    tau_IR = η_IR*k_IR(Z)*Pdot_m_rad_noion(Σ0, Z)/(2*pi*G)

    return tau_IR
end

###############################################################
#######                                                 #######
#######           S T E L L A R   M A S S               #######
#######                                                 #######
###############################################################

function Mstar(M,z,λ,Z,t)
    """ The stellar mass (Msun) formed by time t (yrs)
        in a halo of mass M (Msun) at redshift z,
        with spin parameter λ, and gas metallicity
        Z (relative to Solar). """

    # Surface density at the disk outer edge (for
    # convenience):

    Σ0 = Σdisk0(M, z, λ, t)

    # Critical surface density

    Σcrit = (Pdot_m_rad_noion(Σ0,Z) +  Pdot_m_rest(Σ0, z, Z))/(2*pi*G)

    # tbar (independent of time):

    t_ff  = sqrt(3)*cs(z,Z)*np.sqrt(1 + (Mach(M, z, spin)^2)/3)/(4*G*Σ0)
    tbar  = ϵff*t*yr/t_ff

    # Convenient definitions:

    τ     = τ_IR(Σ0,Z)
    Σbar  = Σ0*(1 + τ)/Σcrit
    A     = (Σbar+1)*tbar*log(1 + 1/Σbar)/((Σbar+1)*tbar - Σbar)
    B     = tbar*log(tbar)/((tbar - 1)*((Σbar + 1)*tbar - Σbar))

    # Only if the disk is Toomre unstable will stars form:

    if Q_Toomre(M,z,λ,Z) < 1
        f_Q = 1.0
    else
        f_Q = 0.0
    end

    # Put it all together:

    Mstar = f_Q*Σbar*Mdisk(M,z,t)*(A-B)/(1 + τ)

    return Mstar
end

