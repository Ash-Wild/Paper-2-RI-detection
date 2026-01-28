'''
 'Python' library for combined sea state and atmospheric boundary layer model

  Developed by Mark A. Bourassa

  Based on Original Bourassa-Vincent-Wood combined sea state and atmospheric
  boundary layer model developed at Purdue University.  Several upgrades and
  options have been added while at the Center for Ocean-Atmospheric
  Prediction Studies (COAPS) at Florida State University.

  Dr. Bourassa is currently (May 2023) an Professor at the Florida
  State Universiy's Department of Earth, Ocean and Atmospheric Science, and is 
  the Associate Director as the Center for Ocean-Atmospheric Prediciton Studies.
  Please send comments, suggestions, and requests for added features to
      Dr. Mark Bourassa
      Center for Ocean-Atmospheric Prediction Studies
      Florida State University
      Tallahassee, FL 32306-2840

      Email: bourassa@coaps.fsu.edu
      Phone: (850) 644-6923

  The bulk of 'in code' documentation is in the subroutine pmix_. Documentation
  is also on the COAPS website: https://www.coaps.fsu.edu/~bourassa/MFT_code/MFT.py

'''
import math as m

### TRUE and FALSE ###
TRUE = 1
FALSE = 0

# # # Globals # # #
#beta-primes are either zero or unity, depending on whether or not
    #their respective type of roughness element exists.  Beta's are the
    #weights for the different types of elements.
betas =  [0,0]
betac =  [0,0]
betag =  [0,0]
u =  [0,0]
u_orb =  [0,0]
u_sfc =  [0,0]
ustar =  [0,0]
wa =  [0,0]
zz =  [0,0]
hsig =  [0,0]
wstar2 = 0
denair = 0

# # # Constants # # #
B = 0.0190000  # Wu's equilibrium const for capillary waves  []
KV = 0.40000  # von Karman's constant []
G = 9.81000  # gravitational acc. (at sea level)  [m/s^2]
H = 1000.0  # approximate height of boundary layer (at low wind speed)  [m]
PI2 = 6.28318  # two PI  []
R = 287.05000  # specific gas constant for dry air  [J/kg]
RW = 441.0000  # specific gas constant for water vapour  [J/kg]
DTR = 3.14159 / 180.0  # degrees to radians conversion factor

stefanBoltzmann = 5.67e-08
Prt = 1.00  # turbulent Prandtl number []
Prair = 0.7
PrWater = 7.0
Sc = 1.00  # turbulent Schmidt number []
thermalDiff = 2e-05  # or air from Gill p. 71 [m^2/s]
molecularDiff = 2.4e-05  # for H2O in air from Gill p. 68 [m^2/s]

#Next value determines water type for solar absorption
    #1 = Type I, 1.1 = Type IA,  1.2 = Type IB, 2.0 = Type II,
    #3.0 = Type III, 0 = no solar absorption
waterType = 1

### CENTRY ###
def pmix_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT, pressure,
            air_moist_prm, air_moist_val_dum, sfc_moist_prm, sfc_moist_val_dum,
            salinity, ss_prm, ss_val_dum, t_air, sst_prm, t_skin, ref_ht_wind,
            ref_ht_tq, astab, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm,
            stable_prm, A_oil, pass_by_ref_params ):
    '''
        'main' subroutine for combined sea state and atmospheric flux model

        Includes:
           BVW sea state parameterization,
           BVW atmospheric boundary layer (roughness length) model (optional),
           Toba's 2/3 wave relation,
           Monin-Obhukov atmoshpheric stability parameterization,
           Godfrey and Beljaar's convective parameterization.

        Options for dynamic   parameters
           Most users will want to input wind speeds.  However, there is a
           growing need to convert scatterometer observations of stress or
           equivalent neutral wind speed to winds and/or stresses.
           dyn_in_prm    Parameter treated as known
           0             Wind speed,
           1             Friction velocity,
           2             Surface stress,
           3             Equivalent neutral wind speed. Note 'equivalent neutral'
                         is NOT equivelent to 'neutral'.

        Options for seastate parameterizations (or parameters):
           There are three possible seastate assumptions: any one of the
           following can be treated as known: wind-wave stability parameter,
           phase speed, or wave age.
           ss_prm     Parameter treated as known
           0          Wind-wave stability parameter (set to 1.0 for local
                      equilibrium),
           1          Phase speed of the dominant waves  [m/s],
           2          Wave age the dominant waves (cp/ustar)  [],
           3          Significant Wave Height  [m],
           4          Significant slope  [],
           5          Period of the dominant waves [s],
           6          Orbital velocity of dominant waves [m/s].

        Options for atmospheric moisture input:
           Choose the moisture parameter that is easiest for you to deal with:
           air_moist_prm  Parameter for moisture of air
           0          Absolute humidity at the reference height of the thermometer
                      and humidity sensor [g vapor / g air],
           1          Relative Humidity [fraction],
           2          Dew point temperature [C],
           3          Wet bulb temperature [C].

        Options for surface moisture input:
           Choose the moisture parameter that is easiest for you to deal with:
           sfc_moist_prm  Parameter for moisture of air
           0          Absolute humidity 'at' (near) the surface [g vapor / g air],
           1          Relative Humidity [fraction],
           2          Dew point temperature [C],
           3          Wet bulb temperature [C].


        Definitions:
        Input parameters  (passed through call):
            CON_P       Convective parameter (unitless).  Recammended value between
                        0.7 and 1.25.  For details see TOGA NOTES #4.
            CRIT        Convergence criterion (unitless fraction).
            air_moist_prm   Moisture parameter index.
            air_moist_val   Value of the parameter corresponding to the above index.
            press       Atmospheric surface pressure (Pa).
            s           Salinity (unitless).
            ss_prm      Seastate parameter index.
            ss_val      Value of the parameter corresponding to the above index.
            ta          Air temperature at the reference height of the thermometer
                        and humidity sensor (degrees C).
            u           Mean wind speed (m/s) at the height (zref) of the annemometer.
            u_sfc       Mean surface current (m/s).
            zref        Height (metres) of the wind observations.
            zrefq       Height (metres) of the humidity observations.  See zreft.
            zreft       Height (metres) of the temperature observations.  Note: in
                        the current version of the code this must equal to height
                        of the humidity observations.
            astab       Option for atmospheric stability condition:
                            0) atmospheric stability is assumed to be neutral,
                            1) assumption of a neutral atmospheric stability is
                               not made.
            theta_ang   Direction of 'dynamic input' (e.g., wind or stress) relative
                        to the direction of propagation of the dominant waves.

        Other parameters:
            bstar       Component of bouyancy flux (-bstar * ustar).
            cmin        Minimun phase speed of water waves - determined through
                        the phase relation   [m/s].
            CP          Specific heat of air [J/kg].
            denair      density of air  [kg/m^3].
            denwat      density of water  [kg/m^3].
            LV          Latent heat of vaporization  [J/kg].
            sfcten      suface tension of air-water interface [N/m].
            ww          Wind-wave stability parameter (unitless).  Equal to one for
                        local equilibrium; greater than one for rising seas; less
                        than one (and greater than zero) for falling seas.
            ww_eql      Value of U(Hs)/cp for local wind-wave equilibrium (0.83).
            zz          wind reference height divided by momentum roughness lentgh [],
            zzq         moisture reference height divided by moisture roughness
                        lentgh  [],
            zzt         temperature reference height divided by temperature roughness
                        lentgh  [],

        Output parameters:
            cp          dominant phase speed of gravity waves   [m/s],
            hsig        significant wave height  [m],
            inv_monin   inverse Monin-Obhukov scale length  [m^-1],
            lhf         latent heat flux  [W/m^2],
            qstar       scaling parameter for moisture  [],
            shf         sensible heat flux  [W/m^2],[W/m^2],
            tau         stress vector   [N/m^2];  tau[0] is parallel the direction
                        of wave propagation, and tau[1] is perpendicular tau[0] in
                        a right handed coordinate system (with the positive
                        vertical axis pointing upward,
            tstar       scaling term for potential temperature  [degrees C]
            ustar       friction velocity  [m/s],
            wa          wave age,  cp/ustar  [],
            zo_m        momentum roughness lentgh (two components: parallel and
                        perpendicular direction of wave propagation) [],
    '''

    global flux_model_dum, no_capw, no_sfcten, smith88, stab_prm
    global TQzo_prm, use_dh, wave_stab_intr, monin_inv
    global betac_prime, betag_prime
    global CP, denair, LV
    global u, u_orb, wa, zz

    #definitions for flux model types
    #   0  BVW
    #   1  Smith 1988

    flux_model_dum = flux_model
    smith88 = FALSE
    no_capw = FALSE
    wave_stab_intr = TRUE
    use_dh = FALSE
    no_sfcten = FALSE

    if flux_model < 0:
        z0_prm = z0_mom_prm
        TQzo_prm = z0_TQ_prm
        stab_prm = stable_prm
    else:
        TQzo_prm = 0
        stab_prm = 0
        if flux_model == 0: #BVW
            z0_prm = 0
        elif flux_model == 1:   #Smith 1988
            z0_prm = 0
            ss_prm = 2
            ss_val_dum = 43.64
            dyn_in_val2 = 0.0
            smith88 = TRUE
            no_capw = TRUE
            wave_stab_intr = FALSE
            stab_prm = 0
        elif flux_model == 2:
            #BVW without interaction between waves and atmospheric stab
            z0_prm = 0
            wave_stab_intr = FALSE
        elif flux_model == 3:
            #BVW with Smith (1988) stability parameterization
            z0_prm = 0
            stab_prm = 1
        elif flux_model == 4:
            #BVW without capillary waves
            z0_prm = 0
            no_capw = TRUE
        elif flux_model == 5:
            #BVW without sfcten in phase relations
            z0_prm = 0
            no_sfcten = TRUE
        elif flux_model == 6:
            #BVW without capillary waves and without sfcten in the phase relation
            z0_prm = 0
            no_sfcten = TRUE
            no_capw = TRUE
        elif flux_model == 7:
            #Taylor and Yelland (2001) roughness length
            z0_prm = 2
            no_sfcten = TRUE
            no_capw = TRUE
        elif flux_model == 8:
            #Taylor and Yelland (2001) roughness length + capillary
            z0_prm = 3
        elif flux_model == 9:
            #Bourassa 2006 friction velocity with CFC zot and zoq
            z0_prm = 1
            TQzo_prm = 1
        elif flux_model == 10:
            #case 9 with displacement height estimated within the code
            z0_prm = 1
            TQzo_prm = 1
            use_dh = TRUE
        elif flux_model == 11:
            #Bourassa 2006 friction velocity with Zilitinkevich et al. zot and zoq
            z0_prm = 1
            TQzo_prm = 2
        elif flux_model == 12:
            #Bourassa 2006 friction velocity with LKB zot and zoq
            z0_prm = 1
            TQzo_prm = 3
        elif flux_model == 13:
            #Bourassa 2006 friction velocity with COARE3.0 zot and zoq
            z0_prm = 1
            TQzo_prm = 4
        elif flux_model == 14:
            #Bourassa 2006 friction velocity with wall theory for moisture
            z0_prm = 1
            TQzo_prm = 0
        elif flux_model == 15:
            #Bourassa 2006 friction velocity with CFC zot and zoq and oil spill z0m and zero LHF
            z0_prm = 4
            TQzo_prm = 1
        else:
            print('Invalid choice of flux_model parameter: '+ flux_model)
            exit(1)

    MS = 58.443000
    # molecular mass of salt  [kg/kmol]
    MW = 18.016000
    # molecular mass of water [kg/kmol]
    betag_prime = 1.0
    betac_prime = 1.0
    if no_capw:
        betac_prime = 0.0
    betas_prime = 0.0 + float(smith88)
    betas[0] = 1.0
    betas[1] = 1.0
    ww_eql = 0.345 / KV

    CRIT = float(CONV_CRIT)
    CON_P = float(CONVECT)
    press = float(pressure)
    air_moist_val = float(air_moist_val_dum)
    sfc_moist_val = float(sfc_moist_val_dum)
    s = float(salinity)
    ss_val = float(ss_val_dum)
    ta = float(t_air)
    tskin = float(t_skin)
    zref = float(ref_ht_wind)
    #/* the height of the temperature sensor must be equal to the height of
    #the humidity sensor.  */
    zrefq = float(ref_ht_tq)
    zreft = zrefq
    if ( z0_prm != 7 ):
       pass_by_ref_params['zo_m'][0] = 0.0001
       pass_by_ref_params['zo_m'][1] = 0.0 
    else:
       zz[0] = zref / pass_by_ref_params['zo_m'][0] 
       zz[1] = 0.0

    # setup the knowns and initial gueses for the seastate.
    wa[0] = 20.0 # initial guess
    wa[1] = 28.0
    cp = 4.0  # initial guess
    fixedcp = FALSE
    fixedwa = FALSE
    hsig[0] = 0.01 # initial guess assumes waves exist
    hsig[1] = 0.01 # initial guess assumes waves exist

    fixedcp = FALSE
    fixedwa = FALSE
    fixed_hsig = FALSE
    ww = 1.0

    if ss_prm == 0:
        # the wind-wave stability parameter is specified
        ww = ss_val
    elif ss_prm == 1:
        # phase speed is known
        fixedcp = TRUE
        cp = ss_val
    elif ss_prm == 2:
        # wave age (ustar/cp) is known
        fixedwa = TRUE
        wa[0] = ss_val
    elif ss_prm == 3:
        # significant wave height is known
        fixed_hsig = TRUE
        hsig[0] = ss_val
    elif ss_prm == 4:
        # significant slope is known
        pass
    elif ss_prm == 5:
        # period of the dominant waves is known
        pass
    elif ss_prm == 6:
        # Orbital velocity is known
        u_orb[0] = ss_val
    else:
        print('Invalid choice of ss_prm: '+ ss_prm)
        exit(1)

    # setup the atmospheric stability parameters
    monin_inv = 0.0
    if astab == 0:
        neutral = TRUE
    elif astab == 1:
        neutral = FALSE
    elif astab == 2: 
        neutral = FALSE
        monin_inv = pass_by_ref_params['z_over_L'] / zref
    else:
        print('Invalid atmospheric stability option: '+ astab)
        print('Allowable options are: 0 (neutral), 1 (non-neutral), and 2 (input z/L)')
        exit(1)

    betag_set = betag_prime
    betac_set = betac_prime
    betas_set = betas_prime

    # determine the surface tension
    sfcten = 7.6100000E-2 - 1.55000E-4 * tskin + 2.77000E-5 * s * MS / MW

    # determine the latent heat of vaporization
    LV = 4186.8 * ( 597.31 - 0.56525 * tskin )

    # use the atmospheric moisture parameter to determine the specific humidity [kg/kg].
    qmixa = find_q_( air_moist_prm, air_moist_val, press, ta )
    if qmixa == -1.0:
        print('Invalid choice of air_moist_prm: '+ air_moist_prm)

    # use the surface moisture parameter to determine the specific humidity [kg/kg].
    qmixw = find_q_( sfc_moist_prm, sfc_moist_val, press, tskin )
    if qmixw == -1.0:
        print('Invalid choice of sfc_moist_prm: '+ sfc_moist_prm)
    if z0_prm == 4:
        qmixw = qmixa

    # determine the heat capacity at constant pressure
    CP = 1004.0 * ( 1.0 + 0.9 * qmixw )

    # determine the density of the moist air
    denair = press / (R * (ta + 273.14) * (1.00000000 + 0.6100000 * qmixa ))

    # determine the density of pure water
    denwat = ((999.8396 + 18.224944 * tskin - 0.007922210 * tskin * tskin ) /
        ( 1.0000000 + 0.018159725 * tskin ))

    # determine the density of saline water
    va = (( 12.97000 + 0.23400000 * tskin - 4.2100000E-3 * tskin * tskin +
        2.8570000E-5 * tskin * tskin * tskin ) * 10.00000E-3 +
        m.sqrt( s * denwat * 10.00000E-3 ) * 2.9820000E-3 - 4.9700000E-5 *
        tskin + 6.0320000E-7 * tskin * tskin)
    denwat = denwat * ( 1.0000000 + s ) / ( 1.000000 + va * denwat * s / MS )

    # determine the viscosity of air
    nu = ( 1.3260000e-5 * (1.00000 + 0.006542 * ta + 8.301000E-6 * ta * ta +
        4.840000E-9 * ta * ta * ta )) / denair

    # setup known inputs for the dynamic input (typically wind speed)
    fixed_ustar = FALSE
    fixed_uen = FALSE

    if dyn_in_prm == 0:   # wind speed is specified
        # split the wind into component parallel and perpendicular to the mean
        #      direction of wave propagation.
        u[0] = float(dyn_in_val)
        u[1] = float(dyn_in_val2)
    elif dyn_in_prm == 1:   # friction velocity is specified
        fixed_ustar = TRUE
        # split the wind into component parallel and perpendicular to the mean
        #      direction of wave propagation.
        ustar[0] = float(dyn_in_val)
        ustar[1] = float(dyn_in_val2)
        if ( m.fabs( ustar[0] ) < 0.00001 ):
            ustar[0] = 0.00001
        if ( m.fabs( ustar[1] ) < 0.00001 ):
            ustar[1] = 0.00001
    elif dyn_in_prm == 2: # surface stress (magnitude) is known
        fixed_ustar = TRUE
        # split the wind into component parallel and perpendicular to the mean
        #     direction of wave propagation.
        ustar[0] = m.sqrt( float(dyn_in_val)  / denair )
        ustar[1] = m.sqrt( float(dyn_in_val2) / denair )
        if ( m.fabs( ustar[0] ) < 0.00001 ):
            ustar[0] = 0.00001
        if ( m.fabs( ustar[1] ) < 0.00001 ):
            ustar[1] = 0.00001
    elif dyn_in_prm == 3:  # equivalent neutral wind speed is known
        fixed_uen = TRUE
        # split the wind into component parallel and perpendicular to the mean
        #     direction of wave propagation.
        u[0] = float(dyn_in_val)
        u[1] = float(dyn_in_val2)
    else:
        print('Invalid choice of dyn_in_prm: '+ dyn_in_prm)
        exit(1)

    # One the first pass assume than waves exist
    hsig[1] = 0.01

    if ( z0_prm != 7 ):
      zz[0] = 10000.0
    zz[1] = 10000.0

    count = solve( fixedcp, fixedwa, fixed_hsig, fixed_uen, fixed_ustar,
        neutral, CON_P, cp, CRIT, denwat, nu, qmixa, qmixw,
        sfcten, ss_prm, ss_val, ta, tskin, Qnet, zref, zreft,
        zrefq, ww, ww_eql, betag_set, betac_set, betas_set, z0_prm, warn, A_oil, astab )

    if ( count <= 1 or ustar[0] == 10.0 or ustar[1] == 10.0 or ( wa[0] == 1.08
         and betag[0] > 0.5 and flux_model != 9 ) or ( wa[0] == 250.0 and betag[0]
         > 0.5 and flux_model != 9 ) or m.fabs(monin_inv) == 10.0 ):
        count = -1
        if ( warn ):
            print( 'non-conv diags: ', count, ustar[0], ustar[1], tstar, qstar,
                wa[0], monin_inv, betag[0])
        if ( not warn and ( ustar[0] == 10.0 or ustar[1] == 10.0 ) ):
            ustar[0] = 0.0  # usually a good approximation
            ustar[1] = 0.0
            print('friction velocity set to zero.')

    if ( count > 1 or not warn ):
        dud = m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
        pass_by_ref_params['tau'] = dict()
        pass_by_ref_params['tau'][0] = float(denair * dud * ustar[0])
        pass_by_ref_params['tau'][1] = float(denair * dud * ustar[1])
        pass_by_ref_params['shf'] = float(-denair * CP * dud * tstar)
        pass_by_ref_params['lhf'] = (1.0 - A_oil) * float(-denair * LV * dud * qstar)
        pass_by_ref_params['u_star'] = dict()
        pass_by_ref_params['u_star'][0] = float(ustar[0])
        pass_by_ref_params['u_star'][1] = float(ustar[1])
        pass_by_ref_params['t_star'] = float(tstar)
        pass_by_ref_params['q_star'] = float(qstar)
        pass_by_ref_params['z_over_L'] = float(not neutral) * float( zref * monin_inv )
        pass_by_ref_params['wave_age'] = float(wa[0])
        pass_by_ref_params['dom_phs_spd'] = float( betag[0] * cp_com )
        pass_by_ref_params['h_sig'] = float( betag[0] * hsig[0] )
        pass_by_ref_params['ww_stab'] = float(ww)
        pass_by_ref_params['zo_m'] = dict()
        pass_by_ref_params['zo_m'][0] = zref / zz[0]
        pass_by_ref_params['zo_m'][1] = zref / zz[1]
    else:
        # printf( "Warning: non-convergence" )
        pass_by_ref_params['tau'] = dict()
        pass_by_ref_params['tau'][0] = -9999.0
        pass_by_ref_params['tau'][1] = -9999.0
        pass_by_ref_params['shf'] = -9999.0
        pass_by_ref_params['lhf'] = -9999.0
        pass_by_ref_params['u_star'] = dict()
        pass_by_ref_params['u_star'][0] = float(ustar[0])
        pass_by_ref_params['u_star'][1] = float(ustar[1])
        pass_by_ref_params['t_star'] = float(tstar)
        pass_by_ref_params['q_star'] = float(qstar)
        pass_by_ref_params['z_over_L'] = float( zref * monin_inv )
        pass_by_ref_params['wave_age'] = float(wa[0])
        pass_by_ref_params['dom_phs_spd'] = float(cp_com)
        pass_by_ref_params['h_sig'] = float(hsig[0])
        pass_by_ref_params['ww_stab'] = float(ww)
        pass_by_ref_params['zo_m'] = dict()
        pass_by_ref_params['zo_m'][0] = -9999.0
        pass_by_ref_params['zo_m'][1] = -9999.0

    return count


def find_q_( moist_prm, moist_val, press, temperature ):

    # printf( "Bq%i %f %f %f\n", moist_prm, moist_val, press, temperature )
    # use the atmospheric moisture parameter to determine the specific humidity [kg/kg].
    if ( moist_prm == 0 ): # input is specific humidity
        spec_hum = moist_val
    elif ( moist_prm == 1 ): # relative humidity (fraction) is used as input
        # determine saturation vapour pressure over a smooth surface of water, at the temperature
        satvp = ( ( 1.000700 + 3.460000E-8 * press ) * 611.2100 *
                m.exp( 17.50200 * temperature / ( 240.97000 + temperature ) ) )
        spec_hum = ( moist_val * satvp * 0.6220 / ( press - 0.3780 * moist_val * satvp ) )
    elif ( moist_prm == 2 ): # dew point temperature is used as input
        # determine the latent heat of vaporization
        satvp = ( ( 1.000700 + 3.460000E-8 * press ) * 611.2100 *
                m.exp( 17.50200 * moist_val / ( 240.97000 + moist_val ) ) )
        spec_hum = ( 0.6220 * satvp / ( press - 0.3780 * satvp ) )
    elif ( moist_prm == 3 ): # wet bulb temperature is used as input
        spec_hum = ( 0.0000066 * ( 1.0 + 0.00115 * moist_val ) ) # dummy use of spec_hum
        satvp = ( (1.000700 + 3.460000E-8 * press) * 611.2100 *
                m.exp( 17.50200 * moist_val / ( 240.97000 + moist_val ) ) -
                spec_hum * press * (temperature - moist_val ) )
        spec_hum = 0.6220 * satvp / ( press - 0.3780 * satvp )
    else:
        print ( 'Invalid choice of moist_prm: ', moist_prm )
        spec_hum = -1.0
    return spec_hum

def ht_adj_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
    pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
    salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
    z_wanted, astab, eqv_neut_prm, Qnet, warn, flux_model, z0_mom_prm,
    z0_TQ_prm, stable_prm, A_oil, pass_by_ref_params):

    global flux_model_dum, no_capw, TQzo_prm, CP, denair
    global ustar, zzq, zzt

    zoN = [None]*2
    TQzo_prm = 0

    # The physics breaks down for zero wind speed.  If the wind speed is zero, make an adjustment
    zero_speed_flag = FALSE
    dyn_in_mag = ( dyn_in_val**2 + dyn_in_val2**2 )**0.5 
    if ( dyn_in_mag <= 0.001 ):
        zero_speed_flag = TRUE
        dyn_in_val = 0.001

    if ( flux_model == 9 or flux_model == 10 ):
        TQzo_prm = 1
    if ( eqv_neut_prm == 2 ):
        eqv_neut = FALSE
    elif ( eqv_neut_prm == 1 ):
        eqv_neut = TRUE
    else:
        eqv_neut = FALSE
        # calculate neutral roughness length
 #       astab_temp = 0
    bvw_flag = pmix_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
        pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
        salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
        astab, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm, A_oil,
        pass_by_ref_params )
    # record neutral roughness length
    zoN[0] = pass_by_ref_params['zo_m'][0]
    zoN[1] = pass_by_ref_params['zo_m'][1]

    # calculate fluxes and height adjustment for non-neutral (or specified) conditions 
 #   if ( eqv_neut != 0 ): 
 #      bvw_flag = pmix_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
 #       pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
 #       salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
 #       astab_temp, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm, A_oil, 
 #       pass_by_ref_params )

    # LV = 4186.8 * ( 597.31 - 0.56525 * t_skin )
    q_sfc = find_q_( sfc_moist_prm, sfc_moist_val, pressure, t_skin )
    # determine the heat capacity at constant pressure
    CP = 1004.0 * ( 1.0 + 0.9 * q_sfc )
    alpha = 1.0/(t_skin+273.15)

#    print( 'ht_adj: zz = ', zz[0], pass_by_ref_params['zo_m'][0]  )

    # correct to the height of z_wanted
    if ( eqv_neut_prm == 0 ):
        phi_u = phi_u_f( z_wanted, z_wanted /  pass_by_ref_params['zo_m'][0] )
        nu = 1.3260000e-5 * (1.00000 + 0.006542 * t_skin +
            8.301000E-6 * t_skin * t_skin + 4.840000E-9 * t_skin * t_skin * t_skin )
        ustar_mag = m.sqrt( ustar[0]*ustar[0] + ustar[1] * ustar[1] )
        zo_mag = m.sqrt(  pass_by_ref_params['zo_m'][0]* pass_by_ref_params['zo_m'][0] +  pass_by_ref_params['zo_m'][1] *  pass_by_ref_params['zo_m'][1] )
        # the following 6 lines seem unnecessary
        viscosity = 1.4e-5
        denair = pressure / (R * ( t_skin + 273.14) * (1.0 + 0.61 * q_sfc ) )
        renewalTime = getRenewalTime( zo_mag, ustar_mag, denair, CP, alpha,
            Qnet, viscosity )
        Stanton = getStanton( ustar_mag, renewalTime)
        Dalton = getDalton( ustar_mag, renewalTime)

        zzt = zzt * z_wanted / ref_ht_tq
        zzq = zzq * z_wanted / ref_ht_tq
        phi_T = phi_T_f( z_wanted, zzt )
        phi_Q = phi_Q_f( z_wanted, zzq )

    elif ( eqv_neut_prm == 1 ):
        # Tang and Liu eqv neut winds
        phi_u = 0.0
        phi_T = 0.0
        phi_Q = 0.0

    elif ( eqv_neut_prm == 2 ):
        # if option was selected, calculate equivalent friction velocities as
        #     defined by Geernaert and Katsaros (1986)
        pass_by_ref_params['ustar'][0] = pow( dyn_in_val / ustar[0] + phi_u / KV -
            m.log( zoN[0] /  pass_by_ref_params['zo_m'][0] ) / KV, -1.0 ) * dyn_in_val
        pass_by_ref_params['ustar'][1] = pow( dyn_in_val2 / ustar[1] + phi_u / KV -
            m.log( zoN[1] /  pass_by_ref_params['zo_m'][1] ) / KV, -1.0 ) * dyn_in_val2
        phi_u = 0.0

    if ( zero_speed_flag == FALSE and pass_by_ref_params['zo_m'][0] > 0.0 ):
        pass_by_ref_params['u_at_z'] = m.sqrt( pow(ustar[0] * ( m.log( z_wanted /  pass_by_ref_params['zo_m'][0] ) + phi_u ) / KV, 2.0) +
        pow(ustar[1] * ( m.log( z_wanted /  pass_by_ref_params['zo_m'][1] ) + phi_u ) / KV, 2.0) )
    else: 
        pass_by_ref_params['u_at_z'] = 0.0
    pass_by_ref_params['t_at_z'] = t_skin + Prt * tstar * ( m.log( zzt ) + phi_T ) / KV
    pass_by_ref_params['q_at_z'] =  q_sfc + Sc * qstar * ( m.log( zzq ) + phi_Q ) / KV

    return( bvw_flag )

###END CENTRY###

def f_bstar( neutral, qmixa, ta, minimum_value ):

    bstar = (float(not neutral) * G * ( ( 1.00 + 0.610 * qmixa ) * tstar + 0.6100 *
            ( ta + 273.15 ) * qstar ) / (  ta + 273.15 * ( 1.0 + 0.610 * qmixa ) ))
    if ( m.fabs( bstar ) < minimum_value ):
        if ( bstar == 0.0 ):
            bstar = minimum_value
        else:
            bstar = minimum_value * bstar / m.fabs( bstar )
    return bstar

# # # # # solve # # # # #
def solve( fixedcp, fixedwa, fixed_hsig, fixed_uen, fixed_ustar, neutral,
    CON_P, cp, CRIT, denwat, nu, qmixa, qmixw, sfcten, ss_prm, ss_val,
    ta, tskin, Qnet, zref, zreft, zrefq, ww, ww_eql, betag_set, betac_set,
    betas_set, z0_prm, warn, A_oil, astab ):

    global betac_prime, betag, betac, betas, betag_prime, betas_prime
    global cmin, monin_inv, qstar, tstar
    global zzq, zzt

    #    convert the temperatures to equivalent potential temperature
    #    estimated through the dry static energy
    ustar_old = [None] * 2
    wa_old = [None] * 2
    zz_old = [None] * 2
    tskin_ept = tskin
    ta_ept = ta + G * zreft / CP

    betag_prime = betag_set
    betag[0] = betag_prime
    betag[1] = betag_prime
    betac_prime = betac_set
    betac[0] = betac_prime
    betac[1] = betac_prime
    betas_prime = betas_set
    betas[0] = betas_prime
    betas[1] = betas_prime

    #    intial guess for value of friction velocity (ustar)
    #    component parallel direction of wave motion
    ustar[0] = ( float(fixed_ustar) * ustar[0] +
        float(not fixed_ustar) * 0.003 )
    #    component perpendicular direction of wave motion
    ustar[1] = ( float(fixed_ustar) * ustar[1] +
        float(not fixed_ustar) * 0.003 )
    tstar = 0.0003 * (ta - tskin)
    if ( m.fabs( tstar ) < 0.001 ):
        tstar = -0.001
    qstar = 0.0003 * (qmixa - qmixw)
    if ( m.fabs( qstar ) < 0.0001 ):
        qstar = -0.0001

    ustar_old[0] = ustar[0] / 2.00
    ustar_old[1] = ustar[1] / 2.00
    tstar_old = tstar / 2.00
    qstar_old = qstar / 2.00
    cmin = m.sqrt( 2.0 * m.sqrt( G * sfcten / denwat ) )
    u_sfc[0] = 0.0
    u_sfc[1] = 0.0
    if ( ss_prm != 6 ):
        u_orb[0] = 0.0
    else:
        u_orb[0] = ss_val
    u_orb[1] = 0.0
    #    initially assume a neutral atmospheric stability (monin_inv = 0) unless z/L is input (astab = 2)
    if ( astab != 2 ):
        monin_inv = 0.0

    unreasonable = TRUE
    while ( unreasonable == TRUE ):
        # iterate over atmospheric stability and friction velocity
        if ( astab != 2 ):
            monin_inv_old2 = -0.001
        else:
            monin_inv_old2 = monin_inv
        ustar_mag_old = 0.0
        ustar_mag = m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
        count_Lu = 0
        while ( ( m.fabs( monin_inv - monin_inv_old2 ) > CRIT*m.fabs( monin_inv_old2 ) or
            m.fabs( ustar_mag - ustar_mag_old ) > CRIT * m.fabs( ustar_mag_old ) or
            count_Lu < 2 ) and count_Lu < 30 and not ( fixed_ustar and count_Lu >= 1 ) ):

            count_Lu = count_Lu + 1
            # determine roughness lengths for heat and moisture profiles
            viscosity = 1.4e-5
            alpha = 1.0/(ta+273.15)
            zo_mag = zref / m.sqrt( zz[0] * zz[0] + zz[1] * zz[1] )
            # zo_mag = zref / zz[0] # test option
            renewalTime = getRenewalTime( zo_mag, ustar_mag, denair, CP, alpha,
                Qnet, viscosity )
            Stanton = getStanton( ustar_mag, renewalTime)
            Dalton = getDalton( ustar_mag, renewalTime)


            zzt = z_o_zt( zref, zreft, nu, zo_mag, Stanton )
            zzq = z_o_zq( zref, zrefq, nu, zo_mag, Dalton )

            ustar_mag_old = ustar_mag
            monin_inv_old2 = monin_inv

            # iterate over solutions to tstar, qstar, and monin_inv
            tstar_old = 0.5 * tstar
            qstar_old = 0.5 * qstar
            monin_inv_old = monin_inv
            count_tqL = 0
            while ( ( m.fabs( tstar - tstar_old ) > CRIT * m.fabs( tstar_old ) or
                m.fabs( qstar - qstar_old ) > CRIT * m.fabs( qstar_old ) or
                m.fabs( monin_inv - monin_inv_old ) > 0.1*CRIT * m.fabs( monin_inv_old )
                or count_tqL < 3 ) and count_tqL < 40 ):

                count_tqL = count_tqL + 1

                tstar_old = tstar
                qstar_old = qstar
                monin_inv_old = monin_inv

                # determine tstar  (assuming the current L is correct)
                # print( "test2a: %f %f %f\n", (float)tstar, (float)qstar, (float)monin_inv );
                # print( "test2c: %f %f %f %f %f\n", (float)delta_theta, (float)KV,
                # (float)( log( zzt+1.0 ), (float)(!neutral) * phi_T ), (float)Prt )
                if ( ta_ept == tskin_ept ):
                    tstar = 0.0
                else:
                    # determine the eqv. pot. temperature stability term
                     phi_T = phi_T_f( zreft, zzt )

                     delta_theta = ta_ept - tskin_ept
                     if ( m.fabs( delta_theta ) < 0.001 ):
                         delta_theta = -0.001
                     tstar = ( KV * delta_theta / ( m.log( zzt+1.0 ) +
                         float(not neutral) * phi_T ) / Prt )

                # determine qstar (assuming current L is correct
                # printf( "test2b: %f %f %f\n", (float)tstar, (float)qstar, (float)monin_inv );

                if ( qmixa == qmixw ):
                     qstar = 0
                else:
                    # determine the moisture stability term  */
                    phi_Q = phi_Q_f( zreft, zzq )

                    delta_q = qmixa - qmixw
                    if ( m.fabs( delta_q ) < 0.0001 ):
                        delta_q = -0.0001
                    qstar = ( KV * ( delta_q ) / ( m.log( zzq+1.0 ) +
                        float(not neutral) * phi_Q ) / Sc )
                        # printf( "%f %f %f\n", delta_q, zzq, phi_Q )

                # determine the bouyancy flux
                bstar = f_bstar( neutral, qmixa, ta, 0.000001 )

                # !!!!This assumes temperature and humidity are measured at THE SAME HEIGHT

                # determine the inverse Monin-Obhukov scale length: 1/L
                if ( astab != 2 ):
                   monin_inv = KV * bstar / ( ustar[0]*ustar[0] + ustar[1] * ustar[1] )
                # insure that 1/L is within reasonable bounds
                max = 8.0
                min = 0.0001
                if ( m.fabs(monin_inv) < min ):
                    monin_inv = min * monin_inv / m.fabs( monin_inv )
                if ( m.fabs( monin_inv ) > max ):
                    monin_inv = max * monin_inv / m.fabs( monin_inv )
                # end of iteration on tstar, qstar, and monin_inv

            # printf( "test3: %f %f %f\n", (float)tstar, (float)qstar, (float)monin_inv )
            done = FALSE
            count = 0
            while ( ( done == FALSE or count < 2 ) and count < 100 ):
                count = count + 1
                done = TRUE
                wa_old[0] = wa[0]
                wa_old[1] = wa[1]
                zz_old[0] = zz[0]
                zz_old[1] = zz[1]
                first = TRUE
                # iteratively solve for u* for non-neutral conditions
                if ( not fixed_ustar ):
                    find_ustar( first, fixedcp, fixedwa, fixed_hsig, fixed_uen,
                        neutral, count, CON_P, cp, CRIT, denwat, nu, qmixa, sfcten, ta,
                        ww, ww_eql, ss_prm, ss_val, zref, zreft, z0_prm, A_oil, astab )
                        # printf( "test A %f %f %f %f\n", (float)ustar[0], (float)ustar[1],
                        #     (float)tstar, (float)qstar )
                else:
                    i = 0
                    if ( fixed_hsig and i == 0 ):
                        hsig_known = TRUE
                    else:
                        hsig_known = FALSE
                    trouble = z0_and_waves( fixedcp, fixedwa, neutral, hsig_known,
                        denwat, nu, sfcten, zref, cp, ww, ww_eql,
                        ss_prm, ss_val, count, i, z0_prm, A_oil )
                    i = 1
                    hsig_known = FALSE
                    trouble = z0_and_waves( fixedcp, fixedwa, neutral, hsig_known,
                        denwat, nu, sfcten, zref, cp, ww, ww_eql,
                        ss_prm, ss_val, count, i, z0_prm, A_oil )

                if ( ( m.fabs( ( wa[0] - wa_old[0]) / wa[0] ) > CRIT or
                    m.fabs( ( wa[1] - wa_old[1] ) / wa[1] ) > CRIT or
                    m.fabs( ( zz[0] - zz_old[0] ) / zz[0] ) > CRIT or
                    m.fabs( ( zz[1] - zz_old[1] ) / zz[1] ) > CRIT ) and count < 60 ):
                        done = FALSE

            ustar_mag = m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
            # end of iteration over friction velocity and atmospheric stability

            # determine if the phase speed of the dominant waves is physically
            # reasonable. If the phase speed is less than the physical minimum,
            # then assume that no capillary waves are pressent. After the second
            # pass, if phase speed is still unacceptable then assume the surface
            # is smooth.

            # if ( smith88 && monin_inv > 0 ) printf( "zo = %f %f\n", zref / zz[0], u[0] )
            # if ( smith88 && monin_inv > 0 ) printf( "1 %f %f %f %f %f %f %f\n", u[0],
            # betas_prime, betac_prime, betag_prime, betas[0], betac[0], betag[0])

        if ( fixedwa or fixedcp or ( betag[0] != 0.0 and
            wa[0] * m.fabs( ustar[0] ) > cmin ) or ( betag[1] != 0.0 and
            wa[1] * m.fabs( ustar[1] ) > cmin ) or betas_prime == 1.0 ):
            unreasonable = FALSE
        else:
            if ( betag[0] != 0.0 and betag[1] != 0.0 and betac_prime != 0.0 ):
                if ( m.fabs( u[0] - u_sfc[0] ) > m.fabs( u[1] - u_sfc[1] ) ):
                    betag[1] = 0.0
                else:
                    betag[0] = 0.0
            elif ( betac_prime > 0.0 ):
                betag[0] = 1.0
                betag[1] = 1.0
                betac_prime = 0.0
            else:
                if ( betag[0] > 0.0 and wa[0] * m.fabs( ustar[0] ) < cmin ):
                    betag[0] = 0.0 + float(smith88)
                if ( betag[1] > 0.0 and wa[1] * m.fabs( ustar[1] ) < cmin ):
                    betag[1] = 0.0 + float(smith88)
                if ( betag[0] == 0.0 and betag[1] == 0.0 ):
                    betag_prime = 0.0 + float(smith88)
                    betas_prime = 1.0
    # if ( smith88 && monin_inv > 0 ) printf( "2 %f %f %f %f %f %f %f\n", u[0], betas_prime,
    # betac_prime, betag_prime, betas[0], betac[0], betag[0])
    return count_Lu + fixed_ustar

# # # # # f_ustar # # # # #
def f_ustar( fixed_uen, neutral, bstar, CON_P, phi_u, delta_u, ww, zref, max_value, i ):
    global wstar2

    # The effect of the surface current is considered in "delta_u":
    #     delta_u = U(zref) - Usfc
    # This slightly reduces the speed of the algorithm.
    wstar2 = ( pow( m.fabs( m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
        * bstar ) * H, 0.6666667) )
    if( monin_inv >= 0.0 ):
        wstar2 = 0.0
    hzzg = zz[i] # in this case hzzg is not related to wave height - it is a dummy variable
    if( hzzg < 0.00001 ):
        hzzg = 0.00001
    if( m.fabs( delta_u ) < 0.000001 ):
        ustar_test = 0.0
    elif( m.log( hzzg + 1.0 ) + float(not fixed_uen) * float(not neutral) * phi_u > 0.01 ):
        delta_u2 = delta_u * delta_u + CON_P * wstar2
        if( m.fabs( delta_u2 ) < 0.000001 ):
            delta_u2 = 0.000001
        delta_u = m.sqrt( delta_u2 ) * delta_u / m.fabs( delta_u )
        #  calculate displacement height, multiplied by zero or one depending on
        #    whether or not it is used.
        d_over_zo = float(use_dh) * hzzg / zref * hsig[0] * 0.8
        ustar_test = ( KV * delta_u / ( m.log( hzzg - d_over_zo + 1.0 ) +
            float(not fixed_uen) * float(not neutral) * phi_u ) )
        # eliminate the problem of excessive stress: limit friction velocity
        if ( pow( KV / m.log( ( 10.0 - d_over_zo ) * zz[i] / zref + 1.0), 2.0 ) > 0.01 ):
            ustar_test = ( 0.01 * pow( m.log( ( 10.0 - d_over_zo ) * zz[i] / zref + 1.0)
            / KV, 2.0 ) )
        if ( ustar_test > 10.0 ):
            ustar_test = 10.0
    else:
        ustar_test = 10.0
    # ustar is not returned because it is a global variable
    return ustar_test

# # # # # phi_u_f # # # # #
def phi_u_f( zref, zzd ):
    a_bh = 0.7
    b_bh = 0.75
    c_bh = 5.0
    d_bh = 0.35

    # unstable case
    if ( monin_inv < 0 ):
        if( stab_prm == 0 or stab_prm == 2):  #Businger-Dyer parameterization  
        # BVW
        # find zeta( z / L ) and zeta( zo / L )
            if zzd < 0.0001:
                zzd = 0.0001
            zeta = pow( 1.0 - 15.0 * zref * monin_inv, 0.25 )
            zeta0 = pow( 1.0 - 15.0 * zref / zzd * monin_inv, 0.25 )
#            phi_u = m.log(( (zeta0 * zeta0 + 1.0 ) * pow( zeta0 + 1.0, 2.0 ) /
#                         ( ( zeta * zeta + 1.0 ) * pow( zeta + 1.0, 2.0 ) ) ) +
#                         2.000000 * ( m.atan(zeta) - m.atan(zeta0) ))
            phi_u = - 2.*m.log((1.+zeta)/(1.+zeta0)) - m.log((1.+zeta**2.)/(1.+zeta0**2.)) + 2.0*m.atan(zeta) - 2.0*m.atan(zeta0)
        elif(stab_prm == 1):   # Smith 88 parameterization  */
            a_bh = m.sqrt( m.sqrt( 1.0 - 16.0 * zref * monin_inv ) )
            phi_u = (-2.0 * m.log( 0.5 * (1.0 + a_bh) ) -
                    m.log( 0.5 * (1.0 + a_bh*a_bh) ) + 2.0 *
                    m.atan( a_bh ) - 3.14159 / 2.0)
        elif(stab_prm == 3):  #COAWST
            zeta = pow( 1.000 - 16.000 * zref * monin_inv, 0.2500 )
            zeta0 = pow( 1.000 - 16.000 * zref / zzd * monin_inv, 0.2500 )
            phi_u = - 2.*m.log((1.+zeta)/(1.+zeta0)) - m.log((1.+zeta**2.)/(1.+zeta0**2.)) + 2.0*m.atan(zeta) - 2.0*m.atan(zeta0)
        elif(stab_prm == 4):  #Paulson 1970
            zeta = pow( 1.000 - 16.000 * zref * monin_inv, 0.2500 )
            phi_u = - 2.*m.log((1.+zeta)/2.) - m.log((1.+zeta**2.)/2.) + 2.0*m.atan(zeta) - 3.14159/2

    else:
    #  stable case  */
        if ( stab_prm == 0):
          # BVW
            phi_u =( a_bh * zref * monin_inv +
                b_bh * ( zref * monin_inv - c_bh / d_bh ) *
                m.exp(-d_bh * zref * monin_inv ) + b_bh * c_bh / d_bh)
        elif (stab_prm == 1 or stab_prm == 2 or stab_prm == 4):   # Hicks parameterization, also used in Smith 88
            phi_u = 5.0 * zref * monin_inv
        elif ( stab_prm == 3):
            phi_u = 5.0 * ( zref * monin_inv - zref / zzd * monin_inv )


    return phi_u

# # # # # phi_Q_f # # # # #
def phi_Q_f( zrefq, zzq ):
    #double a_bh, b_bh, c_bh, d_bh, phi_Q, qlamda, q0lamda;
    a_bh = 1.0
    b_bh = 0.667
    c_bh = 5.0
    d_bh = 0.35

    if ( monin_inv < 0 ):
        if( stab_prm == 0 or stab_prm ==2 ):
        # BVW
            # find zeta( z / L ) and zeta( zo / L )
            qlamda = m.sqrt( 1.000 - 9.000 * zrefq * monin_inv )
            q0lamda = m.sqrt( 1.000 - 9.000 * zrefq / zzq * monin_inv )
            phi_Q = 2.0000 * m.log( ( q0lamda + 1 ) / ( qlamda + 1 ) )
        elif( stab_prm == 1):   # Smith 88 parameterization
            a_bh = m.sqrt( 1.0 - 16.0 * zrefq * monin_inv )
            phi_Q = -2.0 * m.log( ( 1.0 + a_bh ) / 2.0 )
        elif( stab_prm == 3 ):
            qlamda = pow( 1.000 - 16.000 * zrefq * monin_inv, 0.500)
            q0lamda = pow( 1.000 - 16.000 * zrefq / zzq * monin_inv, 0.500)
            phi_Q = -2.0 * m.log( (1.0 + qlamda)/(1.0+q0lamda))
    else:
    #  stable case
        if (stab_prm == 0):  # BVW
             phi_Q = (pow( 1.0 + 0.6667 * a_bh * zrefq * monin_inv, 1.5 ) +
                 b_bh * ( zrefq * monin_inv - c_bh / d_bh ) *
                 m.exp(-d_bh * zrefq * monin_inv ) + b_bh * c_bh / d_bh - 1.0)
        elif(stab_prm == 1 or stab_prm == 2 ):   # Smith 88 parameterization
            phi_Q = 5.0 * zrefq * monin_inv
        elif(stab_prm ==3):
            phi_Q = 5.0 * ( zrefq * monin_inv - zrefq / zzq * monin_inv )
    return phi_Q

# # # # # phi_T_f # # # # #
def phi_T_f( zreft, zzt ):
    a_bh = 1.0
    b_bh = 0.667
    c_bh = 5.0
    d_bh = 0.35

    if ( monin_inv < 0 ):
        if (stab_prm == 0 or stab_prm == 2):  # BVW
        # find zeta( z / L ) and zeta( zo / L )  */
            tlamda = pow( 1.000 - 9.000 * zreft * monin_inv, 0.500)
            t0lamda = pow( 1.000 - 9.000 * zreft / zzt * monin_inv, 0.500)
            phi_T = 2.0000 * m.log( ( t0lamda + 1 ) / ( tlamda + 1 ) )
        elif(stab_prm == 1):   # Smith 88 parameterization
            a_bh = m.sqrt( 1.0 - 16.0 * zreft * monin_inv )
            phi_T = -2.0 * m.log( ( 1.0 + a_bh ) / 2.0 )
        elif( stab_prm == 3 ):
            tlamda = pow( 1.000 - 16.000 * zreft * monin_inv, 0.500)
            t0lamda = pow( 1.000 - 16.000 * zreft / zzt * monin_inv, 0.500)
            phi_T = -2.0 * m.log( (1.0 + tlamda)/(1.0+t0lamda) )
    else:
    #  stable case
        if(stab_prm == 0):  # BVW
            phi_T = (pow( 1.0 + 0.6667 * a_bh * zreft * monin_inv, 1.5 ) +
                    b_bh * ( zreft * monin_inv - c_bh / d_bh ) *
                    m.exp(-d_bh * zreft * monin_inv ) + b_bh * c_bh / d_bh - 1.0)
        elif(stab_prm == 1 or stab_prm == 2 ): # Smith 88 parameterization
            phi_T = 5.0 * zreft * monin_inv
        elif(stab_prm==3):
            phi_T = 5.0 * ( zreft * monin_inv - zreft / zzt * monin_inv )
    return phi_T

# # # # # z_o_zq # # # # #
def z_o_zq( zref, zrefq, nu, zo_mag, Dalton ):
    global zzq

    #  zzq is  zref / z0Q
    ustar_mag = m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
    if ( TQzo_prm == 0 ): # wall theory
        if ( betas[0] == 0.0 ):
            betas_temp = 1.0
        else:
            betas_temp = betas[0]
        zzq = zrefq / ( betas_temp * 0.620 * nu / m.fabs( ustar_mag ) )
    elif ( TQzo_prm == 1 ): # CFC model
        zzq = zrefq / ( zo_mag * m.exp( KV * ( 5.0 - 1.0 / ( Sc * Dalton ) ) ) )
    elif ( TQzo_prm == 2 ): # Zilitinkevich et al. 2001
        zzq = zrefq / ( zo_mag * m.exp( -4.0 * m.sqrt( zo_mag * ustar_mag / nu ) ) )
    elif ( TQzo_prm == 3 ): # LKB 1979
        Rr = zo_mag * ustar_mag / nu
        if ( Rr < 0.11 ):
            const1 = 0.177
            const2 = 0
        elif ( Rr < 0.825 ):
            const1 = 1.376
            const2 = 0.929
        elif ( Rr < 3. ):
            const1 = 1.026
            const2 = -0.599
        elif ( Rr < 10.0 ):
            const1 = 1.625
            const2 = -1.018
        elif ( Rr < 30.0 ):
            const1 = 4.661
            const2 = -1.475
        elif ( Rr < 100.0 ):
            const1 = 34.904
            const2 = -2.067
        else:
            print( "LKB zoq calculation failed due to large Rr: ", Rr)
        # print( "LKB z0 = %f %f %f; LKB ustar = %f; LKB nu = %f\n", zo_mag, zz[0], zz[1], ustar_mag, nu); */
            const1 = 34.904
            const2 = -2.067
        zzq = zrefq * ustar_mag * const1 * pow( Rr, const2 ) / nu
    elif ( TQzo_prm == 4 ): # COARE3.0
        zzq = 0.000055 * pow(  zo_mag * ustar_mag / nu, -0.6 ) # zo_q */
        if ( zzq < 0.00011 ):
            zzq = 0.00011
        zzq = zrefq / zzq
    elif ( TQzo_prm == 5 ): # Modified CFC model
        zzq = zrefq / ( zo_mag * m.exp( KV * ( 0.5 - 1.0 / ( Sc * Dalton ) ) ) )
    else:
        print( "Invalid choice of TQzo_prm: ", TQzo_prm )
        exit(0)
    return zzq

# # # # # z_o_zt # # # # #
def z_o_zt( zref, zreft, nu, zo_mag, Stanton ):
    global zzt

    # double betas_temp, const1, const2, Rr, ustar_mag;
    #print("test:", ustar)
    ustar_mag = m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
    #print("test2:", ustar_mag)
    if ( TQzo_prm == 0 ):
        if ( betas[0] == 0.0 ):
            betas_temp = 1.0
        else:
            betas_temp = betas[0]
        # determine zreft / zT
        zzt = ( zreft / ( betas_temp * 0.40 * nu / m.fabs( ustar_mag ) ) )
    elif ( TQzo_prm == 1 ):
        zzt = ( zreft / ( zo_mag * m.exp( KV * ( 5.0 - 1.0 / ( Prt*Stanton ) ) ) ) )
    elif ( TQzo_prm == 2 ):
        # Zilitinkevich et al. 2001
        zzt = ( zreft / ( zo_mag * m.exp( -4.0 * m.sqrt( zo_mag * ustar_mag / nu ) ) ) )
    elif ( TQzo_prm == 3 ):
      # LKB 1979
        Rr = ( zo_mag * ustar_mag / nu )
        if ( Rr < 0.11 ):
            const1 = 0.292
            const2 = 0
        elif ( Rr < 0.825 ):
            const1 = 1.808
            const2 = 0.826
        elif ( Rr < 3. ):
            const1 = 1.393
            const2 = -0.528
        elif ( Rr < 10.0 ):
            const1 = 1.956
            const2 = -0.870
        elif ( Rr < 30.0 ):
            const1 = 4.994
            const2 = -1.297
        elif ( Rr < 100.0 ):
            const1 = 30.790
            const2 = -1.845
        else:
            print( "LKB zoq calculation failed due to large Rr: \n", Rr )
            const1 = 30.790
            const2 = -1.845
        zzt = ( zreft * ustar_mag * const1 * pow( Rr, const2 ) / nu )
    elif ( TQzo_prm == 4 ):
        # COARE3.0
        zzt = ( 0.000055 * pow(  zo_mag * ustar_mag / nu, -0.6 ) ) # zo t
        if ( zzt < 0.00011 ):
            zzt = 0.00011
        zzt = zreft / zzt
    elif ( TQzo_prm == 5 ): # Modified CFC model
        zzq = zreft / ( zo_mag * m.exp( KV * ( 0.5 - 1.0 / ( Prt * Stanton ) ) ) )
    else:
        print( "Invalid choice of TQzo_prm: \n", TQzo_prm )
        exit(0)
    return zzt

# # # # # z0_and_waves # # # # #
def z0_and_waves( fixedcp, fixedwa, neutral, hsig_known, denwat, nu, sfcten, zref,
    cp, ww, ww_eql, ss_prm, ss_val, count, i, z0_prm, A_oil ):
    global cmin, cp_com, hsig, betag, wa, zz, betas, betac
    global wstar2

    #z0, period_known,  B_Toba, dum, hzzg, p1, p2, p3, p4, oil_mod, phi_u, scale, sig_slope, ueff, Uorb, ustar_temp, wave_len = [None]*17, betac_oil, betac_water
    period = [0,0]
    A_Charnock = 0.035
    oil_mod = 0.25 # crude guess for oil spills - not natural oil
    # A_Charnock = 0.48/ss_val
    test = 0
    ustar_temp = ( m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] ) )
    cmin = m.sqrt( 2.0 * m.sqrt( G * sfcten / denwat ) ) 
    ueff = m.sqrt( ( u[i] - u_sfc[i] + wstar2 ) * ( u[i] - u_sfc[i] + wstar2 ) )

    z0 = zref / zz[0]
    # For the flux models (i.e., only BVW) that calculate sea state
    if ( z0_prm == 0 or z0_prm == 3 ):
        B_Toba = 0.0602
        # the sea state parameter is significant slope
        if ( ss_prm == 4 and i == 0 ):
            sig_slope = ss_val
            cp = ( 2.0 * 3.14159 * m.fabs( ustar[i] ) * B_Toba * B_Toba /
                sig_slope / sig_slope )
            wa[i] = ( cp / m.fabs( ustar[i] ) )
            if ( wa[i] > 12000.0 ):
                wa[i] = 12000.0
            # printf( "%lf %lf %lf\n", cp, wa[i], ustar[i] )
            ueff = ( m.sqrt( u[i] * u[i] + wstar2 ) )
            z0 = ( betas[i] * 0.1100 * nu / m.fabs( ustar[i] ) +
                betac[i] * B * sfcten / ( ustar[i] * ustar[i] * denwat ) +
                betag[i] * 0.016 * ustar[i] * ustar[i] / G )
            phi_u = phi_u_f( betag[i] * hsig[i], hsig[i] / z0 )
            ww = ( m.log( hsig[i] / z0 + 1.0 ) + phi_u ) / ( KV * ww_eql * wa[i] )
            # determine the wave characteristics
            period[i] = ( pow( pow( cp * sig_slope / B_Toba, 2.0 ) / ( G * ustar[i] ),
                0.333333 ) )
            hsig[i] = B_Toba * m.sqrt( G * m.fabs( ustar[i] ) * pow( period[i], 3 ) )
            wave_len = ( cp * period[i] )
        # the sea state parameter is significant wave height OR period
        if ( ( hsig_known or ss_prm == 5 ) and i == 0 ):
            if ( ss_prm == 5 ):
                period_known = TRUE
            else:
                period_known = FALSE

            period[i] = ( float(not period_known) * pow( hsig[i] * hsig[i] / ( G *
                m.fabs( ustar[i] ) * B_Toba * B_Toba ), 0.333333 ) +
                float(period_known) * ss_val )
            hsig[i] = ( float(period_known) * B_Toba * m.sqrt( G * m.fabs( ustar[i] ) *
                pow( period[i], 3 ) ) + float(not period_known) * ss_val )
            p1 = ( G * period[i] / ( 2.0 * 3.14159 ) )
            p2 = pow( cmin, 4.0)
            p3 = ( 0.5 * p2 / p1 + pow( p1, 3.0 ) / 27.0 )
            p4 = ( 0.25 * p2 * p2 / ( p1 * p1 ) + p1 * p1 * p2 / 27.0 )
            cp = ( pow( p3 + m.sqrt( p4 ), 0.333333 ) +
                pow( p3 - m.sqrt( p4 ), 0.333333 ) + p1 / 3.0 )
            wa[i] = ( cp / m.fabs( ustar[i] ) )

            if ( wa[i] > 12000.0 ):
                wa[i] = 12000.0
            # determine the wave characteristics
            wave_len = ( cp * period[i] )
            ueff = m.sqrt( u[i] * u[i] + wstar2 )

            z_over_z0( denwat, nu, sfcten, zref, i, z0_prm, A_oil )

            phi_u = phi_u_f( betag[i] * hsig[i], hsig[i] / z0 )
            ww = ( m.log( hsig[i] / z0 + 1.0 ) + phi_u ) / ( KV * ww_eql * wa[i] )
            # print( "z/z0 T: i = %i; ss_val=%f; period = %f; hsig = %f; wa = %f\n",
                # i, ss_val, period[i], hsig[i], wa[i] )
        else:
            # print( 'test option other' )
            z_over_z0( denwat, nu, sfcten, zref, i, z0_prm, A_oil )

            # determine the wave characteristics
            #print( "betag = ", i, betag[i] )
            if ( betag[i] > 0.1 ):
                dum = ( ustar[i] * ustar[i] * wa[i] * wa[i] )
                #print( "wa test, other: ", i, wa[i], wave_len, period[i], ustar[i], m.sqrt(dum), cmin )
                if ( dum > cmin * cmin ):
                    wave_len = ( 3.14159 * ( dum + ( m.sqrt( dum * dum -
                     float(not no_sfcten) * pow(cmin,4) ) ) ) / G )
                    period[i] = ( pow( hsig[i] * hsig[i] / ( G *
                     m.fabs( ustar[i] ) * B_Toba * B_Toba ), 0.333333 ) )
                else:
                    wave_len = ( PI2 * m.sqrt( sfcten / ( G * denwat ) ) )
                    #hsig[i] = ( float(period_known) * B_Toba * m.sqrt( G * m.fabs( ustar[i] ) *
                    # pow( period[i], 3 ) ) + float(not period_known) * ss_val )

                    period[i] = ( wave_len / ( m.fabs( wa[i] * ustar[i] ) ) )
                    # Toba's relation between wave height, friction velocity and period
                    hsig[i] = ( 0.062 * m.sqrt( G * m.fabs( ustar[i] ) * pow( period[i], 3 ) ) )
                    # limit the height of breaking waves
                    if ( hsig[i] > 0.137 * wave_len ):
                        hsig[i] = ( 0.137 * wave_len )
            else:
                hsig[i] = 0.0
                period[i] = -9999.9
            # printf( "z/z0 Hs: i=%i; hsig=%f; ustar=%f; wa=%f; period=%f\n",
                  # i, hsig[i], ustar[i], wa[i], period[i] );
        # end of block estimating seastate for BVW
        # determine beta-c
    betac[i] = betac_prime
    if ( z0_prm != 3 ):
        scale = 0.8
        hzzg = ( scale * hsig[i] * zz[i] / zref ) # can't be -ve
        phi_u = phi_u_f( betag[i] * scale * hsig[i], hzzg )
        if ( ss_prm == 6 ):
            Uorb = u_orb[i]
        else:
            Uorb = 0.0
        phi_u = phi_u_f( 10.0, hzzg )
        ueff = ( ustar[i] * ( m.log( zz[i] + 1.0 ) + phi_u ) / KV )
        if ( z0_prm == 0 ):
            if ( ueff < 1.8 ):
                betac[i] = 0.0
                betag[i] = 0.0
            else:
                betac[i] =  1.0
                betag[i] = 1.0
        if ( z0_prm == 1 ):
            if ( ueff < 1.0  ):
                betac[i] = 0.0
                betag[i] = 0.0
            else:
                betac[i] = ( m.tanh( pow( 0.4 * ( ueff -1.0 ), 3 ) ) * betac_prime )
                betag[i] = ( m.tanh( pow( 0.3 * ( ueff -1.0 ), 3 ) ) * betag_prime )
        if ( z0_prm == 4 ):
            if ( ueff < 7.0 ):
                betac[i] = 0.0
            else:
                betac[i] = ( m.tanh( pow( 0.4 * ( ueff -7.0 ), 3 ) ) * betac_prime )
                betag[i] = ( m.tanh( pow( 0.3 * ( ueff -7.0 ), 3 ) ) * betag_prime )
        if ( betac[i] < 0.1 and (not no_capw) ):
            betag[i] = betac[i]
        if ( no_capw ):
            betac[i] = 0.0
        betas[i] = ( 1.0 - betac[i] )
        # printf( "betac: hzzg=%f; phi_u=%f; Uorb=%f\n", hzzg, phi_u, Uorb);
        # determine wave age
        phi_u = phi_u_f( betag[i] * hsig[i], hsig[i] * zz[i] / zref )
        hzzg = ( hsig[i] * zz[i] / zref ) # can't be -ve
        if ( hzzg < 0.00001 ):
            hzzg = 0.00001
        if ( ss_prm == 0 or i == 1 ):
            # local equilibrium requires ww = 1
            wa[i] = ( ( m.log( hzzg + 1.0 ) + float(not neutral) * phi_u )
                / ( float( (1.0-i) * ww + float(i) ) * KV * ww_eql ) )
            # if ( wa[i] < 10.0 ): 
            #  wa[i] = 10.0
            cp = ( wa[i] * m.fabs( ustar[i] ) )
            # printf( "test1: %i %f %i %f %f %f %f %f\n", i, hzzg, neutral,
                # phi_u, KV, ww, ww_eql, wa[i] )
            # print( "cp wa test: ", cp, wa[i] )
    # printf( "test2: %i %f %f %f %f\n", i, ueff, betas[i], betac[i], betag[i] )
    if ( not fixedwa ):
        # i must equal 0 to reach this point
        # wave age based on sea state observations (other than wave age)
        # print( " wa test: ",  i, m.fabs( ustar[i] ), cp,  wa[i] )
        if ( m.fabs( ustar[i] ) > 0.0 ):
            wa[i] = ( cp / m.fabs( ustar[i] ) )
        else:
            wa[i] = 999.0
    else:
        cp = ( wa[i] * m.fabs( ustar[i] ) )
        # printf( "z0_and_waves: i=%i; u=%f; wa=%f; ustar=%f; hsig=%f; hzzg=%f;
            # phi_u=%f; cp=%f; bc=%f\n", i, u[i], wa[i], ustar[i], hsig[i],
            # hzzg, phi_u, cp, betac[i] )
        # insure that wa is not horribly unreasonable
    if ( not( fixedwa or hsig_known or ss_prm == 4 ) and wa[i] > 250.0 ):
        wa[i] = 250.0
        test = -1
    else:
        test = ( test - test%2 )
    if ( (not fixedwa) and wa[i] < 1.08 ):
        wa[i] = 1.08
        test = test - 2
    else:
        test = ( test - ( test%4 - test%2 ) )
    if ( not fixedcp ):
        cp_com = ( wa[0] * m.fabs( ustar[0] ) )
    else:
        cp_com = cp
    # printf( "betac & betag: %f %f %f %f\n", betac[0], betag[1], betac[1], betag[1] )
    # printf( "betac & betag: %f %f %f %f\n", betac[i], betag[i], betac_prime, betag_prime );
 
    # determine zref / z0
    z_over_z0( denwat, nu, sfcten, zref, i, z0_prm, A_oil )     
    # insure that the roughness length is not too absurd */
    if ( zref > zz[i] * 25 ):
        zz[i] = zref / 25.0
        test = test - 4
    else:
        test = test - ( test%8 - test%4 - test%2 )
    # printf( "z0_and_waves: i=%i; zz=%f; u*=%f; wa=%f; hs=%f\n", i, zz[i], ustar[i],
        # wa[i], hsig[i] ); */
    return test   # returns an error code


def z_over_z0( denwat, nu, sfcten, zref, i, z0_prm, A_oil ):
    global hsig, betag, wa, zz, betas, betac
    global wstar2

    A_Charnock1 = 0.035
    A_Charnock2 = 0.018
    oil_mod = 0.25 # crude guess for oil spills - not natural oil
    ustar_temp = ( m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] ) )

    if ( z0_prm == 0 ):
        # BVW
        zz[i] = ( zref / ( ( ( betas[i] * 0.1100 * nu / m.fabs( ustar[i] ) ) +
            m.sqrt( pow( betac[i] * B * sfcten / ( ustar_temp * m.fabs( ustar[i] ) *
            denwat ), 2.0 ) + pow( betag[i] * 0.48 / wa[i] * ustar_temp *
            m.fabs( ustar[i] ) / G, 2.0 ) ) ) ) )
        # * exp( -KV * Uorb / fabs( ustar[i] )
    elif ( z0_prm == 1 ):
        # Bourassa (2006)
        zz[i] = ( zref / ( ( ( betas[i] * 0.1100 * nu / m.fabs( ustar[i] ) ) +
            m.sqrt( pow( betac[i] * B * sfcten / ( ustar_temp * m.fabs( ustar[i] ) *
            denwat ), 2.0 ) + pow( betag[i] * A_Charnock2 * ustar_temp *
            m.fabs( ustar[i] ) / G, 2.0 ) ) ) ) )
    elif ( z0_prm == 2 ):
        # BVW with Taylor & Yelland zog
        zz[i] = ( zref / ( ( betas[i] * 0.1100 * nu / m.fabs( ustar[i] ) ) +
            m.sqrt( pow( betac[i] * B * sfcten / ( ustar_temp * m.fabs( ustar[i] ) *
            denwat ), 2.0 ) + pow( betag[i] *
            hsig[i] * 1200.0 * pow( hsig[i] / wave_len, 4.5 ), 2.0 ) ) ) )
    elif ( z0_prm == 3 ):
        # Taylor & Yelland
        zz[i] = ( zref / ( betag[i] *
            hsig[i] * 1200.0 * pow( hsig[i] / wave_len, 4.5 ) ) )
    elif ( z0_prm == 4 ):
        # Bourassa (2006) modified for oil slick
        zz[i] = ( zref / ( ( ( betas[i] * 0.1100 * nu / m.fabs( ustar[i] ) ) +
            A_oil * oil_mod * m.sqrt( pow( betac[i] * B * sfcten / ( ustar_temp *
            m.fabs( ustar[i] ) * denwat ), 2.0 ) + pow( betag[i] * A_Charnock1 *
            ustar_temp * m.fabs( ustar[i] ) / G, 2.0 ) ) ) ) )
    elif ( z0_prm == 5 ):
        # Smooth surface
        zz[i] = ( zref / ( 0.1100 * nu / m.fabs( ustar[i] ) ) )
    elif ( z0_prm == 6 ):
        # 2020 version of the oil slick code
        phi_u = phi_u_f( zref, zz[i] )
        ueff = ustar[0] * ( m.log( zz[0] + 1) + phi_u ) / KV
#        ueff = m.sqrt( u[0]**2 + u[1]**2 )
        if (ueff<1.0):
            betac[i] = 0.0
            betag[i] = 0.0
        else:
            betac[i] = m.tanh( (0.4 * ( ueff - 1.0 ))**3 )
            betag[i] = m.tanh( (0.2 * ( ueff - 1.0 ))**3 )
        if (ueff<7.0):
            betac_oil = 0.0
            betag_oil = 0.0
        else:
            betac_oil = m.tanh( (0.4 * ( ueff - 7.0 ))**3 )
            betag_oil = m.tanh( (0.3 * ( ueff - 7.0 ))**3 )
        weighted_betac = (1.0-A_oil) * betac[i] + A_oil * betac_oil * oil_mod
        weighted_betag = (1.0-A_oil) * betag[i] + A_oil * betag_oil * oil_mod
        betas[i] = 1.0-weighted_betag
        zz[i] = ( zref / ( ( betas[i] * 0.1100 * nu / m.fabs( ustar[i] ) ) +
            m.sqrt( pow( weighted_betac * B * sfcten / ( ustar[i] *
            ustar[i] * denwat ), 2.0 ) + pow( weighted_betag *
            A_Charnock2 * ustar_temp * m.fabs( ustar[i] ) / G, 2.0 ) ) ) )
    elif ( z0_prm == 7 ): 
        # do not recalculate. Use the value that was initialized in the main program, based on
        # a value input to the program. The next line is not used. 
        ueff = m.sqrt( u[0]**2 + u[1]**2 )
    else:
        print( "Unreasonable value of z0_prm: \n", z0_prm )

#    print( "z_over_z0, zz = ", zz[i] )
    return( z0_prm )

# # # # find_ustar # # # #
def find_ustar( first, fixedcp, fixedwa, fixed_hsig, fixed_uen,
    neutral, count, CON_P, cp, CRIT, denwat, nu, qmixa,
    sfcten, ta, ww, ww_eql, ss_prm, ss_val, zref, zreft, z0_prm, A_oil, astab ):

    global monin_inv, ustar

    ustar_old = [None]*2
    ustar_temp = m.sqrt( ustar[0] * ustar[0] + ustar[1] * ustar[1] )
    ucount = 0
    # iteratively solve for u* for non-neutral conditions
    # determine the bouyancy flux
    done_outer = FALSE
    while ( not done_outer and ucount < 30 ):
        ucount = ucount + 1
        i = 0
        while( i < 2 ):
            ustar_old[i] = ustar[i]
            # ustar_temp is ustar squared
            ustar_temp = ustar[0] * ustar[0] + ustar[1] * ustar[1]
            if ( fixed_hsig and i == 0 ):
                hsig_known = TRUE
            else:
                hsig_known = FALSE
            bstar = f_bstar( neutral, qmixa, ta, 0.000001 )

            # !!!!  This assumes temperature and humidity are measured at THE SAME HEIGHT
            # determine the inverse Monin-Obhukov scale length: 1/L
            if ( astab != 2 ): 
               monin_inv = KV * bstar / ustar_temp
            # insure that 1/L is within reasonable bounds
            max = 8.0
            min = 0.0001
            if ( m.fabs(monin_inv) < min ):
                monin_inv = min * monin_inv / m.fabs( monin_inv )
            if ( m.fabs( monin_inv ) > max ):
                monin_inv = max * monin_inv / m.fabs( monin_inv )
            if ( m.fabs( u[i] - u_sfc[i] ) < 0.001 ):
                 ustar[i] = 0.0
            else:
                if ( not first ):
                    done = FALSE
                first = FALSE
                test = z0_and_waves( fixedcp, fixedwa, neutral, hsig_known,
                    denwat, nu, sfcten, zref, cp, ww, ww_eql, ss_prm,
                    ss_val, count, i, z0_prm, A_oil )
                phi_u = phi_u_f( zref, zz[i] )
                # printf( "A %f %f %f %f\n", (float)ustar[0], (float)ustar[1],
                # (float)tstar, (float)qstar )
                ustar[i] = f_ustar( fixed_uen, neutral, bstar, CON_P, phi_u,
                    u[i] - u_sfc[i], ww, zref, 10.0, i )
            i += 1
        if ( ustar[0] <= 0.000001 ):
            done_outer = m.fabs( ustar[0] - ustar_old[0] ) < 0.0001
        else:
            done_outer = m.fabs( ( ustar[0] - ustar_old[0] )/ ustar[0] ) < 0.01
        if ( ustar[1] <= 0.000001 ):
            done_outer = done_outer and m.fabs( ustar[1] - ustar_old[1] ) < 0.0001
        else:
            done_outer = ( done_outer and m.fabs( ( ustar[1] - ustar_old[1] ) /
                ustar[1] ) < 0.01)

# # # getRenewalTime  # # #
def getRenewalTime( zo_mag, ustar_mag, density, cp, alpha, Qnet, viscosity ):
    Rfcr = -2.0e-04
    if(density < 50):
        Rfcr = -2.0e-01
    if(density < 50):
        Cshear = 53.32
    else:
        Cshear = 209 # from Wick thesis
    if(density < 50):
        Cconv = 0.80
    else:
        Cconv = 3.13 # from Wick thesis

    # get surface Richardson number

    Rfo = alpha * G * Qnet * viscosity / (density * cp * pow(ustar_mag, 4.0))
    if(Rfo == m.fabs(Rfo)):
        Rfo *= -1.0
    # now get time rate **/
    #print("DB5 ",ustar_mag,viscosity,zo_mag)
    shearContribution = Cshear * m.sqrt(viscosity * zo_mag / pow(ustar_mag, 3.0))
    convContribution = Cconv * m.sqrt(viscosity * density * cp / (alpha * G * m.fabs(Qnet)))
    renewalTime = (convContribution - shearContribution) * m.exp(-Rfcr / Rfo)
    renewalTime += shearContribution
    return(renewalTime)

# # # # # getSolar # # # # #
def getSolar( solar ):
    # this just determines the amount of solar radiation is
    #     deposited in the skin. Aidong: this isn't needed in new
    #     surface skin model

    if(waterType == 0):
        return(0.0)
    elif(waterType == 1):
        R2 = 0.58
        gamma1 = 0.35
        gamma2 = 23
    elif(waterType == 1.1):
        R2 = 0.62
        gamma1 = 0.60
        gamma2 = 20
    elif(waterType == 1.2):
        R2 = 0.67
        gamma1 = 1.00
        gamma2 = 17
    elif(waterType == 2):
        R2 = 0.77
        gamma1 = 1.50
        gamma2 = 14
    elif(waterType == 3):
        R2 = 0.78
        gamma1 = 1.40
        gamma2 = 7.9

    z = 0.001
    solarAbs = solar - solar*(R2*m.exp(-z/gamma1) + (1-R2)*m.exp(-z/gamma2))
    return(solarAbs)

# # # # # getStanton # # # # #
def getStanton( ustar_mag, renewalTime):
    Stanton = m.sqrt(thermalDiff/(ustar_mag * ustar_mag * renewalTime))
    return( Stanton )

# # # # # getDalton # # # # #
def getDalton( ustar_mag, renewalTime ):
    Dalton = m.sqrt(molecularDiff/(ustar_mag*ustar_mag*renewalTime))
    return(Dalton)

# # # # # getSkin # # # # #
def getSkin( zo_mag, ustar_mag, bulkTemp, rl, Hs, Hl, rs):
    roWater = 1021.7
    CpWater = 4002.0
    roAir = 1.25
    kappawater = 1.4e-07 # from Gill page 71
    thermalExpansion = 3196e-07 # from Gill Table A3.1
    viscosity = 1e-06 # this is kinematic viscosity (Gill pg. 75)
    # get solar radiation absorbed over skin
    Qsolar = getSolar(rs)
    # Get Net Heat Loss from ocean
    Qnet = rl + Hs + Hl + Qsolar
    # get u* for water
    uStarWater = m.sqrt(roAir/roWater)*ustar_mag
    # now get skin temperature
    renewalTime = getRenewalTime(zo_mag,uStarWater, roWater, CpWater,
        thermalExpansion, Qnet, viscosity)

    skinTemperature = m.sqrt(renewalTime)
    skinTemperature *= Qnet/(roWater*CpWater*m.sqrt(kappawater))
    skinTemperature = skinTemperature + bulkTemp

    return(skinTemperature)

