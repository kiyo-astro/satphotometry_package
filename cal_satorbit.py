#--------------------------------------------------------------------------------------------------#
# obs_satorbit.py                                                                                  #
# Developed by Kiyoaki Okudaira * University of Washington / Kyushu University / IAU CPS SatHub    #
#--------------------------------------------------------------------------------------------------#
# Description                                                                                      #
#--------------------------------------------------------------------------------------------------#
# Calculate satellite orbit from TLE with SPICE toolkit                                            #
#--------------------------------------------------------------------------------------------------#
# History                                                                                          #
#--------------------------------------------------------------------------------------------------#
# coding 2025.12.07: 1st coding                                                                    #
#--------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------#
# VERSION                                                                                       #
#-----------------------------------------------------------------------------------------------#
version = "1.0.0"
version_number = 2025120710000

#--------------------------------------------------------------------------------------------------#
# Function                                                                                         #
#--------------------------------------------------------------------------------------------------#
def cal_satorbit(obs_begin,obs_end,obs_step,tle_PATH,obs_params,spice_myfile_PATH):
    # Observatory
    obs_gd_lon_deg = obs_params[0]
    obs_gd_lat_deg = obs_params[1]
    obs_gd_height  = obs_params[2]

    wavelength_m = obs_params[3]
    aperture_m   = obs_params[4]

    #--------------------------------------------------------------------------------------------------#
    # Libraries                                                                                        #
    #--------------------------------------------------------------------------------------------------#
    # Import spiceypy and start up SPICE kernel
    import spiceypy as spice
    spice.furnsh(spice_myfile_PATH)

    # Standard libraries
    import math
    import astropy.units as u
    from astropy.coordinates import Angle
    from astropy.table import Table

    # Satphotometry library
    # satphotometry.satorbit must be imported after starting up SPICE kernel
    from satphotometry import satorbit

    #--------------------------------------------------------------------------------------------------#
    # Constants                                                                                        #
    #--------------------------------------------------------------------------------------------------#
    # Observatory
    obs_gd_lon = math.radians(obs_gd_lon_deg) # [radian]
    obs_gd_lat = math.radians(obs_gd_lat_deg) # [radian]

    # Planetary constants
    # Get the Earth's constants that are needed for SGP4 propagation<br>
    # NAIF integer ID code of the Earth is 399
    earth_constants = satorbit.get_planetconst(399,["J2", "J3", "J4", "KE", "QO", "SO", "ER", "AE"])

    #--------------------------------------------------------------------------------------------------#
    # Main                                                                                             #
    #--------------------------------------------------------------------------------------------------#
    # Observation date and time
    obs_begin_et = spice.utc2et(obs_begin)
    obs_end_et   = spice.utc2et(obs_end)
    obs_ets      = [obs_begin_et+i*obs_step for i in range(0,int((obs_end_et-obs_begin_et)//obs_step+1))]

    # Output list
    output = []

    # Read TLE file
    if type(tle_PATH) is str:
        satname,line1,line2 = satorbit.read_TLEfile(tle_PATH)
    else:
        satname = tle_PATH[0]
        line1 = tle_PATH[1]
        line2 = tle_PATH[2]
    epoch,elems = satorbit.parse_TLE2element(line1,line2)

    epoch_utc = spice.et2utc(epoch, "ISOC", 0)

    if satname is not None:
        if satname[0:2] == "0 ":
            objname = satname[2:]
        objname = objname.rstrip()
    else:
        objname = None
        
    norad_id = str(line1[2:7].strip())
    intldes = line1.split()[2]

    if len(intldes)>8:
        intldes = "0000-000A"
    else:
        if int(intldes[:2]) >= 57:
            intldes = "19" + intldes[:2] + "-" + intldes[2:]
        else:
            intldes = "20" + intldes[:2] + "-" + intldes[2:]

    # SGP4 Propagation
    for et in obs_ets:
        et_utc = spice.et2utc(et, "ISOC", 0)

        state_teme = spice.evsgp4(et, earth_constants, elems)

        obs_itrf = satorbit.geo2itrf(obs_gd_lon,obs_gd_lat,obs_gd_height)
        obs_j2000 = satorbit.itrf2J2000(obs_itrf,et)

        state_j2000 = satorbit.teme2J2000(state_teme,et)

        range_km,ra,dec = satorbit.J20002radec(state_j2000[0:3],obs_j2000)

        ra_hms = Angle(ra*u.rad).to_string(unit=u.hourangle, sep=':', precision=1)
        ra_hms = "+0"+ra_hms if len(ra_hms)<10 else "+"+ra_hms

        dec_dms = Angle(dec*u.rad).to_string(unit=u.deg, sep=':', precision=1)
        dec_dms = "+0"+dec_dms if (dec>0 and len(dec_dms)<10) else ("+"+dec_dms if dec>0 else ("-0"+dec_dms[1:] if len(dec_dms)<11 else dec_dms))
        
        ra_deg = math.degrees(ra)
        dec_deg = math.degrees(dec)

        state_itrf = satorbit.J20002itrf(state_j2000,et)
        _,az,el = satorbit.itrf2azel(state_itrf[0:3],obs_itrf,obs_gd_lon,obs_gd_lat)

        az_deg = math.degrees(az)
        el_deg = math.degrees(el)

        in_umbra = satorbit.check_umbra(state_j2000[0:3],et)
        phase = satorbit.phase_angle(state_j2000[0:3],obs_j2000,et)
        phase_deg = math.degrees(phase)

        apparent_v_km_s = satorbit.apparent_v(state_j2000,obs_itrf,et)
        res_km = range_km * (wavelength_m / aperture_m)
        ex_ms = res_km / apparent_v_km_s

        output.append(
            [
                et_utc,
                in_umbra,
                round(range_km,4),
                ra_hms,
                dec_dms,
                round(ra_deg,4),
                round(dec_deg,4),
                round(phase_deg,4),
                round(az_deg,4),
                round(el_deg,4),
                round(apparent_v_km_s,4),
                round(res_km*1000,4),
                round(1/ex_ms,4)
                ]
            )

    # SPICE toolkit kernel clear
    spice.kclear()

    # Output
    # Save output file
    output = Table(
        list(zip(*output)),
        names=(
            "YYYY-MM-DDThh:mm:ss",
            "umbra",
            "range[km]",
            "ra[hh:mm:ss.s]",
            "dec[dd:mm:ss.s]",
            "ra[deg]",
            "dec[deg]",
            "pha[deg]",
            "az[deg]",
            "el[deg]",
            "v[km/s]",
            "res[m]",
            "ex[1/s]"
            )
        )
    
    return objname,norad_id,intldes,output

#--------------------------------------------------------------------------------------------------#
# Test                                                                                             #
#--------------------------------------------------------------------------------------------------#