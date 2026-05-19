#--------------------------------------------------------------------------------------------------#
# cal_satorbit.py                                                                                  #
# Developed by Kiyoaki Okudaira                                                                    #
#--------------------------------------------------------------------------------------------------#
# Description                                                                                      #
#--------------------------------------------------------------------------------------------------#
# Fast satellite orbit calculation from TLE with SPICE toolkit.                                    #
# Common observer/time calculations are cached in a context.                                       #
#--------------------------------------------------------------------------------------------------#

version = "1.1.0"
version_number = 2026051800000


#--------------------------------------------------------------------------------------------------#
# Libraries                                                                                        #
#--------------------------------------------------------------------------------------------------#
import math
import numpy as np

import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import Table

import spiceypy as spice


#--------------------------------------------------------------------------------------------------#
# Context functions                                                                                 #
#--------------------------------------------------------------------------------------------------#
def init_satorbit_context(obs_begin, obs_end, obs_step, obs_params, spice_myfile_PATH):
    # SPICE kernel load
    spice.furnsh(spice_myfile_PATH)

    # satphotometry.satorbit must be imported after spice.furnsh()
    from satphotometry import satorbit

    # Observatory
    obs_gd_lon_deg = obs_params[0]
    obs_gd_lat_deg = obs_params[1]
    obs_gd_height  = obs_params[2]

    wavelength_m = obs_params[3]
    aperture_m   = obs_params[4]

    obs_gd_lon = math.radians(obs_gd_lon_deg)
    obs_gd_lat = math.radians(obs_gd_lat_deg)

    # Earth constants for SGP4 propagation
    earth_constants = satorbit.get_planetconst(
        399,
        ["J2", "J3", "J4", "KE", "QO", "SO", "ER", "AE"]
    )

    # Observation ETs
    obs_begin_et = spice.utc2et(obs_begin)
    obs_end_et   = spice.utc2et(obs_end)

    n_step = int((obs_end_et - obs_begin_et) // obs_step + 1)

    obs_ets = np.array(
        [obs_begin_et + i * obs_step for i in range(n_step)],
        dtype=float
    )

    obs_utc_list = [
        spice.et2utc(float(et), "ISOC", 0)
        for et in obs_ets
    ]

    # Observatory position
    obs_itrf = satorbit.geo2itrf(
        obs_gd_lon,
        obs_gd_lat,
        obs_gd_height
    )

    # Observatory J2000 positions
    obs_j2000_list = [
        satorbit.itrf2J2000(obs_itrf, float(et))
        for et in obs_ets
    ]

    ctx = {
        "satorbit": satorbit,
        "earth_constants": earth_constants,
        "obs_ets": obs_ets,
        "obs_utc_list": obs_utc_list,
        "obs_itrf": obs_itrf,
        "obs_j2000_list": obs_j2000_list,
        "obs_gd_lon": obs_gd_lon,
        "obs_gd_lat": obs_gd_lat,
        "obs_gd_height": obs_gd_height,
        "wavelength_m": wavelength_m,
        "aperture_m": aperture_m,
        "obs_begin": obs_begin,
        "obs_end": obs_end,
        "obs_step": obs_step,
        "obs_params": obs_params,
        "spice_myfile_PATH": spice_myfile_PATH,
    }

    return ctx


def clear_satorbit_context():
    spice.kclear()


def make_context_with_new_obs_step(ctx, obs_begin, obs_end, obs_step):
    satorbit = ctx["satorbit"]

    obs_begin_et = spice.utc2et(obs_begin)
    obs_end_et   = spice.utc2et(obs_end)

    n_step = int((obs_end_et - obs_begin_et) // obs_step + 1)

    obs_ets = np.array(
        [obs_begin_et + i * obs_step for i in range(n_step)],
        dtype=float
    )

    obs_utc_list = [
        spice.et2utc(float(et), "ISOC", 0)
        for et in obs_ets
    ]

    obs_itrf = ctx["obs_itrf"]

    obs_j2000_list = [
        satorbit.itrf2J2000(obs_itrf, float(et))
        for et in obs_ets
    ]

    ctx_new = ctx.copy()
    ctx_new["obs_ets"] = obs_ets
    ctx_new["obs_utc_list"] = obs_utc_list
    ctx_new["obs_j2000_list"] = obs_j2000_list
    ctx_new["obs_step"] = obs_step
    ctx_new["obs_begin"] = obs_begin
    ctx_new["obs_end"] = obs_end

    return ctx_new


#--------------------------------------------------------------------------------------------------#
# TLE utility                                                                                       #
#--------------------------------------------------------------------------------------------------#
def _read_tle(tle_PATH, satorbit):
    if isinstance(tle_PATH, str):
        satname, line1, line2 = satorbit.read_TLEfile(tle_PATH)
    else:
        satname = tle_PATH[0]
        line1 = tle_PATH[1]
        line2 = tle_PATH[2]

    epoch, elems = satorbit.parse_TLE2element(line1, line2)

    tle_epoch_utc = spice.et2utc(epoch, "ISOC", 0)
    tle_epoch_utc = tle_epoch_utc.replace("T", " ") + " UTC"

    if satname is not None:
        if satname[0:2] == "0 ":
            objname = satname[2:]
        else:
            objname = satname

        objname = objname.rstrip()
    else:
        objname = None

    norad_id = str(line1[2:7].strip())

    intldes = line1.split()[2]

    if len(intldes) > 8:
        intldes = "0000-000A"
    else:
        if int(intldes[:2]) >= 57:
            intldes = "19" + intldes[:2] + "-" + intldes[2:]
        else:
            intldes = "20" + intldes[:2] + "-" + intldes[2:]

    return objname, norad_id, intldes, elems, tle_epoch_utc


#--------------------------------------------------------------------------------------------------#
# Light calculation                                                                                 #
#--------------------------------------------------------------------------------------------------#
def cal_satorbit_radec_only(tle_PATH, ctx):
    satorbit = ctx["satorbit"]
    earth_constants = ctx["earth_constants"]
    obs_ets = ctx["obs_ets"]
    obs_j2000_list = ctx["obs_j2000_list"]

    objname, norad_id, intldes, elems, tle_epoch_utc = _read_tle(
        tle_PATH,
        satorbit
    )

    ra_deg_list = []
    dec_deg_list = []
    umbra_list = []

    for et, obs_j2000 in zip(obs_ets, obs_j2000_list):
        et = float(et)

        state_teme = spice.evsgp4(
            et,
            earth_constants,
            elems
        )

        state_j2000 = satorbit.teme2J2000(
            state_teme,
            et
        )

        _, ra, dec = satorbit.J20002radec(
            state_j2000[0:3],
            obs_j2000
        )

        in_umbra = satorbit.check_umbra(
            state_j2000[0:3],
            et
        )

        ra_deg_list.append(math.degrees(ra) % 360.0)
        dec_deg_list.append(math.degrees(dec))
        umbra_list.append(bool(in_umbra))

    return (
        objname,
        norad_id,
        intldes,
        tle_epoch_utc,
        np.asarray(ra_deg_list, dtype=float),
        np.asarray(dec_deg_list, dtype=float),
        np.asarray(umbra_list, dtype=bool),
    )


#--------------------------------------------------------------------------------------------------#
# Full calculation                                                                                  #
#--------------------------------------------------------------------------------------------------#
def cal_satorbit_full(tle_PATH, ctx):
    satorbit = ctx["satorbit"]
    earth_constants = ctx["earth_constants"]

    obs_ets = ctx["obs_ets"]
    obs_utc_list = ctx["obs_utc_list"]
    obs_itrf = ctx["obs_itrf"]
    obs_j2000_list = ctx["obs_j2000_list"]

    obs_gd_lon = ctx["obs_gd_lon"]
    obs_gd_lat = ctx["obs_gd_lat"]

    wavelength_m = ctx["wavelength_m"]
    aperture_m = ctx["aperture_m"]

    objname, norad_id, intldes, elems, tle_epoch_utc = _read_tle(
        tle_PATH,
        satorbit
    )

    output = []

    for et, et_utc, obs_j2000 in zip(obs_ets, obs_utc_list, obs_j2000_list):
        et = float(et)

        state_teme = spice.evsgp4(
            et,
            earth_constants,
            elems
        )

        state_j2000 = satorbit.teme2J2000(
            state_teme,
            et
        )

        range_km, ra, dec = satorbit.J20002radec(
            state_j2000[0:3],
            obs_j2000
        )

        ra_hms = Angle(ra * u.rad).to_string(
            unit=u.hourangle,
            sep=":",
            precision=1
        )
        ra_hms = "+0" + ra_hms if len(ra_hms) < 10 else "+" + ra_hms

        dec_dms = Angle(dec * u.rad).to_string(
            unit=u.deg,
            sep=":",
            precision=1
        )

        dec_dms = (
            "+0" + dec_dms
            if (dec > 0 and len(dec_dms) < 10)
            else (
                "+" + dec_dms
                if dec > 0
                else (
                    "-0" + dec_dms[1:]
                    if len(dec_dms) < 11
                    else dec_dms
                )
            )
        )

        ra_deg = math.degrees(ra) % 360.0
        dec_deg = math.degrees(dec)

        state_itrf = satorbit.J20002itrf(
            state_j2000,
            et
        )

        _, az, el = satorbit.itrf2azel(
            state_itrf[0:3],
            obs_itrf,
            obs_gd_lon,
            obs_gd_lat
        )

        az_deg = math.degrees(az)
        el_deg = math.degrees(el)

        in_umbra = satorbit.check_umbra(
            state_j2000[0:3],
            et
        )

        phase = satorbit.phase_angle(
            state_j2000[0:3],
            obs_j2000,
            et
        )

        phase_deg = math.degrees(phase)

        apparent_v_km_s = satorbit.apparent_v(
            state_j2000,
            obs_itrf,
            et
        )

        res_km = range_km * (wavelength_m / aperture_m)

        if apparent_v_km_s == 0:
            ex_per_s = np.nan
        else:
            ex_ms = res_km / apparent_v_km_s
            ex_per_s = 1.0 / ex_ms if ex_ms != 0 else np.nan

        output.append(
            [
                et_utc,
                bool(in_umbra),
                round(range_km, 4),
                ra_hms,
                dec_dms,
                round(ra_deg, 4),
                round(dec_deg, 4),
                round(phase_deg, 4),
                round(az_deg, 4),
                round(el_deg, 4),
                round(apparent_v_km_s, 4),
                round(res_km * 1000, 4),
                round(ex_per_s, 4) if np.isfinite(ex_per_s) else np.nan,
            ]
        )

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
            "ex[1/s]",
        )
    )

    return objname, norad_id, intldes, tle_epoch_utc, output


#--------------------------------------------------------------------------------------------------#
# Optional: circular FOV utility                                                                    #
#--------------------------------------------------------------------------------------------------#
def radec_to_unit_vector(ra_deg, dec_deg):
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    return np.stack([x, y, z], axis=-1)


def angular_separation_fast_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    v1 = radec_to_unit_vector(ra1_deg, dec1_deg)
    v2 = radec_to_unit_vector(ra2_deg, dec2_deg)

    dot = np.sum(v1 * v2, axis=-1)
    dot = np.clip(dot, -1.0, 1.0)

    return np.rad2deg(np.arccos(dot))


def min_sep_to_segment_deg(
    center_ra_deg,
    center_dec_deg,
    ra1_deg,
    dec1_deg,
    ra2_deg,
    dec2_deg
):
    c = radec_to_unit_vector(center_ra_deg, center_dec_deg)
    a = radec_to_unit_vector(ra1_deg, dec1_deg)
    b = radec_to_unit_vector(ra2_deg, dec2_deg)

    sep_a = np.arccos(np.clip(np.sum(c * a), -1.0, 1.0))
    sep_b = np.arccos(np.clip(np.sum(c * b), -1.0, 1.0))

    n = np.cross(a, b)
    n_norm = np.linalg.norm(n)

    if n_norm == 0:
        return np.rad2deg(min(sep_a, sep_b))

    n = n / n_norm

    foot = c - np.dot(c, n) * n
    foot_norm = np.linalg.norm(foot)

    if foot_norm == 0:
        return np.rad2deg(min(sep_a, sep_b))

    foot = foot / foot_norm

    total_ab = np.arccos(np.clip(np.sum(a * b), -1.0, 1.0))

    dist_af = np.arccos(np.clip(np.sum(a * foot), -1.0, 1.0))
    dist_fb = np.arccos(np.clip(np.sum(foot * b), -1.0, 1.0))

    foot2 = -foot
    dist_af2 = np.arccos(np.clip(np.sum(a * foot2), -1.0, 1.0))
    dist_f2b = np.arccos(np.clip(np.sum(foot2 * b), -1.0, 1.0))

    eps = 1e-10

    if abs((dist_af + dist_fb) - total_ab) < eps:
        sep_gc = np.arcsin(np.clip(abs(np.dot(c, n)), 0.0, 1.0))
        return np.rad2deg(sep_gc)

    if abs((dist_af2 + dist_f2b) - total_ab) < eps:
        sep_gc = np.arcsin(np.clip(abs(np.dot(c, n)), 0.0, 1.0))
        return np.rad2deg(sep_gc)

    return np.rad2deg(min(sep_a, sep_b))


def in_circular_fov_by_segments(
    ra_deg,
    dec_deg,
    valid_mask,
    fov_center_ra_deg,
    fov_center_dec_deg,
    fov_radius_deg
):
    ra_deg = np.asarray(ra_deg, dtype=float)
    dec_deg = np.asarray(dec_deg, dtype=float)
    valid_mask = np.asarray(valid_mask, dtype=bool)

    if len(ra_deg) == 0:
        return False

    sep_points = angular_separation_fast_deg(
        fov_center_ra_deg,
        fov_center_dec_deg,
        ra_deg,
        dec_deg
    )

    if np.any((sep_points <= fov_radius_deg) & valid_mask):
        return True

    if len(ra_deg) < 2:
        return False

    for j in range(len(ra_deg) - 1):
        if not (valid_mask[j] and valid_mask[j + 1]):
            continue

        min_sep = min_sep_to_segment_deg(
            fov_center_ra_deg,
            fov_center_dec_deg,
            ra_deg[j],
            dec_deg[j],
            ra_deg[j + 1],
            dec_deg[j + 1]
        )

        if min_sep <= fov_radius_deg:
            return True

    return False