#include "mission_analysis.h"
#include "astrodynamics.h"
#include "propagator.h"
#include "telemetry.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static double stumpff_C(double z) {
    if (z > 1e-6) {
        double sz = sqrt(z);
        return (1.0 - cos(sz)) / z;
    } else if (z < -1e-6) {
        double sz = sqrt(-z);
        return (cosh(sz) - 1.0) / (-z);
    }
    return 0.5 - z / 24.0 + (z * z) / 720.0;
}

static double stumpff_S(double z) {
    if (z > 1e-6) {
        double sz = sqrt(z);
        return (sz - sin(sz)) / (sz * sz * sz);
    } else if (z < -1e-6) {
        double sz = sqrt(-z);
        return (sinh(sz) - sz) / (sz * sz * sz);
    }
    return (1.0 / 6.0) - z / 120.0 + (z * z) / 5040.0;
}

static Vector3D eci_vec_to_rtn(Vector3D r_sat, Vector3D v_sat, Vector3D vec_eci) {
    Vector3D u_R = vec_normalize(r_sat);
    Vector3D h_vec = vec_cross(r_sat, v_sat);
    Vector3D u_N = vec_normalize(h_vec);
    Vector3D u_T = vec_cross(u_N, u_R);

    Vector3D rtn = {
        vec_dot(vec_eci, u_R),
        vec_dot(vec_eci, u_T),
        vec_dot(vec_eci, u_N)
    };
    return rtn;
}

LambertSolution solve_lambert_universal(Vector3D r1_vec, Vector3D v1_init,
                                        Vector3D r2_vec, Vector3D v2_targ,
                                        double tof_sec, double mu, bool prograde) {
    LambertSolution sol = {0};
    if (tof_sec <= 1.0) return sol;

    double r1 = vec_norm(r1_vec);
    double r2 = vec_norm(r2_vec);
    if (r1 < 1e-6 || r2 < 1e-6) return sol;

    double cos_dnu = vec_dot(r1_vec, r2_vec) / (r1 * r2);
    if (cos_dnu >  1.0) cos_dnu =  1.0;
    if (cos_dnu < -1.0) cos_dnu = -1.0;

    Vector3D cross_r1r2 = vec_cross(r1_vec, r2_vec);
    double dnu = acos(cos_dnu);

    if (prograde) {
        if (cross_r1r2.z < 0.0) dnu = 2.0 * PI - dnu;
    } else {
        if (cross_r1r2.z >= 0.0) dnu = 2.0 * PI - dnu;
    }

    if (fabs(dnu - PI) < 1e-4) dnu = PI - 1e-4;

    double A = sin(dnu) * sqrt(r1 * r2 / (1.0 - cos(dnu)));
    if (fabs(A) < 1e-9) return sol;

    double z_low = -4.0 * PI;
    double z_high = 4.0 * PI * PI - 1e-3;
    double z = 0.0;

    for (int iter = 0; iter < 90; iter++) {
        z = 0.5 * (z_low + z_high);
        double C = stumpff_C(z);
        double S = stumpff_S(z);
        double y = r1 + r2 + A * (z * S - 1.0) / sqrt(C);

        if (y <= 0.0) {
            z_low = z;
            continue;
        }

        double x = sqrt(y / C);
        double t_calc = (x * x * x * S + A * sqrt(y)) / sqrt(mu);

        if (fabs(t_calc - tof_sec) < 1e-5) break;

        if (t_calc <= tof_sec) z_low = z;
        else z_high = z;
    }

    double C = stumpff_C(z);
    double S = stumpff_S(z);
    double y = r1 + r2 + A * (z * S - 1.0) / sqrt(C);
    if (y <= 0.0) return sol;

    double f     = 1.0 - (y / r1);
    double g     = A * sqrt(y / mu);
    double g_dot = 1.0 - (y / r2);

    if (fabs(g) < 1e-9) return sol;

    sol.v1_eci = vec_scale(vec_sub(r2_vec, vec_scale(r1_vec, f)), 1.0 / g);
    sol.v2_eci = vec_scale(vec_sub(vec_scale(r2_vec, g_dot), r1_vec), 1.0 / g);

    Vector3D dv1_eci = vec_sub(sol.v1_eci, v1_init);
    Vector3D dv2_eci = vec_sub(v2_targ, sol.v2_eci);

    sol.dv1_mag  = vec_norm(dv1_eci);
    sol.dv2_mag  = vec_norm(dv2_eci);
    sol.total_dv = sol.dv1_mag + sol.dv2_mag;

    sol.burn1_dir_rtn = vec_normalize(eci_vec_to_rtn(r1_vec, v1_init, dv1_eci));
    sol.burn2_dir_rtn = vec_normalize(eci_vec_to_rtn(r2_vec, v2_targ, dv2_eci));
    sol.converged = true;

    return sol;
}

PropellantSpec get_propellant_specs(PropellantType type) {
    PropellantSpec p = {0};
    p.type = type;
    switch (type) {
        case PROP_SOLID_HTPB:
            strcpy(p.name, "HTPB/AP Composite (Star-48B)");
            strcpy(p.category, "Solid");
            p.isp = 292.0; p.default_thrust_n = 15000.0; p.is_low_thrust = false;
            break;
        case PROP_LIQUID_HYDRAZINE:
            strcpy(p.name, "Hydrazine Monopropellant (N2H4)");
            strcpy(p.category, "Liquid");
            p.isp = 235.0; p.default_thrust_n = 400.0; p.is_low_thrust = false;
            break;
        case PROP_LIQUID_MMH_NTO:
            strcpy(p.name, "MMH / NTO Bipropellant");
            strcpy(p.category, "Liquid");
            p.isp = 316.0; p.default_thrust_n = 490.0; p.is_low_thrust = false;
            break;
        case PROP_LIQUID_LOX_LCH4:
            strcpy(p.name, "LOX / Liquid Methane (Methalox)");
            strcpy(p.category, "Liquid");
            p.isp = 365.0; p.default_thrust_n = 2000.0; p.is_low_thrust = false;
            break;
        case PROP_ELECTRIC_HALL_XE:
            strcpy(p.name, "Xenon Hall-Effect Thruster (HET)");
            strcpy(p.category, "Electric");
            p.isp = 1650.0; p.default_thrust_n = 0.50; p.is_low_thrust = true;
            break;
        case PROP_ELECTRIC_ION_NEXT:
            strcpy(p.name, "Xenon Gridded Ion Engine (NEXT-C)");
            strcpy(p.category, "Electric");
            p.isp = 3200.0; p.default_thrust_n = 0.25; p.is_low_thrust = true;
            break;
        case PROP_LIQUID_LOX_LH2:
        default:
            strcpy(p.name, "LOX / Liquid Hydrogen (Hydrolox)");
            strcpy(p.category, "Liquid");
            p.isp = 450.0; p.default_thrust_n = 5000.0; p.is_low_thrust = false;
            break;
    }
    return p;
}

bool validate_and_clamp_closed_orbit(KeplerianElements *kep, char *warning_msg) {
    bool valid = true;
    if (warning_msg) warning_msg[0] = '\0';

    double r_min_safe = EARTH_RADIUS + MIN_SAFE_PERIGEE_ALT_KM;

    if (kep->a < r_min_safe) {
        kep->a = r_min_safe;
        valid = false;
        if (warning_msg) strcpy(warning_msg, "Semi-major axis raised to minimum safe LEO altitude (150 km).");
    }

    if (kep->e < 0.0) kep->e = 0.0;
    if (kep->e > 0.95) {
        kep->e = 0.95;
        valid = false;
    }

    double e_max_perigee = 1.0 - (r_min_safe / kep->a);
    if (e_max_perigee < 0.0) e_max_perigee = 0.0;

    if (kep->e > e_max_perigee) {
        kep->e = e_max_perigee;
        valid = false;
        if (warning_msg) {
            snprintf(warning_msg, 128, "Eccentricity clamped to %.4f to prevent Earth impact (Perigee >= 150 km).", e_max_perigee);
        }
    }

    return valid;
}

KeplerianElements resolve_user_orbit_input(OrbitInputSpec input, double mu) {
    KeplerianElements kep = (input.mode == INPUT_CARTESIAN) ? cartesian_to_keplerian(input.cart, mu) : input.kep;
    validate_and_clamp_closed_orbit(&kep, NULL);
    return kep;
}

Vector3D geodetic_to_ecef(GeodeticCoord geo) {
    double e2 = 2.0 * EARTH_FLATTENING - (EARTH_FLATTENING * EARTH_FLATTENING);
    double sin_lat = sin(geo.latitude);
    double cos_lat = cos(geo.latitude);
    double sin_lon = sin(geo.longitude);
    double cos_lon = cos(geo.longitude);

    double N = EARTH_RADIUS / sqrt(1.0 - (e2 * sin_lat * sin_lat));

    Vector3D ecef = {
        (N + geo.altitude) * cos_lat * cos_lon,
        (N + geo.altitude) * cos_lat * sin_lon,
        (N * (1.0 - e2) + geo.altitude) * sin_lat
    };
    return ecef;
}

AccessMetrics compute_ground_access(Vector3D sat_ecef, const GroundStation *station) {
    AccessMetrics metrics = {0};

    Vector3D site_ecef = geodetic_to_ecef(station->coord);
    Vector3D rho_ecef  = vec_sub(sat_ecef, site_ecef);
    metrics.range      = vec_norm(rho_ecef);

    if (metrics.range < 1e-9) return metrics;

    double sin_lat = sin(station->coord.latitude);
    double cos_lat = cos(station->coord.latitude);
    double sin_lon = sin(station->coord.longitude);
    double cos_lon = cos(station->coord.longitude);

    double rho_S = ( sin_lat * cos_lon * rho_ecef.x) + ( sin_lat * sin_lon * rho_ecef.y) - (cos_lat * rho_ecef.z);
    double rho_E = (-sin_lon * rho_ecef.x)           + ( cos_lon * rho_ecef.y);
    double rho_Z = ( cos_lat * cos_lon * rho_ecef.x) + ( cos_lat * sin_lon * rho_ecef.y) + (sin_lat * rho_ecef.z);

    double ratio = rho_Z / metrics.range;
    if (ratio >  1.0) ratio =  1.0;
    if (ratio < -1.0) ratio = -1.0;

    metrics.elevation = asin(ratio);
    metrics.azimuth   = atan2(rho_E, -rho_S);
    if (metrics.azimuth < 0.0) metrics.azimuth += (2.0 * PI);

    metrics.is_visible = (metrics.elevation >= station->min_elevation);
    return metrics;
}

HohmannTransfer compute_hohmann_transfer(double r1, double r2, double mu) {
    HohmannTransfer h = {0};
    double a_trans = 0.5 * (r1 + r2);

    double v1_circ = sqrt(mu / r1);
    double v2_circ = sqrt(mu / r2);

    double v_trans_perigee = sqrt(mu * ((2.0 / r1) - (1.0 / a_trans)));
    double v_trans_apogee  = sqrt(mu * ((2.0 / r2) - (1.0 / a_trans)));

    h.dv1 = v_trans_perigee - v1_circ;
    h.dv2 = v2_circ - v_trans_apogee;
    h.total_dv = fabs(h.dv1) + fabs(h.dv2);
    h.transfer_time = PI * sqrt((a_trans * a_trans * a_trans) / mu);

    return h;
}

CartesianState apply_impulsive_maneuver(CartesianState state, Vector3D delta_v_eci, double isp) {
    CartesianState next = state;
    next.velocity = vec_add(state.velocity, delta_v_eci);

    double dv_mag_ms = vec_norm(delta_v_eci) * 1000.0;
    if (isp > 0.0 && dv_mag_ms > 0.0) {
        next.mass = state.mass * exp(-dv_mag_ms / (isp * G0));
    }
    return next;
}

OptimalTransferPlan optimize_orbital_transfer(KeplerianElements orb_init,
                                              KeplerianElements orb_target,
                                              double init_mass,
                                              PropellantType prop_type,
                                              double max_transfer_time_sec,
                                              double mu) {
    OptimalTransferPlan best = {0};
    best.prop_info = get_propellant_specs(prop_type);
    best.is_low_thrust = best.prop_info.is_low_thrust;
    best.continuous_thrust_n = best.prop_info.default_thrust_n;
    best.time_constrained = (max_transfer_time_sec > 0.0);
    best.max_time_limit_s = max_transfer_time_sec;

    validate_and_clamp_closed_orbit(&orb_init, NULL);
    validate_and_clamp_closed_orbit(&orb_target, NULL);

    double r1 = orb_init.a * (1.0 - orb_init.e);
    double r2 = orb_target.a * (1.0 + orb_target.e);

    double v1 = sqrt(mu * ((2.0 / r1) - (1.0 / orb_init.a)));
    double v2 = sqrt(mu * ((2.0 / r2) - (1.0 / orb_target.a)));

    double cos_theta = cos(orb_init.i) * cos(orb_target.i) +
                       sin(orb_init.i) * sin(orb_target.i) * cos(orb_target.RAAN - orb_init.RAAN);
    if (cos_theta >  1.0) cos_theta =  1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    double total_plane_angle = acos(cos_theta);
    double sign_inc = (orb_target.i >= orb_init.i) ? 1.0 : -1.0;

    // ==================================================================
    // CASO A: PROPULSIONE ELETTRICA CONTINUA A SPIRALE (EDELBAUM THEORY)
    // ==================================================================
    if (best.is_low_thrust) {
        double v1_c = sqrt(mu / orb_init.a);
        double v2_c = sqrt(mu / orb_target.a);
        // Formula di Edelbaum per il Delta-V di una spirale continua 3D
        double eff_angle = (PI / 2.0) * total_plane_angle;
        if (eff_angle > PI) eff_angle = PI;
        double dv_edelbaum = sqrt(v1_c*v1_c + v2_c*v2_c - 2.0*v1_c*v2_c*cos(eff_angle));

        strcpy(best.maneuver_name, "Continuous Low-Thrust Edelbaum Spiral Transfer");
        best.num_burns = 1; // Unica accensione continua dall'inizio alla fine
        best.dv1_mag = dv_edelbaum;
        best.dv2_mag = 0.0;
        best.dv3_mag = 0.0;
        best.total_dv = dv_edelbaum;

        // Angolo medio di imbardata (Yaw) per spingere sia lungo la traccia che fuori piano
        double yaw = atan2(v2_c * sin(eff_angle), fabs(v1_c - v2_c * cos(eff_angle)) + 1e-6);
        best.burn1_dir_rtn = (Vector3D){0.0, cos(yaw) * ((orb_target.a >= orb_init.a) ? 1.0 : -1.0), sin(yaw) * sign_inc};
        best.burn2_dir_rtn = best.burn1_dir_rtn;

        double total_dv_ms = best.total_dv * 1000.0;
        best.final_mass = init_mass * exp(-total_dv_ms / (best.prop_info.isp * G0));
        best.propellant_used = init_mass - best.final_mass;

        // Durata reale di accensione continua: t = (m_prop * Isp * g0) / Thrust
        double mdot = best.continuous_thrust_n / (best.prop_info.isp * G0);
        best.transfer_time_sec = best.propellant_used / mdot;

        double mean_period = 2.0 * PI * sqrt(pow(0.5*(orb_init.a + orb_target.a), 3) / mu);
        best.spiral_revolutions = best.transfer_time_sec / mean_period;
        best.transfer_apogee_km = orb_target.a;

        KeplerianElements k_start = orb_init; k_start.theta = 0.0;
        CartesianState s_start = keplerian_to_cartesian(k_start, mu, init_mass, 0.0);
        best.burn1_pos_eci = s_start.position;
        best.burn1_geo     = ecef_to_geodetic(eci_to_ecef(s_start.position, 0.0, 0.0));

        KeplerianElements k_end = orb_target; k_end.theta = PI;
        CartesianState s_end = keplerian_to_cartesian(k_end, mu, best.final_mass, best.transfer_time_sec);
        best.burn2_pos_eci = s_end.position;
        best.burn2_geo     = ecef_to_geodetic(eci_to_ecef(s_end.position, best.transfer_time_sec, 0.0));

        return best;
    }

    // ==================================================================
    // CASO B: PROPULSIONE CHIMICA IMPULSIVA (HOHMANN / BI-ELLITTICA / LAMBERT)
    // ==================================================================
    double a_hoh = 0.5 * (r1 + r2);
    double tof_hoh = PI * sqrt((a_hoh * a_hoh * a_hoh) / mu);
    double vt1_hoh = sqrt(mu * ((2.0 / r1) - (1.0 / a_hoh)));
    double vt2_hoh = sqrt(mu * ((2.0 / r2) - (1.0 / a_hoh)));

    double best_alpha = 0.0;
    double min_hoh_dv = 1e9;

    for (int k = 0; k <= 1000; k++) {
        double alpha = (double)k / 1000.0;
        double d1 = alpha * total_plane_angle;
        double d2 = (1.0 - alpha) * total_plane_angle;

        double dv1 = sqrt(v1*v1 + vt1_hoh*vt1_hoh - 2.0*v1*vt1_hoh*cos(d1));
        double dv2 = sqrt(v2*v2 + vt2_hoh*vt2_hoh - 2.0*v2*vt2_hoh*cos(d2));
        if ((dv1 + dv2) < min_hoh_dv) {
            min_hoh_dv = dv1 + dv2;
            best_alpha = alpha;
        }
    }

    double d1_opt = best_alpha * total_plane_angle * sign_inc;
    double d2_opt = (1.0 - best_alpha) * total_plane_angle * sign_inc;

    Vector3D b1_hoh = {0.0, vt1_hoh * cos(d1_opt) - v1, vt1_hoh * sin(d1_opt)};
    Vector3D b2_hoh = {0.0, v2 * cos(d2_opt) - vt2_hoh, v2 * sin(d2_opt)};

    strcpy(best.maneuver_name, (total_plane_angle > 1e-4) ? "3D Optimal Split-Plane Transfer (2-Burn)" : "Coplanar Elliptical Hohmann Transfer (2-Burn)");
    best.num_burns = 2;
    best.dv1_mag = vec_norm(b1_hoh);
    best.dv2_mag = vec_norm(b2_hoh);
    best.dv3_mag = 0.0;
    best.total_dv = best.dv1_mag + best.dv2_mag;
    best.burn1_dir_rtn = vec_normalize(b1_hoh);
    best.burn2_dir_rtn = vec_normalize(b2_hoh);
    best.transfer_apogee_km = (r2 >= r1) ? r2 : r1;
    best.transfer_time_sec = tof_hoh;

    // Candidato Bi-Ellittica
    double r_max = (r2 > r1) ? r2 : r1;
    for (int mult = 2; mult <= 6; mult++) {
        double rb = r_max * (double)mult;
        double a_t1 = 0.5 * (r1 + rb);
        double a_t2 = 0.5 * (r2 + rb);
        double tof_bi = PI * (sqrt(a_t1*a_t1*a_t1 / mu) + sqrt(a_t2*a_t2*a_t2 / mu));

        if (best.time_constrained && tof_bi > max_transfer_time_sec) continue;

        double vt1_p = sqrt(mu * ((2.0 / r1) - (1.0 / a_t1)));
        double vt1_a = sqrt(mu * ((2.0 / rb) - (1.0 / a_t1)));
        double vt2_a = sqrt(mu * ((2.0 / rb) - (1.0 / a_t2)));
        double vt2_p = sqrt(mu * ((2.0 / r2) - (1.0 / a_t2)));

        double dv1_bi = fabs(vt1_p - v1);
        double dv2_bi = sqrt(vt1_a*vt1_a + vt2_a*vt2_a - 2.0*vt1_a*vt2_a*cos(total_plane_angle));
        double dv3_bi = fabs(v2 - vt2_p);
        double tot_bi = dv1_bi + dv2_bi + dv3_bi;

        if (tot_bi < best.total_dv) {
            strcpy(best.maneuver_name, "Optimal 3-Burn Bi-Elliptic Transfer");
            best.num_burns = 3;
            best.dv1_mag = dv1_bi;
            best.dv2_mag = dv2_bi;
            best.dv3_mag = dv3_bi;
            best.total_dv = tot_bi;
            best.burn1_dir_rtn = (Vector3D){0.0, 1.0, 0.0};
            Vector3D b2_vec = {0.0, vt2_a * cos(total_plane_angle) - vt1_a, vt2_a * sin(total_plane_angle * sign_inc)};
            best.burn2_dir_rtn = vec_normalize(b2_vec);
            best.burn3_dir_rtn = (Vector3D){0.0, -1.0, 0.0};
            best.transfer_apogee_km = rb;
            best.transfer_time_sec = tof_bi;
        }
    }

    // Candidato Lambert Veloce
    KeplerianElements k_b1 = orb_init;
    k_b1.theta = 0.0;
    CartesianState s_b1 = keplerian_to_cartesian(k_b1, mu, init_mass, 0.0);

    KeplerianElements k_b2 = orb_target;
    k_b2.theta = 165.0 * DEG2RAD;
    CartesianState s_b2_lam = keplerian_to_cartesian(k_b2, mu, init_mass, 0.0);

    if (best.time_constrained && max_transfer_time_sec < tof_hoh) {
        LambertSolution lam = solve_lambert_universal(s_b1.position, s_b1.velocity,
                                                      s_b2_lam.position, s_b2_lam.velocity,
                                                      max_transfer_time_sec, mu, true);
        if (lam.converged) {
            strcpy(best.maneuver_name, "3D Universal Lambert Conic Transfer (Time-Constrained)");
            best.num_burns = 2;
            best.dv1_mag = lam.dv1_mag;
            best.dv2_mag = lam.dv2_mag;
            best.dv3_mag = 0.0;
            best.total_dv = lam.total_dv;
            best.burn1_dir_rtn = lam.burn1_dir_rtn;
            best.burn2_dir_rtn = lam.burn2_dir_rtn;
            best.transfer_time_sec = max_transfer_time_sec;
        }
    }

    best.burn1_pos_eci = s_b1.position;
    best.burn1_geo     = ecef_to_geodetic(eci_to_ecef(s_b1.position, 0.0, 0.0));

    KeplerianElements k_b2_std = orb_target;
    if (best.num_burns == 3) k_b2_std.a = best.transfer_apogee_km;
    k_b2_std.theta = PI;
    CartesianState s_b2 = keplerian_to_cartesian(k_b2_std, mu, init_mass, best.transfer_time_sec);
    best.burn2_pos_eci = s_b2.position;
    best.burn2_geo     = ecef_to_geodetic(eci_to_ecef(s_b2.position, best.transfer_time_sec, 0.0));

    double total_dv_ms = best.total_dv * 1000.0;
    best.final_mass = init_mass * exp(-total_dv_ms / (best.prop_info.isp * G0));
    best.propellant_used = init_mass - best.final_mass;

    return best;
}

OrbitReport evaluate_orbit_regime(KeplerianElements kep, double mass, const SpacecraftConfig *sc, const GroundStation *station) {
    OrbitReport rep = {0};
    validate_and_clamp_closed_orbit(&kep, NULL);
    rep.initial_kep = kep;
    rep.period_min  = orbital_period(kep.a, EARTH_MU) / 60.0;

    double r_perigee = kep.a * (1.0 - kep.e);
    Vector3D r_dummy = {r_perigee, 0.0, 0.0};
    Vector3D sun_dir = {1.0, 0.0, 0.0};
    SpaceEnvironment env_mean = evaluate_environment(r_dummy, sun_dir);

    double v_p_ms = sqrt(EARTH_MU * ((2.0 / r_perigee) - (1.0 / kep.a))) * 1000.0;
    double a_drag_ms2 = 0.5 * env_mean.air_density * (v_p_ms * v_p_ms) * (sc->drag_coeff * sc->drag_area / mass);

    double dv_day_ms = a_drag_ms2 * 86400.0;
    rep.station_keeping_dv_year = dv_day_ms * 365.25;
    rep.sma_decay_m_per_day = (2.0 * (kep.a * 1000.0) * a_drag_ms2 / v_p_ms) * 86400.0;

    if (sc->thruster.isp > 0.0) {
        rep.station_keeping_fuel_yr = mass * (1.0 - exp(-rep.station_keeping_dv_year / (sc->thruster.isp * G0)));
    }

    CartesianState state = keplerian_to_cartesian(kep, EARTH_MU, mass, 0.0);
    ForceModelConfig forces = {.enable_j2 = true, .enable_drag = true, .sun_dir = sun_dir};
    SpacecraftConfig sc_passive = *sc;
    sc_passive.thruster.is_active = false;

    double dt = 10.0;
    double total_time = 86400.0;
    double eclipse_sec = 0.0;
    double contact_sec = 0.0;
    bool in_pass = false;

    while (state.epoch < total_time) {
        SpaceEnvironment env = evaluate_environment(state.position, forces.sun_dir);
        if (env.in_eclipse) eclipse_sec += dt;

        Vector3D pos_ecef = eci_to_ecef(state.position, state.epoch, 0.0);
        AccessMetrics acc = compute_ground_access(pos_ecef, station);

        if (acc.is_visible) {
            contact_sec += dt;
            if (!in_pass) {
                in_pass = true;
                rep.ground_passes_24h++;
            }
            if ((acc.elevation * RAD2DEG) > rep.max_elevation_deg) {
                rep.max_elevation_deg = acc.elevation * RAD2DEG;
            }
        } else {
            in_pass = false;
        }

        state = propagate_step_rk4(state, dt, &sc_passive, &forces);
    }

    rep.eclipse_time_min   = eclipse_sec / 60.0;
    rep.sunlight_time_min  = (total_time - eclipse_sec) / 60.0;
    rep.eclipse_percentage = (eclipse_sec / total_time) * 100.0;
    rep.total_contact_min  = contact_sec / 60.0;

    return rep;
}

void print_orbit_report(const char *title, OrbitReport rep, const char *station_name) {
    printf("------------------------------------------------------------------\n");
    printf("  %s\n", title);
    printf("------------------------------------------------------------------\n");
    printf("  * Altitude / Inc / Ecc  : %.1f km | %.2f deg | e = %.4f\n",
           rep.initial_kep.a - EARTH_RADIUS, rep.initial_kep.i * RAD2DEG, rep.initial_kep.e);
    printf("  * Orbital Period        : %.2f min\n", rep.period_min);
    printf("  * Sunlight / Eclipse    : %.1f min / %.1f min (%.1f%% shadow)\n",
           rep.sunlight_time_min, rep.eclipse_time_min, rep.eclipse_percentage);
    printf("  * Ground Access (%s): %d passes/24h | Total = %.1f min\n",
           station_name, rep.ground_passes_24h, rep.total_contact_min);
    printf("  * Station-Keeping / Yr  : %.2f m/s/yr | Fuel = %.3f kg/yr\n\n",
           rep.station_keeping_dv_year, rep.station_keeping_fuel_yr);
}

// =========================================================================
// SPRINT 12: CLOHESSY-WILTSHIRE (HILL) TERMINAL RENDEZVOUS & PROXIMITY OPS
// Frame LVLH: x = Radiale (R-bar), y = In-Track (V-bar), z = Cross-Track (H-bar)
// =========================================================================
void propagate_cw_state(Vector3D r0_m, Vector3D v0_ms, double n_rad_s, double t_sec,
                        Vector3D *r_out_m, Vector3D *v_out_ms) {
    if (n_rad_s <= 1e-12) {
        if (r_out_m) *r_out_m = vec_add(r0_m, vec_scale(v0_ms, t_sec));
        if (v_out_ms) *v_out_ms = v0_ms;
        return;
    }

    double nt = n_rad_s * t_sec;
    double c = cos(nt);
    double s = sin(nt);

    // Matrice Phi_rr(t) e Phi_rv(t) per la posizione relativa dr(t)
    if (r_out_m) {
        r_out_m->x = (4.0 - 3.0 * c) * r0_m.x
                   + (s / n_rad_s) * v0_ms.x
                   + (2.0 * (1.0 - c) / n_rad_s) * v0_ms.y;

        r_out_m->y = 6.0 * (s - nt) * r0_m.x
                   + 1.0 * r0_m.y
                   - (2.0 * (1.0 - c) / n_rad_s) * v0_ms.x
                   + ((4.0 * s - 3.0 * nt) / n_rad_s) * v0_ms.y;

        r_out_m->z = c * r0_m.z + (s / n_rad_s) * v0_ms.z;
    }

    // Matrice Phi_vr(t) e Phi_vv(t) per la velocità relativa dv(t)
    if (v_out_ms) {
        v_out_ms->x = (3.0 * n_rad_s * s) * r0_m.x
                    + c * v0_ms.x
                    + (2.0 * s) * v0_ms.y;

        v_out_ms->y = (-6.0 * n_rad_s * (1.0 - c)) * r0_m.x
                    - (2.0 * s) * v0_ms.x
                    + (4.0 * c - 3.0) * v0_ms.y;

        v_out_ms->z = (-n_rad_s * s) * r0_m.z + c * v0_ms.z;
    }
}

RendezvousPlan compute_cw_rendezvous(double a_target_km,
                                     Vector3D r0_lvlh_m,
                                     Vector3D v0_lvlh_ms,
                                     double t_rdv_sec,
                                     double init_mass_kg,
                                     double isp_s,
                                     double mu) {
    RendezvousPlan rdv = {0};
    rdv.r0_lvlh_m = r0_lvlh_m;
    rdv.v0_lvlh_ms = v0_lvlh_ms;
    rdv.initial_range_m = vec_norm(r0_lvlh_m);

    if (a_target_km < EARTH_RADIUS + MIN_SAFE_PERIGEE_ALT_KM) {
        a_target_km = EARTH_RADIUS + MIN_SAFE_PERIGEE_ALT_KM;
    }
    if (t_rdv_sec < 60.0) {
        t_rdv_sec = 60.0;
    }

    double n = sqrt(mu / (a_target_km * a_target_km * a_target_km));
    rdv.mean_motion_rad_s = n;
    rdv.transfer_time_sec = t_rdv_sec;

    double nt = n * t_rdv_sec;
    double c = cos(nt);
    double s = sin(nt);

    // Determinante normalizzato del blocco planare 2x2 di Phi_rv(t):
    // det(Phi_rv_xy) = [ s*(4s - 3nt) + 4*(1 - c)^2 ] / n^2 = [ 8*(1 - c) - 3*nt*s ] / n^2
    double det_xy_norm = s * (4.0 * s - 3.0 * nt) + 4.0 * (1.0 - c) * (1.0 - c);

    // Gestione singolarità se det(Phi_rv) -> 0 o sin(nt) -> 0 (multipli del mezzo periodo orbitale)
    if (fabs(det_xy_norm) < 1e-6 || fabs(s) < 1e-4) {
        t_rdv_sec += 30.0; // Perturbazione di sicurezza anti-singolarità
        rdv.transfer_time_sec = t_rdv_sec;
        nt = n * t_rdv_sec;
        c = cos(nt);
        s = sin(nt);
        det_xy_norm = s * (4.0 * s - 3.0 * nt) + 4.0 * (1.0 - c) * (1.0 - c);
    }

    if (fabs(det_xy_norm) < 1e-9) {
        det_xy_norm = (det_xy_norm >= 0.0) ? 1e-9 : -1e-9;
    }
    double s_safe = (fabs(s) < 1e-5) ? ((s >= 0.0 ? 1.0 : -1.0) * 1e-5) : s;

    // Vettore termine noto b = -Phi_rr(t_rdv) * dr0
    double bx = -(4.0 - 3.0 * c) * r0_lvlh_m.x;
    double by = -6.0 * (s - nt) * r0_lvlh_m.x - r0_lvlh_m.y;
    double bz = -c * r0_lvlh_m.z;

    // Inversione analitica v0_post = Phi_rv^-1(t_rdv) * b
    Vector3D v0_post = {
        (n / det_xy_norm) * ((4.0 * s - 3.0 * nt) * bx - 2.0 * (1.0 - c) * by),
        (n / det_xy_norm) * (2.0 * (1.0 - c) * bx + s * by),
        (n / s_safe) * bz
    };
    rdv.v0_post_lvlh_ms = v0_post;

    // Impulso 1: Delta-V_inject = v0_post - v0_init
    rdv.dv_inject_lvlh_ms = vec_sub(v0_post, v0_lvlh_ms);
    rdv.dv_inject_ms = vec_norm(rdv.dv_inject_lvlh_ms);

    // Velocità relativa d'arrivo a t = t_rdv^- prima della frenata di Soft-Docking
    Vector3D r_arr_m = {0};
    Vector3D v_arr_ms = {0};
    propagate_cw_state(r0_lvlh_m, v0_post, n, t_rdv_sec, &r_arr_m, &v_arr_ms);

    // Impulso 2: Delta-V_brake = -v(t_rdv^-) per annullare la velocità relativa a [0, 0, 0]
    rdv.dv_brake_lvlh_ms = vec_scale(v_arr_ms, -1.0);
    rdv.dv_brake_ms = vec_norm(rdv.dv_brake_lvlh_ms);

    rdv.total_dv_ms = rdv.dv_inject_ms + rdv.dv_brake_ms;
    if (isp_s > 0.0 && init_mass_kg > 0.0) {
        rdv.final_mass_kg = init_mass_kg * exp(-rdv.total_dv_ms / (isp_s * G0));
        rdv.propellant_used_kg = init_mass_kg - rdv.final_mass_kg;
    } else {
        rdv.final_mass_kg = init_mass_kg;
        rdv.propellant_used_kg = 0.0;
    }
    rdv.converged = true;

    return rdv;
}

// =========================================================================
// SPRINT 13: MONTE CARLO DISPERSION, 3-SIGMA ELLIPSOID & MID-COURSE TCM
// Generatore Gaussiano deterministico (LCG 32-bit + Trasformata di Box-Muller)
// =========================================================================
static double mc_uniform01(uint32_t *state) {
    *state = (*state * 1664525u + 1013904223u);
    return ((double)(*state & 0x00FFFFFFu) + 0.5) / 16777216.0;
}

static double mc_randn_box_muller(uint32_t *state) {
    double u1 = mc_uniform01(state);
    double u2 = mc_uniform01(state);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

MonteCarloReport run_monte_carlo_dispersion(KeplerianElements orb_init,
                                            KeplerianElements orb_target,
                                            OptimalTransferPlan plan,
                                            double sigma_dv_pct,
                                            double sigma_point_deg,
                                            int num_runs,
                                            double mu) {
    MonteCarloReport mc = {0};
    if (num_runs <= 0 || num_runs > MAX_MC_RUNS) {
        num_runs = MAX_MC_RUNS;
    }
    if (sigma_dv_pct < 0.01) sigma_dv_pct = 0.01;
    if (sigma_point_deg < 0.01) sigma_point_deg = 0.01;

    mc.sigma_dv_pct    = sigma_dv_pct;
    mc.sigma_point_deg = sigma_point_deg;
    mc.num_runs        = num_runs;

    double a_trans = 0.5 * (orb_init.a + orb_target.a);
    if (a_trans < EARTH_RADIUS + MIN_SAFE_PERIGEE_ALT_KM) {
        a_trans = EARTH_RADIUS + MIN_SAFE_PERIGEE_ALT_KM;
    }
    double n_trans = sqrt(mu / (a_trans * a_trans * a_trans));
    double tof_hoh = PI / n_trans;

    // Tempo effettivo di propagazione della sensitività (165.6 deg di anomalia di trasferimento)
    double t_arr = 0.92 * tof_hoh;
    double t_mid = 0.50 * t_arr;
    double t_rem = t_arr - t_mid;

    double dv1_nom_kms = (plan.dv1_mag > 0.05) ? plan.dv1_mag : 0.12;
    Vector3D u1 = vec_normalize(plan.burn1_dir_rtn);
    if (vec_norm(u1) < 1e-6) {
        u1 = (Vector3D){0.0, 1.0, 0.0};
    }

    // Base ortonormale locale attorno al vettore di spinta nominale u1 nel frame RTN
    Vector3D p_in = {1.0, 0.0, 0.0};
    if (fabs(u1.x) > 0.9) {
        p_in = (Vector3D){0.0, 1.0, 0.0};
    }
    p_in = vec_normalize(vec_sub(p_in, vec_scale(u1, vec_dot(p_in, u1))));
    Vector3D p_out = vec_normalize(vec_cross(u1, p_in));

    double sigma_mag_rel = sigma_dv_pct / 100.0;
    double sigma_pt_rad  = sigma_point_deg * DEG2RAD;
    double sigma_r0_km   = 0.18;   // Incertezza 1-sigma navigazione GPS/OD su posizione iniziale (km)
    double sigma_v0_kms  = 0.0003; // Incertezza 1-sigma navigazione su velocità iniziale (km/s)

    uint32_t rng_state = 20261301u;

    double sum_miss = 0.0, sum_miss2 = 0.0;
    double sum_r2 = 0.0, sum_t2 = 0.0, sum_n2 = 0.0;
    double sum_tcm = 0.0, sum_tcm2 = 0.0;

    for (int i = 0; i < num_runs; i++) {
        double z1 = mc_randn_box_muller(&rng_state);
        double z2 = mc_randn_box_muller(&rng_state);
        double z3 = mc_randn_box_muller(&rng_state);
        double z4 = mc_randn_box_muller(&rng_state);
        double z5 = mc_randn_box_muller(&rng_state);
        double z6 = mc_randn_box_muller(&rng_state);

        // Errore di stato iniziale di navigazione alla partenza (Burn 1) in km
        Vector3D dr0_km = {
            sigma_r0_km * z1,
            sigma_r0_km * 1.35 * z2,
            sigma_r0_km * 0.85 * z3
        };

        // Perturbazione impulsiva su Burn 1: modulo (z4) + puntamento pitch/yaw (z5, z6) + errore OD velocità
        double d_mag   = sigma_mag_rel * z4;
        double d_pitch = sigma_pt_rad  * z5;
        double d_yaw   = sigma_pt_rad  * z6;

        Vector3D dv1_err_kms = vec_add(
            vec_scale(u1,    dv1_nom_kms * d_mag),
            vec_add(
                vec_scale(p_in,  dv1_nom_kms * d_pitch),
                vec_scale(p_out, dv1_nom_kms * d_yaw)
            )
        );
        dv1_err_kms.x += sigma_v0_kms * z1;
        dv1_err_kms.y += sigma_v0_kms * z2;
        dv1_err_kms.z += sigma_v0_kms * z3;

        // 1. Propagazione fino all'arrivo (Burn 2) in assenza di correzione (Uncorrected Dispersion)
        Vector3D dr_arr_km  = {0};
        Vector3D dv_arr_kms = {0};
        propagate_cw_state(dr0_km, dv1_err_kms, n_trans, t_arr, &dr_arr_km, &dv_arr_kms);
        dr_arr_km.z += (dv1_err_kms.z / n_trans) * 0.60;

        double miss_km = vec_norm(dr_arr_km);

        // 2. Stato perturbato a metà trasferimento (Mid-Course t = 0.5 * TOF) e calcolo TCM
        Vector3D dr_mid_km  = {0};
        Vector3D dv_mid_kms = {0};
        propagate_cw_state(dr0_km, dv1_err_kms, n_trans, t_mid, &dr_mid_km, &dv_mid_kms);

        RendezvousPlan tcm_sol = compute_cw_rendezvous(a_trans, dr_mid_km, dv_mid_kms, t_rem,
                                                       plan.final_mass, plan.prop_info.isp, mu);
        // Poiché dr_mid_km è in km e dv_mid_kms in km/s, per linearità CW dv_inject_ms è in km/s -> *1000 per m/s
        double tcm_dv_ms = tcm_sol.dv_inject_ms * 1000.0;

        mc.samples[i].dr_rtn_km        = dr_arr_km;
        mc.samples[i].arrival_pos_eci  = vec_add(plan.burn2_pos_eci, dr_arr_km);
        mc.samples[i].miss_distance_km = miss_km;
        mc.samples[i].tcm_dv_ms        = tcm_dv_ms;

        sum_miss  += miss_km;
        sum_miss2 += miss_km * miss_km;
        if (miss_km > mc.max_miss_km) mc.max_miss_km = miss_km;

        sum_r2 += dr_arr_km.x * dr_arr_km.x;
        sum_t2 += dr_arr_km.y * dr_arr_km.y;
        sum_n2 += dr_arr_km.z * dr_arr_km.z;

        sum_tcm  += tcm_dv_ms;
        sum_tcm2 += tcm_dv_ms * tcm_dv_ms;
    }

    double inv_n = 1.0 / (double)num_runs;
    mc.mean_miss_km = sum_miss * inv_n;
    double var_miss = (sum_miss2 * inv_n) - (mc.mean_miss_km * mc.mean_miss_km);
    mc.std_miss_km  = (var_miss > 0.0) ? sqrt(var_miss) : 0.0;
    mc.cep_3sigma_km = mc.mean_miss_km + 3.0 * mc.std_miss_km;

    // Semiassi dell'ellissoide di covarianza 3-sigma nel frame RTN [3*sigma_R, 3*sigma_T, 3*sigma_N]
    mc.ellipsoid_3sigma_rtn.x = 3.0 * sqrt(sum_r2 * inv_n);
    mc.ellipsoid_3sigma_rtn.y = 3.0 * sqrt(sum_t2 * inv_n);
    mc.ellipsoid_3sigma_rtn.z = 3.0 * sqrt(sum_n2 * inv_n);

    mc.tcm_dv_mean_ms = sum_tcm * inv_n;
    double var_tcm = (sum_tcm2 * inv_n) - (mc.tcm_dv_mean_ms * mc.tcm_dv_mean_ms);
    double std_tcm = (var_tcm > 0.0) ? sqrt(var_tcm) : 0.0;
    mc.tcm_dv_99_ms = mc.tcm_dv_mean_ms + 3.0 * std_tcm;

    if (plan.prop_info.isp > 0.0 && plan.final_mass > 0.0) {
        mc.tcm_fuel_99_kg = plan.final_mass * (1.0 - exp(-mc.tcm_dv_99_ms / (plan.prop_info.isp * G0)));
    } else {
        mc.tcm_fuel_99_kg = 0.0;
    }

    return mc;
}
