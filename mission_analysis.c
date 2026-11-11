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
        if (C <= 1e-12) {
            z_high = z;
            continue;
        }
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
    if (C <= 1e-12) return sol;
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
    double r_max_hill = EARTH_RADIUS + MAX_HILL_APOGEE_ALT_KM;

    if (isnan(kep->a) || isinf(kep->a) || kep->a < r_min_safe) {
        kep->a = r_min_safe;
        valid = false;
        if (warning_msg) strcpy(warning_msg, "Semi-major axis raised to minimum safe LEO altitude (150 km).");
    }
    if (kep->a > r_max_hill) {
        kep->a = r_max_hill;
        valid = false;
        if (warning_msg) strcpy(warning_msg, "Semi-major axis clamped to Earth Hill sphere limit (900,000 km).");
    }

    if (isnan(kep->e) || kep->e < 0.0) {
        kep->e = 0.0;
        valid = false;
    }
    if (kep->e >= 0.95) {
        kep->e = 0.9499;
        valid = false;
        if (warning_msg && warning_msg[0] == '\0') {
            strcpy(warning_msg, "Eccentricity clamped below 0.95 (closed-orbit safety guardrail).");
        }
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

    if (kep->a * (1.0 + kep->e) > r_max_hill) {
        double e_max_apogee = (r_max_hill / kep->a) - 1.0;
        if (e_max_apogee < 0.0) e_max_apogee = 0.0;
        if (kep->e > e_max_apogee) kep->e = e_max_apogee;
        valid = false;
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
    double sign_inc = (orb_target.i > orb_init.i || (fabs(orb_target.i - orb_init.i) < 1e-9 && orb_target.RAAN >= orb_init.RAAN)) ? 1.0 : -1.0;

    // ==================================================================
    // CASO A: PROPULSIONE ELETTRICA CONTINUA A SPIRALE (EDELBAUM THEORY)
    // ==================================================================
    if (best.is_low_thrust) {
        if (best.continuous_thrust_n <= 1e-6) best.continuous_thrust_n = 0.50;
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

// =========================================================================
// SPRINT 14: 6-DOF ATTITUDE DYNAMICS (ADCS), QUATERNION SLEW & REACTION WHEELS
// =========================================================================
static double adcs_clamp(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static Quaternion quat_normalize(Quaternion q) {
    double n = sqrt(q.q0*q.q0 + q.q1*q.q1 + q.q2*q.q2 + q.q3*q.q3);
    if (n < 1e-12) {
        Quaternion id = {1.0, 0.0, 0.0, 0.0};
        return id;
    }
    Quaternion res = {q.q0 / n, q.q1 / n, q.q2 / n, q.q3 / n};
    return res;
}

static Quaternion euler_deg_to_quat(double roll_deg, double pitch_deg, double yaw_deg) {
    double phi   = 0.5 * roll_deg  * DEG2RAD;
    double theta = 0.5 * pitch_deg * DEG2RAD;
    double psi   = 0.5 * yaw_deg   * DEG2RAD;

    double cp = cos(phi),   sp = sin(phi);
    double ct = cos(theta), st = sin(theta);
    double cy = cos(psi),   sy = sin(psi);

    Quaternion q = {
        cp * ct * cy + sp * st * sy,
        sp * ct * cy - cp * st * sy,
        cp * st * cy + sp * ct * sy,
        cp * ct * sy - sp * st * cy
    };
    return quat_normalize(q);
}

static void quat_to_euler_deg(Quaternion q, double *roll_deg, double *pitch_deg, double *yaw_deg) {
    q = quat_normalize(q);
    double sinp = adcs_clamp(2.0 * (q.q0 * q.q2 - q.q3 * q.q1), -1.0, 1.0);
    if (roll_deg)  *roll_deg  = atan2(2.0 * (q.q0 * q.q1 + q.q2 * q.q3), 1.0 - 2.0 * (q.q1 * q.q1 + q.q2 * q.q2)) * RAD2DEG;
    if (pitch_deg) *pitch_deg = asin(sinp) * RAD2DEG;
    if (yaw_deg)   *yaw_deg   = atan2(2.0 * (q.q0 * q.q3 + q.q1 * q.q2), 1.0 - 2.0 * (q.q2 * q.q2 + q.q3 * q.q3)) * RAD2DEG;
}

// Errore quaternione q_e = q_cmd^* (x) q_cur con scelta del cammino geodetico minimo (q_e0 >= 0)
static Quaternion quat_error_shortest(Quaternion q_cmd, Quaternion q_cur) {
    Quaternion qe = {
         q_cmd.q0 * q_cur.q0 + q_cmd.q1 * q_cur.q1 + q_cmd.q2 * q_cur.q2 + q_cmd.q3 * q_cur.q3,
         q_cmd.q0 * q_cur.q1 - q_cmd.q1 * q_cur.q0 - q_cmd.q2 * q_cur.q3 + q_cmd.q3 * q_cur.q2,
         q_cmd.q0 * q_cur.q2 + q_cmd.q1 * q_cur.q3 - q_cmd.q2 * q_cur.q0 - q_cmd.q3 * q_cur.q1,
         q_cmd.q0 * q_cur.q3 - q_cmd.q1 * q_cur.q2 + q_cmd.q2 * q_cur.q1 - q_cmd.q3 * q_cur.q0
    };
    if (qe.q0 < 0.0) {
        qe.q0 = -qe.q0;
        qe.q1 = -qe.q1;
        qe.q2 = -qe.q2;
        qe.q3 = -qe.q3;
    }
    return quat_normalize(qe);
}

// Ruota un vettore dal frame orbitale LVLH al frame solidale di corpo (Body) tramite q^* (x) v (x) q
static Vector3D quat_rotate_lvlh_to_body(Quaternion q, Vector3D v) {
    Vector3D u = {-q.q1, -q.q2, -q.q3};
    Vector3D t = vec_scale(vec_cross(u, v), 2.0);
    return vec_add(v, vec_add(vec_scale(t, q.q0), vec_cross(u, t)));
}

// Cinematica differenziale dei quaternioni: dq/dt = 0.5 * Omega(omega) * q
static Quaternion quat_kinematics_rhs(Quaternion q, Vector3D w) {
    Quaternion dq = {
        0.5 * (-w.x * q.q1 - w.y * q.q2 - w.z * q.q3),
        0.5 * ( w.x * q.q0 + w.z * q.q2 - w.y * q.q3),
        0.5 * ( w.y * q.q0 - w.z * q.q1 + w.x * q.q3),
        0.5 * ( w.z * q.q0 + w.y * q.q1 - w.x * q.q2)
    };
    return dq;
}

ADCSReport simulate_adcs_profile(KeplerianElements orb_init,
                                 KeplerianElements orb_target,
                                 OptimalTransferPlan plan,
                                 Vector3D inertia_kgm2,
                                 double rw_h_max_nms,
                                 double mu) {
    ADCSReport rep = {0};

    // Validazione e default del tensore principale d'inerzia J = diag(Jxx, Jyy, Jzz) [kg*m^2]
    if (inertia_kgm2.x < 10.0) inertia_kgm2.x = 120.0;
    if (inertia_kgm2.y < 10.0) inertia_kgm2.y = 110.0;
    if (inertia_kgm2.z < 10.0) inertia_kgm2.z = 85.0;
    if (rw_h_max_nms < 1.0)    rw_h_max_nms   = 10.0;

    rep.inertia_kgm2 = inertia_kgm2;
    rep.rw_h_max_nms = rw_h_max_nms;

    double j_mean = (inertia_kgm2.x + inertia_kgm2.y + inertia_kgm2.z) / 3.0;
    double j_max  = inertia_kgm2.x;
    if (inertia_kgm2.y > j_max) j_max = inertia_kgm2.y;
    if (inertia_kgm2.z > j_max) j_max = inertia_kgm2.z;
    double gamma_j = j_mean / 105.0; // Fattore di scala rispetto a bus SmallSat 120/110/85 kg*m^2

    // Durata complessiva delle 3 fasi della missione per allineare l'asse temporale t_min
    double p_init_s = orbital_period(orb_init.a, mu);
    double p_targ_s = orbital_period(orb_target.a, mu);
    double t_phase1_s = 2.0 * p_init_s;
    double t_phase2_s = plan.is_low_thrust ? (6.0 * 5600.0) : plan.transfer_time_sec;
    if (t_phase2_s < 600.0) t_phase2_s = 2700.0;
    double t_phase3_s = 2.0 * p_targ_s;
    double total_mission_min = (t_phase1_s + t_phase2_s + t_phase3_s) / 60.0;

    // Direzioni di spinta RTN per Burn 1 e Burn 2 -> Angoli di assetto target (Roll, Pitch, Yaw)
    Vector3D u1 = vec_normalize(plan.burn1_dir_rtn);
    if (vec_norm(u1) < 1e-6) u1 = (Vector3D){0.0, 1.0, 0.0};
    Vector3D u2 = vec_normalize(plan.burn2_dir_rtn);
    if (vec_norm(u2) < 1e-6) u2 = (Vector3D){0.0, 1.0, 0.0};

    double psi1_rtn   = atan2(u1.z, u1.y) * RAD2DEG;
    double theta1_rtn = asin(adcs_clamp(u1.x, -1.0, 1.0)) * RAD2DEG;
    double b1_yaw     = psi1_rtn + ((fabs(psi1_rtn) < 12.0) ? 32.0 : 0.0);
    double b1_pitch   = theta1_rtn - 22.0;
    double b1_roll    = -12.0 + 10.0 * u1.z;

    double psi2_rtn   = atan2(u2.z, u2.y) * RAD2DEG;
    double theta2_rtn = asin(adcs_clamp(u2.x, -1.0, 1.0)) * RAD2DEG;
    double b2_yaw     = plan.is_low_thrust ? (-b1_yaw * 0.65) : (psi2_rtn - ((fabs(psi2_rtn - b1_yaw) < 18.0) ? 48.0 : 0.0));
    double b2_pitch   = theta2_rtn + 26.0;
    double b2_roll    = 15.0 - 10.0 * u2.z;

    Quaternion q_nadir  = euler_deg_to_quat(0.0, 0.0, 0.0);
    Quaternion q_sun    = euler_deg_to_quat(14.0, -16.0, -20.0);
    Quaternion q_burn1  = euler_deg_to_quat(b1_roll, b1_pitch, b1_yaw);
    Quaternion q_burn2  = euler_deg_to_quat(b2_roll, b2_pitch, b2_yaw);
    Quaternion q_target = euler_deg_to_quat(8.0, -11.0, 16.0);

    Quaternion qe_slew1 = quat_error_shortest(q_burn1, q_sun);
    Quaternion qe_slew2 = quat_error_shortest(q_burn2, q_burn1);
    rep.slew_angle_b1_deg = 2.0 * acos(adcs_clamp(qe_slew1.q0, 0.0, 1.0)) * RAD2DEG;
    rep.slew_angle_b2_deg = 2.0 * acos(adcs_clamp(qe_slew2.q0, 0.0, 1.0)) * RAD2DEG;

    // Parametri legge di controllo PD sui quaternioni e limiti fisici Reaction Wheels
    double tau_max_axis_nm = 0.115 * pow(gamma_j, 0.25); // Scala coppia max asse RW [N*m]
    if (tau_max_axis_nm > 0.240) tau_max_axis_nm = 0.240;
    double wn = 0.048 / pow(gamma_j, 0.32); // Pulsazione naturale ad anello chiuso [rad/s]
    double zeta = 0.86;                     // Smorzamento critico per evitare overshoot
    Vector3D kp = {
        2.0 * wn * wn * inertia_kgm2.x,
        2.0 * wn * wn * inertia_kgm2.y,
        2.0 * wn * wn * inertia_kgm2.z
    };
    Vector3D kd = {
        2.0 * zeta * wn * inertia_kgm2.x,
        2.0 * zeta * wn * inertia_kgm2.y,
        2.0 * zeta * wn * inertia_kgm2.z
    };

    // Tempo fisico di Pre-Burn Slew (accelerazione + decelerazione + settling fine)
    double dtheta_rad = rep.slew_angle_b1_deg * DEG2RAD;
    rep.slew_time_s = 2.0 * sqrt((dtheta_rad * j_max) / 0.180) + (1.8 / wn);

    // Stima analitica e numerica delle coppie disturbanti ambientali [mN*m]
    double r_mean_km = 0.5 * (orb_init.a + orb_target.a);
    double mu_r3 = mu / (r_mean_km * r_mean_km * r_mean_km); // [s^-2]
    double dj_eff = fabs(inertia_kgm2.x - inertia_kgm2.z);
    if (dj_eff < 0.18 * j_mean) dj_eff = 0.18 * j_mean;
    rep.gg_torque_mnm = 1.5 * mu_r3 * dj_eff * 1000.0; // 3*mu/(2*r^3)*|Jxx-Jzz|*sin(2*theta) in mN*m

    double alt_init_km = orb_init.a - EARTH_RADIUS;
    double rho_init = (alt_init_km < 1000.0) ? (1.0e-11 * exp(-(alt_init_km - 400.0) / 60.0)) : 1.0e-15;
    double v_init_ms = sqrt(mu / orb_init.a) * 1000.0;
    double q_dyn = 0.5 * rho_init * v_init_ms * v_init_ms;
    double l_cp_m = 0.24 * pow(gamma_j, 0.33);
    rep.aero_torque_mnm = q_dyn * 2.2 * 4.0 * l_cp_m * 1000.0;
    if (rep.aero_torque_mnm < 0.02) rep.aero_torque_mnm = 0.02;

    double p_srp = 4.56e-6; // Pressione di radiazione solare a 1 AU [N/m^2]
    double area_srp = 6.5 * pow(gamma_j, 0.33);
    double l_sp_m = 0.35 * pow(gamma_j, 0.33);
    rep.srp_torque_mnm = p_srp * area_srp * 1.4 * l_sp_m * 1000.0;

    // Stato iniziale 6-DOF (Nadir Pointing con piccolo transitorio di acquisizione iniziale)
    Quaternion q_cur = euler_deg_to_quat(2.5, -3.0, 4.0);
    Vector3D omega_rad_s = {0.0, 0.0, 0.0};
    Vector3D h_rw_nms = {0.15, -0.10, 0.12};

    int num_pts = 151; // 151 campioni da i = 0 a i = 150
    rep.num_samples = num_pts;

    double sum_steady_err2 = 0.0;
    int steady_count = 0;
    double sum_power_w = 0.0;
    bool desat_active_prev = false;

    for (int i = 0; i < num_pts; i++) {
        double tau = (double)i / (double)(num_pts - 1);
        double t_min = tau * total_mission_min;

        // Selezione della Pointing Mode e del Quaternione Comandato q_cmd lungo le fasi della missione
        AttitudeMode mode = MODE_NADIR_POINTING;
        Quaternion q_cmd = q_nadir;
        if (tau < 0.18) {
            mode  = MODE_NADIR_POINTING;
            q_cmd = q_nadir;
        } else if (tau < 0.34) {
            mode  = MODE_SUN_POINTING;
            q_cmd = q_sun;
        } else if (tau < 0.56) {
            mode  = MODE_THRUST_VECTOR_ALIGN;
            q_cmd = q_burn1;
        } else if (tau < 0.76) {
            mode  = MODE_THRUST_VECTOR_ALIGN;
            q_cmd = q_burn2;
        } else if (tau < 0.88) {
            mode  = MODE_NADIR_POINTING;
            q_cmd = q_nadir;
        } else {
            mode  = MODE_TARGET_TRACKING;
            q_cmd = q_target;
        }

        // Calcolo stato corrente ed esportazione nel campione i
        double roll_d = 0.0, pitch_d = 0.0, yaw_d = 0.0;
        quat_to_euler_deg(q_cur, &roll_d, &pitch_d, &yaw_d);

        Quaternion qe_now = quat_error_shortest(q_cmd, q_cur);
        double point_err_deg = 2.0 * acos(adcs_clamp(qe_now.q0, 0.0, 1.0)) * RAD2DEG;
        // Jitter residuo sensore stellare (Star Tracker + IMU gyro noise floor ~0.015 - 0.030 deg)
        double star_tracker_jitter = (0.016 + 0.008 * fabs(sin(tau * 37.0))) * pow(gamma_j, 0.12);
        if (point_err_deg < star_tracker_jitter) {
            point_err_deg = star_tracker_jitter;
        }

        // Legge di controllo PD istantanea per valutare la coppia delle Reaction Wheels
        Vector3D j_omega = {
            inertia_kgm2.x * omega_rad_s.x,
            inertia_kgm2.y * omega_rad_s.y,
            inertia_kgm2.z * omega_rad_s.z
        };
        Vector3D h_tot = vec_add(j_omega, h_rw_nms);
        Vector3D gyro_term = vec_cross(omega_rad_s, h_tot);

        Vector3D tau_cmd = {
            -kp.x * qe_now.q1 - kd.x * omega_rad_s.x + gyro_term.x,
            -kp.y * qe_now.q2 - kd.y * omega_rad_s.y + gyro_term.y,
            -kp.z * qe_now.q3 - kd.z * omega_rad_s.z + gyro_term.z
        };
        Vector3D tau_rw = {
            adcs_clamp(-tau_cmd.x, -tau_max_axis_nm, tau_max_axis_nm),
            adcs_clamp(-tau_cmd.y, -tau_max_axis_nm, tau_max_axis_nm),
            adcs_clamp(-tau_cmd.z, -tau_max_axis_nm, tau_max_axis_nm)
        };

        // Aggiunta del ripple di coppia ambientale in regime stazionario
        double env_ripple_mnm = (rep.gg_torque_mnm + rep.aero_torque_mnm + rep.srp_torque_mnm) * (0.85 + 0.30 * sin(tau * 18.0));
        double tau_rw_mnm = vec_norm(tau_rw) * 1000.0 + env_ripple_mnm;

        double h_rw_mag = vec_norm(h_rw_nms);
        double sat_pct  = (h_rw_mag / rw_h_max_nms) * 100.0;
        if (sat_pct > 99.5) sat_pct = 99.5;

        Vector3D omega_deg_s = vec_scale(omega_rad_s, RAD2DEG);
        double omega_norm_deg_s = vec_norm(omega_deg_s);

        // Potenza elettrica assorbita dall'ADCS: avionica base (22 W) + motori Reaction Wheels + Magnetorquers
        double p_adcs_w = 22.0 + 0.28 * tau_rw_mnm + 2.4 * h_rw_mag;

        rep.samples[i].t_min            = t_min;
        rep.samples[i].q                = q_cur;
        rep.samples[i].roll_deg         = roll_d;
        rep.samples[i].pitch_deg        = pitch_d;
        rep.samples[i].yaw_deg          = yaw_d;
        rep.samples[i].omega_deg_s      = omega_deg_s;
        rep.samples[i].omega_norm_deg_s = omega_norm_deg_s;
        rep.samples[i].point_err_deg    = point_err_deg;
        rep.samples[i].rw_torque_mnm    = tau_rw_mnm;
        rep.samples[i].rw_mom_nms       = h_rw_mag;
        rep.samples[i].rw_sat_pct       = sat_pct;
        rep.samples[i].adcs_power_w     = p_adcs_w;
        rep.samples[i].mode             = mode;

        if (tau_rw_mnm > rep.peak_rw_torque_mnm) rep.peak_rw_torque_mnm = tau_rw_mnm;
        if (h_rw_mag   > rep.max_rw_mom_nms)     rep.max_rw_mom_nms     = h_rw_mag;
        if (sat_pct    > rep.max_rw_sat_pct)     rep.max_rw_sat_pct     = sat_pct;
        if (p_adcs_w   > rep.peak_power_w)       rep.peak_power_w       = p_adcs_w;
        sum_power_w += p_adcs_w;

        if (point_err_deg < 0.50) {
            sum_steady_err2 += point_err_deg * point_err_deg;
            steady_count++;
        }

        // Integrazione numerica (8 sotto-passi con dt_sub = 2.5 s) delle equazioni dinamiche di Eulero e cinematica quaternioni
        double dt_sub = 2.5;
        for (int sub = 0; sub < 8; sub++) {
            Quaternion qe = quat_error_shortest(q_cmd, q_cur);

            Vector3D jw = {
                inertia_kgm2.x * omega_rad_s.x,
                inertia_kgm2.y * omega_rad_s.y,
                inertia_kgm2.z * omega_rad_s.z
            };
            Vector3D h_total = vec_add(jw, h_rw_nms);
            Vector3D wx_h    = vec_cross(omega_rad_s, h_total);

            Vector3D t_cmd = {
                -kp.x * qe.q1 - kd.x * omega_rad_s.x + wx_h.x,
                -kp.y * qe.q2 - kd.y * omega_rad_s.y + wx_h.y,
                -kp.z * qe.q3 - kd.z * omega_rad_s.z + wx_h.z
            };
            Vector3D t_rw = {
                adcs_clamp(-t_cmd.x, -tau_max_axis_nm, tau_max_axis_nm),
                adcs_clamp(-t_cmd.y, -tau_max_axis_nm, tau_max_axis_nm),
                adcs_clamp(-t_cmd.z, -tau_max_axis_nm, tau_max_axis_nm)
            };

            // Coppia di gradiente di gravità nel frame di corpo: tau_gg = (3*mu/r^3) * r_b x (J * r_b)
            Vector3D r_lvlh = vec_normalize((Vector3D){0.15 * sin(4.0 * PI * tau), 0.12 * cos(4.0 * PI * tau), 1.0});
            Vector3D r_b    = quat_rotate_lvlh_to_body(q_cur, r_lvlh);
            Vector3D j_rb   = {inertia_kgm2.x * r_b.x, inertia_kgm2.y * r_b.y, inertia_kgm2.z * r_b.z};
            Vector3D t_gg   = vec_scale(vec_cross(r_b, j_rb), 3.0 * mu_r3);

            // Equazione di Eulero: J * d(omega)/dt = tau_dist - tau_rw - omega x (J*omega + h_rw)
            Vector3D net_torque = vec_sub(vec_sub(t_gg, t_rw), wx_h);
            Vector3D d_omega = {
                net_torque.x / inertia_kgm2.x,
                net_torque.y / inertia_kgm2.y,
                net_torque.z / inertia_kgm2.z
            };

            omega_rad_s = vec_add(omega_rad_s, vec_scale(d_omega, dt_sub));

            // Limite fisico sulla velocità angolare di slew (max 2.2 deg/s)
            double max_w_rad = 2.2 * DEG2RAD;
            double w_mag = vec_norm(omega_rad_s);
            if (w_mag > max_w_rad) {
                omega_rad_s = vec_scale(omega_rad_s, max_w_rad / w_mag);
            }

            // Cinematica dei quaternioni: q(t + dt) = normalize(q + 0.5 * Omega(omega)*q * dt)
            Quaternion dq = quat_kinematics_rhs(q_cur, omega_rad_s);
            q_cur.q0 += dq.q0 * dt_sub;
            q_cur.q1 += dq.q1 * dt_sub;
            q_cur.q2 += dq.q2 * dt_sub;
            q_cur.q3 += dq.q3 * dt_sub;
            q_cur = quat_normalize(q_cur);

            // Accumulo di momento angolare sulle Reaction Wheels (coppie secolari ambientali + scambio dinamico di slew)
            double sec_scale = pow(gamma_j, 0.35);
            Vector3D tau_sec = {
                0.00095 * sec_scale * (1.0 + 0.5 * sin(2.0 * PI * tau)),
                0.00130 * sec_scale * ((mode == MODE_THRUST_VECTOR_ALIGN) ? 2.1 : 1.0),
                0.00080 * sec_scale * cos(2.0 * PI * tau)
            };
            h_rw_nms = vec_add(h_rw_nms, vec_scale(vec_add(vec_scale(t_rw, 0.25), tau_sec), dt_sub));

            // Desaturazione tramite Magnetorquers (Momentum Dumping: tau_mtq = m x B_earth)
            double cur_sat = vec_norm(h_rw_nms) / rw_h_max_nms;
            bool desat_window = (tau >= 0.76 && tau <= 0.87) || (cur_sat >= 0.88);
            if (desat_window && cur_sat > 0.22) {
                h_rw_nms = vec_scale(h_rw_nms, 1.0 - 0.018 * dt_sub);
                if (!desat_active_prev && sub == 0) {
                    rep.desat_cycles++;
                    desat_active_prev = true;
                }
            } else if (!desat_window) {
                desat_active_prev = false;
            }
        }
    }

    if (rep.desat_cycles == 0) rep.desat_cycles = 1;
    rep.steady_point_err_deg = (steady_count > 0) ? sqrt(sum_steady_err2 / (double)steady_count) : 0.022;
    rep.steady_point_err_arcsec = rep.steady_point_err_deg * 3600.0;
    rep.mean_power_w = sum_power_w / (double)num_pts;

    return rep;
}

// =========================================================================
// SPRINT 15: THERMAL CONTROL SYSTEM (TCS) 2-NODE MODEL & RADIATION (VAN ALLEN / SAA)
// m*c_p * dT/dt = Q_sun + Q_albedo + Q_IR + Q_int + Q_heater - sigma*eps*A*T^4
// =========================================================================
SurfaceOpticalProps get_surface_optical_props(ThermalCoatingType type) {
    SurfaceOpticalProps p = {0};
    p.type = type;
    switch (type) {
        case COATING_OSR_MIRROR:
            strcpy(p.name, "Optical Solar Reflector (OSR)");
            p.alpha_s = 0.10;
            p.eps_ir  = 0.82;
            break;
        case COATING_MLI_GOLD_KAPTON:
            strcpy(p.name, "MLI Gold Kapton");
            p.alpha_s = 0.35;
            p.eps_ir  = 0.55;
            break;
        case COATING_BLACK_ANODIZED:
            strcpy(p.name, "Black Anodized Al");
            p.alpha_s = 0.88;
            p.eps_ir  = 0.86;
            break;
        case COATING_WHITE_PAINT_SG121:
        default:
            strcpy(p.name, "White Paint SG121");
            p.alpha_s = 0.22;
            p.eps_ir  = 0.86;
            break;
    }
    return p;
}

ThermalRadReport simulate_thermal_and_radiation(KeplerianElements orb_init,
                                                KeplerianElements orb_target,
                                                OptimalTransferPlan plan,
                                                SurfaceOpticalProps coating,
                                                double al_shield_mm,
                                                double t_min_op_c,
                                                double t_max_op_c,
                                                double mu) {
    ThermalRadReport rep = {0};

    double alpha_solar = adcs_clamp(coating.alpha_s, 0.05, 0.98);
    double eps_ir      = adcs_clamp(coating.eps_ir,  0.05, 0.98);
    al_shield_mm       = adcs_clamp(al_shield_mm,    0.50, 30.0);
    if (t_max_op_c <= t_min_op_c + 5.0) {
        t_min_op_c = -10.0;
        t_max_op_c =  45.0;
    }

    coating.alpha_s = alpha_solar;
    coating.eps_ir  = eps_ir;
    if (coating.name[0] == '\0') {
        if (fabs(alpha_solar - 0.10) < 0.04)      strcpy(coating.name, "Optical Solar Reflector (OSR)");
        else if (fabs(alpha_solar - 0.35) < 0.05) strcpy(coating.name, "MLI Gold Kapton");
        else if (alpha_solar > 0.70)              strcpy(coating.name, "Black Anodized Al");
        else                                      strcpy(coating.name, "White Paint SG121");
    }

    rep.coating              = coating;
    rep.alpha_solar          = alpha_solar;
    rep.eps_ir               = eps_ir;
    rep.al_shield_mm         = al_shield_mm;
    rep.thermal_mass_j_k     = 14500.0; // Capacità termica nodo Bus Interno m*c_p [J/K]
    rep.ext_thermal_mass_j_k = 3400.0;  // Capacità termica nodo Pannelli/Superficie Esterna [J/K]
    rep.rad_area_m2          = 3.85;    // Area radiante efficace dei radiatori del Bus [m^2]
    rep.t_min_op_c           = t_min_op_c;
    rep.t_max_op_c           = t_max_op_c;

    double a_proj_m2   = 2.60; // Area efficace esposta al flusso solare diretto [m^2]
    double a_nadir_m2  = 2.30; // Area efficace esposta ad Albedo e IR terrestre [m^2]
    double a_rad_ext   = 3.15; // Area radiante efficace nodo esterno [m^2]
    double g_link_w_k  = 8.50; // Conduttanza termica d'accoppiamento tra nodo esterno e bus [W/K]

    double p_init_s = orbital_period(orb_init.a, mu);
    double p_targ_s = orbital_period(orb_target.a, mu);
    double t_phase1_s = 2.0 * p_init_s;
    double t_phase2_s = plan.is_low_thrust ? (6.0 * 5600.0) : plan.transfer_time_sec;
    if (t_phase2_s < 600.0) t_phase2_s = 2700.0;
    double t_phase3_s = 2.0 * p_targ_s;
    double total_mission_min = (t_phase1_s + t_phase2_s + t_phase3_s) / 60.0;

    // Temperature teoriche di equilibrio stazionario (Hot Case diurno vs Cold Case eclissi)
    double r_mean_km = 0.5 * (orb_init.a + orb_target.a);
    double vf_mean   = (EARTH_RADIUS / r_mean_km) * (EARTH_RADIUS / r_mean_km);
    double q_sol_max = alpha_solar * SOLAR_CONSTANT_W_M2 * a_proj_m2;
    double q_alb_max = alpha_solar * SOLAR_CONSTANT_W_M2 * EARTH_ALBEDO_COEFF * vf_mean * a_nadir_m2 * 0.70;
    double q_ir_mean = eps_ir * EARTH_IR_FLUX_W_M2 * vf_mean * a_nadir_m2;
    double q_int_nom = 165.0;

    double denom_rad_bus = STEFAN_BOLTZMANN * eps_ir * rep.rad_area_m2;
    double denom_rad_ext = STEFAN_BOLTZMANN * eps_ir * a_rad_ext;
    rep.eq_sun_temp_c = pow((q_sol_max + q_alb_max + q_ir_mean + 240.0) / denom_rad_bus, 0.25) - 273.15;
    rep.eq_ecl_temp_c = pow((q_ir_mean + q_int_nom) / denom_rad_bus, 0.25) - 273.15;

    // Fattore di attenuazione schermatura in Alluminio rispetto al riferimento di 4.0 mm
    double f_shield = exp(-0.36 * (al_shield_mm - 4.0)) * pow(4.0 / al_shield_mm, 0.60);

    // Temperature iniziali dei 2 nodi (T_bus = +18.0 °C, T_ext = +21.0 °C)
    double t_bus_k = 291.15;
    double t_ext_k = 294.15;
    rep.min_temp_c     =  1e9;
    rep.max_temp_c     = -1e9;
    rep.ext_min_temp_c =  1e9;
    rep.ext_max_temp_c = -1e9;

    int num_pts = 151;
    rep.num_samples = num_pts;

    double sum_temp_c = 0.0, sum_sol = 0.0, sum_alb = 0.0, sum_ir = 0.0, sum_rad = 0.0;
    double sum_tid_krad_yr = 0.0;
    int heater_active_count = 0;

    // Passo fisico di integrazione termica tra campioni (dt_sample = 65 s diviso in 10 sotto-passi da 6.5 s)
    double dt_sub = 6.5;
    double dt_sample_s = (total_mission_min * 60.0) / (double)(num_pts - 1);

    for (int i = 0; i < num_pts; i++) {
        double tau = (double)i / (double)(num_pts - 1);
        double t_min = tau * total_mission_min;

        // Raggio e inclinazione interpolati lungo le 3 fasi e fattore di vista terrestre F_E = (R_E / r)^2
        double r_cur_km = orb_init.a;
        double inc_cur_rad = orb_init.i;
        if (tau >= 0.34 && tau <= 0.76) {
            double f_tr = (tau - 0.34) / 0.42;
            r_cur_km    = orb_init.a + f_tr * (orb_target.a - orb_init.a);
            inc_cur_rad = orb_init.i + f_tr * (orb_target.i - orb_init.i);
        } else if (tau > 0.76) {
            r_cur_km    = orb_target.a;
            inc_cur_rad = orb_target.i;
        }
        double alt_cur_km = r_cur_km - EARTH_RADIUS;
        if (alt_cur_km < 150.0) alt_cur_km = 150.0;
        double vf_earth = (EARTH_RADIUS / r_cur_km) * (EARTH_RADIUS / r_cur_km);

        // Fase angolare orbitale (4.5 rivoluzioni complete nella finestra di missione) e cono d'ombra (eclissi)
        double orb_angle = 4.5 * 2.0 * PI * tau;
        double cos_sun   = cos(orb_angle);
        int in_eclipse   = (cos_sun < -0.38) ? 1 : 0;

        // 1. Flusso Solare Diretto: Q_sun = (1 - eclipse) * G_s * A_proj * alpha_s
        double q_sun_w = in_eclipse ? 0.0 : (alpha_solar * SOLAR_CONSTANT_W_M2 * a_proj_m2 * (0.88 + 0.12 * fabs(cos_sun)));

        // 2. Flusso di Albedo Terrestre: Q_albedo = (1 - eclipse) * G_s * a_alb * F_E * A_nadir * alpha_s * cos(theta_sun)
        double cos_theta_sun = (cos_sun > 0.0) ? (0.35 + 0.65 * cos_sun) : 0.20;
        double q_alb_w = in_eclipse ? 0.0 : (alpha_solar * SOLAR_CONSTANT_W_M2 * EARTH_ALBEDO_COEFF * vf_earth * a_nadir_m2 * cos_theta_sun);

        // 3. Emissione Infrarossa Terrestre (Earth IR): Q_IR = q_IR * F_E * A_nadir * eps_IR
        double q_ir_w = eps_ir * EARTH_IR_FLUX_W_M2 * vf_earth * a_nadir_m2;
        double q_ir_alb_w = q_alb_w + q_ir_w;

        // 4. Dissipazione Elettronica Interna avionica + propulsione durante i burn / spinta elettrica Q_int [W]
        double q_base_int_w = 210.0;
        if (tau >= 0.34 && tau <= 0.76) {
            q_base_int_w = plan.is_low_thrust ? 360.0 : 290.0;
        }

        // 5. Potenza dei Riscaldatori termostatati di sopravvivenza Q_heater [W]
        double t_bus_c = t_bus_k - 273.15;
        double t_ext_c = t_ext_k - 273.15;
        double q_heater_w = 0.0;
        if (t_bus_c < (t_min_op_c + 12.0)) {
            double frac_h = adcs_clamp(((t_min_op_c + 12.0) - t_bus_c) / 12.0, 0.0, 1.0);
            q_heater_w = 320.0 * frac_h;
        }
        if (q_heater_w > 5.0) heater_active_count++;
        if (q_heater_w > rep.peak_heater_w) rep.peak_heater_w = q_heater_w;
        rep.heater_energy_wh += q_heater_w * (dt_sample_s / 3600.0);

        double q_int_total_w = q_base_int_w + q_heater_w;
        double q_rad_w = denom_rad_bus * (t_bus_k * t_bus_k * t_bus_k * t_bus_k - 81.0);

        // 6. Modello di Radiazione Orbitale (Van Allen Proton/Electron Belts + South Atlantic Anomaly SAA)
        double alt_cap_leo = (alt_cur_km < 1200.0) ? alt_cur_km : 1200.0;
        double d_leo   = 0.72 * exp((alt_cap_leo - 400.0) / 310.0);
        double d_inner = (alt_cur_km > 750.0) ? (195.0 * exp(-pow((alt_cur_km - 3200.0) / 1050.0, 2.0))) : 0.0;
        double d_outer = (alt_cur_km > 7000.0) ? (135.0 * exp(-pow((alt_cur_km - 19000.0) / 5200.0, 2.0))) : 0.0;
        double s_inc   = fabs(sin(inc_cur_rad));
        double f_inc   = 0.78 + 0.42 * s_inc + 0.18 * (s_inc * s_inc);

        double tid_local_krad_yr = (d_leo + d_inner + d_outer) * f_inc * f_shield;
        sum_tid_krad_yr += tid_local_krad_yr;

        // Traccia sub-satellite e passaggio nella South Atlantic Anomaly (SAA: Lat ~ -25 deg, Lon ~ -45 deg)
        double lat_deg = (inc_cur_rad * RAD2DEG) * sin(orb_angle - 0.35);
        double lon_raw = -15.0 + (orb_angle * 0.55 - 2.2 * PI * tau) * RAD2DEG;
        double lon_deg = fmod(lon_raw + 540.0, 360.0) - 180.0;

        double dlat_saa = lat_deg - (-25.0);
        double dlon_saa = fmod((lon_deg - (-45.0)) + 540.0, 360.0) - 180.0;
        double saa_gauss = exp(-(dlat_saa * dlat_saa) / (2.0 * 18.0 * 18.0) - (dlon_saa * dlon_saa) / (2.0 * 34.0 * 34.0));

        // Finestre periodiche di attraversamento del core protonico SAA nell'emisfero australe per inclinazioni >= 18 deg
        double saa_pass_wave = 0.0;
        if ((inc_cur_rad * RAD2DEG) >= 18.0 && alt_cur_km <= 3500.0) {
            double p1 = exp(-pow((tau - 0.19) / 0.032, 2.0));
            double p2 = exp(-pow((tau - 0.62) / 0.035, 2.0));
            double p3 = exp(-pow((tau - 0.86) / 0.030, 2.0));
            saa_pass_wave = p1 + 1.15 * p2 + 0.95 * p3;
        }
        double saa_intensity = (saa_gauss > saa_pass_wave) ? saa_gauss : saa_pass_wave;
        int in_saa = (saa_intensity > 0.35) ? 1 : 0;
        if (in_saa) rep.saa_encountered = true;

        // Rateo di dose istantanea [urad/h] con picchi nell'attraversamento SAA e corni polari
        double polar_horn = (fabs(lat_deg) > 45.0) ? (0.35 * pow((fabs(lat_deg) - 45.0) / 45.0, 2.0)) : 0.0;
        double inst_dose_urad_h = tid_local_krad_yr * 114.08 * (0.65 + 4.80 * saa_intensity + polar_horn);
        if (inst_dose_urad_h > rep.peak_dose_urad_h) {
            rep.peak_dose_urad_h = inst_dose_urad_h;
        }

        rep.samples[i].t_min            = t_min;
        rep.samples[i].t_ext_c          = t_ext_c;
        rep.samples[i].t_bus_c          = t_bus_c;
        rep.samples[i].temp_c           = t_bus_c;
        rep.samples[i].temp_k           = t_bus_k;
        rep.samples[i].q_sun_w          = q_sun_w;
        rep.samples[i].q_sol_w          = q_sun_w;
        rep.samples[i].q_alb_w          = q_alb_w;
        rep.samples[i].q_ir_w           = q_ir_w;
        rep.samples[i].q_ir_alb_w       = q_ir_alb_w;
        rep.samples[i].q_int_w          = q_int_total_w;
        rep.samples[i].q_heater_w       = q_heater_w;
        rep.samples[i].q_rad_w          = q_rad_w;
        rep.samples[i].inst_dose_urad_h = inst_dose_urad_h;
        rep.samples[i].in_eclipse       = in_eclipse;
        rep.samples[i].in_saa           = in_saa;

        if (t_bus_c < rep.min_temp_c)     rep.min_temp_c     = t_bus_c;
        if (t_bus_c > rep.max_temp_c)     rep.max_temp_c     = t_bus_c;
        if (t_ext_c < rep.ext_min_temp_c) rep.ext_min_temp_c = t_ext_c;
        if (t_ext_c > rep.ext_max_temp_c) rep.ext_max_temp_c = t_ext_c;

        sum_temp_c += t_bus_c;
        sum_sol    += q_sun_w;
        sum_alb    += q_alb_w;
        sum_ir     += q_ir_w;
        sum_rad    += q_rad_w;

        // Integrazione numerica accoppiata dei 2 nodi (Pannelli Esterni T_ext e Bus Interno T_bus)
        for (int sub = 0; sub < 10; sub++) {
            double tc_bus_sub = t_bus_k - 273.15;
            double q_h_sub = 0.0;
            if (tc_bus_sub < (t_min_op_c + 12.0)) {
                q_h_sub = 320.0 * adcs_clamp(((t_min_op_c + 12.0) - tc_bus_sub) / 12.0, 0.0, 1.0);
            }

            // Nodo 1: Superficie Esterna / Pannelli Solari (bassa inerzia termica, forte escursione Sole/Eclissi)
            double q_cond_bus_to_ext = g_link_w_k * (t_bus_k - t_ext_k);
            double q_in_ext  = 1.18 * q_sun_w + q_alb_w + q_ir_w + q_cond_bus_to_ext;
            double q_out_ext = denom_rad_ext * (t_ext_k * t_ext_k * t_ext_k * t_ext_k - 81.0);
            double dt_ext_k  = ((q_in_ext - q_out_ext) / rep.ext_thermal_mass_j_k) * dt_sub;

            // Nodo 2: Bus Interno / Batterie Li-Ion (alta inerzia termica, dissipazione avionica + survival heaters)
            double q_in_bus  = q_sun_w + q_alb_w + q_ir_w + q_base_int_w + q_h_sub - 0.22 * q_cond_bus_to_ext;
            double q_out_bus = denom_rad_bus * (t_bus_k * t_bus_k * t_bus_k * t_bus_k - 81.0);
            double dt_bus_k  = ((q_in_bus - q_out_bus) / rep.thermal_mass_j_k) * dt_sub;

            t_ext_k = adcs_clamp(t_ext_k + dt_ext_k, 150.0, 500.0);
            t_bus_k = adcs_clamp(t_bus_k + dt_bus_k, 150.0, 500.0);
        }
    }

    double inv_n = 1.0 / (double)num_pts;
    rep.mean_temp_c     = sum_temp_c * inv_n;
    rep.mean_q_sol_w    = sum_sol * inv_n;
    rep.mean_q_alb_w    = sum_alb * inv_n;
    rep.mean_q_ir_w     = sum_ir  * inv_n;
    rep.mean_q_rad_w    = sum_rad * inv_n;
    rep.heater_duty_pct = ((double)heater_active_count / (double)num_pts) * 100.0;

    // Metriche annuali di radiazione ionizzante (TID, SEU/day e Vita Utile Elettronica)
    rep.tid_krad_yr = sum_tid_krad_yr * inv_n * (rep.saa_encountered ? 1.28 : 1.0);
    if (rep.tid_krad_yr < 0.05) rep.tid_krad_yr = 0.05;
    rep.seu_rate_day = (0.42 * rep.tid_krad_yr + (rep.saa_encountered ? 0.65 : 0.15)) * pow(f_shield, 0.45);
    rep.avionics_life_yr = 25.0 / rep.tid_krad_yr; // Tolleranza avionica Rad-Tolerant = 25 krad(Si)
    if (rep.avionics_life_yr > 25.0) rep.avionics_life_yr = 25.0;
    if (rep.avionics_life_yr < 0.2)  rep.avionics_life_yr = 0.2;

    double h_mean_km = 0.5 * ((orb_init.a - EARTH_RADIUS) + (orb_target.a - EARTH_RADIUS));
    if (h_mean_km >= 1000.0 && h_mean_km <= 6500.0) {
        strcpy(rep.saa_status, "INNER VAN ALLEN PROTON BELT + SAA CORE");
    } else if (h_mean_km > 6500.0) {
        strcpy(rep.saa_status, "OUTER VAN ALLEN ELECTRON BELT TRANSIT");
    } else if (rep.saa_encountered) {
        strcpy(rep.saa_status, "SAA PROTON FLUX ENCOUNTER (-25S, -45W)");
    } else {
        strcpy(rep.saa_status, "NOMINAL LOW-INCLINATION LEO FLUX");
    }

    return rep;
}

ThermalReport simulate_thermal_profile(KeplerianElements orb_init,
                                       KeplerianElements orb_target,
                                       OptimalTransferPlan plan,
                                       double alpha_solar,
                                       double eps_ir,
                                       double t_min_op_c,
                                       double t_max_op_c,
                                       double mu) {
    SurfaceOpticalProps coating = get_surface_optical_props(COATING_WHITE_PAINT_SG121);
    coating.alpha_s = alpha_solar;
    coating.eps_ir  = eps_ir;
    return simulate_thermal_and_radiation(orb_init, orb_target, plan, coating, 4.0, t_min_op_c, t_max_op_c, mu);
}

// ============================================================================
// SPRINT 16: RF TELECOMMUNICATION LINK BUDGET, Eb/N0 MARGIN & DATA DOWNLINK
// ============================================================================
TelecomReport evaluate_rf_link_budget(KeplerianElements orb_init,
                                      KeplerianElements orb_target,
                                      OptimalTransferPlan plan,
                                      OrbitReport rep_target,
                                      RFBandType band,
                                      double tx_power_w,
                                      const GroundStation *station,
                                      double mu) {
    TelecomReport rep = {0};
    rep.band = band;
    rep.tx_power_w = adcs_clamp(tx_power_w, 0.5, 250.0);
    rep.tx_power_dbw = 10.0 * log10(rep.tx_power_w);
    rep.rx_dish_m = 3.5; // Parabola da 3.5m presso la Ground Station di Milano-Bovisa
    rep.lock_threshold_db = 4.5; // Soglia di decodifica QPSK/8PSK + LDPC (BER < 10^-6)
    rep.req_ebn0_db = 7.5;       // 4.5 dB Lock Threshold + 3.0 dB Safety Margin

    switch (band) {
        case BAND_S_2_2GHZ:
            strcpy(rep.band_name, "S-Band (2.2 GHz TT&C)");
            strcpy(rep.modulation_scheme, "QPSK + LDPC (1/2 - 3/4) ACM");
            rep.freq_ghz           = 2.20;
            rep.nom_bitrate_mbps   = 2.0;
            rep.tx_gain_dbi        = 6.5;   // Antenna TT&C emisferica/patch a largo fascio
            rep.rx_gain_dbi        = 36.2;  // Parabola 3.5m a 2.2 GHz (eff = 58%)
            rep.sys_noise_temp_k   = 135.0;
            rep.zenith_atm_loss_db = 0.08;
            rep.pointing_loss_db   = 0.25;
            break;

        case BAND_KA_26GHZ:
            strcpy(rep.band_name, "Ka-Band (26.5 GHz High-Rate)");
            strcpy(rep.modulation_scheme, "16APSK / 32APSK + DVB-S2X ACM");
            rep.freq_ghz           = 26.50;
            rep.nom_bitrate_mbps   = 250.0;
            rep.tx_gain_dbi        = 25.0;  // Antenna direttiva Ka-band HGA gimballed
            rep.rx_gain_dbi        = 57.6;  // Parabola 3.5m a 26.5 GHz
            rep.sys_noise_temp_k   = 255.0;
            rep.zenith_atm_loss_db = 0.78;  // Sensibile a vapore acqueo/pioggia a bassa elevazione
            rep.pointing_loss_db   = 0.75;
            break;

        case BAND_X_8_4GHZ:
        default:
            rep.band = BAND_X_8_4GHZ;
            strcpy(rep.band_name, "X-Band (8.4 GHz Science)");
            strcpy(rep.modulation_scheme, "8PSK / 16APSK + LDPC ACM");
            rep.freq_ghz           = 8.40;
            rep.nom_bitrate_mbps   = 50.0;
            rep.tx_gain_dbi        = 15.5;  // Antenna X-Band Horn / Phased Array
            rep.rx_gain_dbi        = 47.8;  // Parabola 3.5m a 8.4 GHz
            rep.sys_noise_temp_k   = 165.0;
            rep.zenith_atm_loss_db = 0.22;
            rep.pointing_loss_db   = 0.45;
            break;
    }

    double line_loss_db = 0.80; // Perdite di cavo/diplexer RF di bordo
    rep.eirp_dbw  = rep.tx_power_dbw + rep.tx_gain_dbi - line_loss_db;
    rep.n0_dbw_hz = BOLTZMANN_DBW_K_HZ + 10.0 * log10(rep.sys_noise_temp_k);

    double a_init = (orb_init.a > EARTH_RADIUS + 100.0) ? orb_init.a : (EARTH_RADIUS + 400.0);
    double a_targ = (orb_target.a > EARTH_RADIUS + 100.0) ? orb_target.a : (EARTH_RADIUS + 650.0);
    double p_init_s = 2.0 * PI * sqrt((a_init * a_init * a_init) / mu);
    double p_targ_s = 2.0 * PI * sqrt((a_targ * a_targ * a_targ) / mu);
    double t_phase1_s = 2.0 * p_init_s;
    double t_phase2_s = plan.is_low_thrust ? (6.0 * 5600.0) : ((plan.transfer_time_sec > 600.0) ? plan.transfer_time_sec : 2700.0);
    double t_phase3_s = 2.0 * p_targ_s;
    double total_mission_min = (t_phase1_s + t_phase2_s + t_phase3_s) / 60.0;

    double min_el_deg = (station != NULL && station->min_elevation > 0.0)
                        ? (station->min_elevation * RAD2DEG)
                        : 10.0;
    double st_lat_deg = (station != NULL) ? (station->coord.latitude * RAD2DEG) : 45.502;

    // Centri temporali normalizzati dei 4 passaggi radar sopra Milano-Bovisa lungo la missione
    const double pass_centers[4] = {0.14, 0.41, 0.67, 0.89};
    const double pass_half_w[4]  = {0.055, 0.062, 0.058, 0.065};
    const double pass_max_el_base[4] = {52.0, 68.0, 44.0, 74.0};

    double inc_deg = orb_target.i * RAD2DEG;
    double inc_factor = 1.0 - 0.22 * fabs(inc_deg - st_lat_deg) / 90.0;
    inc_factor = adcs_clamp(inc_factor, 0.68, 1.08);

    int num_pts = 151;
    rep.num_samples = num_pts;
    double dt_sample_s = (total_mission_min * 60.0) / (double)(num_pts - 1);

    rep.mean_fspl_db      = 0.0;
    rep.max_fspl_db       = 0.0;
    rep.peak_ebn0_db      = -100.0;
    rep.min_margin_db     = 100.0;
    rep.peak_margin_db    = -100.0;
    rep.max_achieved_mbps = 0.0;
    rep.pass_windows      = 0;
    rep.contact_time_min  = 0.0;

    double sum_fspl_pass = 0.0;
    int    pass_sample_count = 0;
    double cum_data_gb = 0.0;
    bool   in_pass_prev = false;
    double log10_nom_bps = 10.0 * log10(rep.nom_bitrate_mbps * 1.0e6);

    for (int i = 0; i < num_pts; i++) {
        double tau   = (double)i / (double)(num_pts - 1);
        double t_min = tau * total_mission_min;

        double r_cur_km = a_init;
        if (tau >= 0.34 && tau <= 0.76) {
            double f_tr = (tau - 0.34) / 0.42;
            r_cur_km = a_init + f_tr * (a_targ - a_init);
        } else if (tau > 0.76) {
            r_cur_km = a_targ;
        }
        double alt_cur_km = r_cur_km - EARTH_RADIUS;
        if (alt_cur_km < 150.0) alt_cur_km = 150.0;

        // Calcolo dell'angolo di elevazione El(t) rispetto alla Ground Station
        double el_deg = -18.0 + 6.0 * sin(2.0 * PI * 4.5 * tau);
        for (int p = 0; p < 4; p++) {
            double u = (tau - pass_centers[p]) / pass_half_w[p];
            if (fabs(u) < 1.35) {
                double arch = pass_max_el_base[p] * inc_factor * cos(0.5 * PI * (u / 1.05));
                if (u * u > 1.05 * 1.05) {
                    arch = -12.0 * (fabs(u) - 1.05) / 0.30;
                }
                if (arch > el_deg) el_deg = arch;
            }
        }

        // Distanza obliqua (Slant Range) rho(t) [km] dalla geometria sferica Terra-Satellite
        double slant_range_km = 0.0;
        if (el_deg > 0.0) {
            double sin_el = sin(el_deg * DEG2RAD);
            slant_range_km = -EARTH_RADIUS * sin_el +
                             sqrt(EARTH_RADIUS * EARTH_RADIUS * sin_el * sin_el +
                                  r_cur_km * r_cur_km - EARTH_RADIUS * EARTH_RADIUS);
        } else {
            double horiz_km = sqrt(fmax(100.0, r_cur_km * r_cur_km - EARTH_RADIUS * EARTH_RADIUS));
            slant_range_km = horiz_km * (1.0 + 0.12 * fabs(el_deg) / 20.0);
        }
        if (slant_range_km < alt_cur_km) slant_range_km = alt_cur_km;

        // Free-Space Path Loss: FSPL(t) = 20*log10(rho_km) + 20*log10(f_GHz) + 92.45 [dB]
        double fspl_db = 20.0 * log10(slant_range_km) + 20.0 * log10(rep.freq_ghz) + 92.45;
        if (fspl_db > rep.max_fspl_db) rep.max_fspl_db = fspl_db;

        // Attenuazione atmosferica dipendente dall'elevazione: L_atm(t) = L_zenith / sin(max(El, 5 deg)) + L_point
        double el_clamped_deg = (el_deg > 5.0) ? el_deg : 5.0;
        double atm_loss_db = (rep.zenith_atm_loss_db / sin(el_clamped_deg * DEG2RAD)) + rep.pointing_loss_db;
        if (el_deg < min_el_deg) {
            // Forte attenuazione di multippercorso/ostruzione orizzonte sotto l'elevazione minima di stazione
            double below_deg = min_el_deg - el_deg;
            atm_loss_db += 0.55 * below_deg;
        }

        // Potenza ricevuta P_rx(t) [dBW], C/N0(t) [dB-Hz], Eb/N0(t) [dB] e Link Margin(t) [dB]
        double prx_dbw  = rep.eirp_dbw - fspl_db - atm_loss_db + rep.rx_gain_dbi;
        double cn0_dbhz = prx_dbw - rep.n0_dbw_hz;
        double ebn0_db  = cn0_dbhz - log10_nom_bps;
        double link_margin_db = ebn0_db - rep.req_ebn0_db;

        // Adaptive Coding & Modulation (ACM) e integrazione del volume dati durante i passaggi validi
        double achieved_mbps = 0.0;
        int in_pass = 0;
        if (el_deg >= min_el_deg && link_margin_db > 0.0) {
            in_pass = 1;
            if (!in_pass_prev) {
                rep.pass_windows++;
                in_pass_prev = true;
            }
            rep.contact_time_min += dt_sample_s / 60.0;
            sum_fspl_pass += fspl_db;
            pass_sample_count++;

            if (link_margin_db >= 9.0) {
                achieved_mbps = rep.nom_bitrate_mbps * 1.60; // 32APSK 4/5 High-Order ACM
            } else if (link_margin_db >= 5.0) {
                achieved_mbps = rep.nom_bitrate_mbps * 1.25; // 16APSK 3/4 ACM
            } else if (link_margin_db >= 1.5) {
                achieved_mbps = rep.nom_bitrate_mbps * 1.00; // 8PSK 2/3 Nominal Rate
            } else {
                achieved_mbps = rep.nom_bitrate_mbps * 0.65; // QPSK 1/2 Robust Rate
            }

            // Conversione Mbps * dt_sec in Gigabyte (1 Byte = 8 bit, 1 GB = 8000 Mbit)
            cum_data_gb += (achieved_mbps * dt_sample_s) / 8000.0;

            if (achieved_mbps > rep.max_achieved_mbps) rep.max_achieved_mbps = achieved_mbps;
            if (link_margin_db < rep.min_margin_db)    rep.min_margin_db     = link_margin_db;
        } else {
            in_pass_prev = false;
        }

        if (ebn0_db > rep.peak_ebn0_db)          rep.peak_ebn0_db   = ebn0_db;
        if (link_margin_db > rep.peak_margin_db) rep.peak_margin_db = link_margin_db;

        rep.samples[i].t_min          = t_min;
        rep.samples[i].slant_range_km = slant_range_km;
        rep.samples[i].elevation_deg  = el_deg;
        rep.samples[i].fspl_db        = fspl_db;
        rep.samples[i].atm_loss_db    = atm_loss_db;
        rep.samples[i].prx_dbw        = prx_dbw;
        rep.samples[i].cn0_dbhz       = cn0_dbhz;
        rep.samples[i].ebn0_db        = ebn0_db;
        rep.samples[i].link_margin_db = link_margin_db;
        rep.samples[i].achieved_mbps  = achieved_mbps;
        rep.samples[i].cum_data_gb    = cum_data_gb;
        rep.samples[i].in_pass        = in_pass;
    }

    if (pass_sample_count > 0) {
        rep.mean_fspl_db = sum_fspl_pass / (double)pass_sample_count;
    } else {
        rep.mean_fspl_db = rep.max_fspl_db - 4.5;
        rep.min_margin_db = rep.peak_margin_db;
    }

    rep.mission_data_gb = cum_data_gb;

    // Proiezione volume dati scaricato sulle 24 ore (combinando i passaggi della missione e il contatto giornaliero)
    double daily_contact_min = (rep_target.total_contact_min > 18.0) ? rep_target.total_contact_min : 38.0;
    if (rep.contact_time_min > 1.0 && cum_data_gb > 0.0) {
        double mean_rate_gb_per_min = cum_data_gb / rep.contact_time_min;
        rep.total_data_gb_24h = fmax(cum_data_gb, mean_rate_gb_per_min * daily_contact_min);
    } else {
        rep.total_data_gb_24h = 0.0;
    }

    return rep;
}

// ============================================================================
// SPRINT 17: PAYLOAD SENSOR SWATH, WALKER DELTA CONSTELLATION & MASS/COST MODEL
// ============================================================================
SystemEngineeringReport evaluate_payload_constellation_and_cost(KeplerianElements orb_target,
                                                                OrbitReport rep_target,
                                                                OptimalTransferPlan plan,
                                                                RendezvousPlan rdv,
                                                                MonteCarloReport mc,
                                                                PayloadAndConstellationSpec spec,
                                                                const GroundStation *station,
                                                                double mu) {
    SystemEngineeringReport sys = {0};

    spec.half_fov_deg      = adcs_clamp((spec.half_fov_deg > 0.5) ? spec.half_fov_deg : 25.0, 1.0, 68.0);
    spec.opt_aperture_m    = adcs_clamp((spec.opt_aperture_m > 0.02) ? spec.opt_aperture_m : 0.35, 0.05, 4.0);
    spec.wavelength_nm     = adcs_clamp((spec.wavelength_nm > 100.0) ? spec.wavelength_nm : 550.0, 250.0, 15000.0);
    spec.walker_total_sats = (spec.walker_total_sats >= 1) ? spec.walker_total_sats : 24;
    if (spec.walker_total_sats > MAX_WALKER_SATS) spec.walker_total_sats = MAX_WALKER_SATS;
    spec.walker_planes     = (spec.walker_planes >= 1) ? spec.walker_planes : 3;
    if (spec.walker_planes > spec.walker_total_sats) spec.walker_planes = spec.walker_total_sats;
    spec.walker_phasing    = (spec.walker_phasing >= 0) ? (spec.walker_phasing % spec.walker_planes) : 1;
    sys.spec = spec;

    // 1. Geometria Sensore Ottico/SAR, Swath Width a Terra [km] e Ground Sample Distance [m]
    double a_targ_km = (orb_target.a > EARTH_RADIUS + 100.0) ? orb_target.a : (EARTH_RADIUS + 650.0);
    double h_km      = fmax(150.0, a_targ_km - EARTH_RADIUS);
    double r_km      = EARTH_RADIUS + h_km;
    sys.target_alt_km = h_km;

    double fov_half_rad = spec.half_fov_deg * DEG2RAD;
    double sin_arg = (r_km / EARTH_RADIUS) * sin(fov_half_rad);
    double lambda_rad = 0.0;
    if (sin_arg < 0.9998) {
        lambda_rad = asin(sin_arg) - fov_half_rad;
    } else {
        // Limite geometrico dell'orizzonte terrestre
        lambda_rad = acos(EARTH_RADIUS / r_km);
    }
    if (lambda_rad < 1e-5) lambda_rad = 1e-5;

    sys.earth_central_angle_deg = lambda_rad * RAD2DEG;
    sys.swath_width_km          = 2.0 * EARTH_RADIUS * lambda_rad;
    sys.footprint_radius_km     = 0.5 * sys.swath_width_km;

    // Risoluzione spaziale al suolo (Limite di diffrazione di Rayleigh: GSD = 1.22 * lambda * h / D_opt)
    sys.gsd_nadir_m = 1.22 * ((spec.wavelength_nm * 1.0e-9) * (h_km * 1000.0)) / spec.opt_aperture_m;

    double v_orb_kms = sqrt(mu / r_km);
    double v_gnd_kms = v_orb_kms * (EARTH_RADIUS / r_km);
    sys.area_coverage_rate_km2_s = sys.swath_width_km * v_gnd_kms;

    // 2. Costellazione Simmetrica di Walker Delta (i : T / P / F)
    int T = spec.walker_total_sats;
    int P = spec.walker_planes;
    int F = spec.walker_phasing;
    sys.inc_deg        = orb_target.i * RAD2DEG;
    sys.sats_per_plane = T / P;
    if (sys.sats_per_plane < 1) sys.sats_per_plane = 1;

    double raan0_deg = orb_target.RAAN * RAD2DEG;
    double inc_rad   = orb_target.i;
    int sat_idx = 0;

    for (int p = 0; p < P && sat_idx < MAX_WALKER_SATS; p++) {
        int sats_in_p = (T / P) + ((p < (T % P)) ? 1 : 0);
        double raan_p_deg = fmod(raan0_deg + (double)p * (360.0 / (double)P), 360.0);
        if (raan_p_deg < 0.0) raan_p_deg += 360.0;

        for (int m = 0; m < sats_in_p && sat_idx < MAX_WALKER_SATS; m++) {
            double u_deg = fmod((double)m * (360.0 / (double)sats_in_p) +
                                (double)(p * F) * (360.0 / (double)T), 360.0);
            double u_rad = u_deg * DEG2RAD;

            double sin_lat = sin(inc_rad) * sin(u_rad);
            sin_lat = adcs_clamp(sin_lat, -1.0, 1.0);
            double lat_deg = asin(sin_lat) * RAD2DEG;

            double dlon_rad = atan2(cos(inc_rad) * sin(u_rad), cos(u_rad));
            double lon_raw  = raan_p_deg + dlon_rad * RAD2DEG;
            double lon_deg  = fmod(lon_raw + 540.0, 360.0) - 180.0;

            sys.walker_subpoints[sat_idx].lat_deg     = lat_deg;
            sys.walker_subpoints[sat_idx].lon_deg     = lon_deg;
            sys.walker_subpoints[sat_idx].plane_id    = p;
            sys.walker_subpoints[sat_idx].sat_id      = sat_idx;
            sys.walker_subpoints[sat_idx].raan_deg    = raan_p_deg;
            sys.walker_subpoints[sat_idx].arg_lat_deg = u_deg;
            sat_idx++;
        }
    }
    sys.num_walker_sats = sat_idx;

    // Copertura globale istantanea (calotta sferica di semi-apertura Lambda per ciascun satellite)
    double f_cap_single = 0.5 * (1.0 - cos(lambda_rad));
    sys.global_coverage_pct = 100.0 * (1.0 - pow(1.0 - f_cap_single, (double)T));
    if (sys.global_coverage_pct > 99.9) sys.global_coverage_pct = 99.9;

    // Numero medio di satelliti in vista sopra El_min e PDOP
    double min_el_rad = (station != NULL && station->min_elevation > 0.0) ? station->min_elevation : (10.0 * DEG2RAD);
    double st_lat_deg = (station != NULL) ? (station->coord.latitude * RAD2DEG) : 45.502;
    double cos_el = cos(min_el_rad);
    double arg_acc = adcs_clamp((EARTH_RADIUS / r_km) * cos_el, -1.0, 1.0);
    double lambda_acc_rad = acos(arg_acc) - min_el_rad;
    if (lambda_acc_rad < lambda_rad) lambda_acc_rad = lambda_rad;
    double f_acc_single = 0.5 * (1.0 - cos(lambda_acc_rad));

    double inc_match = 1.05 - 0.30 * fabs(sys.inc_deg - fabs(st_lat_deg)) / 90.0;
    inc_match = adcs_clamp(inc_match, 0.65, 1.20);
    sys.mean_vis_sats = (double)T * f_acc_single * inc_match * 1.35;
    sys.pdop_geom     = 4.60 / sqrt(fmax(0.85, sys.mean_vis_sats + 0.35 * (double)P));

    // Mean Revisit Time [min] sulla latitudine della stazione di terra
    double period_min = (rep_target.period_min > 60.0) ? rep_target.period_min : (2.0 * PI * sqrt(r_km * r_km * r_km / mu) / 60.0);
    double lat_circ_km = 2.0 * PI * EARTH_RADIUS * fmax(0.25, cos(st_lat_deg * DEG2RAD));
    double swath_eff_km = fmax(80.0, sys.swath_width_km);
    double revs_for_full_band = lat_circ_km / swath_eff_km;
    double plane_diversity = 0.72 + 0.10 * fmin(8.0, (double)P);
    sys.mean_revisit_min = (period_min * revs_for_full_band) / ((double)T * plane_diversity * inc_match);
    sys.mean_revisit_min = adcs_clamp(sys.mean_revisit_min, 1.5, 1440.0);

    // 3. Bilancio di Massa (SMAD Dry Mass Breakdown + Propellant Margin)
    double wet_mass_kg = plan.final_mass + plan.propellant_used;
    if (wet_mass_kg < 20.0) wet_mass_kg = 500.0;
    sys.wet_mass_kg            = wet_mass_kg;
    sys.propellant_transfer_kg = plan.propellant_used;
    sys.propellant_margin_kg   = rdv.propellant_used_kg + mc.tcm_fuel_99_kg + rep_target.station_keeping_fuel_yr;
    sys.total_propellant_kg    = sys.propellant_transfer_kg + sys.propellant_margin_kg;
    if (sys.total_propellant_kg > 0.82 * wet_mass_kg) {
        sys.total_propellant_kg = 0.82 * wet_mass_kg;
    }
    sys.dry_mass_kg       = wet_mass_kg - sys.total_propellant_kg;
    sys.dry_wet_ratio_pct = (sys.dry_mass_kg / wet_mass_kg) * 100.0;

    // Ripartizione SMAD della Massa a Secco (28%, 22%, 19%, 11%, 9%, 6%, 5%)
    sys.mass_payload_kg  = 0.28 * sys.dry_mass_kg;
    sys.mass_struct_kg   = 0.22 * sys.dry_mass_kg;
    sys.mass_eps_kg      = 0.19 * sys.dry_mass_kg;
    sys.mass_prop_dry_kg = 0.11 * sys.dry_mass_kg;
    sys.mass_adcs_kg     = 0.09 * sys.dry_mass_kg;
    sys.mass_ttc_obdh_kg = 0.06 * sys.dry_mass_kg;
    sys.mass_tcs_kg      = 0.05 * sys.dry_mass_kg;

    // 4. Stima Parametrica dei Costi di Missione (SMAD / NASA NICM in M€)
    double opt_factor = pow(spec.opt_aperture_m / 0.35, 1.15);
    sys.cost_platform_meur      = 6.50 + 0.082 * pow(sys.dry_mass_kg, 0.85);
    sys.cost_payload_meur       = (4.20 + 0.145 * pow(sys.mass_payload_kg, 0.88)) * opt_factor;
    sys.cost_spacecraft_meur    = sys.cost_platform_meur + sys.cost_payload_meur;

    double launch_rate_meur_kg  = (h_km > 1500.0) ? 0.0120 : 0.0085; // Tariffa specifica LEO / MEO-GTO [M€/kg]
    sys.cost_launch_meur        = wet_mass_kg * launch_rate_meur_kg;
    sys.cost_ops_5yr_meur       = 5.0 * (1.35 + 0.035 * sys.cost_spacecraft_meur);
    sys.total_mission_cost_meur = sys.cost_spacecraft_meur + sys.cost_launch_meur + sys.cost_ops_5yr_meur;

    // Costo della costellazione completa di T satelliti con fattore di apprendimento industriale di Wright
    sys.constellation_cost_meur = sys.cost_spacecraft_meur * pow((double)T, 0.78) +
                                  (double)T * sys.cost_launch_meur * 0.72 +
                                  sys.cost_ops_5yr_meur * (1.0 + 0.08 * (double)T);

    return sys;
}
