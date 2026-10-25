#include "mission_analysis.h"
#include "astrodynamics.h"
#include "propagator.h"
#include "telemetry.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// =========================================================================
// FUNZIONI DI STUMPFF C(z) E S(z) PER VARIABILI UNIVERSALI DI LAMBERT
// =========================================================================
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

// Converte un vettore Delta-V dal sistema inerziale ECI al frame locale RTN del satellite
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

    // Se l'angolo è quasi esattamente 180 gradi (trasferimento antipodale),
    // applichiamo una perturbazione millimetrica per evitare la singolarità di piano di Lambert
    if (fabs(dnu - PI) < 1e-4) {
        dnu = PI - 1e-4;
    }

    double A = sin(dnu) * sqrt(r1 * r2 / (1.0 - cos(dnu)));
    if (fabs(A) < 1e-9) return sol;

    // Ricerca robusta di z tramite bisezione + raffinamento su [-4*PI, 4*PI^2]
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

        if (t_calc <= tof_sec) {
            z_low = z;
        } else {
            z_high = z;
        }
    }

    double C = stumpff_C(z);
    double S = stumpff_S(z);
    double y = r1 + r2 + A * (z * S - 1.0) / sqrt(C);
    if (y <= 0.0) return sol;

    // Coefficienti di Lagrange f, g, g_dot
    double f     = 1.0 - (y / r1);
    double g     = A * sqrt(y / mu);
    double g_dot = 1.0 - (y / r2);

    if (fabs(g) < 1e-9) return sol;

    // v1 = (r2 - f * r1) / g
    sol.v1_eci = vec_scale(vec_sub(r2_vec, vec_scale(r1_vec, f)), 1.0 / g);
    // v2 = (g_dot * r2 - r1) / g
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
    PropellantSpec p;
    p.type = type;
    switch (type) {
        case PROP_SOLID_HTPB:
            strcpy(p.name, "HTPB/AP Composite (Star-48B)");
            strcpy(p.category, "Solid");
            p.isp = 292.0;
            break;
        case PROP_LIQUID_HYDRAZINE:
            strcpy(p.name, "Hydrazine Monopropellant (N2H4)");
            strcpy(p.category, "Liquid");
            p.isp = 235.0;
            break;
        case PROP_LIQUID_MMH_NTO:
            strcpy(p.name, "MMH / NTO Bipropellant");
            strcpy(p.category, "Liquid");
            p.isp = 316.0;
            break;
        case PROP_LIQUID_LOX_LCH4:
            strcpy(p.name, "LOX / Liquid Methane (Methalox)");
            strcpy(p.category, "Liquid");
            p.isp = 365.0;
            break;
        case PROP_LIQUID_LOX_LH2:
        default:
            strcpy(p.name, "LOX / Liquid Hydrogen (Hydrolox)");
            strcpy(p.category, "Liquid");
            p.isp = 450.0;
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

    // ------------------------------------------------------------------
    // CANDIDATO 1: 3D OPTIMAL SPLIT-PLANE HOHMANN
    // ------------------------------------------------------------------
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

    // ------------------------------------------------------------------
    // CANDIDATO 2: BI-ELLITTICA 3D (3 Impulsi)
    // ------------------------------------------------------------------
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

    // ------------------------------------------------------------------
    // CANDIDATO 3: RISOLUTORE UNIVERSALE DI LAMBERT 3D (Se c'è vincolo tempo)
    // ------------------------------------------------------------------
    KeplerianElements k_b1 = orb_init;
    k_b1.theta = 0.0;
    CartesianState s_b1 = keplerian_to_cartesian(k_b1, mu, init_mass, 0.0);

    KeplerianElements k_b2 = orb_target;
    k_b2.theta = 165.0 * DEG2RAD; // Target non-antipodale per trasferimento veloce di Lambert
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

    // Coordinate spaziali ECI e Geodetiche dei Burn Points
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