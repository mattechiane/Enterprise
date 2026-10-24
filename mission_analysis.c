#include "mission_analysis.h"
#include "astrodynamics.h"
#include "propagator.h"
#include "telemetry.h"
#include <stdio.h>
#include <math.h>

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

OptimalTransferPlan plan_optimal_transfer(KeplerianElements orb_init, KeplerianElements orb_target, double init_mass, double isp, double mu) {
    OptimalTransferPlan plan = {0};
    double r1 = orb_init.a;
    double r2 = orb_target.a;
    double delta_inc = orb_target.i - orb_init.i;

    double a_trans = 0.5 * (r1 + r2);
    double v1 = sqrt(mu / r1);
    double v2 = sqrt(mu / r2);

    double vt_1 = sqrt(mu * ((2.0 / r1) - (1.0 / a_trans)));
    double vt_2 = sqrt(mu * ((2.0 / r2) - (1.0 / a_trans)));

    // Impulso 1: puramente tangenziale nel piano iniziale (al perigeo)
    double dv1_tangential = vt_1 - v1;
    plan.dv1_mag = fabs(dv1_tangential);
    plan.burn1_dir_rtn = (Vector3D){0.0, (dv1_tangential >= 0.0) ? 1.0 : -1.0, 0.0};

    // Impulso 2: combinato (circolarizzazione + cambio piano orbitale all'apogeo)
    // Componenti nel sistema locale RTN (Radiale, Tangenziale/In-Track, Normale/Cross-Track)
    double dv2_T = v2 * cos(delta_inc) - vt_2;
    double dv2_N = v2 * sin(delta_inc);
    Vector3D dv2_vec_rtn = {0.0, dv2_T, dv2_N};

    plan.dv2_mag = vec_norm(dv2_vec_rtn);
    plan.burn2_dir_rtn = vec_normalize(dv2_vec_rtn);

    plan.total_dv = plan.dv1_mag + plan.dv2_mag;
    plan.transfer_time_sec = PI * sqrt((a_trans * a_trans * a_trans) / mu);

    // Calcolo consumo propellente (Equazione di Tsiolkovsky)
    double total_dv_ms = plan.total_dv * 1000.0;
    plan.final_mass = init_mass * exp(-total_dv_ms / (isp * G0));
    plan.propellant_used = init_mass - plan.final_mass;

    return plan;
}

OrbitReport evaluate_orbit_regime(KeplerianElements kep, double mass, const SpacecraftConfig *sc, const GroundStation *station) {
    OrbitReport rep = {0};
    rep.initial_kep = kep;
    rep.period_min  = orbital_period(kep.a, EARTH_MU) / 60.0;

    // 1. Stima analitica ad alta precisione del decadimento per Drag (Station-Keeping)
    Vector3D r_dummy = {kep.a, 0.0, 0.0};
    Vector3D sun_dir = {1.0, 0.0, 0.0};
    SpaceEnvironment env_mean = evaluate_environment(r_dummy, sun_dir);

    double v_circ_ms = sqrt(EARTH_MU / kep.a) * 1000.0; // m/s
    // Accelerazione di drag in m/s^2: a_d = 0.5 * rho * v^2 * Cd * A / m
    double a_drag_ms2 = 0.5 * env_mean.air_density * (v_circ_ms * v_circ_ms) * (sc->drag_coeff * sc->drag_area / mass);

    // Delta-V perso al giorno e all'anno per colpa dell'attrito atmosferico
    double dv_day_ms = a_drag_ms2 * 86400.0;
    rep.station_keeping_dv_year = dv_day_ms * 365.25;

    // Decadimento del semiasse maggiore: da/dt = -2 * a * a_drag / v
    rep.sma_decay_m_per_day = (2.0 * (kep.a * 1000.0) * a_drag_ms2 / v_circ_ms) * 86400.0;

    // Carburante annuo per compensare il drag
    if (sc->thruster.isp > 0.0) {
        rep.station_keeping_fuel_yr = mass * (1.0 - exp(-rep.station_keeping_dv_year / (sc->thruster.isp * G0)));
    }

    // 2. Simulazione RK4 su 24 ore per Eclissi e Finestre di Visibilità Stazione
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
    printf("  * Parametri geometrici  : Quota = %.1f km | Inc = %.2f deg | e = %.4f\n",
           rep.initial_kep.a - EARTH_RADIUS, rep.initial_kep.i * RAD2DEG, rep.initial_kep.e);
    printf("  * Periodo Orbitale      : %.2f minuti (%.1f orbite/giorno)\n",
           rep.period_min, 1440.0 / rep.period_min);
    printf("  * Bilancio Luce / Ombra : Luce = %.1f min | Eclissi = %.1f min (%.1f%% ombra)\n",
           rep.sunlight_time_min, rep.eclipse_time_min, rep.eclipse_percentage);
    printf("  * Visibilita' (%s): %d passaggi/24h | Contatto tot = %.1f min | Max El = %.1f deg\n",
           station_name, rep.ground_passes_24h, rep.total_contact_min, rep.max_elevation_deg);
    printf("  * Decadimento Atmosfer. : -%.2f metri/giorno\n", rep.sma_decay_m_per_day);
    printf("  * Station-Keeping Annuo : Delta-V = %.2f m/s/anno | Carburante = %.3f kg/anno\n\n",
           rep.station_keeping_dv_year, rep.station_keeping_fuel_yr);
}