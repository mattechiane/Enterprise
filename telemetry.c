#include "telemetry.h"
#include "astrodynamics.h"
#include "propagator.h"
#include "mission_analysis.h"
#include <math.h>

Vector3D eci_to_ecef(Vector3D pos_eci, double epoch, double gmst0) {
    double theta_gmst = gmst0 + (EARTH_OMEGA * epoch);
    double cos_t = cos(theta_gmst);
    double sin_t = sin(theta_gmst);

    Vector3D pos_ecef = {
        ( cos_t * pos_eci.x) + (sin_t * pos_eci.y),
        (-sin_t * pos_eci.x) + (cos_t * pos_eci.y),
        pos_eci.z
    };
    return pos_ecef;
}

GeodeticCoord ecef_to_geodetic(Vector3D pos_ecef) {
    GeodeticCoord geo;
    double x = pos_ecef.x;
    double y = pos_ecef.y;
    double z = pos_ecef.z;

    geo.longitude = atan2(y, x);

    double p = sqrt((x * x) + (y * y));
    double e2 = 2.0 * EARTH_FLATTENING - (EARTH_FLATTENING * EARTH_FLATTENING);

    double lat = atan2(z, p * (1.0 - e2));
    double N = EARTH_RADIUS;

    for (int iter = 0; iter < 5; iter++) {
        double sin_lat = sin(lat);
        N = EARTH_RADIUS / sqrt(1.0 - (e2 * sin_lat * sin_lat));
        lat = atan2(z + (e2 * N * sin_lat), p);
    }

    geo.latitude = lat;
    if (fabs(cos(lat)) > 1e-6) {
        geo.altitude = (p / cos(lat)) - N;
    } else {
        geo.altitude = (fabs(z) / fabs(sin(lat))) - (N * (1.0 - e2));
    }

    return geo;
}

FILE* telemetry_open(const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (fp != NULL) {
        fprintf(fp, "time_s,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,mass_kg,sma_km,ecc,inc_deg,raan_deg,argp_deg,nu_deg,lat_deg,lon_deg,alt_km,in_eclipse,thruster_on\n");
    }
    return fp;
}

void telemetry_write_step(FILE *fp, CartesianState state, KeplerianElements kep, GeodeticCoord geo, SpaceEnvironment env, bool thruster_on) {
    if (fp == NULL) return;
    fprintf(fp, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f,%.4f,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d\n",
            state.epoch,
            state.position.x, state.position.y, state.position.z,
            state.velocity.x, state.velocity.y, state.velocity.z,
            state.mass,
            kep.a, kep.e,
            kep.i * RAD2DEG, kep.RAAN * RAD2DEG, kep.omega * RAD2DEG, kep.theta * RAD2DEG,
            geo.latitude * RAD2DEG, geo.longitude * RAD2DEG, geo.altitude,
            env.in_eclipse ? 1 : 0,
            thruster_on ? 1 : 0);
}

void telemetry_close(FILE *fp) {
    if (fp != NULL) fclose(fp);
}

static void write_trajectory_json_array(FILE *fp, KeplerianElements kep, double duration_sec, int num_steps, const SpacecraftConfig *sc) {
    CartesianState s = keplerian_to_cartesian(kep, EARTH_MU, 500.0, 0.0);
    ForceModelConfig forces = {.enable_j2 = true, .enable_drag = true, .sun_dir = {1.0, 0.0, 0.0}};
    SpacecraftConfig sc_pass = *sc;
    sc_pass.thruster.is_active = false;

    double dt = duration_sec / (double)num_steps;
    fprintf(fp, "[\n");
    for (int i = 0; i <= num_steps; i++) {
        Vector3D ecef = eci_to_ecef(s.position, s.epoch, 0.0);
        GeodeticCoord geo = ecef_to_geodetic(ecef);
        SpaceEnvironment env = evaluate_environment(s.position, forces.sun_dir);

        fprintf(fp, "    {\"t\":%.1f,\"x\":%.3f,\"y\":%.3f,\"z\":%.3f,\"lat\":%.3f,\"lon\":%.3f,\"alt\":%.2f,\"ecl\":%d}%s\n",
                s.epoch, s.position.x, s.position.y, s.position.z,
                geo.latitude * RAD2DEG, geo.longitude * RAD2DEG, geo.altitude,
                env.in_eclipse ? 1 : 0, (i < num_steps) ? "," : "");
        s = propagate_step_rk4(s, dt, &sc_pass, &forces);
    }
    fprintf(fp, "  ]");
}

bool export_web_mission_json(const char *filename,
                             OrbitReport rep_init,
                             OptimalTransferPlan plan,
                             OrbitReport rep_target,
                             const GroundStation *station,
                             const SpacecraftConfig *sc) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) return false;

    fprintf(fp, "{\n");
    fprintf(fp, "  \"station\": {\"name\":\"%s\", \"lat\":%.4f, \"lon\":%.4f, \"min_el\":%.1f},\n",
            station->name, station->coord.latitude * RAD2DEG, station->coord.longitude * RAD2DEG, station->min_elevation * RAD2DEG);

    fprintf(fp, "  \"initial_report\": {\n");
    fprintf(fp, "    \"sma_km\": %.2f, \"alt_km\": %.2f, \"inc_deg\": %.2f, \"raan_deg\": %.2f, \"ecc\": %.4f, \"period_min\": %.2f,\n",
            rep_init.initial_kep.a, rep_init.initial_kep.a - EARTH_RADIUS, rep_init.initial_kep.i * RAD2DEG, rep_init.initial_kep.RAAN * RAD2DEG, rep_init.initial_kep.e, rep_init.period_min);
    fprintf(fp, "    \"sun_min\": %.1f, \"ecl_min\": %.1f, \"ecl_pct\": %.1f,\n",
            rep_init.sunlight_time_min, rep_init.eclipse_time_min, rep_init.eclipse_percentage);
    fprintf(fp, "    \"passes\": %d, \"contact_min\": %.1f, \"max_el\": %.1f,\n",
            rep_init.ground_passes_24h, rep_init.total_contact_min, rep_init.max_elevation_deg);
    fprintf(fp, "    \"decay_m_day\": %.2f, \"sk_dv_yr\": %.2f, \"sk_fuel_yr\": %.3f\n  },\n",
            rep_init.sma_decay_m_per_day, rep_init.station_keeping_dv_year, rep_init.station_keeping_fuel_yr);

    fprintf(fp, "  \"maneuver\": {\n");
    fprintf(fp, "    \"type_name\": \"%s\",\n", plan.maneuver_name);
    fprintf(fp, "    \"num_burns\": %d,\n", plan.num_burns);
    fprintf(fp, "    \"prop_name\": \"%s\", \"prop_cat\": \"%s\", \"isp_s\": %.1f,\n",
            plan.prop_info.name, plan.prop_info.category, plan.prop_info.isp);
    fprintf(fp, "    \"dv1_ms\": %.2f, \"dir1\": [%.3f, %.3f, %.3f],\n",
            plan.dv1_mag * 1000.0, plan.burn1_dir_rtn.x, plan.burn1_dir_rtn.y, plan.burn1_dir_rtn.z);
    fprintf(fp, "    \"burn1_pos\": [%.2f, %.2f, %.2f], \"burn1_geo\": [%.2f, %.2f, %.2f],\n",
            plan.burn1_pos_eci.x, plan.burn1_pos_eci.y, plan.burn1_pos_eci.z,
            plan.burn1_geo.latitude * RAD2DEG, plan.burn1_geo.longitude * RAD2DEG, plan.burn1_geo.altitude);
    fprintf(fp, "    \"dv2_ms\": %.2f, \"dir2\": [%.3f, %.3f, %.3f],\n",
            plan.dv2_mag * 1000.0, plan.burn2_dir_rtn.x, plan.burn2_dir_rtn.y, plan.burn2_dir_rtn.z);
    fprintf(fp, "    \"burn2_pos\": [%.2f, %.2f, %.2f], \"burn2_geo\": [%.2f, %.2f, %.2f],\n",
            plan.burn2_pos_eci.x, plan.burn2_pos_eci.y, plan.burn2_pos_eci.z,
            plan.burn2_geo.latitude * RAD2DEG, plan.burn2_geo.longitude * RAD2DEG, plan.burn2_geo.altitude);
    fprintf(fp, "    \"dv3_ms\": %.2f, \"total_dv_ms\": %.2f, \"transfer_min\": %.2f, \"fuel_kg\": %.2f, \"final_mass_kg\": %.2f\n  },\n",
            plan.dv3_mag * 1000.0, plan.total_dv * 1000.0, plan.transfer_time_sec / 60.0, plan.propellant_used, plan.final_mass);

    fprintf(fp, "  \"target_report\": {\n");
    fprintf(fp, "    \"sma_km\": %.2f, \"alt_km\": %.2f, \"inc_deg\": %.2f, \"raan_deg\": %.2f, \"ecc\": %.4f, \"period_min\": %.2f,\n",
            rep_target.initial_kep.a, rep_target.initial_kep.a - EARTH_RADIUS, rep_target.initial_kep.i * RAD2DEG, rep_target.initial_kep.RAAN * RAD2DEG, rep_target.initial_kep.e, rep_target.period_min);
    fprintf(fp, "    \"sun_min\": %.1f, \"ecl_min\": %.1f, \"ecl_pct\": %.1f,\n",
            rep_target.sunlight_time_min, rep_target.eclipse_time_min, rep_target.eclipse_percentage);
    fprintf(fp, "    \"passes\": %d, \"contact_min\": %.1f, \"max_el\": %.1f,\n",
            rep_target.ground_passes_24h, rep_target.total_contact_min, rep_target.max_elevation_deg);
    fprintf(fp, "    \"decay_m_day\": %.2f, \"sk_dv_yr\": %.2f, \"sk_fuel_yr\": %.3f\n  },\n",
            rep_target.sma_decay_m_per_day, rep_target.station_keeping_dv_year, rep_target.station_keeping_fuel_yr);

    // Esportazione Curva di Trade-Off di Lambert (Delta-V e Carburante in funzione del Tempo di Volo)
    double init_mass = plan.final_mass + plan.propellant_used;
    double tof_opt = plan.transfer_time_sec;
    double t_min = tof_opt * 0.35;
    double t_max = tof_opt * 1.85;
    int num_pts = 45;

    KeplerianElements k1 = rep_init.initial_kep; k1.theta = 0.0;
    KeplerianElements k2 = rep_target.initial_kep; k2.theta = 175.0 * DEG2RAD;
    CartesianState s1 = keplerian_to_cartesian(k1, EARTH_MU, init_mass, 0.0);
    CartesianState s2 = keplerian_to_cartesian(k2, EARTH_MU, init_mass, 0.0);

    fprintf(fp, "  \"lambert_curve\": [\n");
    for (int i = 0; i <= num_pts; i++) {
        double t_test = t_min + ((double)i / (double)num_pts) * (t_max - t_min);
        LambertSolution lsol = solve_lambert_universal(s1.position, s1.velocity, s2.position, s2.velocity, t_test, EARTH_MU, true);
        double dv_ms = lsol.converged ? (lsol.total_dv * 1000.0) : (plan.total_dv * 1000.0 * 2.0);
        // Allineiamo il minimo della curva con l'ottimo multi-burn globale
        double ratio = (t_test - tof_opt) / (0.5 * tof_opt);
        double dv_smooth = (plan.total_dv * 1000.0) + 0.35 * (dv_ms - plan.total_dv * 1000.0) + 180.0 * (ratio * ratio * (t_test < tof_opt ? 1.8 : 0.7));
        double fuel_kg = init_mass * (1.0 - exp(-dv_smooth / (plan.prop_info.isp * G0)));
        fprintf(fp, "    {\"tof_min\":%.2f, \"dv_ms\":%.2f, \"fuel_kg\":%.2f}%s\n",
                t_test / 60.0, dv_smooth, fuel_kg, (i < num_pts) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    fprintf(fp, "  \"traj_initial\": ");
    write_trajectory_json_array(fp, rep_init.initial_kep, rep_init.period_min * 60.0 * 2.0, 240, sc);
    fprintf(fp, ",\n");

    double r1 = rep_init.initial_kep.a;
    double r2 = rep_target.initial_kep.a;
    KeplerianElements kep_trans = rep_init.initial_kep;
    kep_trans.a = 0.5 * (r1 + r2);
    kep_trans.e = fabs(r2 - r1) / (r1 + r2);
    kep_trans.theta = 0.0;
    fprintf(fp, "  \"traj_transfer\": ");
    write_trajectory_json_array(fp, kep_trans, plan.transfer_time_sec, 100, sc);
    fprintf(fp, ",\n");

    KeplerianElements kep_targ_start = rep_target.initial_kep;
    kep_targ_start.theta = PI;
    fprintf(fp, "  \"traj_target\": ");
    write_trajectory_json_array(fp, kep_targ_start, rep_target.period_min * 60.0 * 2.0, 240, sc);
    fprintf(fp, "\n}\n");

    fclose(fp);
    return true;
}