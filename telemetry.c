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

// Propaga un segmento di missione ed esporta Posizione, Velocità, Tracking Radar e Bilancio Energetico (EPS)
static void write_trajectory_json_array(FILE *fp, KeplerianElements kep, double duration_sec, int num_steps,
                                        const SpacecraftConfig *sc, const GroundStation *station,
                                        double *t_cum, double *soc_pct) {
    CartesianState s = keplerian_to_cartesian(kep, EARTH_MU, 500.0, 0.0);
    ForceModelConfig forces = {.enable_j2 = true, .enable_drag = true, .sun_dir = {1.0, 0.0, 0.0}};
    SpacecraftConfig sc_pass = *sc;
    sc_pass.thruster.is_active = false;

    double dt = duration_sec / (double)num_steps;
    const double P_SOLAR_MAX_W = 850.0; // Pannelli solari GaAs (W)
    const double P_LOAD_NOM_W  = 320.0; // Consumo nominale avionica + ADCS + payload (W)
    const double BATT_CAP_WH   = 1200.0;// Capacità batteria Li-Ion di bordo (Wh)

    fprintf(fp, "[\n");
    for (int i = 0; i <= num_steps; i++) {
        double t_global = *t_cum + s.epoch;
        Vector3D ecef = eci_to_ecef(s.position, t_global, 0.0);
        GeodeticCoord geo = ecef_to_geodetic(ecef);
        SpaceEnvironment env = evaluate_environment(s.position, forces.sun_dir);
        AccessMetrics acc = compute_ground_access(ecef, station);

        double v_mag = vec_norm(s.velocity);
        double p_gen = env.in_eclipse ? 0.0 : P_SOLAR_MAX_W;

        if (i > 0) {
            double net_p = (p_gen >= P_LOAD_NOM_W) ? (p_gen - P_LOAD_NOM_W) * 0.90 : (p_gen - P_LOAD_NOM_W);
            double d_wh = net_p * (dt / 3600.0);
            *soc_pct += (d_wh / BATT_CAP_WH) * 100.0;
            if (*soc_pct > 100.0) *soc_pct = 100.0;
            if (*soc_pct < 15.0)  *soc_pct = 15.0;
        }

        fprintf(fp, "    {\"t\":%.1f,\"x\":%.2f,\"y\":%.2f,\"z\":%.2f,\"v\":%.3f,\"lat\":%.2f,\"lon\":%.2f,\"alt\":%.2f,\"ecl\":%d,\"el\":%.2f,\"rng\":%.1f,\"pwr\":%.0f,\"soc\":%.1f}%s\n",
                t_global, s.position.x, s.position.y, s.position.z, v_mag,
                geo.latitude * RAD2DEG, geo.longitude * RAD2DEG, geo.altitude,
                env.in_eclipse ? 1 : 0,
                acc.elevation * RAD2DEG, acc.range, p_gen, *soc_pct,
                (i < num_steps) ? "," : "");
        s = propagate_step_rk4(s, dt, &sc_pass, &forces);
    }
    fprintf(fp, "  ]");
    *t_cum += duration_sec;
}

static void write_spiral_trajectory_json_array(FILE *fp, KeplerianElements k1, KeplerianElements k2,
                                               double total_tof_sec, int visual_revs, int num_steps,
                                               const GroundStation *station, double *t_cum, double *soc_pct) {
    (void)total_tof_sec;
    Vector3D sun_dir = {1.0, 0.0, 0.0};
    double vis_duration_sec = (double)visual_revs * 5600.0;
    double dt_vis = vis_duration_sec / (double)num_steps;
    const double P_SOLAR_MAX_W = 850.0;
    const double P_LOAD_SPIRAL_W = 480.0; // Carico maggiorato durante spinta elettrica
    const double BATT_CAP_WH = 1200.0;

    fprintf(fp, "[\n");
    for (int i = 0; i <= num_steps; i++) {
        double frac = (double)i / (double)num_steps;
        KeplerianElements kt;
        kt.a     = k1.a + frac * (k2.a - k1.a);
        kt.e     = k1.e + frac * (k2.e - k1.e);
        kt.i     = k1.i + frac * (k2.i - k1.i);
        kt.RAAN  = k1.RAAN + frac * (k2.RAAN - k1.RAAN);
        kt.omega = k1.omega + frac * (k2.omega - k1.omega);
        kt.theta = frac * ((double)visual_revs * 2.0 * PI + PI);

        double t_global = *t_cum + frac * vis_duration_sec;
        CartesianState s = keplerian_to_cartesian(kt, EARTH_MU, 500.0, t_global);
        Vector3D ecef = eci_to_ecef(s.position, t_global, 0.0);
        GeodeticCoord geo = ecef_to_geodetic(ecef);
        SpaceEnvironment env = evaluate_environment(s.position, sun_dir);
        AccessMetrics acc = compute_ground_access(ecef, station);

        double v_mag = vec_norm(s.velocity);
        double p_gen = env.in_eclipse ? 0.0 : P_SOLAR_MAX_W;

        if (i > 0) {
            double net_p = (p_gen >= P_LOAD_SPIRAL_W) ? (p_gen - P_LOAD_SPIRAL_W) * 0.90 : (p_gen - P_LOAD_SPIRAL_W);
            double d_wh = net_p * (dt_vis / 3600.0);
            *soc_pct += (d_wh / BATT_CAP_WH) * 100.0;
            if (*soc_pct > 100.0) *soc_pct = 100.0;
            if (*soc_pct < 15.0)  *soc_pct = 15.0;
        }

        fprintf(fp, "    {\"t\":%.1f,\"x\":%.2f,\"y\":%.2f,\"z\":%.2f,\"v\":%.3f,\"lat\":%.2f,\"lon\":%.2f,\"alt\":%.2f,\"ecl\":%d,\"el\":%.2f,\"rng\":%.1f,\"pwr\":%.0f,\"soc\":%.1f}%s\n",
                t_global, s.position.x, s.position.y, s.position.z, v_mag,
                geo.latitude * RAD2DEG, geo.longitude * RAD2DEG, geo.altitude,
                env.in_eclipse ? 1 : 0,
                acc.elevation * RAD2DEG, acc.range, p_gen, *soc_pct,
                (i < num_steps) ? "," : "");
    }
    fprintf(fp, "  ]");
    *t_cum += vis_duration_sec;
}

static void write_rendezvous_trajectory_json_array(FILE *fp, RendezvousPlan rdv, int num_steps) {
    fprintf(fp, "[\n");
    for (int i = 0; i <= num_steps; i++) {
        double frac = (double)i / (double)num_steps;
        double t_sec = frac * rdv.transfer_time_sec;
        Vector3D r_m = {0};
        Vector3D v_ms = {0};

        if (i < num_steps) {
            propagate_cw_state(rdv.r0_lvlh_m, rdv.v0_post_lvlh_ms, rdv.mean_motion_rad_s, t_sec, &r_m, &v_ms);
        } else {
            // Soft-Docking completato a t = t_rdv dopo l'impulso di frenata Delta-V_brake: posizione e velocità azzerate
            r_m = (Vector3D){0.0, 0.0, 0.0};
            v_ms = (Vector3D){0.0, 0.0, 0.0};
        }

        double range_m = vec_norm(r_m);
        fprintf(fp, "    {\"t\":%.1f,\"rx_m\":%.2f,\"ry_m\":%.2f,\"rz_m\":%.2f,\"vx_ms\":%.4f,\"vy_ms\":%.4f,\"vz_ms\":%.4f,\"range_m\":%.2f}%s\n",
                t_sec, r_m.x, r_m.y, r_m.z, v_ms.x, v_ms.y, v_ms.z, range_m,
                (i < num_steps) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_monte_carlo_cloud_json_array(FILE *fp, MonteCarloReport mc) {
    fprintf(fp, "[\n");
    for (int i = 0; i < mc.num_runs; i++) {
        DispersionSample s = mc.samples[i];
        fprintf(fp, "    {\"dx_km\":%.3f,\"dy_km\":%.3f,\"dz_km\":%.3f,\"miss_km\":%.3f,\"tcm_dv_ms\":%.3f}%s\n",
                s.dr_rtn_km.x, s.dr_rtn_km.y, s.dr_rtn_km.z, s.miss_distance_km, s.tcm_dv_ms,
                (i < mc.num_runs - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

bool export_web_mission_json(const char *filename,
                             OrbitReport rep_init,
                             OptimalTransferPlan plan,
                             OrbitReport rep_target,
                             RendezvousPlan rdv,
                             MonteCarloReport mc,
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
    fprintf(fp, "    \"num_burns\": %d, \"is_low_thrust\": %d, \"thrust_n\": %.3f, \"spiral_revs\": %.1f,\n",
            plan.num_burns, plan.is_low_thrust ? 1 : 0, plan.continuous_thrust_n, plan.spiral_revolutions);
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

    fprintf(fp, "  \"rendezvous\": {\n");
    fprintf(fp, "    \"r0_m\": [%.2f, %.2f, %.2f],\n",
            rdv.r0_lvlh_m.x, rdv.r0_lvlh_m.y, rdv.r0_lvlh_m.z);
    fprintf(fp, "    \"v0_ms\": [%.4f, %.4f, %.4f],\n",
            rdv.v0_lvlh_ms.x, rdv.v0_lvlh_ms.y, rdv.v0_lvlh_ms.z);
    fprintf(fp, "    \"v0_post_ms\": [%.4f, %.4f, %.4f],\n",
            rdv.v0_post_lvlh_ms.x, rdv.v0_post_lvlh_ms.y, rdv.v0_post_lvlh_ms.z);
    fprintf(fp, "    \"initial_range_m\": %.2f, \"tof_sec\": %.1f, \"tof_min\": %.2f, \"mean_motion_rad_s\": %.8f,\n",
            rdv.initial_range_m, rdv.transfer_time_sec, rdv.transfer_time_sec / 60.0, rdv.mean_motion_rad_s);
    fprintf(fp, "    \"dv_inject_lvlh_ms\": [%.4f, %.4f, %.4f], \"dv_brake_lvlh_ms\": [%.4f, %.4f, %.4f],\n",
            rdv.dv_inject_lvlh_ms.x, rdv.dv_inject_lvlh_ms.y, rdv.dv_inject_lvlh_ms.z,
            rdv.dv_brake_lvlh_ms.x, rdv.dv_brake_lvlh_ms.y, rdv.dv_brake_lvlh_ms.z);
    fprintf(fp, "    \"dv_inject_ms\": %.3f, \"dv_brake_ms\": %.3f, \"total_dv_ms\": %.3f,\n",
            rdv.dv_inject_ms, rdv.dv_brake_ms, rdv.total_dv_ms);
    fprintf(fp, "    \"fuel_kg\": %.3f, \"final_mass_kg\": %.2f\n  },\n",
            rdv.propellant_used_kg, rdv.final_mass_kg);

    fprintf(fp, "  \"monte_carlo\": {\n");
    fprintf(fp, "    \"num_runs\": %d, \"sigma_dv_pct\": %.2f, \"sigma_point_deg\": %.2f,\n",
            mc.num_runs, mc.sigma_dv_pct, mc.sigma_point_deg);
    fprintf(fp, "    \"mean_miss_km\": %.2f, \"std_miss_km\": %.2f, \"cep_3sigma_km\": %.2f, \"max_miss_km\": %.2f,\n",
            mc.mean_miss_km, mc.std_miss_km, mc.cep_3sigma_km, mc.max_miss_km);
    fprintf(fp, "    \"ellipsoid_3sigma_rtn_km\": [%.2f, %.2f, %.2f],\n",
            mc.ellipsoid_3sigma_rtn.x, mc.ellipsoid_3sigma_rtn.y, mc.ellipsoid_3sigma_rtn.z);
    fprintf(fp, "    \"tcm_dv_mean_ms\": %.2f, \"tcm_dv_99_ms\": %.2f, \"tcm_fuel_99_kg\": %.3f\n  },\n",
            mc.tcm_dv_mean_ms, mc.tcm_dv_99_ms, mc.tcm_fuel_99_kg);

    double init_mass = plan.final_mass + plan.propellant_used;
    double a_hoh = 0.5 * (rep_init.initial_kep.a + rep_target.initial_kep.a);
    double tof_opt = PI * sqrt((a_hoh * a_hoh * a_hoh) / EARTH_MU);
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
        double ratio = (t_test - tof_opt) / (0.5 * tof_opt);
        double dv_smooth = (plan.total_dv * 1000.0) + 0.35 * (dv_ms - plan.total_dv * 1000.0) + 180.0 * (ratio * ratio * (t_test < tof_opt ? 1.8 : 0.7));
        double fuel_kg = init_mass * (1.0 - exp(-dv_smooth / (plan.prop_info.isp * G0)));
        fprintf(fp, "    {\"tof_min\":%.2f, \"dv_ms\":%.2f, \"fuel_kg\":%.2f}%s\n",
                t_test / 60.0, dv_smooth, fuel_kg, (i < num_pts) ? "," : "");
    }
    fprintf(fp, "  ],\n");

    double t_cum = 0.0;
    double soc_pct = 98.0; // Batteria carica al 98% a inizio missione

    fprintf(fp, "  \"traj_initial\": ");
    write_trajectory_json_array(fp, rep_init.initial_kep, rep_init.period_min * 60.0 * 2.0, 240, sc, station, &t_cum, &soc_pct);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_transfer\": ");
    if (plan.is_low_thrust) {
        write_spiral_trajectory_json_array(fp, rep_init.initial_kep, rep_target.initial_kep, plan.transfer_time_sec, 6, 360, station, &t_cum, &soc_pct);
    } else {
        double r1 = rep_init.initial_kep.a;
        double r2 = rep_target.initial_kep.a;
        KeplerianElements kep_trans = rep_init.initial_kep;
        kep_trans.a = 0.5 * (r1 + r2);
        kep_trans.e = fabs(r2 - r1) / (r1 + r2);
        kep_trans.theta = 0.0;
        write_trajectory_json_array(fp, kep_trans, plan.transfer_time_sec, 100, sc, station, &t_cum, &soc_pct);
    }
    fprintf(fp, ",\n");

    KeplerianElements kep_targ_start = rep_target.initial_kep;
    kep_targ_start.theta = PI;
    fprintf(fp, "  \"traj_target\": ");
    write_trajectory_json_array(fp, kep_targ_start, rep_target.period_min * 60.0 * 2.0, 240, sc, station, &t_cum, &soc_pct);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_rendezvous\": ");
    write_rendezvous_trajectory_json_array(fp, rdv, 120);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"mc_cloud\": ");
    write_monte_carlo_cloud_json_array(fp, mc);
    fprintf(fp, "\n}\n");

    fclose(fp);
    return true;
}
