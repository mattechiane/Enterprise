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

static void write_adcs_trajectory_json_array(FILE *fp, ADCSReport adcs) {
    fprintf(fp, "[\n");
    for (int i = 0; i < adcs.num_samples; i++) {
        AttitudeSample s = adcs.samples[i];
        fprintf(fp, "    {\"t_min\":%.2f,\"roll_deg\":%.2f,\"pitch_deg\":%.2f,\"yaw_deg\":%.2f,\"point_err_deg\":%.3f,\"omega_deg_s\":%.3f,\"rw_torque_mnm\":%.2f,\"rw_mom_nms\":%.3f,\"rw_sat_pct\":%.2f,\"mode_id\":%d}%s\n",
                s.t_min, s.roll_deg, s.pitch_deg, s.yaw_deg, s.point_err_deg,
                s.omega_norm_deg_s, s.rw_torque_mnm, s.rw_mom_nms, s.rw_sat_pct, (int)s.mode,
                (i < adcs.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_thermal_trajectory_json_array(FILE *fp, ThermalRadReport tcs) {
    fprintf(fp, "[\n");
    for (int i = 0; i < tcs.num_samples; i++) {
        ThermalSample s = tcs.samples[i];
        fprintf(fp, "    {\"t_min\":%.2f,\"t_ext_c\":%.2f,\"t_bus_c\":%.2f,\"q_sun_w\":%.1f,\"q_ir_alb_w\":%.1f,\"q_heater_w\":%.1f,\"inst_dose_urad_h\":%.1f,\"temp_c\":%.2f,\"q_sol_w\":%.1f,\"q_alb_w\":%.1f,\"q_ir_w\":%.1f,\"q_int_w\":%.1f,\"q_rad_w\":%.1f,\"ecl\":%d,\"saa\":%d}%s\n",
                s.t_min, s.t_ext_c, s.t_bus_c, s.q_sun_w, s.q_ir_alb_w, s.q_heater_w, s.inst_dose_urad_h,
                s.temp_c, s.q_sol_w, s.q_alb_w, s.q_ir_w, s.q_int_w, s.q_rad_w, s.in_eclipse, s.in_saa,
                (i < tcs.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_telecom_trajectory_json_array(FILE *fp, TelecomReport telecom) {
    fprintf(fp, "[\n");
    for (int i = 0; i < telecom.num_samples; i++) {
        LinkSample s = telecom.samples[i];
        fprintf(fp, "    {\"t_min\":%.2f,\"fspl_db\":%.2f,\"ebn0_db\":%.2f,\"link_margin_db\":%.2f,\"achieved_mbps\":%.2f,\"cum_data_gb\":%.3f,\"el_deg\":%.2f,\"range_km\":%.1f,\"in_pass\":%d}%s\n",
                s.t_min, s.fspl_db, s.ebn0_db, s.link_margin_db, s.achieved_mbps, s.cum_data_gb,
                s.elevation_deg, s.slant_range_km, s.in_pass,
                (i < telecom.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_ekf_trajectory_json_array(FILE *fp, NavigationReport nav) {
    fprintf(fp, "[\n");
    for (int i = 0; i < nav.num_samples; i++) {
        EKFSample s = nav.samples[i];
        fprintf(fp, "    {\"t_min\":%.2f,\"raw_err_m\":%.2f,\"ekf_err_m\":%.2f,\"cov_3sigma_m\":%.2f,\"vel_err_ms\":%.4f,\"cov_vel_3sigma_ms\":%.4f,\"kalman_gain\":%.4f}%s\n",
                s.t_min, s.raw_err_m, s.ekf_err_m, s.cov_3sigma_m,
                s.vel_err_ms, s.cov_vel_3sigma_ms, s.kalman_gain_norm,
                (i < nav.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_cam_tradeoff_json_array(FILE *fp, CAMReport cam) {
    fprintf(fp, "[\n");
    for (int i = 0; i < cam.num_samples; i++) {
        BPlaneSample s = cam.samples[i];
        fprintf(fp, "    {\"dv_cam_cms\":%.2f,\"dv_cam_ms\":%.4f,\"xi_m\":%.2f,\"zeta_m\":%.2f,\"miss_dist_m\":%.2f,\"pc\":%.4e,\"log10_pc\":%.3f}%s\n",
                s.dv_cam_cms, s.dv_cam_ms, s.xi_m, s.zeta_m, s.miss_dist_m, s.pc, s.log10_pc,
                (i < cam.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_reentry_trajectory_json_array(FILE *fp, EndOfLifeReport eol) {
    fprintf(fp, "[\n");
    for (int i = 0; i < eol.num_samples; i++) {
        ReentrySample s = eol.samples[i];
        fprintf(fp, "    {\"t_sec\":%.1f,\"alt_km\":%.2f,\"vel_kms\":%.3f,\"fpa_deg\":%.2f,\"decel_g\":%.2f,\"heat_flux_wcm2\":%.2f,\"qdyn_kpa\":%.2f,\"downrange_km\":%.1f}%s\n",
                s.t_sec, s.alt_km, s.vel_kms, s.fpa_deg, s.decel_g, s.heat_flux_wcm2, s.qdyn_kpa, s.downrange_km,
                (i < eol.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_cr3bp_trajectory_json_array(FILE *fp, const CR3BPSample *samples, int num_samples) {
    fprintf(fp, "[\n");
    for (int i = 0; i < num_samples; i++) {
        CR3BPSample s = samples[i];
        fprintf(fp, "    {\"t_days\":%.4f,\"x_syn_km\":%.2f,\"y_syn_km\":%.2f,\"z_syn_km\":%.2f,\"jacobi_c\":%.10f,\"dist_moon_km\":%.2f}%s\n",
                s.t_days, s.x_syn_km, s.y_syn_km, s.z_syn_km, s.jacobi_c, s.dist_moon_km,
                (i < num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_csi_trajectory_json_array(FILE *fp, StructuralReport csi) {
    fprintf(fp, "[\n");
    for (int i = 0; i < csi.num_samples; i++) {
        SloshFlexSample s = csi.samples[i];
        fprintf(fp, "    {\"t_sec\":%.2f,\"tip_defl_mm\":%.3f,\"slosh_force_n\":%.3f,\"slosh_torque_mnm\":%.3f,\"jitter_arcsec\":%.2f,\"bode_freq_hz\":%.4f,\"bode_mag_db\":%.2f,\"bode_unfilt_db\":%.2f}%s\n",
                s.t_sec, s.tip_defl_mm, s.slosh_force_n, s.slosh_torque_mnm, s.jitter_arcsec,
                s.bode_freq_hz, s.bode_mag_db, s.bode_unfilt_db,
                (i < csi.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void write_fdir_trajectory_json_array(FILE *fp, FDIRReport fdir) {
    fprintf(fp, "[\n");
    for (int i = 0; i < fdir.num_samples; i++) {
        FDIRSample s = fdir.samples[i];
        fprintf(fp, "    {\"t_sec\":%.2f,\"sys_health_pct\":%.2f,\"att_err_deg\":%.4f,\"bus_voltage_v\":%.3f,\"fdir_state_id\":%d,\"primary_unit_pct\":%.1f,\"redundant_unit_pct\":%.1f}%s\n",
                s.t_sec, s.sys_health_pct, s.att_err_deg, s.bus_voltage_v, s.fdir_state_id,
                s.primary_unit_pct, s.redundant_unit_pct,
                (i < fdir.num_samples - 1) ? "," : "");
    }
    fprintf(fp, "  ]");
}

static void format_oem_utc_timestamp(double elapsed_sec, char *buf, size_t buf_size) {
    // Base Epoch: 2026-11-12T00:00:00.000 UTC
    int total_ms = (int)llround(fmax(0.0, elapsed_sec) * 1000.0);
    int ms = total_ms % 1000;
    int total_s = total_ms / 1000;
    int sec = total_s % 60;
    int total_m = total_s / 60;
    int min = total_m % 60;
    int total_h = total_m / 60;
    int hour = total_h % 24;
    int day_offset = total_h / 24;
    int day = 12 + day_offset;
    snprintf(buf, buf_size, "2026-11-%02dT%02d:%02d:%02d.%03d", day, hour, min, sec, ms);
}

bool export_ccsds_oem_file(const char *filename,
                           OrbitReport rep_init,
                           OptimalTransferPlan plan,
                           OrbitReport rep_target,
                           const SpacecraftConfig *sc) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) return false;

    double dur_init  = rep_init.period_min * 60.0;
    double dur_trans = plan.transfer_time_sec;
    double dur_targ  = rep_target.period_min * 60.0;
    double total_dur = dur_init + dur_trans + dur_targ;

    char start_utc[40], stop_utc[40];
    format_oem_utc_timestamp(0.0, start_utc, sizeof(start_utc));
    format_oem_utc_timestamp(total_dur, stop_utc, sizeof(stop_utc));

    fprintf(fp, "CCSDS_OEM_VERS = 2.0\n");
    fprintf(fp, "COMMENT  GENERATED BY PROJECT ENTERPRISE C11 ASTRODYNAMICS & FDIR SIMULATOR\n");
    fprintf(fp, "COMMENT  STANDARD CCSDS 502.0-B-2 ORBIT EPHEMERIS MESSAGE (GMAT / STK / ESOC READY)\n");
    fprintf(fp, "CREATION_DATE  = 2026-11-12T12:00:00.000\n");
    fprintf(fp, "ORIGINATOR     = PROJECT_ENTERPRISE_FD_GNC\n\n");

    fprintf(fp, "META_START\n");
    fprintf(fp, "OBJECT_NAME          = ENTERPRISE-1\n");
    fprintf(fp, "OBJECT_ID            = 2026-001A\n");
    fprintf(fp, "CENTER_NAME          = EARTH\n");
    fprintf(fp, "REF_FRAME            = EME2000\n");
    fprintf(fp, "REF_FRAME_EPOCH      = 2000-01-01T12:00:00.000\n");
    fprintf(fp, "TIME_SYSTEM          = UTC\n");
    fprintf(fp, "START_TIME           = %s\n", start_utc);
    fprintf(fp, "USEABLE_START_TIME   = %s\n", start_utc);
    fprintf(fp, "USEABLE_STOP_TIME    = %s\n", stop_utc);
    fprintf(fp, "STOP_TIME            = %s\n", stop_utc);
    fprintf(fp, "INTERPOLATION        = HERMITE\n");
    fprintf(fp, "INTERPOLATION_DEGREE = 7\n");
    fprintf(fp, "META_STOP\n\n");

    fprintf(fp, "COMMENT  Columns: Epoch [UTC]  X [km]  Y [km]  Z [km]  VX [km/s]  VY [km/s]  VZ [km/s]\n");
    fprintf(fp, "COMMENT  --- PHASE 1: INITIAL PARKING ORBIT (h = %.1f km, i = %.2f deg) ---\n",
            rep_init.initial_kep.a - EARTH_RADIUS, rep_init.initial_kep.i * RAD2DEG);

    ForceModelConfig forces = {.enable_j2 = true, .enable_drag = true, .sun_dir = {1.0, 0.0, 0.0}};
    SpacecraftConfig sc_pass = *sc;
    sc_pass.thruster.is_active = false;
    double init_mass = plan.final_mass + plan.propellant_used;
    CartesianState state = keplerian_to_cartesian(rep_init.initial_kep, EARTH_MU, init_mass, 0.0);
    double t_elapsed = 0.0;

    int n_init = 60, n_trans = 60, n_targ = 60;
    double dt_init = dur_init / (double)n_init;
    for (int i = 0; i <= n_init; i++) {
        char t_str[40];
        format_oem_utc_timestamp(t_elapsed, t_str, sizeof(t_str));
        fprintf(fp, "%s  %14.6f  %14.6f  %14.6f  %12.7f  %12.7f  %12.7f\n",
                t_str,
                state.position.x, state.position.y, state.position.z,
                state.velocity.x, state.velocity.y, state.velocity.z);
        if (i < n_init) {
            state = propagate_step_rk4(state, dt_init, &sc_pass, &forces);
            t_elapsed += dt_init;
        }
    }

    fprintf(fp, "\nCOMMENT  --- PHASE 2: ORBITAL TRANSFER SEGMENT (%s, dV = %.2f m/s) ---\n",
            plan.maneuver_name, plan.total_dv * 1000.0);
    double r1 = rep_init.initial_kep.a;
    double r2 = rep_target.initial_kep.a;
    KeplerianElements kep_trans = rep_init.initial_kep;
    kep_trans.a = 0.5 * (r1 + r2);
    kep_trans.e = fabs(r2 - r1) / (r1 + r2);
    kep_trans.omega = (r2 >= r1) ? rep_init.initial_kep.omega : (rep_init.initial_kep.omega + PI);
    kep_trans.theta = (r2 >= r1) ? 0.0 : PI;
    state = keplerian_to_cartesian(kep_trans, EARTH_MU, init_mass - 0.5 * plan.propellant_used, t_elapsed);

    double dt_trans = dur_trans / (double)n_trans;
    for (int i = 1; i <= n_trans; i++) {
        state = propagate_step_rk4(state, dt_trans, &sc_pass, &forces);
        t_elapsed += dt_trans;
        char t_str[40];
        format_oem_utc_timestamp(t_elapsed, t_str, sizeof(t_str));
        fprintf(fp, "%s  %14.6f  %14.6f  %14.6f  %12.7f  %12.7f  %12.7f\n",
                t_str,
                state.position.x, state.position.y, state.position.z,
                state.velocity.x, state.velocity.y, state.velocity.z);
    }

    fprintf(fp, "\nCOMMENT  --- PHASE 3: TARGET OPERATIONAL ORBIT (h = %.1f km, i = %.2f deg) ---\n",
            rep_target.initial_kep.a - EARTH_RADIUS, rep_target.initial_kep.i * RAD2DEG);
    KeplerianElements kep_targ_start = rep_target.initial_kep;
    kep_targ_start.theta = PI;
    state = keplerian_to_cartesian(kep_targ_start, EARTH_MU, plan.final_mass, t_elapsed);

    double dt_targ = dur_targ / (double)n_targ;
    for (int i = 1; i <= n_targ; i++) {
        state = propagate_step_rk4(state, dt_targ, &sc_pass, &forces);
        t_elapsed += dt_targ;
        char t_str[40];
        format_oem_utc_timestamp(t_elapsed, t_str, sizeof(t_str));
        fprintf(fp, "%s  %14.6f  %14.6f  %14.6f  %12.7f  %12.7f  %12.7f\n",
                t_str,
                state.position.x, state.position.y, state.position.z,
                state.velocity.x, state.velocity.y, state.velocity.z);
    }

    fclose(fp);
    return true;
}

bool export_web_mission_json(const char *filename,
                             OrbitReport rep_init,
                             OptimalTransferPlan plan,
                             OrbitReport rep_target,
                             RendezvousPlan rdv,
                             MonteCarloReport mc,
                             ADCSReport adcs,
                             ThermalReport tcs,
                             TelecomReport telecom,
                             SystemEngineeringReport sys,
                             NavigationReport nav,
                             CAMReport cam,
                             EndOfLifeReport eol,
                             CislunarReport cislunar,
                             StructuralReport csi,
                             FDIRReport fdir,
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

    fprintf(fp, "  \"adcs_summary\": {\n");
    fprintf(fp, "    \"inertia_kgm2\": [%.1f, %.1f, %.1f], \"rw_h_max_nms\": %.2f,\n",
            adcs.inertia_kgm2.x, adcs.inertia_kgm2.y, adcs.inertia_kgm2.z, adcs.rw_h_max_nms);
    fprintf(fp, "    \"slew_time_s\": %.1f, \"slew_angle_b1_deg\": %.1f, \"slew_angle_b2_deg\": %.1f,\n",
            adcs.slew_time_s, adcs.slew_angle_b1_deg, adcs.slew_angle_b2_deg);
    fprintf(fp, "    \"steady_point_err_deg\": %.4f, \"steady_point_err_arcsec\": %.1f,\n",
            adcs.steady_point_err_deg, adcs.steady_point_err_arcsec);
    fprintf(fp, "    \"gg_torque_mnm\": %.3f, \"aero_torque_mnm\": %.3f, \"srp_torque_mnm\": %.3f,\n",
            adcs.gg_torque_mnm, adcs.aero_torque_mnm, adcs.srp_torque_mnm);
    fprintf(fp, "    \"peak_rw_torque_mnm\": %.2f, \"max_rw_mom_nms\": %.3f, \"max_rw_sat_pct\": %.1f,\n",
            adcs.peak_rw_torque_mnm, adcs.max_rw_mom_nms, adcs.max_rw_sat_pct);
    fprintf(fp, "    \"desat_cycles\": %d, \"mean_power_w\": %.1f, \"peak_power_w\": %.1f\n  },\n",
            adcs.desat_cycles, adcs.mean_power_w, adcs.peak_power_w);

    fprintf(fp, "  \"thermal_rad_summary\": {\n");
    fprintf(fp, "    \"coating_name\": \"%s\", \"alpha_solar\": %.2f, \"eps_ir\": %.2f, \"al_shield_mm\": %.1f,\n",
            tcs.coating.name, tcs.alpha_solar, tcs.eps_ir, tcs.al_shield_mm);
    fprintf(fp, "    \"thermal_mass_j_k\": %.1f, \"ext_thermal_mass_j_k\": %.1f, \"rad_area_m2\": %.2f,\n",
            tcs.thermal_mass_j_k, tcs.ext_thermal_mass_j_k, tcs.rad_area_m2);
    fprintf(fp, "    \"t_min_op_c\": %.1f, \"t_max_op_c\": %.1f,\n",
            tcs.t_min_op_c, tcs.t_max_op_c);
    fprintf(fp, "    \"t_min_c\": %.2f, \"t_max_c\": %.2f, \"min_temp_c\": %.2f, \"max_temp_c\": %.2f, \"mean_temp_c\": %.2f,\n",
            tcs.min_temp_c, tcs.max_temp_c, tcs.min_temp_c, tcs.max_temp_c, tcs.mean_temp_c);
    fprintf(fp, "    \"ext_min_temp_c\": %.2f, \"ext_max_temp_c\": %.2f, \"eq_sun_temp_c\": %.2f, \"eq_ecl_temp_c\": %.2f,\n",
            tcs.ext_min_temp_c, tcs.ext_max_temp_c, tcs.eq_sun_temp_c, tcs.eq_ecl_temp_c);
    fprintf(fp, "    \"mean_q_sol_w\": %.1f, \"mean_q_alb_w\": %.1f, \"mean_q_ir_w\": %.1f, \"mean_q_rad_w\": %.1f,\n",
            tcs.mean_q_sol_w, tcs.mean_q_alb_w, tcs.mean_q_ir_w, tcs.mean_q_rad_w);
    fprintf(fp, "    \"peak_heater_w\": %.1f, \"heater_duty_pct\": %.1f, \"heater_energy_wh\": %.2f,\n",
            tcs.peak_heater_w, tcs.heater_duty_pct, tcs.heater_energy_wh);
    fprintf(fp, "    \"tid_krad_yr\": %.2f, \"peak_dose_urad_h\": %.1f, \"seu_rate_day\": %.2f, \"avionics_life_yr\": %.1f,\n",
            tcs.tid_krad_yr, tcs.peak_dose_urad_h, tcs.seu_rate_day, tcs.avionics_life_yr);
    fprintf(fp, "    \"saa_encountered\": %d, \"saa_status\": \"%s\"\n  },\n",
            tcs.saa_encountered ? 1 : 0, tcs.saa_status);

    fprintf(fp, "  \"thermal_summary\": {\n");
    fprintf(fp, "    \"coating_name\": \"%s\", \"alpha_solar\": %.2f, \"eps_ir\": %.2f, \"al_shield_mm\": %.1f,\n",
            tcs.coating.name, tcs.alpha_solar, tcs.eps_ir, tcs.al_shield_mm);
    fprintf(fp, "    \"thermal_mass_j_k\": %.1f, \"ext_thermal_mass_j_k\": %.1f, \"rad_area_m2\": %.2f,\n",
            tcs.thermal_mass_j_k, tcs.ext_thermal_mass_j_k, tcs.rad_area_m2);
    fprintf(fp, "    \"t_min_op_c\": %.1f, \"t_max_op_c\": %.1f,\n",
            tcs.t_min_op_c, tcs.t_max_op_c);
    fprintf(fp, "    \"t_min_c\": %.2f, \"t_max_c\": %.2f, \"min_temp_c\": %.2f, \"max_temp_c\": %.2f, \"mean_temp_c\": %.2f,\n",
            tcs.min_temp_c, tcs.max_temp_c, tcs.min_temp_c, tcs.max_temp_c, tcs.mean_temp_c);
    fprintf(fp, "    \"ext_min_temp_c\": %.2f, \"ext_max_temp_c\": %.2f, \"eq_sun_temp_c\": %.2f, \"eq_ecl_temp_c\": %.2f,\n",
            tcs.ext_min_temp_c, tcs.ext_max_temp_c, tcs.eq_sun_temp_c, tcs.eq_ecl_temp_c);
    fprintf(fp, "    \"mean_q_sol_w\": %.1f, \"mean_q_alb_w\": %.1f, \"mean_q_ir_w\": %.1f, \"mean_q_rad_w\": %.1f,\n",
            tcs.mean_q_sol_w, tcs.mean_q_alb_w, tcs.mean_q_ir_w, tcs.mean_q_rad_w);
    fprintf(fp, "    \"peak_heater_w\": %.1f, \"heater_duty_pct\": %.1f, \"heater_energy_wh\": %.2f,\n",
            tcs.peak_heater_w, tcs.heater_duty_pct, tcs.heater_energy_wh);
    fprintf(fp, "    \"tid_krad_yr\": %.2f, \"peak_dose_urad_h\": %.1f, \"seu_rate_day\": %.2f, \"avionics_life_yr\": %.1f,\n",
            tcs.tid_krad_yr, tcs.peak_dose_urad_h, tcs.seu_rate_day, tcs.avionics_life_yr);
    fprintf(fp, "    \"saa_encountered\": %d, \"saa_status\": \"%s\"\n  },\n",
            tcs.saa_encountered ? 1 : 0, tcs.saa_status);

    fprintf(fp, "  \"telecom_summary\": {\n");
    fprintf(fp, "    \"band_id\": %d, \"band_name\": \"%s\", \"modulation_scheme\": \"%s\",\n",
            (int)telecom.band, telecom.band_name, telecom.modulation_scheme);
    fprintf(fp, "    \"freq_ghz\": %.2f, \"tx_power_w\": %.1f, \"tx_power_dbw\": %.2f,\n",
            telecom.freq_ghz, telecom.tx_power_w, telecom.tx_power_dbw);
    fprintf(fp, "    \"tx_gain_dbi\": %.1f, \"rx_dish_m\": %.1f, \"rx_gain_dbi\": %.1f, \"eirp_dbw\": %.2f,\n",
            telecom.tx_gain_dbi, telecom.rx_dish_m, telecom.rx_gain_dbi, telecom.eirp_dbw);
    fprintf(fp, "    \"sys_noise_temp_k\": %.1f, \"n0_dbw_hz\": %.2f, \"zenith_atm_loss_db\": %.2f, \"pointing_loss_db\": %.2f,\n",
            telecom.sys_noise_temp_k, telecom.n0_dbw_hz, telecom.zenith_atm_loss_db, telecom.pointing_loss_db);
    fprintf(fp, "    \"lock_threshold_db\": %.2f, \"req_ebn0_db\": %.2f,\n",
            telecom.lock_threshold_db, telecom.req_ebn0_db);
    fprintf(fp, "    \"nom_bitrate_mbps\": %.1f, \"max_achieved_mbps\": %.2f,\n",
            telecom.nom_bitrate_mbps, telecom.max_achieved_mbps);
    fprintf(fp, "    \"mean_fspl_db\": %.2f, \"max_fspl_db\": %.2f, \"peak_ebn0_db\": %.2f,\n",
            telecom.mean_fspl_db, telecom.max_fspl_db, telecom.peak_ebn0_db);
    fprintf(fp, "    \"min_margin_db\": %.2f, \"peak_margin_db\": %.2f,\n",
            telecom.min_margin_db, telecom.peak_margin_db);
    fprintf(fp, "    \"pass_windows\": %d, \"contact_time_min\": %.2f, \"mission_data_gb\": %.3f, \"total_data_gb_24h\": %.2f\n  },\n",
            telecom.pass_windows, telecom.contact_time_min, telecom.mission_data_gb, telecom.total_data_gb_24h);

    fprintf(fp, "  \"systems_engineering\": {\n");
    fprintf(fp, "    \"half_fov_deg\": %.2f, \"opt_aperture_m\": %.2f, \"wavelength_nm\": %.1f,\n",
            sys.spec.half_fov_deg, sys.spec.opt_aperture_m, sys.spec.wavelength_nm);
    fprintf(fp, "    \"target_alt_km\": %.2f, \"earth_central_angle_deg\": %.3f,\n",
            sys.target_alt_km, sys.earth_central_angle_deg);
    fprintf(fp, "    \"swath_width_km\": %.2f, \"footprint_radius_km\": %.2f, \"gsd_nadir_m\": %.2f, \"area_coverage_rate_km2_s\": %.1f,\n",
            sys.swath_width_km, sys.footprint_radius_km, sys.gsd_nadir_m, sys.area_coverage_rate_km2_s);
    fprintf(fp, "    \"walker_total_sats\": %d, \"walker_planes\": %d, \"walker_phasing\": %d, \"sats_per_plane\": %d, \"inc_deg\": %.2f,\n",
            sys.spec.walker_total_sats, sys.spec.walker_planes, sys.spec.walker_phasing, sys.sats_per_plane, sys.inc_deg);
    fprintf(fp, "    \"mean_revisit_min\": %.2f, \"global_coverage_pct\": %.2f, \"mean_vis_sats\": %.2f, \"pdop_geom\": %.2f,\n",
            sys.mean_revisit_min, sys.global_coverage_pct, sys.mean_vis_sats, sys.pdop_geom);
    fprintf(fp, "    \"wet_mass_kg\": %.2f, \"dry_mass_kg\": %.2f, \"dry_wet_ratio_pct\": %.1f,\n",
            sys.wet_mass_kg, sys.dry_mass_kg, sys.dry_wet_ratio_pct);
    fprintf(fp, "    \"propellant_transfer_kg\": %.2f, \"propellant_margin_kg\": %.2f, \"total_propellant_kg\": %.2f,\n",
            sys.propellant_transfer_kg, sys.propellant_margin_kg, sys.total_propellant_kg);
    fprintf(fp, "    \"mass_breakdown\": [\n");
    fprintf(fp, "      {\"subsystem\":\"Payload (Optical/SAR)\", \"mass_kg\":%.2f, \"pct_dry\":28.0},\n", sys.mass_payload_kg);
    fprintf(fp, "      {\"subsystem\":\"Structure & Mechanisms\", \"mass_kg\":%.2f, \"pct_dry\":22.0},\n", sys.mass_struct_kg);
    fprintf(fp, "      {\"subsystem\":\"EPS Solar & Li-Ion\", \"mass_kg\":%.2f, \"pct_dry\":19.0},\n", sys.mass_eps_kg);
    fprintf(fp, "      {\"subsystem\":\"Propulsion (Dry)\", \"mass_kg\":%.2f, \"pct_dry\":11.0},\n", sys.mass_prop_dry_kg);
    fprintf(fp, "      {\"subsystem\":\"ADCS (RW & Star Track)\", \"mass_kg\":%.2f, \"pct_dry\":9.0},\n", sys.mass_adcs_kg);
    fprintf(fp, "      {\"subsystem\":\"TT&C & OBDH Avionics\", \"mass_kg\":%.2f, \"pct_dry\":6.0},\n", sys.mass_ttc_obdh_kg);
    fprintf(fp, "      {\"subsystem\":\"Thermal Control (TCS)\", \"mass_kg\":%.2f, \"pct_dry\":5.0},\n", sys.mass_tcs_kg);
    fprintf(fp, "      {\"subsystem\":\"Transfer Propellant\", \"mass_kg\":%.2f, \"pct_dry\":0.0},\n", sys.propellant_transfer_kg);
    fprintf(fp, "      {\"subsystem\":\"GNC TCM & SK Reserve\", \"mass_kg\":%.2f, \"pct_dry\":0.0}\n", sys.propellant_margin_kg);
    fprintf(fp, "    ],\n");
    fprintf(fp, "    \"cost_platform_meur\": %.2f, \"cost_payload_meur\": %.2f, \"cost_spacecraft_meur\": %.2f,\n",
            sys.cost_platform_meur, sys.cost_payload_meur, sys.cost_spacecraft_meur);
    fprintf(fp, "    \"cost_launch_meur\": %.2f, \"cost_ops_5yr_meur\": %.2f,\n",
            sys.cost_launch_meur, sys.cost_ops_5yr_meur);
    fprintf(fp, "    \"total_mission_cost_meur\": %.2f, \"constellation_cost_meur\": %.2f,\n",
            sys.total_mission_cost_meur, sys.constellation_cost_meur);
    fprintf(fp, "    \"walker_subpoints\": [\n");
    for (int i = 0; i < sys.num_walker_sats; i++) {
        WalkerSatPoint w = sys.walker_subpoints[i];
        fprintf(fp, "      {\"sat_id\":%d,\"plane_id\":%d,\"lat\":%.2f,\"lon\":%.2f,\"raan_deg\":%.1f,\"u_deg\":%.1f}%s\n",
                w.sat_id, w.plane_id, w.lat_deg, w.lon_deg, w.raan_deg, w.arg_lat_deg,
                (i < sys.num_walker_sats - 1) ? "," : "");
    }
    fprintf(fp, "    ]\n  },\n");

    fprintf(fp, "  \"ekf_summary\": {\n");
    fprintf(fp, "    \"gps_sigma_pos_m\": %.2f, \"gps_sigma_vel_ms\": %.3f,\n",
            nav.spec.gps_sigma_pos_m, nav.spec.gps_sigma_vel_ms);
    fprintf(fp, "    \"init_sigma_pos_m\": %.1f, \"init_sigma_vel_ms\": %.2f, \"proc_noise_accel_ms2\": %.5f,\n",
            nav.spec.init_sigma_pos_m, nav.spec.init_sigma_vel_ms, nav.spec.proc_noise_accel_ms2);
    fprintf(fp, "    \"raw_pos_rmse_m\": %.2f, \"steady_pos_rmse_m\": %.2f,\n",
            nav.raw_pos_rmse_m, nav.steady_pos_rmse_m);
    fprintf(fp, "    \"steady_vel_rmse_ms\": %.4f, \"steady_vel_rmse_cms\": %.2f,\n",
            nav.steady_vel_rmse_ms, nav.steady_vel_rmse_cms);
    fprintf(fp, "    \"final_cov_3sigma_m\": %.2f, \"final_cov_vel_3sigma_ms\": %.4f,\n",
            nav.final_cov_3sigma_m, nav.final_cov_vel_3sigma_ms);
    fprintf(fp, "    \"convergence_time_min\": %.2f, \"noise_reduction_pct\": %.1f,\n",
            nav.convergence_time_min, nav.noise_reduction_pct);
    fprintf(fp, "    \"consistency_pct\": %.1f, \"consistency_3sigma_pct\": %.1f, \"mean_kalman_gain\": %.4f, \"final_kalman_gain_norm\": %.4f, \"filter_converged\": 1\n  },\n",
            nav.consistency_pct, nav.consistency_pct, nav.mean_kalman_gain, nav.mean_kalman_gain);

    fprintf(fp, "  \"ssa_cam_summary\": {\n");
    fprintf(fp, "    \"debris_id\": \"%s\",\n", cam.event.debris_id);
    fprintf(fp, "    \"rel_vel_kms\": %.2f, \"hbr_m\": %.2f,\n", cam.event.rel_vel_kms, cam.event.hbr_m);
    fprintf(fp, "    \"sigma_xi_m\": %.2f, \"sigma_zeta_m\": %.2f,\n", cam.event.sigma_xi_m, cam.event.sigma_zeta_m);
    fprintf(fp, "    \"lead_orbits\": %.2f, \"tca_lead_time_min\": %.2f,\n", cam.event.lead_orbits, cam.tca_lead_time_min);
    fprintf(fp, "    \"pre_xi_m\": %.2f, \"pre_zeta_m\": %.2f, \"pre_miss_m\": %.2f,\n", cam.pre_xi_m, cam.pre_zeta_m, cam.pre_miss_m);
    fprintf(fp, "    \"pre_pc\": %.4e, \"pre_log10_pc\": %.3f, \"pre_mahalanobis_d\": %.2f,\n", cam.pre_pc, cam.pre_log10_pc, cam.pre_mahalanobis_d);
    fprintf(fp, "    \"req_dv_cam_ms\": %.4f, \"req_dv_cam_cms\": %.2f, \"cam_fuel_g\": %.2f, \"delta_a_m\": %.2f,\n",
            cam.req_dv_cam_ms, cam.req_dv_cam_cms, cam.cam_fuel_g, cam.delta_a_m);
    fprintf(fp, "    \"post_xi_m\": %.2f, \"post_zeta_m\": %.2f, \"post_miss_m\": %.2f,\n", cam.post_xi_m, cam.post_zeta_m, cam.post_miss_m);
    fprintf(fp, "    \"post_pc\": %.4e, \"post_log10_pc\": %.3f, \"post_mahalanobis_d\": %.2f,\n", cam.post_pc, cam.post_log10_pc, cam.post_mahalanobis_d);
    fprintf(fp, "    \"alert_status\": \"%s\"\n  },\n", cam.alert_status);

    fprintf(fp, "  \"eol_reentry_summary\": {\n");
    fprintf(fp, "    \"entry_alt_km\": %.1f, \"entry_fpa_deg\": %.2f, \"nose_radius_m\": %.2f,\n",
            eol.entry_alt_km, eol.entry_fpa_deg, eol.nose_radius_m);
    fprintf(fp, "    \"ballistic_coeff_kgm2\": %.2f, \"deorbit_perigee_km\": %.1f,\n",
            eol.ballistic_coeff_kgm2, eol.deorbit_perigee_km);
    fprintf(fp, "    \"deorbit_dv_ms\": %.2f, \"deorbit_fuel_kg\": %.2f, \"passive_decay_years\": %.1f,\n",
            eol.deorbit_dv_ms, eol.deorbit_fuel_kg, eol.passive_decay_years);
    fprintf(fp, "    \"compliance_status\": \"%s\",\n", eol.compliance_status);
    fprintf(fp, "    \"peak_heat_flux_wcm2\": %.2f, \"peak_heat_alt_km\": %.2f, \"total_heat_load_jcm2\": %.1f,\n",
            eol.peak_heat_flux_wcm2, eol.peak_heat_alt_km, eol.total_heat_load_jcm2);
    fprintf(fp, "    \"peak_decel_g\": %.2f, \"peak_decel_alt_km\": %.2f, \"peak_qdyn_kpa\": %.2f,\n",
            eol.peak_decel_g, eol.peak_decel_alt_km, eol.peak_qdyn_kpa);
    fprintf(fp, "    \"breakup_alt_km\": %.2f, \"breakup_vel_kms\": %.2f,\n",
            eol.breakup_alt_km, eol.breakup_vel_kms);
    fprintf(fp, "    \"total_entry_time_sec\": %.1f, \"total_downrange_km\": %.1f\n  },\n",
            eol.total_entry_time_sec, eol.total_downrange_km);

    fprintf(fp, "  \"cislunar_summary\": {\n");
    fprintf(fp, "    \"target_regime\": %d, \"regime_name\": \"%s\",\n",
            cislunar.spec.target_regime, cislunar.regime_name);
    fprintf(fp, "    \"mu_em\": %.10f, \"l_star_km\": %.1f, \"t_star_sec\": %.2f, \"v_star_kms\": %.6f,\n",
            cislunar.mu_em, cislunar.l_star_km, cislunar.t_star_sec, cislunar.v_star_kms);
    fprintf(fp, "    \"earth_x_km\": %.2f, \"moon_x_km\": %.2f,\n",
            cislunar.earth_x_km, cislunar.moon_x_km);
    fprintf(fp, "    \"l1_x_km\": %.2f, \"l2_x_km\": %.2f, \"l3_x_km\": %.2f,\n",
            cislunar.l1_x_km, cislunar.l2_x_km, cislunar.l3_x_km);
    fprintf(fp, "    \"l4_x_km\": %.2f, \"l4_y_km\": %.2f, \"l5_x_km\": %.2f, \"l5_y_km\": %.2f,\n",
            cislunar.l4_x_km, cislunar.l4_y_km, cislunar.l5_x_km, cislunar.l5_y_km);
    fprintf(fp, "    \"cj_l1\": %.6f, \"cj_l2\": %.6f, \"cj_l3\": %.6f, \"cj_l4_l5\": %.6f,\n",
            cislunar.cj_l1, cislunar.cj_l2, cislunar.cj_l3, cislunar.cj_l4_l5);
    fprintf(fp, "    \"dv_tli_ms\": %.2f, \"dv_loi_ms\": %.2f, \"total_cislunar_dv_ms\": %.2f, \"cislunar_fuel_kg\": %.2f,\n",
            cislunar.dv_tli_ms, cislunar.dv_loi_ms, cislunar.total_cislunar_dv_ms, cislunar.cislunar_fuel_kg);
    fprintf(fp, "    \"tof_days\": %.3f, \"halo_period_days\": %.3f,\n",
            cislunar.tof_days, cislunar.halo_period_days);
    fprintf(fp, "    \"halo_ax_km\": %.1f, \"halo_ay_km\": %.1f, \"halo_az_km\": %.1f, \"perilune_alt_km\": %.1f,\n",
            cislunar.halo_ax_km, cislunar.halo_ay_km, cislunar.halo_az_km, cislunar.perilune_alt_km);
    fprintf(fp, "    \"jacobi_c\": %.8f, \"jacobi_halo_c\": %.8f, \"delta_cj\": %.3e, \"delta_cj_rel\": %.3e\n  },\n",
            cislunar.jacobi_c, cislunar.jacobi_halo_c, cislunar.delta_cj_rel * cislunar.jacobi_c, cislunar.delta_cj_rel);

    fprintf(fp, "  \"structural_csi_summary\": {\n");
    fprintf(fp, "    \"flex_freq_hz\": %.3f, \"flex_freq_mode2_hz\": %.3f, \"flex_damping_zeta\": %.4f,\n",
            csi.flex_freq_hz, csi.flex_freq_mode2_hz, csi.spec.flex_damping_zeta);
    fprintf(fp, "    \"slosh_freq_hz\": %.3f, \"slosh_omega_rads\": %.3f, \"slosh_damping_zeta\": %.4f,\n",
            csi.slosh_freq_hz, csi.slosh_omega_rads, csi.spec.slosh_damping_zeta);
    fprintf(fp, "    \"total_fuel_mass_kg\": %.2f, \"slosh_mass_kg\": %.2f, \"fixed_fuel_mass_kg\": %.2f,\n",
            csi.total_fuel_mass_kg, csi.slosh_mass_kg, csi.fixed_fuel_mass_kg);
    fprintf(fp, "    \"tank_radius_m\": %.2f, \"tank_fill_ratio\": %.2f, \"slosh_pendulum_len_m\": %.3f, \"eff_accel_ms2\": %.3f,\n",
            csi.spec.tank_radius_m, csi.spec.tank_fill_ratio, csi.slosh_pendulum_len_m, csi.eff_accel_ms2);
    fprintf(fp, "    \"array_span_m\": %.2f, \"modal_coupling_gamma\": %.2f,\n",
            csi.spec.array_span_m, csi.modal_coupling_gamma);
    fprintf(fp, "    \"peak_tip_defl_mm\": %.2f, \"peak_slosh_disp_mm\": %.2f, \"peak_slosh_force_n\": %.2f, \"peak_slosh_torque_mnm\": %.2f,\n",
            csi.peak_tip_defl_mm, csi.peak_slosh_disp_mm, csi.peak_slosh_force_n, csi.peak_slosh_torque_mnm);
    fprintf(fp, "    \"rms_jitter_arcsec\": %.2f, \"peak_jitter_arcsec\": %.2f, \"settling_time_sec\": %.1f,\n",
            csi.rms_jitter_arcsec, csi.peak_jitter_arcsec, csi.settling_time_sec);
    fprintf(fp, "    \"adcs_bandwidth_hz\": %.2f, \"notch_center_hz\": %.3f, \"notch_attenuation_db\": %.2f,\n",
            csi.adcs_bandwidth_hz, csi.notch_center_hz, csi.notch_attenuation_db);
    fprintf(fp, "    \"peak_res_unfilt_db\": %.2f, \"peak_res_filt_db\": %.2f, \"csi_gain_margin_db\": %.2f,\n",
            csi.peak_res_unfilt_db, csi.peak_res_filt_db, csi.csi_gain_margin_db);
    fprintf(fp, "    \"notch_status\": \"%s\",\n", csi.notch_status);
    fprintf(fp, "    \"launch_axial_freq_hz\": %.1f, \"launch_lateral_freq_hz\": %.1f,\n",
            csi.launch_axial_freq_hz, csi.launch_lateral_freq_hz);
    fprintf(fp, "    \"qsl_axial_g\": %.1f, \"qsl_lateral_g\": %.1f,\n",
            csi.qsl_axial_g, csi.qsl_lateral_g);
    fprintf(fp, "    \"max_von_mises_mpa\": %.1f, \"allowable_stress_mpa\": %.1f, \"launch_mos_pct\": %.1f,\n",
            csi.max_von_mises_mpa, csi.allowable_stress_mpa, csi.launch_mos_pct);
    fprintf(fp, "    \"launch_compliance\": \"%s\"\n  },\n", csi.launch_compliance);

    fprintf(fp, "  \"fdir_summary\": {\n");
    fprintf(fp, "    \"fault_scenario\": %d, \"fault_name\": \"%s\", \"fault_subsystem\": \"%s\",\n",
            fdir.spec.fault_scenario, fdir.fault_name, fdir.fault_subsystem);
    fprintf(fp, "    \"primary_unit_name\": \"%s\", \"redundant_unit_name\": \"%s\",\n",
            fdir.primary_unit_name, fdir.redundant_unit_name);
    fprintf(fp, "    \"fault_inject_time_s\": %.2f, \"detection_latency_s\": %.2f, \"isolation_switchover_s\": %.2f,\n",
            fdir.fault_inject_time_s, fdir.detection_latency_s, fdir.isolation_switchover_s);
    fprintf(fp, "    \"safe_mode_duration_s\": %.2f, \"total_recovery_time_s\": %.2f,\n",
            fdir.safe_mode_duration_s, fdir.total_recovery_time_s);
    fprintf(fp, "    \"peak_att_err_deg\": %.3f, \"steady_post_err_deg\": %.4f,\n",
            fdir.peak_att_err_deg, fdir.steady_post_err_deg);
    fprintf(fp, "    \"nominal_bus_voltage_v\": %.2f, \"min_bus_voltage_v\": %.2f, \"min_sys_health_pct\": %.1f,\n",
            fdir.nominal_bus_voltage_v, fdir.min_bus_voltage_v, fdir.min_sys_health_pct);
    fprintf(fp, "    \"redundancy_status\": \"%s\",\n", fdir.redundancy_status);
    fprintf(fp, "    \"recovery_outcome\": \"%s\",\n", fdir.recovery_outcome);
    fprintf(fp, "    \"ccsds_oem_filename\": \"%s\", \"ccsds_oem_status\": \"%s\", \"oem_state_vectors\": %d,\n",
            fdir.ccsds_oem_filename, fdir.ccsds_oem_status, fdir.oem_state_vectors);
    fprintf(fp, "    \"tle_norad_id\": %d, \"tle_sat_name\": \"%s\"\n  },\n",
            fdir.tle_norad_id, fdir.tle_sat_name);

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
        kep_trans.omega = (r2 >= r1) ? rep_init.initial_kep.omega : (rep_init.initial_kep.omega + PI);
        kep_trans.theta = (r2 >= r1) ? 0.0 : PI;
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
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_adcs\": ");
    write_adcs_trajectory_json_array(fp, adcs);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_thermal\": ");
    write_thermal_trajectory_json_array(fp, tcs);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_telecom\": ");
    write_telecom_trajectory_json_array(fp, telecom);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_ekf\": ");
    write_ekf_trajectory_json_array(fp, nav);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"cam_tradeoff\": ");
    write_cam_tradeoff_json_array(fp, cam);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_reentry\": ");
    write_reentry_trajectory_json_array(fp, eol);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_cislunar_transfer\": ");
    write_cr3bp_trajectory_json_array(fp, cislunar.traj_cislunar_transfer, cislunar.num_transfer_samples);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_halo_orbit\": ");
    write_cr3bp_trajectory_json_array(fp, cislunar.traj_halo_orbit, cislunar.num_halo_samples);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_csi\": ");
    write_csi_trajectory_json_array(fp, csi);
    fprintf(fp, ",\n");

    fprintf(fp, "  \"traj_fdir\": ");
    write_fdir_trajectory_json_array(fp, fdir);
    fprintf(fp, "\n}\n");

    fclose(fp);
    return true;
}
