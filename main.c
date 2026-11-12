#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "mission_analysis.h"
#include "telemetry.h"

int main(void) {
    printf("======================================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 21: CISLUNAR CR3BP & LUNAR HALO ORBIT  \n");
    printf("======================================================================\n\n");

    // 1. INPUT UTENTE: Massa, Propellente e Vincolo di Tempo Massimo
    double initial_mass = 500.0; // kg
    PropellantType selected_prop = PROP_LIQUID_MMH_NTO; // Liquido Bipropellente (Isp = 316 s)
    double max_time_limit_sec = 0.0; // 0 = Nessun limite (cerca il minimo Delta-V assoluto)

    // 2. INPUT ORBITA INIZIALE (Scelta tra INPUT_KEPLERIAN o INPUT_CARTESIAN)
    OrbitInputSpec input_init = {
        .mode = INPUT_KEPLERIAN,
        .kep  = {.a = EARTH_RADIUS + 400.0, .e = 0.001, .i = 51.6 * DEG2RAD, .RAAN = 30.0 * DEG2RAD, .omega = 0.0, .theta = 0.0}
    };

    // 3. INPUT ORBITA FINALE TARGET
    OrbitInputSpec input_target = {
        .mode = INPUT_KEPLERIAN,
        .kep  = {.a = EARTH_RADIUS + 650.0, .e = 0.001, .i = 53.0 * DEG2RAD, .RAAN = 30.0 * DEG2RAD, .omega = 0.0, .theta = 0.0}
    };

    // Risoluzione automatica coordinate
    KeplerianElements orb_initial = resolve_user_orbit_input(input_init, EARTH_MU);
    KeplerianElements orb_target  = resolve_user_orbit_input(input_target, EARTH_MU);

    PropellantSpec prop = get_propellant_specs(selected_prop);
    SpacecraftConfig sc = {
        .drag_coeff = 2.2,
        .drag_area  = 4.0,
        .thruster   = {.is_active = false, .thrust = 50.0, .isp = prop.isp}
    };

    GroundStation milano = {
        .name = "Milano-Bovisa",
        .coord = {.latitude = 45.502 * DEG2RAD, .longitude = 9.156 * DEG2RAD, .altitude = 0.140},
        .min_elevation = 10.0 * DEG2RAD
    };

    // Calcolo Report ed Esecuzione Ottimizzatore Globale di Trasferimento
    OrbitReport rep_init = evaluate_orbit_regime(orb_initial, initial_mass, &sc, &milano);
    OptimalTransferPlan plan = optimize_orbital_transfer(orb_initial, orb_target, initial_mass, selected_prop, max_time_limit_sec, EARTH_MU);
    OrbitReport rep_target = evaluate_orbit_regime(orb_target, plan.final_mass, &sc, &milano);

    // 4. SPRINT 12: TERMINAL RENDEZVOUS & DOCKING (CLOHESSY-WILTSHIRE / HILL LVLH)
    Vector3D r0_lvlh_m  = {-250.0, -2000.0, 100.0};
    Vector3D v0_lvlh_ms = {0.0, 0.0, 0.0};
    double t_rdv_sec    = 45.0 * 60.0; // 45 minuti di approccio finale

    RendezvousPlan rdv = compute_cw_rendezvous(orb_target.a, r0_lvlh_m, v0_lvlh_ms, t_rdv_sec,
                                               plan.final_mass, plan.prop_info.isp, EARTH_MU);

    // 5. SPRINT 13: ANALISI DI DISPERSIONE MONTE CARLO (N=100) E CALCOLO MID-COURSE TCM
    double sigma_dv_pct    = 1.0; // 1.0% (1-sigma) Thrust Magnitude Error
    double sigma_point_deg = 0.5; // 0.5 deg (1-sigma) Nozzle Attitude Pointing Error
    MonteCarloReport mc = run_monte_carlo_dispersion(orb_initial, orb_target, plan,
                                                     sigma_dv_pct, sigma_point_deg, MAX_MC_RUNS, EARTH_MU);

    // 6. SPRINT 14: DINAMICA D'ASSETTO 6-DOF (ADCS), SLEW QUATERNIONI & REACTION WHEELS
    Vector3D sc_inertia_kgm2 = {120.0, 110.0, 85.0}; // Tensore principale d'inerzia J = diag(Jxx, Jyy, Jzz) [kg*m^2]
    double rw_h_max_nms      = 10.0;                 // Capacità momento angolare Reaction Wheels [N*m*s]
    ADCSReport adcs = simulate_adcs_profile(orb_initial, orb_target, plan, sc_inertia_kgm2, rw_h_max_nms, EARTH_MU);

    // 7. SPRINT 15: THERMAL CONTROL SYSTEM (TCS) 2-NODE MODEL & RADIATION ENVIRONMENT (VAN ALLEN / SAA)
    SurfaceOpticalProps coating = get_surface_optical_props(COATING_WHITE_PAINT_SG121);
    double al_shield_mm = 4.0;  // Spessore schermatura anti-radiazioni in Alluminio [mm]
    double t_min_op_c   = -10.0;
    double t_max_op_c   = 45.0;
    ThermalRadReport tcs = simulate_thermal_and_radiation(orb_initial, orb_target, plan,
                                                          coating, al_shield_mm, t_min_op_c, t_max_op_c, EARTH_MU);

    // 8. SPRINT 16: RF TELECOMMUNICATION LINK BUDGET, Eb/N0 MARGIN & GROUND STATION DATA DOWNLINK
    RFBandType rf_band = BAND_X_8_4GHZ; // X-Band (8.4 GHz Science Downlink, 50 Mbps nominal)
    double tx_power_w  = 15.0;          // Potenza trasmettitore di bordo P_tx = 15.0 W (11.76 dBW)
    TelecomReport telecom = evaluate_rf_link_budget(orb_initial, orb_target, plan, rep_target,
                                                    rf_band, tx_power_w, &milano, EARTH_MU);

    // 9. SPRINT 17: PAYLOAD SENSOR SWATH, WALKER DELTA CONSTELLATION & MASS/COST BUDGET
    PayloadAndConstellationSpec sys_spec = {
        .half_fov_deg      = 25.0,  // Semiapertura angolare sensore theta_(FOV/2) = 25.0 deg
        .opt_aperture_m    = 0.35,  // Diametro apertura ottica D_opt = 0.35 m
        .wavelength_nm     = 550.0, // Lunghezza d'onda visibile lambda = 550 nm
        .walker_total_sats = 24,    // Costellazione Walker T = 24 satelliti
        .walker_planes     = 3,     // P = 3 piani orbitali (8 satelliti per piano)
        .walker_phasing    = 1      // F = 1 fasatura inter-piano
    };
    SystemEngineeringReport sys = evaluate_payload_constellation_and_cost(orb_target, rep_target, plan, rdv, mc,
                                                                          sys_spec, &milano, EARTH_MU);

    // 10. SPRINT 18: AUTONOMOUS ONBOARD NAVIGATION & REAL-TIME ORBIT DETERMINATION VIA EKF
    NavigationEKFSpec ekf_spec = {
        .gps_sigma_pos_m      = 15.0,   // Rumore di misura GNSS 1-sigma sulla posizione [m]
        .gps_sigma_vel_ms     = 0.15,   // Rumore di misura GNSS 1-sigma sulla velocità [m/s]
        .init_sigma_pos_m     = 500.0,  // Incertezza iniziale a priori 1-sigma sulla posizione [m]
        .init_sigma_vel_ms    = 2.50,   // Incertezza iniziale a priori 1-sigma sulla velocità [m/s]
        .proc_noise_accel_ms2 = 0.0038  // Rumore di processo dinamico q_a [m/s^2]
    };
    NavigationReport nav = run_onboard_ekf_navigation(orb_initial, orb_target, plan, ekf_spec, &sc, EARTH_MU);

    // 11. SPRINT 19: SPACE SITUATIONAL AWARENESS (SSA), B-PLANE CONJUNCTION & CAM
    ConjunctionEvent deb_event = {
        .debris_id      = "COSMOS-2251 DEB [NORAD #34454]",
        .nominal_miss_m = 180.0,  // Nominal Miss Distance pre-CAM al TCA [m]
        .rel_vel_kms    = 0.0,    // Calcolato dinamicamente dalla velocità orbitale LEO (~12.4 km/s)
        .hbr_m          = 5.0,    // Hard-Body Radius combinato R_sat + R_debris [m]
        .sigma_xi_m     = 85.0,   // Semiasse 1-sigma B-Plane cross-track [m]
        .sigma_zeta_m   = 175.0,  // Semiasse 1-sigma B-Plane along-track projected [m]
        .lead_orbits    = 1.5     // Anticipo manovra di evitamento prima del TCA [orbite]
    };
    CAMReport cam = evaluate_debris_conjunction_and_cam(orb_target, rep_target, plan, deb_event, EARTH_MU);

    // 12. SPRINT 20: END-OF-LIFE (EoL) DEORBIT, 3-DOF RE-ENTRY & SUTTON-GRAVES HEATING
    double entry_fpa_deg = -1.5; // Angolo di rampa all'Entry Interface (h_EI = 120 km) [deg]
    double nose_radius_m = 0.50; // Raggio di curvatura effettivo R_n [m]
    EndOfLifeReport eol = simulate_eol_deorbit_and_reentry(orb_target, rep_target, plan,
                                                           entry_fpa_deg, nose_radius_m, &sc, EARTH_MU);

    // 13. SPRINT 21: CISLUNAR ASTRODYNAMICS, CR3BP, LAGRANGE POINTS & LUNAR HALO ORBIT
    CR3BPSpec cr3bp_spec = {
        .target_regime      = CISLUNAR_L2_HALO, // L2 Far-Side Southern Halo (Lunar Gateway)
        .halo_az_km         = 15000.0,          // Ampiezza verticale Halo A_z = 15000 km
        .leo_parking_alt_km = 400.0             // Quota LEO pre-TLI [km]
    };
    CislunarReport cislunar = simulate_cislunar_cr3bp(orb_initial, plan, cr3bp_spec);

    printf("[OPTIMAL MANEUVER SELECTED]: %s\n", plan.maneuver_name);
    printf("  * Propellant            : %s [%s] (Isp = %.1f s)\n", plan.prop_info.name, plan.prop_info.category, plan.prop_info.isp);
    printf("  * Burn 1 Point (ECI)    : [ %+9.2f, %+9.2f, %+9.2f ] km (Lat %.1f deg, Lon %.1f deg)\n",
           plan.burn1_pos_eci.x, plan.burn1_pos_eci.y, plan.burn1_pos_eci.z,
           plan.burn1_geo.latitude * RAD2DEG, plan.burn1_geo.longitude * RAD2DEG);
    printf("  * Burn 1 Cost & Dir RTN : %7.2f m/s | Dir: [ %+5.3f, %+5.3f, %+5.3f ]\n",
           plan.dv1_mag * 1000.0, plan.burn1_dir_rtn.x, plan.burn1_dir_rtn.y, plan.burn1_dir_rtn.z);
    printf("  * Burn 2 Point (ECI)    : [ %+9.2f, %+9.2f, %+9.2f ] km (Lat %.1f deg, Lon %.1f deg)\n",
           plan.burn2_pos_eci.x, plan.burn2_pos_eci.y, plan.burn2_pos_eci.z,
           plan.burn2_geo.latitude * RAD2DEG, plan.burn2_geo.longitude * RAD2DEG);
    printf("  * Burn 2 Cost & Dir RTN : %7.2f m/s | Dir: [ %+5.3f, %+5.3f, %+5.3f ]\n",
           plan.dv2_mag * 1000.0, plan.burn2_dir_rtn.x, plan.burn2_dir_rtn.y, plan.burn2_dir_rtn.z);
    printf("  * Total Delta-V / Fuel  : %7.2f m/s | Fuel = %.2f kg | TOF = %.2f min\n\n",
           plan.total_dv * 1000.0, plan.propellant_used, plan.transfer_time_sec / 60.0);

    printf("[TERMINAL RENDEZVOUS & PROXIMITY OPS - CLOHESSY-WILTSHIRE LVLH]\n");
    printf("  * Initial LVLH Offset   : dX(R-bar) = %+.1f m, dY(V-bar) = %+.1f m, dZ(H-bar) = %+.1f m (Range = %.1f m)\n",
           rdv.r0_lvlh_m.x, rdv.r0_lvlh_m.y, rdv.r0_lvlh_m.z, rdv.initial_range_m);
    printf("  * CW Targeting Burn 1   : %6.3f m/s | LVLH: [ %+6.3f, %+6.3f, %+6.3f ] m/s\n",
           rdv.dv_inject_ms, rdv.dv_inject_lvlh_ms.x, rdv.dv_inject_lvlh_ms.y, rdv.dv_inject_lvlh_ms.z);
    printf("  * Soft-Docking Burn 2   : %6.3f m/s | LVLH: [ %+6.3f, %+6.3f, %+6.3f ] m/s\n",
           rdv.dv_brake_ms, rdv.dv_brake_lvlh_ms.x, rdv.dv_brake_lvlh_ms.y, rdv.dv_brake_lvlh_ms.z);
    printf("  * Total Rendezvous Cost : %6.3f m/s | Fuel = %.3f kg | Approach TOF = %.1f min\n\n",
           rdv.total_dv_ms, rdv.propellant_used_kg, rdv.transfer_time_sec / 60.0);

    printf("[MONTE CARLO DISPERSION ANALYSIS (%d RUNS) & MID-COURSE TCM BUDGET]\n", mc.num_runs);
    printf("  * GNC Execution 1-Sigma : Thrust Error = %.2f%% | Attitude Pointing = %.2f deg\n",
           mc.sigma_dv_pct, mc.sigma_point_deg);
    printf("  * Uncorrected Miss Dist : Mean = %.2f km | 3-Sigma (99.7%% CEP) = %.2f km | Max = %.2f km\n",
           mc.mean_miss_km, mc.cep_3sigma_km, mc.max_miss_km);
    printf("  * 3-Sigma Ellipsoid RTN : [ 3s_R = %.2f km, 3s_T = %.2f km, 3s_N = %.2f km ]\n",
           mc.ellipsoid_3sigma_rtn.x, mc.ellipsoid_3sigma_rtn.y, mc.ellipsoid_3sigma_rtn.z);
    printf("  * Mid-Course TCM Budget : Mean = %.2f m/s | 99.7%% (3-Sigma) = %.2f m/s | Reserve Fuel = %.3f kg\n\n",
           mc.tcm_dv_mean_ms, mc.tcm_dv_99_ms, mc.tcm_fuel_99_kg);

    printf("[6-DOF ATTITUDE DYNAMICS (ADCS), QUATERNION SLEW & REACTION WHEELS]\n");
    printf("  * Principal Inertia J   : diag(%.1f, %.1f, %.1f) kg*m^2 | RW Capacity = %.1f N*m*s\n",
           adcs.inertia_kgm2.x, adcs.inertia_kgm2.y, adcs.inertia_kgm2.z, adcs.rw_h_max_nms);
    printf("  * Pre-Burn Slew Agile   : Duration = %.1f s | Burn 1 Slew = %.1f deg | Burn 2 Slew = %.1f deg\n",
           adcs.slew_time_s, adcs.slew_angle_b1_deg, adcs.slew_angle_b2_deg);
    printf("  * Steady Pointing RMS   : %.4f deg (%.1f arcsec) | Env Torques (GG/Aero/SRP) = %.2f / %.2f / %.2f mN*m\n",
           adcs.steady_point_err_deg, adcs.steady_point_err_arcsec,
           adcs.gg_torque_mnm, adcs.aero_torque_mnm, adcs.srp_torque_mnm);
    printf("  * Reaction Wheel Load   : Peak Torque = %.1f mN*m | Max H_rw = %.2f N*m*s (%.1f%% Sat) | MTQ Desats = %d\n\n",
           adcs.peak_rw_torque_mnm, adcs.max_rw_mom_nms, adcs.max_rw_sat_pct, adcs.desat_cycles);

    printf("[SPACECRAFT THERMAL CONTROL (TCS) & IONIZING RADIATION ENVIRONMENT (VAN ALLEN / SAA)]\n");
    printf("  * Surface Coating       : %s (alpha_s = %.2f, eps_IR = %.2f, alpha/eps = %.2f)\n",
           tcs.coating.name, tcs.alpha_solar, tcs.eps_ir, tcs.alpha_solar / tcs.eps_ir);
    printf("  * 2-Node Temp Excursion : Bus = [ %+.2f, %+.2f ] C | Ext Array = [ %+.2f, %+.2f ] C (Li-Ion Op: [%.0f, +%.0f] C)\n",
           tcs.min_temp_c, tcs.max_temp_c, tcs.ext_min_temp_c, tcs.ext_max_temp_c, tcs.t_min_op_c, tcs.t_max_op_c);
    printf("  * Mean Heat Fluxes      : Solar = %.1f W | Albedo = %.1f W | Earth IR = %.1f W | Radiated = %.1f W\n",
           tcs.mean_q_sol_w, tcs.mean_q_alb_w, tcs.mean_q_ir_w, tcs.mean_q_rad_w);
    printf("  * Survival Heaters      : Peak = %.1f W | Duty Cycle = %.1f%% | Energy = %.2f Wh\n",
           tcs.peak_heater_w, tcs.heater_duty_pct, tcs.heater_energy_wh);
    printf("  * Radiation & SAA (%.1fmm): Annual TID = %.2f krad/yr | Peak = %.1f urad/h | SEU = %.2f /day | Life = %.1f yrs\n",
           tcs.al_shield_mm, tcs.tid_krad_yr, tcs.peak_dose_urad_h, tcs.seu_rate_day, tcs.avionics_life_yr);
    printf("  * Magnetic Belt Status  : %s\n\n", tcs.saa_status);

    printf("[RF TELECOMMUNICATION LINK BUDGET, Eb/N0 MARGIN & GROUND STATION DOWNLINK]\n");
    printf("  * RF Band & Modulation  : %s | %s\n", telecom.band_name, telecom.modulation_scheme);
    printf("  * Onboard Transmitter   : P_tx = %.1f W (%.2f dBW) | G_tx = %.1f dBi -> EIRP = %.2f dBW\n",
           telecom.tx_power_w, telecom.tx_power_dbw, telecom.tx_gain_dbi, telecom.eirp_dbw);
    printf("  * Ground Station (%.1fm) : %s | G_rx = %.1f dBi | T_sys = %.0f K (N0 = %.2f dBW/Hz)\n",
           telecom.rx_dish_m, milano.name, telecom.rx_gain_dbi, telecom.sys_noise_temp_k, telecom.n0_dbw_hz);
    printf("  * Path Loss & Link Eb/N0: Mean FSPL = %.2f dB | Peak Eb/N0 = %.2f dB | Link Margin = [ %+.2f, %+.2f ] dB\n",
           telecom.mean_fspl_db, telecom.peak_ebn0_db, telecom.min_margin_db, telecom.peak_margin_db);
    printf("  * ACM Data Downlink     : Peak Rate = %.1f Mbps | Mission Vol = %.2f GB | 24h Volume = %.2f GB/day\n\n",
           telecom.max_achieved_mbps, telecom.mission_data_gb, telecom.total_data_gb_24h);

    printf("[PAYLOAD SENSOR SWATH, WALKER DELTA CONSTELLATION & MISSION MASS/COST BUDGET]\n");
    printf("  * Sensor Optics & Swath : Half-FOV = %.1f deg | Swath Width = %.1f km | GSD (Nadir) = %.2f m | Rate = %.0f km^2/s\n",
           sys.spec.half_fov_deg, sys.swath_width_km, sys.gsd_nadir_m, sys.area_coverage_rate_km2_s);
    printf("  * Walker Delta (i:T/P/F): %.1f deg : %d / %d / %d (%d sats/plane) | Mean Revisit = %.1f min | Global Cov = %.1f%%\n",
           sys.inc_deg, sys.spec.walker_total_sats, sys.spec.walker_planes, sys.spec.walker_phasing,
           sys.sats_per_plane, sys.mean_revisit_min, sys.global_coverage_pct);
    printf("  * Mass Budget (SMAD)    : Wet = %.1f kg | Dry = %.1f kg (%.1f%%) | Propellant = %.2f kg (Trans %.2f + Res %.2f kg)\n",
           sys.wet_mass_kg, sys.dry_mass_kg, sys.dry_wet_ratio_pct,
           sys.total_propellant_kg, sys.propellant_transfer_kg, sys.propellant_margin_kg);
    printf("  * Lifecycle Cost (NICM) : S/C = %.2f M EUR | Launch = %.2f M EUR | Ops (5yr) = %.2f M EUR -> Total = %.2f M EUR\n\n",
           sys.cost_spacecraft_meur, sys.cost_launch_meur, sys.cost_ops_5yr_meur, sys.total_mission_cost_meur);

    printf("[AUTONOMOUS ONBOARD NAVIGATION & EXTENDED KALMAN FILTER (6-STATE EKF)]\n");
    printf("  * GNSS Sensor & Prior   : GPS Noise 1s = %.1f m / %.2f m/s | Prior Uncertainty = %.1f m / %.2f m/s\n",
           nav.spec.gps_sigma_pos_m, nav.spec.gps_sigma_vel_ms, nav.spec.init_sigma_pos_m, nav.spec.init_sigma_vel_ms);
    printf("  * Steady-State EKF RMSE : Position = %.2f m (Raw GNSS = %.2f m) | Velocity = %.2f cm/s (%.4f m/s)\n",
           nav.steady_pos_rmse_m, nav.raw_pos_rmse_m, nav.steady_vel_rmse_cms, nav.steady_vel_rmse_ms);
    printf("  * Filter Performance    : Noise Reduction = %.1f%% | Convergence = %.1f min | 3-Sigma Bound = %.1f%% (Final 3s = %.2f m)\n\n",
           nav.noise_reduction_pct, nav.convergence_time_min, nav.consistency_pct, nav.final_cov_3sigma_m);

    printf("[SPACE SITUATIONAL AWARENESS (SSA), B-PLANE CONJUNCTION & CAM]\n");
    printf("  * Cataloged Secondary   : %s | v_rel = %.2f km/s | R_HBR = %.1f m\n",
           cam.event.debris_id, cam.event.rel_vel_kms, cam.event.hbr_m);
    printf("  * Pre-CAM B-Plane State : Miss = %.1f m (xi = %+.1f m, zeta = %+.1f m) | P_c = %.2e (log10 = %.2f)\n",
           cam.pre_miss_m, cam.pre_xi_m, cam.pre_zeta_m, cam.pre_pc, cam.pre_log10_pc);
    printf("  * ESA/NASA Alert Status : %s (Mahalanobis = %.2f sigma)\n",
           cam.alert_status, cam.pre_mahalanobis_d);
    printf("  * Optimal Phasing CAM   : dV_CAM = %.2f cm/s (%.4f m/s) at -%.1f orbits (-%.1f min) | Fuel = %.2f g\n",
           cam.req_dv_cam_cms, cam.req_dv_cam_ms, cam.event.lead_orbits, cam.tca_lead_time_min, cam.cam_fuel_g);
    printf("  * Post-CAM Safe State   : Separation = %.1f m (%.1f sigma) | Post-CAM P_c = %.2e (< 1e-7 SAFE)\n\n",
           cam.post_miss_m, cam.post_mahalanobis_d, cam.post_pc);

    printf("[END-OF-LIFE (EoL) DEORBIT, ATMOSPHERIC RE-ENTRY & SUTTON-GRAVES HEATING]\n");
    printf("  * Deorbit Burn Budget   : dV_deorbit = %.2f m/s (to h_p = %.0f km) | Fuel = %.2f kg | Passive Decay = %.1f yrs\n",
           eol.deorbit_dv_ms, eol.deorbit_perigee_km, eol.deorbit_fuel_kg, eol.passive_decay_years);
    printf("  * Sustainability Status : %s\n", eol.compliance_status);
    printf("  * Entry Interface (EI)  : h_EI = %.0f km | gamma_EI = %.2f deg | R_n = %.2f m | beta = %.1f kg/m^2\n",
           eol.entry_alt_km, eol.entry_fpa_deg, eol.nose_radius_m, eol.ballistic_coeff_kgm2);
    printf("  * Aerothermal & G-Loads : Peak Sutton-Graves = %.1f W/cm^2 (at %.1f km) | Peak Decel = %.2f g (at %.1f km)\n",
           eol.peak_heat_flux_wcm2, eol.peak_heat_alt_km, eol.peak_decel_g, eol.peak_decel_alt_km);
    printf("  * Structural Demise     : Breakup Altitude = %.1f km (V = %.2f km/s) | Entry Time = %.0f s | Downrange = %.0f km\n\n",
           eol.breakup_alt_km, eol.breakup_vel_kms, eol.total_entry_time_sec, eol.total_downrange_km);

    printf("[CISLUNAR CR3BP ASTRODYNAMICS, LAGRANGE POINTS (L1-L5) & LUNAR HALO ORBIT]\n");
    printf("  * Cislunar Regime       : %s (mu_EM = %.8f, L* = %.0f km)\n",
           cislunar.regime_name, cislunar.mu_em, cislunar.l_star_km);
    printf("  * Collinear L-Points    : L1 = %+.1f km (CJ=%.5f) | L2 = %+.1f km (CJ=%.5f) | L3 = %+.1f km\n",
           cislunar.l1_x_km, cislunar.cj_l1, cislunar.l2_x_km, cislunar.cj_l2, cislunar.l3_x_km);
    printf("  * Triangular L-Points   : L4/L5 = [ %+.1f, +/-%.1f ] km (CJ_min = %.5f)\n",
           cislunar.l4_x_km, cislunar.l4_y_km, cislunar.cj_l4_l5);
    printf("  * Propulsive Budget     : dV_TLI = %.2f m/s | dV_LOI (Halo) = %.2f m/s | Total = %.2f m/s | TOF = %.2f days\n",
           cislunar.dv_tli_ms, cislunar.dv_loi_ms, cislunar.total_cislunar_dv_ms, cislunar.tof_days);
    printf("  * Halo Orbit & Jacobi   : A_z = %.0f km (Period = %.2f d) | Jacobi C_J = %.6f | |dC_J/C_J0| = %.2e (< 1e-8 OK)\n\n",
           cislunar.halo_az_km, cislunar.halo_period_days, cislunar.jacobi_c, cislunar.delta_cj_rel);

    if (export_web_mission_json("../mission_report.json", rep_init, plan, rep_target, rdv, mc, adcs, tcs, telecom, sys, nav, cam, eol, cislunar, &milano, &sc)) {
        printf("[SUCCESS] 'mission_report.json' updated with Cislunar CR3BP, Lagrange Points & All Subsystems!\n");
    }
    export_web_mission_json("mission_report.json", rep_init, plan, rep_target, rdv, mc, adcs, tcs, telecom, sys, nav, cam, eol, cislunar, &milano, &sc);

    printf("======================================================================\n");
    return 0;
}
