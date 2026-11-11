#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "mission_analysis.h"
#include "telemetry.h"

int main(void) {
    printf("======================================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 16: TCS, RAD & RF TELECOM LINK BUDGET  \n");
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

    if (export_web_mission_json("../mission_report.json", rep_init, plan, rep_target, rdv, mc, adcs, tcs, telecom, &milano, &sc)) {
        printf("[SUCCESS] 'mission_report.json' updated with RF Telecom Link Budget, Thermal TCS, 6-DOF ADCS & GNC Data!\n");
    }
    export_web_mission_json("mission_report.json", rep_init, plan, rep_target, rdv, mc, adcs, tcs, telecom, &milano, &sc);

    printf("======================================================================\n");
    return 0;
}
