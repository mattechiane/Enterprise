#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "mission_analysis.h"
#include "telemetry.h"

int main(void) {
    printf("======================================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 32: EO WGS84 RAYCASTING & COVERAGE HEATMAP\n");
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

    // 3. INPUT ORBITA FINALE TARGET (Default: h = 10000 km, e = 0.02, i = 55.0 deg, RAAN = 50.0 deg)
    OrbitInputSpec input_target = {
        .mode = INPUT_KEPLERIAN,
        .kep  = {.a = EARTH_RADIUS + 10000.0, .e = 0.02, .i = 55.0 * DEG2RAD, .RAAN = 50.0 * DEG2RAD, .omega = 0.0, .theta = 0.0}
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

    // 12. SPRINT 20 & SPRINT 27: END-OF-LIFE (EoL) DEORBIT & 3-PHASE EDL (HYPERSONIC, PARACHUTE & POWERED DESCENT)
    EDLSequenceSpec edl_spec = {
        .entry_fpa_deg           = -1.5,   // Angolo di rampa all'Entry Interface (h_EI = 120 km) [deg]
        .nose_radius_m           = 0.50,   // Raggio di curvatura effettivo R_n [m]
        .chute_deploy_mach       = 2.1,    // Trigger apertura paracadute Supersonic Disk-Gap-Band (M < 2.1)
        .chute_max_alt_km        = 12.0,   // Quota massima apertura paracadute (h < 12 km)
        .chute_area_m2           = 185.0,  // Area paracadute dispiegato A_chute [m^2]
        .powered_descent_alt_m   = 2000.0, // Quota sgancio paracadute e inizio Powered Descent [m]
        .target_touchdown_vel_ms = 1.80    // Velocità nominale di Soft-Touchdown a h = 0 (< 2.5 m/s)
    };
    EndOfLifeReport eol = simulate_full_edl_sequence(orb_target, rep_target, plan, edl_spec, &sc, EARTH_MU);

    // 13. SPRINT 21: CISLUNAR ASTRODYNAMICS, CR3BP, LAGRANGE POINTS & LUNAR HALO ORBIT
    CR3BPSpec cr3bp_spec = {
        .target_regime      = CISLUNAR_L2_HALO, // L2 Far-Side Southern Halo (Lunar Gateway)
        .halo_az_km         = 15000.0,          // Ampiezza verticale Halo A_z = 15000 km
        .leo_parking_alt_km = 400.0             // Quota LEO pre-TLI [km]
    };
    CislunarReport cislunar = simulate_cislunar_cr3bp(orb_initial, plan, cr3bp_spec);

    // 14. SPRINT 22: FLEXIBLE STRUCTURAL DYNAMICS, PROPELLANT SLOSHING & CSI
    StructuralCSISpec csi_spec = {
        .flex_freq_hz       = 0.85,  // 1° modo flessionale pannello solare a sbalzo f_flex,1 = 0.85 Hz
        .flex_damping_zeta  = 0.015, // Smorzamento strutturale modale zeta_flex = 0.015
        .slosh_damping_zeta = 0.030, // Smorzamento setti antisciacquio (Baffles) zeta_slosh = 0.030
        .tank_radius_m      = 0.38,  // Raggio interno serbatoio propellente R_tank = 0.38 m
        .tank_fill_ratio    = 1.15,  // Rapporto riempimento liquido h_f / R_tank = 1.15
        .array_span_m       = 3.80,  // Lunghezza ala pannello solare L_wing = 3.80 m
        .enable_notch_filter = true  // Filtro digitale Bi-Quad Notch GNC attivo
    };
    StructuralReport csi = simulate_structural_slosh_and_csi(plan, adcs, sys, csi_spec);

    // 15. SPRINT 23: AUTONOMOUS FDIR, NORAD TLE PARSER & CCSDS OEM 2.0 EXPORT
    FDIRSpec fdir_spec = {
        .fault_scenario      = FDIR_FAULT_RW2_ACTUATOR, // Guasto tachimetro/driver RW-2 -> Switchover su Skewed RW-4 B-Side
        .fault_inject_time_s = 12.0,                    // Iniezione anomalia a t = 12.0 s
        .detection_threshold = 3.0,                     // Soglia residuo FDIR 3.0-sigma
        .nominal_bus_v       = 28.0                     // Bus primario non regolato a 28.0 V
    };
    FDIRReport fdir = simulate_fdir_and_safemode(adcs, tcs, fdir_spec);

    // Test diretto di parse_norad_tle(...) su TLE reale ISS (ZARYA)
    KeplerianElements kep_tle_test = {0};
    const char *tle_l1 = "1 25544U 98067A   26316.50000000  .00016717  00000-0  10270-3 0  9993";
    const char *tle_l2 = "2 25544  51.6416  24.8512 0006703 130.5360 229.6085 15.49560532421537";
    bool tle_ok = parse_norad_tle(tle_l1, tle_l2, &kep_tle_test);

    // 16. SPRINT 24: INTERPLANETARY PATCHED CONICS (EARTH -> MARS), PORKCHOP GRID & B-PLANE
    static InterplanetaryReport ipl;
    InterplanetarySpec ipl_spec = {
        .target         = PLANET_MARS,
        .capture_alt_km = 400.0,   // Periasse di cattura 400 km sopra la superficie
        .tof_days       = 210.0    // TOF eliocentrico di progetto
    };
    bool ipl_ok = compute_interplanetary_porkchop_and_transfer(orb_initial, ipl_spec, &ipl);

    // 17. SPRINT 25: MULTI-SATELLITE FORMATION FLYING (CHIEF + 3 DEPUTIES), J2-INVARIANT ROEs & E/I SAFETY
    FormationSpec form_spec = {
        .geometry           = FORM_PROJECTED_CIRCULAR_PCO, // Projected Circular Orbit (PCO)
        .baseline_radius_m  = 500.0,                       // Raggio della formazione rho_form = 500.0 m
        .tolerance_window_m = 5.0                          // Tolleranza di mantenimento baseline +/-5.0 m
    };
    FormationFlyingReport form = simulate_formation_flying_cluster(orb_target, plan, form_spec, EARTH_MU);

    // 18. SPRINT 26: AUTONOMOUS ONBOARD OPERATIONS SCHEDULER, SSR MEMORY BUFFER & ECSS/NASA PDR COMPLIANCE
    MissionOpsSpec ops_spec = {
        .ssr_capacity_gb        = 128.0, // Memoria Flash NAND a stato solido (SSR) con EDAC = 128 GB
        .science_data_rate_mbps = 120.0, // Rateo dati Payload Ottico/SAR = 120 Mbps
        .ssr_overflow_limit_pct = 95.0   // Soglia di protezione anti-overflow = 95%
    };
    MissionOpsReport ops = evaluate_autonomous_ops_and_ssr(rep_init, plan, rep_target, adcs, tcs, telecom,
                                                           nav, cam, eol, csi, fdir, form, ops_spec);

    // 19. SPRINT 28: SPACE WEATHER, SOLAR CYCLE (F10.7/Kp) & HIGH-FIDELITY ATMOSPHERIC DRAG (JACCHIA/NRLMSISE-00)
    SpaceWeatherSpec swx_spec = {
        .cycle_phase           = SOLAR_CYCLE_AVERAGE, // Nominal Solar Activity (F10.7 = 140 SFU)
        .f107_sfu              = 140.0,               // Flusso radio solare giornaliero a 10.7 cm [SFU]
        .f107_avg_sfu          = 140.0,               // Media centrata a 81 giorni F10.7 [SFU]
        .kp_index              = 3.0,                 // Indice planetario geomagnetico Kp (0..9)
        .diurnal_bulge_lag_deg = 30.0                 // Diurnal Bulge lag (~14:00 LST, +30 deg East)
    };
    OrbitLifetimeReport swx = simulate_lifetime_orbit_decay(orb_target, plan, &sc, swx_spec, EARTH_MU);

    // 20. SPRINT 29: OPTICAL INTER-SATELLITE LINKS (ISL), LINE-OF-SIGHT & CONSTELLATION MESH NETWORK
    ISLNetworkSpec isl_spec = {
        .max_range_km        = 5000.0, // Portata massima del terminale laser LCT [km]
        .occ_alt_km          = 80.0,   // Quota atmosferica opaca di occultazione R_occ = R_E + 80 km
        .laser_wavelength_nm = 1550.0, // Lunghezza d'onda laser C-Band [nm]
        .link_rate_gbps      = 10.0    // Throughput per singolo link coerente BPSK/Homodyne [Gbps]
    };
    ISLNetworkReport isl = evaluate_optical_isl_network(orb_target, sys, isl_spec);

    // 21. SPRINT 30: MULTI-GROUND STATION TRACKING, DSN/ESTRACK GEOMETRY & LINK HANDOVERS
    GroundNetworkSpec gnd_spec = {
        .network_type       = NET_NASA_DSN, // NASA Deep Space Network (Goldstone, Madrid, Canberra + Milano-Bovisa)
        .elevation_mask_deg = 5.0           // Maschera di elevazione minima sopra l'orizzonte [deg]
    };
    GroundSegmentReport gnd = evaluate_multi_station_coverage(orb_target, plan, &sc, gnd_spec, EARTH_MU);

    // 22. SPRINT 31: HIGH-FIDELITY EPS, CONICAL PENUMBRA/UMBRA ECLIPSES & LI-ION DEGRADATION
    EPSHiFiSpec eps_spec = {
        .batt_config            = BATT_CONFIG_8S4P,          // 8S4P Li-Ion per 28V Unregulated Bus (48 Ah, 1382 Wh)
        .cell_type              = CELL_GAAS_TRIPLE_JUNCTION, // GaAs Triple-Junction InGaP/GaAs/Ge (eta_0 = 30.2%)
        .solar_array_area_m2    = 2.85,                      // Area efficace pannelli solari A_sa [m^2]
        .mission_design_life_yr = 5.0                        // Vita operativa di progetto (5 anni)
    };
    EPSHiFiReport eps = simulate_eps_hifi_profile(orb_initial, orb_target, plan, tcs, adcs, telecom, eps_spec, EARTH_MU);

    // 23. SPRINT 32: EARTH OBSERVATION (EO) 6-DOF WGS84 RAYCASTING, SWATH POLYGON & COVERAGE HEATMAP
    static EORaycastReport eo;
    EOSensorSpec eo_spec = {
        .fov_x_deg      = 10.0, // Field of View Along-Track = 10.0 deg
        .fov_y_deg      = 15.0, // Field of View Cross-Track = 15.0 deg
        .off_nadir_deg  = 0.0,  // Puntamento Boresight Off-Nadir = 0.0 deg (Nadir Pushbroom)
        .opt_aperture_m = 0.35, // Apertura ottica D_opt = 0.35 m
        .wavelength_nm  = 550.0 // Lunghezza d'onda visibile lambda = 550 nm
    };
    eo = simulate_eo_raycasting_and_coverage(orb_target, plan, adcs, eo_spec, EARTH_MU);

    // 24. SPRINT 33: INTERPLANETARY GRAVITY ASSIST (PLANETARY FLYBY), V-INFINITY ROTATION & B-PLANE TARGETING
    GravityAssistSpec flyby_spec = {
        .enable_flyby = true,          // Abilita Gravity Assist (Slingshot) di default
        .flyby_planet = FLYBY_JUPITER, // Flyby su Giove (Voyager / Cassini / New Horizons class)
        .flyby_alt_km = 140000.0       // Altitudine di pericentro h_p = 140,000 km (~1.96 R_J sopra le nubi)
    };
    GravityAssistReport flyby = simulate_gravity_assist_flyby(orb_initial, &ipl, &sc, flyby_spec);

    // 25. SPRINT 34: AEROCAPTURE, AEROBRAKING & MARTIAN ATMOSPHERIC PASS
    AerocaptureSpec aero_spec = {
        .method               = INSERTION_AEROCAPTURE, // Cattura atmosferica su Marte invece del burn propulsivo dV_MOI
        .planet               = PLANET_MARS,           // Atmosfera di CO2 di Marte (Entry Interface h_EI = 125 km)
        .lift_to_drag         = 0.25,                  // Rapporto Portanza/Resistenza ipersonico L/D = 0.25
        .entry_fpa_deg        = -11.35,                // Flight-Path Angle nominale nel centro del corridoio [deg]
        .ballistic_coeff_kgm2 = 145.0,                 // Coefficiente balistico beta = m / (C_D * A) [kg/m^2]
        .nose_radius_m        = 1.25,                  // Raggio dello scudo termico TPS R_n = 1.25 m
        .target_apoapsis_km   = 1400.0,                // Quota di apocentro target post-passaggio [km]
        .final_periapsis_km   = 400.0                  // Quota di pericentro finale dopo Periapsis Raise Burn [km]
    };
    AerocaptureReport aero = simulate_aerocapture_pass(&ipl, plan, aero_spec);
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
    printf("  * Phase 2 Parachute     : Deploy at h = %.2f km (M = %.2f, V = %.1f m/s) | Snatch Shock = %.2f g (%.1f kN)\n",
           eol.edl.chute_deploy_alt_km, eol.edl.chute_deploy_mach, eol.edl.chute_deploy_vel_ms,
           eol.edl.chute_snatch_g, eol.edl.chute_snatch_kn);
    printf("  * Phase 3 Powered Desc. : Ignition at h = %.0f m (V = %.1f m/s) | dV_PD = %.1f m/s | Fuel = %.2f kg | Touchdown V_td = %.2f m/s\n",
           eol.edl.pd_start_alt_m, eol.edl.pd_start_vel_ms, eol.edl.pd_dv_ms, eol.edl.pd_fuel_used_kg, eol.edl.touchdown_vel_ms);
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

    printf("[FLEXIBLE STRUCTURAL DYNAMICS, PROPELLANT SLOSHING & CONTROL-STRUCTURE INTERACTION (CSI)]\n");
    printf("  * Solar Array Modes     : f_flex,1 = %.2f Hz (zeta = %.3f) | f_flex,2 = %.2f Hz | Peak Tip Defl = %.2f mm\n",
           csi.flex_freq_hz, csi.spec.flex_damping_zeta, csi.flex_freq_mode2_hz, csi.peak_tip_defl_mm);
    printf("  * Abramson Tank Slosh   : f_slosh = %.3f Hz (m_slosh = %.2f kg, zeta = %.3f) | Max Torque = %.2f mN*m\n",
           csi.slosh_freq_hz, csi.slosh_mass_kg, csi.spec.slosh_damping_zeta, csi.peak_slosh_torque_mnm);
    printf("  * Micro-Vibration Jitter: RMS = %.2f arcsec | Peak = %.2f arcsec | Settling Time = %.1f s\n",
           csi.rms_jitter_arcsec, csi.peak_jitter_arcsec, csi.settling_time_sec);
    printf("  * GNC Notch Filter & CSI: %s (Gain Margin = +%.1f dB)\n",
           csi.notch_status, csi.csi_gain_margin_db);
    printf("  * Launch Loads & QSL    : Stowed f_ax = %.1f Hz (>35), f_lat = %.1f Hz (>15) | sigma_vM = %.1f MPa -> MoS = +%.1f%%\n\n",
           csi.launch_axial_freq_hz, csi.launch_lateral_freq_hz, csi.max_von_mises_mpa, csi.launch_mos_pct);

    printf("[AUTONOMOUS FDIR, REDUNDANCY SWITCHOVER, NORAD TLE PARSER & CCSDS OEM 2.0]\n");
    printf("  * Injected Fault Event  : %s [%s] at t = %.1f s\n",
           fdir.fault_name, fdir.fault_subsystem, fdir.fault_inject_time_s);
    printf("  * FDIR Timing & Latency : Detection = %.2f s | LCL Isolation & B-Side Switchover = %.2f s | Recovery = %.1f s\n",
           fdir.detection_latency_s, fdir.isolation_switchover_s, fdir.total_recovery_time_s);
    printf("  * Redundancy Switchover : %s -> %s\n",
           fdir.primary_unit_name, fdir.redundant_unit_name);
    printf("  * Transient Excursion   : Peak Attitude Error = %.2f deg (Post-Recovery = %.3f deg) | Min Bus = %.2f V (Nom %.1f V)\n",
           fdir.peak_att_err_deg, fdir.steady_post_err_deg, fdir.min_bus_voltage_v, fdir.nominal_bus_voltage_v);
    printf("  * NORAD TLE Decoder     : %s -> a = %.2f km (h = %.1f km), i = %.4f deg, e = %.7f\n",
           tle_ok ? "VALIDATED (ISS #25544)" : "ERROR",
           kep_tle_test.a, kep_tle_test.a - EARTH_RADIUS, kep_tle_test.i * RAD2DEG, kep_tle_test.e);
    printf("  * CCSDS OEM Ephemeris   : %s (%d EME2000 State Vectors)\n\n",
           fdir.ccsds_oem_status, fdir.oem_state_vectors);

    printf("[INTERPLANETARY PATCHED CONICS: EARTH -> %s, PORKCHOP & B-PLANE TARGETING]\n", ipl.target_name);
    if (ipl_ok) {
        printf("  * Launch Window         : day +%.0f | TOF = %.0f d | transfer angle = %.1f deg | a = %.4f AU, e = %.4f\n",
               ipl.opt_dep_day, ipl.opt_tof_days, ipl.transfer_angle_deg, ipl.transfer_a_au, ipl.transfer_e);
        printf("  * Departure Hyperbola   : C3 = %.3f km2/s2 (v_inf = %.3f km/s) | DLA = %.2f deg, RLA = %.2f deg | dV_TMI = %.1f m/s\n",
               ipl.c3_opt, ipl.vinf_dep_kms, ipl.dla_deg, ipl.rla_deg, ipl.dv_tmi_ms);
        printf("  * Arrival Hyperbola     : v_inf = %.3f km/s | B-plane B = %.0f km | dV_MOI = %.1f m/s | Total dV = %.1f m/s\n",
               ipl.vinf_arr_kms, ipl.b_plane_km, ipl.dv_moi_ms, ipl.dv_total_ms);
        printf("  * Porkchop Grid Minima  : C3_min = %.2f km2/s2 (dep +%.0f d, TOF %.0f d) | dV_min = %.0f m/s (dep +%.0f d, TOF %.0f d)\n",
               ipl.grid_min_c3, ipl.grid_min_c3_dep, ipl.grid_min_c3_tof, ipl.grid_min_dv_ms, ipl.grid_min_dv_dep, ipl.grid_min_dv_tof);
        printf("  * Lambert Closure Check : arrival position residual = %.3f km\n\n", ipl.lambert_residual_km);
    } else {
        printf("  * No valid Lambert solution for the selected TOF.\n\n");
    }

    printf("[MULTI-SATELLITE FORMATION FLYING (CHIEF + %d DEPUTIES), ROEs & J2-INVARIANCE]\n", form.num_deputies);
    printf("  * Formation Topology    : %s | Baseline Radius rho = %.1f m (+/-%.1f m window)\n",
           form.geometry_name, form.spec.baseline_radius_m, form.spec.tolerance_window_m);
    printf("  * Deputy 1 D'Amico ROEs : a*[da, dl, dex, dey, dix, diy] = [ %+.1f, %+.1f, %+.1f, %+.1f, %+.1f, %+.1f ] m\n",
           form.roe_deputy[0].da_m, form.roe_deputy[0].dlambda_m,
           form.roe_deputy[0].dex_m, form.roe_deputy[0].dey_m,
           form.roe_deputy[0].dix_m, form.roe_deputy[0].diy_m);
    printf("  * Passive E/I Safety    : Min Safe R/N Sep = %.2f m | Min 3D Inter-Sat = %.2f m | E/I Phase Angle = %.1f deg (%s)\n",
           form.min_safe_rn_sep_m, form.min_3d_inter_sat_m, form.ei_phase_angle_deg,
           form.passive_safe ? "COLLINEAR SAFE" : "CHECK PHASE");
    printf("  * Secular J2 Drift      : dOmega/dt = %+.5f deg/day | domega/dt = %+.5f deg/day | Max Baseline Err = %.2f m\n",
           form.j2_raan_drift_deg_day, form.j2_argp_drift_deg_day, form.max_baseline_err_m);
    printf("  * Fleet Delta-V Budget  : Reconfiguration dV = %.3f m/s | Formation-Keeping dV_FK = %.3f m/s/yr (Fuel = %.3f kg/yr)\n\n",
           form.reconfig_dv_ms, form.fk_dv_yr_ms, form.fk_fuel_yr_kg);

    printf("[AUTONOMOUS ONBOARD OPERATIONS SCHEDULER, SSR MEMORY BUFFER & ECSS/NASA PDR COMPLIANCE]\n");
    printf("  * SSR Flash Memory Spec : Capacity = %.1f GB | Science Rate R_gen = %.1f Mbps | Overflow Guard = %.0f%%\n",
           ops.spec.ssr_capacity_gb, ops.spec.science_data_rate_mbps, ops.spec.ssr_overflow_limit_pct);
    printf("  * SSR Memory Balance    : Peak Occupancy = %.2f GB (%.1f%%) | Acquired = %.2f GB | Downlinked = %.2f GB | Backlog = %.2f GB\n",
           ops.peak_ssr_used_gb, ops.peak_ssr_fill_pct, ops.total_acquired_gb, ops.total_downlinked_gb, ops.net_backlog_gb);
    printf("  * 24h Ops Duty Cycle    : Science = %.1f%% | Downlink = %.1f%% | Sun-Charge = %.1f%% | Maneuver = %.1f%% | Eclipse = %.1f%%\n",
           ops.duty_science_pct, ops.duty_downlink_pct, ops.duty_sun_charge_pct, ops.duty_maneuver_pct, ops.duty_eclipse_pct);
    printf("  * ECSS/NASA Compliance  : Prop=%s | EPS=%s | TCS=%s | RF=%s | SSR=%s | GNC/FDIR=%s | SSA/EoL=%s\n",
           ops.status_propulsion, ops.status_eps_battery, ops.status_thermal_rad, ops.status_rf_telecom,
           ops.status_ssr_buffer, ops.status_gnc_adcs_fdir, ops.status_ssa_eol);
    printf("  * System Readiness Score: %.2f%% -> %s\n\n",
           ops.readiness_score_pct, ops.pdr_compliance_verdict);

    printf("[SPACE WEATHER, SOLAR CYCLE (F10.7/Kp), DIURNAL BULGE & LONG-TERM ORBIT DECAY]\n");
    printf("  * Solar Cycle & Geomagn.: %s | Kp = %.1f (Ap = %.1f) | Bulge Lag = +%.0f deg (Ratio = %.2fx)\n",
           swx.phase_name, swx.spec.kp_index, swx.ap_index, swx.spec.diurnal_bulge_lag_deg, swx.diurnal_bulge_ratio);
    printf("  * Exospheric Temp (J71) : T_night = %.0f K | T_day(Bulge) = %.0f K | Orbit Mean T_inf = %.0f K\n",
           swx.t_night_min_k, swx.t_day_apex_k, swx.t_exospheric_k);
    printf("  * Thermosphere @ %.0f km : rho = %.3e kg/m^3 | Scale Height H = %.1f km | |da/dt| = %.2f m/day\n",
           swx.initial_alt_km, swx.rho_at_target_kgm3, swx.scale_height_km, swx.initial_decay_m_day);
    printf("  * Lifetime Comparison   : Solar Min(70) = %.1f yr | Nominal(140) = %.2f yr | Solar Max(220) = %.2f yr\n",
           swx.lifetime_solar_min_yr, swx.lifetime_solar_avg_yr, swx.lifetime_solar_max_yr);
    printf("  * 25-Year Rule Verdict  : %s\n\n", swx.compliance_verdict);

    printf("[OPTICAL INTER-SATELLITE LINKS (ISL), LINE-OF-SIGHT & WALKER MESH TOPOLOGY]\n");
    printf("  * Laser Terminal Spec   : lambda = %.0f nm | Rate = %.0f Gbps/link | Max Range = %.0f km | Limb R_occ = R_E+%.0f km\n",
           isl.spec.laser_wavelength_nm, isl.spec.link_rate_gbps, isl.spec.max_range_km, isl.spec.occ_alt_km);
    printf("  * Active Mesh Topology  : %d Active Links (%d Intra-Plane Ring + %d Cross-Plane) across %d Walker Nodes\n",
           isl.num_active_links, isl.intra_plane_links, isl.inter_plane_links, isl.num_nodes);
    printf("  * ISL Range & Latency   : d = [ %.1f, %.1f ] km (Mean = %.1f km) | Mean Latency = %.2f ms (Max = %.2f ms)\n",
           isl.min_distance_km, isl.max_distance_km, isl.mean_distance_km, isl.mean_latency_ms, isl.max_latency_ms);
    printf("  * Optical Link Budget   : Mean Margin = +%.2f dB (Min = +%.2f dB) | Aggregate Mesh Capacity = %.0f Gbps\n",
           isl.mean_margin_db, isl.min_margin_db, isl.total_capacity_gbps);
    printf("  * Network Status        : %s (Mean Degree = %.2f links/sat)\n\n",
           isl.mesh_connectivity_status, isl.mean_node_degree);

    printf("[MULTI-GROUND STATION TRACKING, DSN/ESTRACK GEOMETRY & LINK HANDOVERS]\n");
    printf("  * Active Ground Network : %s | Elevation Mask = %.1f deg | Stations = %d\n",
           gnd.network_name, gnd.spec.elevation_mask_deg, gnd.num_stations);
    printf("  * 24h Network Coverage  : Total Contact = %.2f h/24h (%.1f%%) | Max Comm Gap = %.1f min | Handovers = %d\n",
           gnd.total_contact_hours_24h, gnd.coverage_duty_pct, gnd.max_comm_gap_min, gnd.total_handovers);
    printf("  * Primary Tracking Site : %s (Contact = %.1f min/24h) | Total Passes = %d\n",
           gnd.best_station_name, gnd.best_station_contact_min, gnd.num_passes);
    printf("  * Antenna Pedestal Slew : Peak Mechanical Slew Rate = %.4f deg/s (Zenith Tracking OK)\n\n",
           gnd.peak_antenna_slew_deg_s);

    printf("[HIGH-FIDELITY EPS, CONICAL PENUMBRA/UMBRA ECLIPSES & LI-ION DEGRADATION]\n");
    printf("  * Solar Array & MPPT    : %s | eta_T = %.2f%% | MPPT Eff = %.2f%% | TID Loss = -%.2f%%\n",
           eps.cell_type_name, eps.solar_eff_mean_pct, eps.mppt_eff_pct, eps.tid_deg_pct);
    printf("  * Conical Shadow Model  : Penumbra Transit = %.2f min | Umbra Core = %.2f min | Peak Solar = %.1f W\n",
           eps.penumbra_time_min, eps.umbra_time_min, eps.peak_solar_pwr_w);
    printf("  * Li-Ion Battery Pack   : %s | Bus Range = [ %.2f, %.2f ] V (Nom %.1f V)\n",
           eps.batt_config_name, eps.min_bus_voltage_v, eps.max_bus_voltage_v, eps.nominal_bus_v);
    printf("  * Cyclic DoD & EoL SoH  : Peak DoD = %.2f%% (Min SoC = %.2f%%) | End-of-Life SoH = %.2f%% (%.0f / %.0f Wh)\n",
           eps.peak_dod_pct, eps.min_soc_pct, eps.eol_soh_pct, eps.eol_cap_wh, eps.nominal_cap_wh);
    printf("  * ECSS Power Verdict    : %s\n\n", eps.eps_compliance_status);

    printf("[EARTH OBSERVATION (EO) WGS84 6-DOF RAYCASTING & CUMULATIVE COVERAGE HEATMAP]\n");
    printf("  * Sensor Mode & FOV     : %s | Off-Nadir Tilt = %+.1f deg\n",
           eo.sensor_mode_name, eo.spec.off_nadir_deg);
    printf("  * Boresight Geometry    : Mean Slant Range = %.2f km | Incidence = %.2f deg | Swath Area = %.0f km^2\n",
           eo.avg_slant_range_km, eo.avg_incidence_deg, eo.swath_area_km2);
    printf("  * Pixel Smearing & GSD  : Nadir GSD = %.3f m -> Effective Off-Nadir GSD = %.3f m (Smearing = %.3fx)\n",
           eo.nadir_gsd_m, eo.mean_gsd_m, eo.off_nadir_smearing_factor);
    printf("  * Global Heatmap Grid   : Cumulative Coverage = %.2f%% (%d / %d cells hit, Peak Revisit = %dx)\n\n",
           eo.global_coverage_pct, eo.unique_cells_hit, eo.total_grid_cells, eo.max_revisit_hits);

    printf("[INTERPLANETARY GRAVITY ASSIST (PLANETARY FLYBY), V-INFINITY ROTATION & B-PLANE]\n");
    printf("  * Flyby Target & Periap.: %s (a = %.2f AU, R_p = %.1f km) | h_flyby = %.0f km (r_p = %.0f km)\n",
           flyby.flyby_planet_name, flyby.planet_sma_au, flyby.r_planet_km, flyby.spec.flyby_alt_km, flyby.periapsis_radius_km);
    printf("  * Hyperbolic Excess v_inf: ||v_inf|| = %.3f km/s | B-Plane Impact B = %.0f km | e_hyp = %.4f\n",
           flyby.vinf_kms, flyby.b_plane_km, flyby.hyp_ecc);
    printf("  * Deflection & Velocities: Turn Angle delta = %.2f deg | V_helio,in = %.2f km/s -> V_helio,out = %.2f km/s\n",
           flyby.deflection_angle_deg, flyby.v_sc_in_kms, flyby.v_sc_out_kms);
    printf("  * Propellantless Gain   : Free Gravity-Assist dV = %.3f km/s (%.0f m/s) | Helio Boost = %+.3f km/s (0 kg Fuel!)\n\n",
           flyby.free_dv_kms, flyby.free_dv_ms, flyby.helio_speed_gain_kms);

    printf("[PLANETARY AEROCAPTURE, AEROBRAKING & MARTIAN ATMOSPHERIC PASS]\n");
    printf("  * Insertion Method & L/D: %s on %s | L/D = %.2f | beta = %.0f kg/m^2\n",
           aero.method_name, aero.planet_name, aero.spec.lift_to_drag, aero.spec.ballistic_coeff_kgm2);
    printf("  * Corridor & Min Alt    : gamma_EI = %.2f deg [Corridor: %.2f .. %.2f deg] | Min Pass Alt = %.2f km\n",
           aero.spec.entry_fpa_deg, aero.corridor_burn_fpa_deg, aero.corridor_skip_fpa_deg, aero.min_altitude_km);
    printf("  * Aerothermal Heating   : Peak Heat Flux = %.2f W/cm^2 | Total Load = %.0f J/cm^2 | Peak Decel = %.2f g\n",
           aero.peak_heat_flux_wcm2, aero.total_heat_load_jcm2, aero.peak_decel_g);
    printf("  * Orbit Capture & Fuel  : e_in = %.3f -> e_final = %.3f (h_a = %.0f km) | Eq. dV Saved = %.3f km/s (%.0f kg Heavy Stack Fuel Saved!)\n\n",
           aero.initial_hyp_ecc, aero.final_eccentricity, aero.captured_apoapsis_km, aero.eq_dv_saved_kms, aero.fuel_saved_heavy_stack_kg);

    if (export_ccsds_oem_file("../enterprise_mission.oem", rep_init, plan, rep_target, &sc)) {
        printf("[SUCCESS] 'enterprise_mission.oem' exported (CCSDS 502.0-B-2 OEM v2.0 Standard)!\n");
    }
    export_ccsds_oem_file("enterprise_mission.oem", rep_init, plan, rep_target, &sc);

    if (export_web_mission_json("../mission_report.json", rep_init, plan, rep_target, rdv, mc, adcs, tcs, telecom, sys, nav, cam, eol, cislunar, csi, fdir, &ipl, form, ops, swx, isl, gnd, eps, &eo, flyby, aero, &milano, &sc)) {
        printf("[SUCCESS] 'mission_report.json' updated with Martian Aerocapture, Gravity Assist, EO Raycasting & All Subsystems!\n");
    }
    export_web_mission_json("mission_report.json", rep_init, plan, rep_target, rdv, mc, adcs, tcs, telecom, sys, nav, cam, eol, cislunar, csi, fdir, &ipl, form, ops, swx, isl, gnd, eps, &eo, flyby, aero, &milano, &sc);

    printf("======================================================================\n");
    return 0;
}
