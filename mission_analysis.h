#ifndef MISSION_ANALYSIS_H
#define MISSION_ANALYSIS_H

#include "orbit_types.h"

#define MIN_SAFE_PERIGEE_ALT_KM 150.0   // Limite minimo atmosferico per orbita chiusa stabile
#define MAX_HILL_APOGEE_ALT_KM 900000.0 // Limite massimo Sfera di Hill terrestre

// Risultato del Risolutore Universale di Lambert 3D
typedef struct {
    bool converged;          // true se l'algoritmo numerico è converso
    Vector3D v1_eci;         // Velocità richiesta all'inizio dell'arco di trasferimento (km/s)
    Vector3D v2_eci;         // Velocità d'arrivo alla fine dell'arco di trasferimento (km/s)
    double dv1_mag;          // Costo impulsivo Burn 1 (km/s)
    double dv2_mag;          // Costo impulsivo Burn 2 (km/s)
    double total_dv;         // Delta-V totale di Lambert (km/s)
    Vector3D burn1_dir_rtn;  // Versore direzione spinta 1 nel frame RTN
    Vector3D burn2_dir_rtn;  // Versore direzione spinta 2 nel frame RTN
} LambertSolution;

Vector3D geodetic_to_ecef(GeodeticCoord geo);
AccessMetrics compute_ground_access(Vector3D sat_ecef, const GroundStation *station);
HohmannTransfer compute_hohmann_transfer(double r1, double r2, double mu);
CartesianState apply_impulsive_maneuver(CartesianState state, Vector3D delta_v_eci, double isp);

PropellantSpec get_propellant_specs(PropellantType type);
KeplerianElements resolve_user_orbit_input(OrbitInputSpec input, double mu);
bool validate_and_clamp_closed_orbit(KeplerianElements *kep, char *warning_msg);

// Risolutore Universale di Lambert (Variabili Universali + Funzioni di Stumpff)
LambertSolution solve_lambert_universal(Vector3D r1, Vector3D v1_init,
                                        Vector3D r2, Vector3D v2_targ,
                                        double tof_sec, double mu, bool prograde);

// Ottimizzatore Globale 3D (include Hohmann Split-Plane, Bi-Ellittica e Risolutore di Lambert)
OptimalTransferPlan optimize_orbital_transfer(KeplerianElements orb_init,
                                              KeplerianElements orb_target,
                                              double init_mass,
                                              PropellantType prop_type,
                                              double max_transfer_time_sec,
                                              double mu);

OrbitReport evaluate_orbit_regime(KeplerianElements kep, double mass, const SpacecraftConfig *sc, const GroundStation *station);
void print_orbit_report(const char *title, OrbitReport rep, const char *station_name);

// Propagazione analitica chiusa tramite Matrice di Transizione di Stato (STM) di Clohessy-Wiltshire (LVLH)
void propagate_cw_state(Vector3D r0_m, Vector3D v0_ms, double n_rad_s, double t_sec,
                        Vector3D *r_out_m, Vector3D *v_out_ms);

// Risolutore Two-Impulse Clohessy-Wiltshire per Terminal Rendezvous & Soft-Docking a [0, 0, 0]
RendezvousPlan compute_cw_rendezvous(double a_target_km,
                                     Vector3D r0_lvlh_m,
                                     Vector3D v0_lvlh_ms,
                                     double t_rdv_sec,
                                     double init_mass_kg,
                                     double isp_s,
                                     double mu);

// Sprint 13: Analisi di Dispersione Monte Carlo (Box-Muller), Ellissoide 3-Sigma e Calcolo Mid-Course TCM
MonteCarloReport run_monte_carlo_dispersion(KeplerianElements orb_init,
                                            KeplerianElements orb_target,
                                            OptimalTransferPlan plan,
                                            double sigma_dv_pct,
                                            double sigma_point_deg,
                                            int num_runs,
                                            double mu);

// Sprint 14: Dinamica d'Assetto 6-DOF (Eulero + Quaternioni), Slew Control PD e Gestione Reaction Wheels
ADCSReport simulate_adcs_profile(KeplerianElements orb_init,
                                 KeplerianElements orb_target,
                                 OptimalTransferPlan plan,
                                 Vector3D inertia_kgm2,
                                 double rw_h_max_nms,
                                 double mu);

// Sprint 15: Controllo Termico (TCS) a 2 Nodi (Bus + Array) e Ambiente Radiativo (Van Allen TID & SAA)
SurfaceOpticalProps get_surface_optical_props(ThermalCoatingType type);

ThermalRadReport simulate_thermal_and_radiation(KeplerianElements orb_init,
                                                KeplerianElements orb_target,
                                                OptimalTransferPlan plan,
                                                SurfaceOpticalProps coating,
                                                double al_shield_mm,
                                                double t_min_op_c,
                                                double t_max_op_c,
                                                double mu);

ThermalReport simulate_thermal_profile(KeplerianElements orb_init,
                                       KeplerianElements orb_target,
                                       OptimalTransferPlan plan,
                                       double alpha_solar,
                                       double eps_ir,
                                       double t_min_op_c,
                                       double t_max_op_c,
                                       double mu);

// Sprint 16: RF Telecommunication Link Budget, Eb/N0 Margin & Ground Station Data Downlink Volume
TelecomReport evaluate_rf_link_budget(KeplerianElements orb_init,
                                      KeplerianElements orb_target,
                                      OptimalTransferPlan plan,
                                      OrbitReport rep_target,
                                      RFBandType band,
                                      double tx_power_w,
                                      const GroundStation *station,
                                      double mu);

// Sprint 17: Payload Sensor Swath, Walker Delta Constellation Coverage & Parametric Mission Cost/Mass Model
SystemEngineeringReport evaluate_payload_constellation_and_cost(KeplerianElements orb_target,
                                                                OrbitReport rep_target,
                                                                OptimalTransferPlan plan,
                                                                RendezvousPlan rdv,
                                                                MonteCarloReport mc,
                                                                PayloadAndConstellationSpec spec,
                                                                const GroundStation *station,
                                                                double mu);

// Sprint 18: Autonomous Onboard Navigation & Real-Time Orbit Determination via Extended Kalman Filter (EKF)
NavigationReport run_onboard_ekf_navigation(KeplerianElements orb_init,
                                            KeplerianElements orb_target,
                                            OptimalTransferPlan plan,
                                            NavigationEKFSpec spec,
                                            const SpacecraftConfig *sc,
                                            double mu);

// Sprint 19: Space Situational Awareness (SSA), B-Plane Conjunction Analysis & Collision Avoidance Maneuver (CAM)
CAMReport evaluate_debris_conjunction_and_cam(KeplerianElements orb_target,
                                              OrbitReport rep_target,
                                              OptimalTransferPlan plan,
                                              ConjunctionEvent event,
                                              double mu);

// Sprint 20 & Sprint 27: End-of-Life (EoL) Deorbit Disposal, 3-Phase EDL (Hypersonic Entry,
// Supersonic Parachute Deployment & Powered Descent Gravity-Turn/Skycrane) & Sutton-Graves Heating
EndOfLifeReport simulate_full_edl_sequence(KeplerianElements orb_target,
                                           OrbitReport rep_target,
                                           OptimalTransferPlan plan,
                                           EDLSequenceSpec edl_spec,
                                           const SpacecraftConfig *sc,
                                           double mu);

EndOfLifeReport simulate_eol_deorbit_and_reentry(KeplerianElements orb_target,
                                                 OrbitReport rep_target,
                                                 OptimalTransferPlan plan,
                                                 double entry_fpa_deg,
                                                 double nose_radius_m,
                                                 const SpacecraftConfig *sc,
                                                 double mu);

// Sprint 21: Cislunar Astrodynamics, Circular Restricted Three-Body Problem (CR3BP), Lagrange Points & Lunar Halo
CislunarReport simulate_cislunar_cr3bp(KeplerianElements orb_init,
                                       OptimalTransferPlan plan,
                                       CR3BPSpec spec);

// Sprint 22: Flexible Structural Dynamics, Propellant Sloshing Pendulum & Control-Structure Interaction (CSI)
StructuralReport simulate_structural_slosh_and_csi(OptimalTransferPlan plan,
                                                   ADCSReport adcs,
                                                   SystemEngineeringReport sys,
                                                   StructuralCSISpec spec);

// Sprint 23: Autonomous FDIR (Fault Detection, Isolation & Recovery), Redundancy Switchover & Safe-Mode
FDIRReport simulate_fdir_and_safemode(ADCSReport adcs,
                                      ThermalReport tcs,
                                      FDIRSpec spec);

// Sprint 24: Interplanetary Patched Conics (Earth -> Mars/Venus), Porkchop Grid (C3 / V_inf) & B-Plane Targeting
// Ritorna false se nessuna soluzione Lambert valida esiste per il TOF di progetto.
bool compute_interplanetary_porkchop_and_transfer(KeplerianElements orb_init,
                                                  InterplanetarySpec spec,
                                                  InterplanetaryReport *out);

// Sprint 25: Multi-Satellite Formation Flying (Chief + 3 Deputies), J2-Invariant ROEs (D'Amico / Schaub-Alfriend)
// & Passive E/I-Vector Separation Safety
FormationFlyingReport simulate_formation_flying_cluster(KeplerianElements orb_chief,
                                                        OptimalTransferPlan plan,
                                                        FormationSpec spec,
                                                        double mu);

// Sprint 26: Autonomous Onboard Operations Scheduler, Solid-State Recorder (SSR) Memory Buffer
// & ECSS/NASA Preliminary Design Review (PDR) Compliance Matrix
MissionOpsReport evaluate_autonomous_ops_and_ssr(OrbitReport rep_init,
                                                 OptimalTransferPlan plan,
                                                 OrbitReport rep_target,
                                                 ADCSReport adcs,
                                                 ThermalReport tcs,
                                                 TelecomReport telecom,
                                                 NavigationReport nav,
                                                 CAMReport cam,
                                                 EndOfLifeReport eol,
                                                 StructuralReport csi,
                                                 FDIRReport fdir,
                                                 FormationFlyingReport form,
                                                 MissionOpsSpec spec);

// Sprint 28: Space Weather, Solar Cycle (F10.7 / Kp), Diurnal Bulge & High-Fidelity Atmospheric Drag
// (Jacchia-Bowman / NRLMSISE-00 Equivalent) for Long-Term Orbital Decay & 25-Year Rule Evaluation
OrbitLifetimeReport simulate_lifetime_orbit_decay(KeplerianElements orb_target,
                                                  OptimalTransferPlan plan,
                                                  const SpacecraftConfig *sc,
                                                  SpaceWeatherSpec spec,
                                                  double mu);

// Sprint 29: Optical Inter-Satellite Links (ISL), 3D Line-of-Sight (LoS) Limb Occultation
// & Walker Constellation Laser Mesh Network Topology (10 Gbps LCT Link Budget)
ISLNetworkReport evaluate_optical_isl_network(KeplerianElements orb_target,
                                              SystemEngineeringReport sys,
                                              ISLNetworkSpec spec);

// Sprint 30: Multi-Ground Station Tracking (NASA DSN / ESA ESTRACK / Polar LEO),
// Topocentric SEZ Azimuth-Elevation-Range (AER), Link Handovers & Antenna Slew Rate
GroundSegmentReport evaluate_multi_station_coverage(KeplerianElements orb_target,
                                                    OptimalTransferPlan plan,
                                                    const SpacecraftConfig *sc,
                                                    GroundNetworkSpec spec,
                                                    double mu);

// Sprint 31: High-Fidelity Electrical Power System (EPS), Conical Umbra/Penumbra Eclipse Model,
// Temperature/TID Solar Array MPPT & Nonlinear Li-Ion Battery Discharge / Cyclic DoD SoH Degradation
double compute_conical_eclipse_fraction(Vector3D sat_pos_eci, Vector3D sun_dir_hat);

EPSHiFiReport simulate_eps_hifi_profile(KeplerianElements orb_init,
                                        KeplerianElements orb_target,
                                        OptimalTransferPlan plan,
                                        ThermalRadReport tcs,
                                        ADCSReport adcs,
                                        TelecomReport telecom,
                                        EPSHiFiSpec spec,
                                        double mu);

#endif // MISSION_ANALYSIS_H
