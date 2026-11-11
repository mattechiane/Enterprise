#ifndef ORBIT_TYPES_H
#define ORBIT_TYPES_H

#include <stdbool.h>

// ==========================================
// 1. COSTANTI FISICHE E PLANETARIE (TERRA - WGS84)
// ==========================================
#define EARTH_MU 398600.4418                    // Costante gravitazionale terrestre (km^3/s^2)
#define EARTH_RADIUS 6378.137                   // Raggio equatoriale terrestre WGS84 (km)
#define EARTH_FLATTENING (1.0 / 298.257223563)  // Schiacciamento terrestre WGS84 (f)
#define EARTH_J2 0.00108262668                  // Coefficiente armonico zonale J2
#define EARTH_OMEGA 7.2921158553e-5             // Velocità angolare di rotazione terrestre (rad/s)
#define G0 9.80665                              // Accelerazione di gravità standard (m/s^2)

// ==========================================
// 2. MATEMATICA E GEOMETRIA
// ==========================================
typedef struct {
    double x;
    double y;
    double z;
} Vector3D;

// Quaternione unitario d'assetto q = [q0 (scalare), q1, q2, q3 (vettoriale)], ||q|| = 1
typedef struct {
    double q0; // Parte scalare: cos(theta / 2)
    double q1; // Componente vettoriale x: e_x * sin(theta / 2)
    double q2; // Componente vettoriale y: e_y * sin(theta / 2)
    double q3; // Componente vettoriale z: e_z * sin(theta / 2)
} Quaternion;

// ==========================================
// 3. STATO DEL VEICOLO SPAZIALE E INPUT FLESSIBILE
// ==========================================
typedef struct {
    Vector3D position; // Posizione spaziale ECI in km (X, Y, Z)
    Vector3D velocity; // Velocità spaziale ECI in km/s (Vx, Vy, Vz)
    double mass;       // Massa attuale del satellite in kg
    double epoch;      // Tempo di simulazione in secondi
} CartesianState;

typedef struct {
    double a;       // Semiasse maggiore (km)
    double e;       // Eccentricità (adimensionale)
    double i;       // Inclinazione (radianti)
    double omega;   // Argomento del pericentro (radianti)
    double RAAN;    // Ascensione Retta del Nodo Ascendente (radianti)
    double theta;   // Anomalia vera (radianti)
} KeplerianElements;

typedef struct {
    double latitude;  // Latitudine geodetica WGS84 (radianti)
    double longitude; // Longitudine terrestre (radianti)
    double altitude;  // Quota reale sopra l'ellissoide WGS84 (km)
} GeodeticCoord;

typedef enum {
    INPUT_KEPLERIAN = 0,
    INPUT_CARTESIAN = 1
} OrbitInputMode;

typedef struct {
    OrbitInputMode mode;
    KeplerianElements kep;
    CartesianState cart;
} OrbitInputSpec;

// ==========================================
// 4. DATABASE PROPELLENTI (SOLIDI, LIQUIDI ED ELETTRICI LOW-THRUST)
// ==========================================
typedef enum {
    PROP_SOLID_HTPB = 0,          // Solido: HTPB/AP (Star-48B) -> Isp = 292 s
    PROP_LIQUID_HYDRAZINE = 1,    // Liquido Monopropellente: Idrazina -> Isp = 235 s
    PROP_LIQUID_MMH_NTO = 2,      // Liquido Bipropellente: MMH/NTO -> Isp = 316 s
    PROP_LIQUID_LOX_LCH4 = 3,     // Liquido Green: LOX/Metano -> Isp = 365 s
    PROP_LIQUID_LOX_LH2 = 4,      // Liquido Criogenico: LOX/Idrogeno -> Isp = 450 s
    PROP_ELECTRIC_HALL_XE = 5,    // Elettrico Hall-Effect: Xenon HET -> Isp = 1650 s, Spinta continua
    PROP_ELECTRIC_ION_NEXT = 6    // Elettrico Gridded Ion: Xenon Ion -> Isp = 3200 s, Spinta continua
} PropellantType;

typedef struct {
    PropellantType type;
    char name[64];
    char category[32];      // "Solid", "Liquid", oppure "Electric"
    double isp;             // Impulso specifico (s)
    double default_thrust_n;// Spinta nominale del motore in Newton (N)
    bool is_low_thrust;     // true se è propulsione elettrica a spirale continua
} PropellantSpec;

// ==========================================
// 5. STAZIONI DI TERRA E OTTIMIZZAZIONE MANOVRE
// ==========================================
typedef struct {
    char name[64];
    GeodeticCoord coord;
    double min_elevation;
} GroundStation;

typedef struct {
    bool is_visible;
    double azimuth;
    double elevation;
    double range;
} AccessMetrics;

typedef struct {
    double dv1;
    double dv2;
    double total_dv;
    double transfer_time;
} HohmannTransfer;

// Piano di trasferimento (supporta sia Impulsivo Chimico sia Low-Thrust Elettrico a Spirale)
typedef struct {
    char maneuver_name[80];
    int num_burns;            // 2 o 3 per chimico, 1 (accensione continua) per elettrico
    bool is_low_thrust;       // true se trasferimento a spirale continua
    double continuous_thrust_n;// Spinta continua applicata (N)
    double spiral_revolutions;// Numero di rivoluzioni complessive della spirale
    bool time_constrained;
    double max_time_limit_s;

    double dv1_mag;
    double dv2_mag;
    double dv3_mag;
    double total_dv;

    Vector3D burn1_dir_rtn;
    Vector3D burn2_dir_rtn;
    Vector3D burn3_dir_rtn;

    Vector3D burn1_pos_eci;
    GeodeticCoord burn1_geo;
    Vector3D burn2_pos_eci;
    GeodeticCoord burn2_geo;
    Vector3D burn3_pos_eci;
    GeodeticCoord burn3_geo;

    double transfer_apogee_km;
    double transfer_time_sec;
    double propellant_used;
    double final_mass;
    PropellantSpec prop_info;
} OptimalTransferPlan;

// Piano di Rendezvous Terminale e Proximity Operations (Clohessy-Wiltshire / Hill Frame LVLH)
typedef struct {
    Vector3D r0_lvlh_m;         // Posizione relativa iniziale LVLH [dx0 (R-bar), dy0 (V-bar), dz0 (H-bar)] (m)
    Vector3D v0_lvlh_ms;        // Velocità relativa iniziale pre-manovra nel frame LVLH (m/s)
    Vector3D v0_post_lvlh_ms;   // Velocità relativa iniziale post-impulso di iniezione CW (m/s)
    double transfer_time_sec;   // Tempo di trasferimento relativo di approccio t_rdv (s)
    double mean_motion_rad_s;   // Moto medio dell'orbita target n = sqrt(mu / a_target^3) (rad/s)
    double initial_range_m;     // Distanza relativa iniziale ||dr0|| (m)

    Vector3D dv_inject_lvlh_ms; // Impulso 1: Iniezione su traiettoria di intercetto CW (m/s)
    Vector3D dv_brake_lvlh_ms;  // Impulso 2: Frenata all'arrivo per Soft-Docking a velocità zero (m/s)
    double dv_inject_ms;        // Modulo del primo impulso (m/s)
    double dv_brake_ms;         // Modulo del secondo impulso (m/s)
    double total_dv_ms;         // Costo totale Delta-V di Rendezvous (m/s)
    double propellant_used_kg;  // Propellente consumato per il Rendezvous (kg)
    double final_mass_kg;       // Massa finale del Chaser dopo il Docking (kg)
    bool converged;             // true se il sistema lineare CW è ben condizionato
} RendezvousPlan;

// ==========================================
// SPRINT 13: MONTE CARLO DISPERSION & TCM BUDGET
// ==========================================
#define MAX_MC_RUNS 100

typedef struct {
    Vector3D arrival_pos_eci; // Coordinate ECI di arrivo perturbate [x, y, z] (km)
    Vector3D dr_rtn_km;       // Scostamento rispetto al punto nominale Burn 2 nel frame RTN [dR, dT, dN] (km)
    double miss_distance_km;  // Errore scalare di posizione all'arrivo dr_miss (km)
    double tcm_dv_ms;         // Delta-V correttivo Mid-Course TCM richiesto a t = 0.5 * TOF (m/s)
} DispersionSample;

typedef struct {
    double sigma_dv_pct;           // Errore 1-sigma sul modulo di spinta Delta-V1 (%)
    double sigma_point_deg;        // Errore 1-sigma angolare di puntamento dell'ugello (deg)
    int num_runs;                  // Numero di traiettorie Monte Carlo campionate (100)

    double mean_miss_km;           // Errore medio di posizione all'arrivo senza correzione (km)
    double std_miss_km;            // Deviazione standard dell'errore di posizione (km)
    double cep_3sigma_km;          // 99.7% 3-Sigma Miss Radius all'arrivo (km)
    double max_miss_km;            // Errore massimo osservato nella nuvola Monte Carlo (km)
    Vector3D ellipsoid_3sigma_rtn; // Semiassi dell'ellissoide di incertezza 3-sigma (3*sigma_R, 3*sigma_T, 3*sigma_N) [km]

    double tcm_dv_mean_ms;         // Delta-V medio per la Mid-Course TCM a t = 0.5 * TOF (m/s)
    double tcm_dv_99_ms;           // Delta-V 99% (3-sigma) TCM Margin per il budget GNC (m/s)
    double tcm_fuel_99_kg;         // Massa di propellente stocastico di riserva per la TCM (kg)

    DispersionSample samples[MAX_MC_RUNS];
} MonteCarloReport;

// ==========================================
// SPRINT 14: 6-DOF ATTITUDE DYNAMICS (ADCS), QUATERNION SLEW & REACTION WHEELS
// ==========================================
#define MAX_ADCS_SAMPLES 160

typedef enum {
    MODE_NADIR_POINTING = 0,      // Asse +Z_body verso Terra (-R), +X_body lungo +T (Osservazione/Comm)
    MODE_SUN_POINTING = 1,        // Normale pannelli solari allineata al vettore Sole s_hat (Carica Batteria)
    MODE_THRUST_VECTOR_ALIGN = 2, // Ugello motore allineato al vettore Delta-V_RTN (Pre-Burn Slew / Spirale)
    MODE_TARGET_TRACKING = 3      // Puntamento verso il Target (Rendezvous) o Ground Station (Radar Pass)
} AttitudeMode;

typedef struct {
    double t_min;             // Tempo di missione cumulativo (min)
    Quaternion q;             // Quaternione unitario di assetto corrente rispetto a LVLH
    double roll_deg;          // Angolo di Roll (phi) rispetto a LVLH (deg)
    double pitch_deg;         // Angolo di Pitch (theta) rispetto a LVLH (deg)
    double yaw_deg;           // Angolo di Yaw (psi) rispetto a LVLH (deg)
    Vector3D omega_deg_s;     // Velocità angolare di corpo [wx, wy, wz] (deg/s)
    double omega_norm_deg_s;  // Modulo velocità angolare ||omega|| (deg/s)
    double point_err_deg;     // Errore istantaneo di puntamento theta_err rispetto a q_cmd (deg)
    double rw_torque_mnm;     // Coppia di controllo erogata dalle Reaction Wheels ||tau_rw|| (mN*m)
    double rw_mom_nms;        // Momento angolare accumulato nelle Reaction Wheels ||h_rw|| (N*m*s)
    double rw_sat_pct;        // Livello di saturazione percentuale rispetto a h_max (%)
    double adcs_power_w;      // Potenza elettrica istantanea assorbita dall'ADCS (W)
    AttitudeMode mode;        // Modalità operativa di puntamento corrente (0..3)
} AttitudeSample;

typedef struct {
    Vector3D inertia_kgm2;          // Tensore principale d'inerzia diag(Jxx, Jyy, Jzz) [kg*m^2]
    double rw_h_max_nms;            // Capacità massima di momento angolare delle Reaction Wheels [N*m*s]
    double slew_time_s;             // Durata della manovra di Pre-Burn Slew per allineamento ugello [s]
    double slew_angle_b1_deg;       // Ampiezza rotazione Pre-Burn 1 Slew [deg]
    double slew_angle_b2_deg;       // Ampiezza rotazione Pre-Burn 2 Slew [deg]
    double steady_point_err_deg;    // Errore di puntamento RMS in regime stazionario [deg]
    double steady_point_err_arcsec; // Errore di puntamento RMS in regime stazionario [arcsec]
    double gg_torque_mnm;           // Coppia disturbante di Gradiente di Gravità max [mN*m]
    double aero_torque_mnm;         // Coppia disturbante Aerodinamica max [mN*m]
    double srp_torque_mnm;          // Coppia disturbante di Pressione di Radiazione Solare [mN*m]
    double peak_rw_torque_mnm;      // Coppia di controllo massima erogata dalle Reaction Wheels [mN*m]
    double max_rw_mom_nms;          // Momento angolare massimo accumulato ||h_rw||_max [N*m*s]
    double max_rw_sat_pct;          // Percentuale massima di saturazione delle ruote [%]
    int    desat_cycles;            // Numero di eventi di desaturazione magnetica (Momentum Dumping)
    double mean_power_w;            // Potenza elettrica media assorbita dal sottosistema ADCS [W]
    double peak_power_w;            // Potenza elettrica di picco durante le manovre di Slew [W]
    int    num_samples;             // Numero di campioni temporali nel profilo ADCS
    AttitudeSample samples[MAX_ADCS_SAMPLES];
} ADCSReport;

// ==========================================
// SPRINT 15: THERMAL CONTROL SYSTEM (TCS) & IONIZING RADIATION (VAN ALLEN / SAA)
// ==========================================
#define MAX_THERMAL_SAMPLES 160
#define STEFAN_BOLTZMANN 5.670374419e-8 // Costante di Stefan-Boltzmann [W/(m^2*K^4)]
#define SOLAR_CONSTANT_W_M2 1361.0      // Flusso solare medio a 1 AU [W/m^2]
#define EARTH_IR_FLUX_W_M2 237.0        // Emissione IR terrestre media OLR [W/m^2]
#define EARTH_ALBEDO_COEFF 0.30         // Coefficiente medio di Albedo terrestre

typedef enum {
    COATING_WHITE_PAINT_SG121 = 0, // White Paint SG121 (alpha_s = 0.20, eps_IR = 0.88)
    COATING_OSR_MIRROR        = 1, // Optical Solar Reflector OSR (alpha_s = 0.10, eps_IR = 0.82)
    COATING_MLI_GOLD_KAPTON   = 2, // MLI Gold Kapton (alpha_s = 0.35, eps_IR = 0.55)
    COATING_BLACK_ANODIZED    = 3, // Black Anodized Al (alpha_s = 0.88, eps_IR = 0.86)
    COATING_CUSTOM            = 4  // Custom Thermo-Optical Properties
} ThermalCoatingType;

typedef struct {
    ThermalCoatingType type;
    char name[64];
    double alpha_s; // Assorbivita solare (alpha_s)
    double eps_ir;  // Emissivita infrarossa (epsilon_IR)
} SurfaceOpticalProps;

typedef struct {
    double t_min;            // Tempo di missione cumulativo (min)
    double t_ext_c;          // Temperatura nodo esterno: Pannelli/Superficie Esterna (°C)
    double t_bus_c;          // Temperatura nodo interno: Bus / Batterie Li-Ion (°C)
    double temp_c;           // Alias temperatura principale del bus (°C)
    double temp_k;           // Temperatura assoluta del bus (K)
    double q_sun_w;          // Flusso solare diretto assorbito Q_sun [W]
    double q_sol_w;          // Alias Q_sun [W]
    double q_alb_w;          // Flusso di Albedo terrestre assorbito Q_albedo [W]
    double q_ir_w;           // Flusso Infrarosso terrestre assorbito Q_IR [W]
    double q_ir_alb_w;       // Somma flussi planetari Q_albedo + Q_IR [W]
    double q_int_w;          // Dissipazione elettronica interna + Riscaldatori [W]
    double q_heater_w;       // Potenza erogata dai Riscaldatori termostatati di sopravvivenza [W]
    double q_rad_w;          // Potenza termica irradiata nello spazio (sigma*eps*A*T^4) [W]
    double inst_dose_urad_h; // Rateo di dose ionizzante istantanea dietro schermatura Al [urad/h]
    int    in_eclipse;       // 1 se il satellite è nel cono d'ombra terrestre, 0 se illuminato
    int    in_saa;           // 1 se il satellite attraversa la South Atlantic Anomaly (SAA)
} ThermalSample;

typedef struct {
    SurfaceOpticalProps coating;     // Proprietà termo-ottiche del rivestimento selezionato
    double alpha_solar;              // Assorbività solare del rivestimento (alpha_s)
    double eps_ir;                   // Emissività infrarossa dei radiatori (epsilon_IR)
    double al_shield_mm;             // Spessore schermatura in Alluminio [mm] (default 4.0 mm)
    double thermal_mass_j_k;         // Capacità termica del nodo Bus interno m*c_p [J/K]
    double ext_thermal_mass_j_k;     // Capacità termica del nodo Pannelli/Guscio esterno [J/K]
    double rad_area_m2;              // Area radiante efficace del satellite [m^2]
    double t_min_op_c;               // Limite operativo minimo batterie Li-Ion [°C] (-10 °C)
    double t_max_op_c;               // Limite operativo massimo batterie Li-Ion [°C] (+45 °C)
    double min_temp_c;               // Temperatura minima raggiunta dal Bus interno [°C]
    double max_temp_c;               // Temperatura massima raggiunta dal Bus interno [°C]
    double mean_temp_c;              // Temperatura media orbitale del Bus interno [°C]
    double ext_min_temp_c;           // Temperatura minima raggiunta dalla superficie esterna [°C]
    double ext_max_temp_c;           // Temperatura massima raggiunta dalla superficie esterna [°C]
    double eq_sun_temp_c;            // Temperatura teorica di equilibrio diurno (Hot Case) [°C]
    double eq_ecl_temp_c;            // Temperatura teorica di equilibrio in eclissi (Cold Case) [°C]
    double mean_q_sol_w;             // Flusso solare medio assorbito [W]
    double mean_q_alb_w;             // Flusso di albedo medio assorbito [W]
    double mean_q_ir_w;              // Flusso IR terrestre medio assorbito [W]
    double mean_q_rad_w;             // Potenza media radiata nello spazio [W]
    double peak_heater_w;            // Potenza massima erogata dai riscaldatori (Heaters) [W]
    double heater_duty_pct;          // Percentuale di tempo con riscaldatori attivi [%]
    double heater_energy_wh;         // Energia elettrica totale consumata dagli heater in missione [Wh]
    double tid_krad_yr;              // Total Ionizing Dose annuale dietro schermatura Al [krad/yr]
    double peak_dose_urad_h;         // Picco di dose radiativa istantanea (SAA / Van Allen) [urad/h]
    double seu_rate_day;             // Rateo stimato di Single Event Upsets [SEU/day]
    double avionics_life_yr;         // Vita utile stimata dell'elettronica di bordo [anni]
    bool   saa_encountered;          // true se l'orbita interseca l'Anomalia del Sud Atlantico (SAA)
    char   saa_status[64];           // Stringa di stato attraversamento SAA / Fasce di Van Allen
    int    num_samples;              // Numero di campioni temporali nel profilo termico/radiativo
    ThermalSample samples[MAX_THERMAL_SAMPLES];
} ThermalRadReport;

typedef ThermalRadReport ThermalReport;

typedef struct {
    KeplerianElements initial_kep;
    double period_min;
    double sunlight_time_min;
    double eclipse_time_min;
    double eclipse_percentage;
    int    ground_passes_24h;
    double total_contact_min;
    double max_elevation_deg;
    double sma_decay_m_per_day;
    double station_keeping_dv_year;
    double station_keeping_fuel_yr;
} OrbitReport;

// ==========================================
// 6. SOTTOSISTEMI E AMBIENTE
// ==========================================
typedef struct {
    bool is_active;
    double thrust;
    double isp;
    Vector3D thrust_dir;
} Thruster;

typedef struct {
    double drag_coeff;
    double drag_area;
    Thruster thruster;
} SpacecraftConfig;

typedef struct {
    bool enable_j2;
    bool enable_drag;
    Vector3D sun_dir;
} ForceModelConfig;

typedef struct {
    bool in_eclipse;
    double air_density;
} SpaceEnvironment;

#endif // ORBIT_TYPES_H
