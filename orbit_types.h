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

// ==========================================
// 5E. SPRINT 16: RF TELECOMMUNICATION LINK BUDGET, Eb/N0 MARGIN & DATA DOWNLINK
// ==========================================
#define MAX_TELECOM_SAMPLES 160
#define BOLTZMANN_DBW_K_HZ  (-228.60) // Costante di Boltzmann k_B in [dBW/K/Hz]

typedef enum {
    BAND_S_2_2GHZ  = 0, // S-Band (f = 2.2 GHz, TT&C robusto, bassa attenuazione atm, 2 Mbps nom)
    BAND_X_8_4GHZ  = 1, // X-Band (f = 8.4 GHz, Science Downlink ad alta velocità, 50 Mbps nom)
    BAND_KA_26GHZ  = 2  // Ka-Band (f = 26.5 GHz, High-Throughput Downlink, 250 Mbps nom, sensibile a pioggia/El)
} RFBandType;

typedef struct {
    double t_min;             // Tempo missione [minuti]
    double slant_range_km;    // Distanza obliqua satellite-stazione di terra rho(t) [km]
    double elevation_deg;     // Angolo di elevazione sull'orizzonte El(t) [deg]
    double fspl_db;           // Free-Space Path Loss FSPL(t) [dB]
    double atm_loss_db;       // Attenuazione atmosferica e perdite di puntamento L_atm(t) + L_point [dB]
    double prx_dbw;           // Potenza ricevuta a terra P_rx(t) [dBW]
    double cn0_dbhz;          // Rapporto Portante/Densità Spettrale di Rumore C/N0(t) [dB-Hz]
    double ebn0_db;           // Rapporto Energia per Bit / Rumore Eb/N0(t) [dB]
    double link_margin_db;    // Margine di collegamento RF Margin(t) = Eb/N0(t) - (Eb/N0)_req [dB]
    double achieved_mbps;     // Throughput istantaneo con Adaptive Coding & Modulation (ACM) [Mbps]
    double cum_data_gb;       // Volume cumulativo di dati scaricati a terra [GB]
    int    in_pass;           // 1 = Finestra di visibilità attiva (El >= El_min e Link Locked), 0 = Fuori vista
} LinkSample;

typedef struct {
    RFBandType band;                 // Banda di frequenza RF selezionata (S / X / Ka)
    char   band_name[48];            // Nome descrittivo della banda (es. "X-Band (8.4 GHz Science)")
    char   modulation_scheme[48];    // Schema di modulazione e codifica (es. "8PSK / 16APSK + LDPC ACM")
    double freq_ghz;                 // Frequenza portante di downlink f [GHz]
    double tx_power_w;               // Potenza trasmettitore di bordo P_tx [W]
    double tx_power_dbw;             // Potenza trasmettitore convertita in [dBW] = 10*log10(P_tx)
    double tx_gain_dbi;              // Guadagno antenna di bordo G_tx [dBi]
    double rx_dish_m;                // Diametro antenna parabolica stazione di terra [m] (3.5 m Milano-Bovisa)
    double rx_gain_dbi;              // Guadagno antenna stazione di terra G_rx [dBi]
    double eirp_dbw;                 // Equivalent Isotropically Radiated Power EIRP = P_tx + G_tx - L_line [dBW]
    double sys_noise_temp_k;         // Temperatura equivalente di rumore di sistema T_sys [K]
    double n0_dbw_hz;                // Densità spettrale di potenza di rumore N0 = k_B + 10*log10(T_sys) [dBW/Hz]
    double zenith_atm_loss_db;       // Attenuazione atmosferica allo zenit L_zenith [dB]
    double pointing_loss_db;         // Perdita di puntamento d'assetto ADCS L_point [dB]
    double lock_threshold_db;        // Soglia minima di decodifica QPSK/8PSK + LDPC (4.5 dB)
    double req_ebn0_db;              // Eb/N0 richiesto incluso margine di sicurezza (+3.0 dB => 7.5 dB)
    double nom_bitrate_mbps;         // Data rate nominale della banda [Mbps] (2 / 50 / 250 Mbps)
    double max_achieved_mbps;        // Data rate di picco raggiunto con ACM [Mbps]
    double mean_fspl_db;             // Free-Space Path Loss medio durante i passaggi [dB]
    double max_fspl_db;              // Free-Space Path Loss massimo lungo l'orbita [dB]
    double peak_ebn0_db;             // Picco di Eb/N0 al culminare del passaggio (Max Elevation) [dB]
    double min_margin_db;            // Margine di link minimo registrato durante i passaggi utili [dB]
    double peak_margin_db;           // Margine di link massimo al culminare del passaggio [dB]
    int    pass_windows;             // Numero di finestre di visibilità/downlink nella timeline
    double contact_time_min;         // Tempo totale di contatto radar in finestra di missione [min]
    double mission_data_gb;          // Volume dati scaricato lungo la finestra di missione simulata [GB]
    double total_data_gb_24h;        // Volume totale di dati scaricabili a terra in 24 ore [GB/24h]
    int    num_samples;              // Numero di campioni temporali nel profilo telecom
    LinkSample samples[MAX_TELECOM_SAMPLES];
} TelecomReport;

// ==========================================
// 5F. SPRINT 17: PAYLOAD SENSOR SWATH, WALKER DELTA CONSTELLATION & MASS/COST MODEL
// ==========================================
#define MAX_WALKER_SATS 128

typedef struct {
    double lat_deg;      // Latitudine sub-satellitare a t = 0 [deg]
    double lon_deg;      // Longitudine sub-satellitare a t = 0 [deg]
    int    plane_id;     // Indice del piano orbitale p in [0, P-1]
    int    sat_id;       // Indice globale del satellite in [0, T-1]
    double raan_deg;     // Ascensione retta del nodo ascendente del piano [deg]
    double arg_lat_deg;  // Argomento di latitudine u = omega + theta [deg]
} WalkerSatPoint;

typedef struct {
    double half_fov_deg;       // Semiapertura angolare del sensore theta_(FOV/2) [deg] (default 25.0 deg)
    double opt_aperture_m;     // Diametro apertura ottica D_opt [m] (default 0.35 m)
    double wavelength_nm;      // Lunghezza d'onda operativa lambda [nm] (default 550.0 nm)
    int    walker_total_sats;  // Numero totale di satelliti della costellazione T (default 24)
    int    walker_planes;      // Numero di piani orbitali P (default 3)
    int    walker_phasing;     // Parametro di fasatura inter-piano F (default 1)
} PayloadAndConstellationSpec;

typedef struct {
    PayloadAndConstellationSpec spec;

    // 1. Geometria Sensore Ottico/SAR, Swath Width e Ground Sample Distance (GSD)
    double target_alt_km;            // Quota orbitale operativa media h = a_target - R_E [km]
    double earth_central_angle_deg;  // Angolo centrale terrestre di semi-swath Lambda [deg]
    double swath_width_km;           // Larghezza striscia osservata a terra S_swath [km]
    double footprint_radius_km;      // Raggio dell'impronta circolare al suolo [km]
    double gsd_nadir_m;              // Risoluzione spaziale al nadir (Ground Sample Distance) [m]
    double area_coverage_rate_km2_s; // Tasso di copertura areale istantaneo [km^2/s]

    // 2. Costellazione Simmetrica di Walker Delta (i : T / P / F)
    double inc_deg;                  // Inclinazione dei piani della costellazione [deg]
    int    sats_per_plane;           // Numero di satelliti per piano S = T / P
    double mean_revisit_min;         // Tempo medio di rivisita sulla latitudine target [min]
    double global_coverage_pct;      // Percentuale di superficie terrestre coperta istantaneamente [%]
    double mean_vis_sats;            // Numero medio di satelliti simultaneamente in vista N_vis
    double pdop_geom;                // Geometric Position Dilution of Precision (PDOP)
    int    num_walker_sats;          // Numero effettivo di sub-points esportati
    WalkerSatPoint walker_subpoints[MAX_WALKER_SATS];

    // 3. Bilancio di Massa (SMAD Dry Mass Breakdown + Propellant Margin)
    double wet_mass_kg;              // Massa totale al lancio (Wet Mass m_0) [kg]
    double propellant_transfer_kg;   // Propellente consumato per il trasferimento orbitale [kg]
    double propellant_margin_kg;     // Riserva propellente (Rendezvous CW + TCM 3-Sigma + SK 1 anno) [kg]
    double total_propellant_kg;      // Massa totale di propellente imbarcato m_prop [kg]
    double dry_mass_kg;              // Massa a secco del satellite m_dry = m_wet - m_prop [kg]
    double dry_wet_ratio_pct;        // Rapporto Dry/Wet Mass [%]
    double mass_payload_kg;          // Payload Ottico/SAR (28% m_dry) [kg]
    double mass_struct_kg;           // Struttura Primaria e Meccanismi (22% m_dry) [kg]
    double mass_eps_kg;              // EPS Pannelli Solari, Li-Ion & PCDU (19% m_dry) [kg]
    double mass_prop_dry_kg;         // Serbatoi, Tubazioni e Propulsori a secco (11% m_dry) [kg]
    double mass_adcs_kg;             // ADCS Reaction Wheels, Star Trackers & MTQ (9% m_dry) [kg]
    double mass_ttc_obdh_kg;         // TT&C Transponder RF & On-Board Computer (6% m_dry) [kg]
    double mass_tcs_kg;              // Controllo Termico MLI, Radiatori & Heaters (5% m_dry) [kg]

    // 4. Stima Parametrica dei Costi di Missione (SMAD / NASA NICM in M€)
    double cost_platform_meur;       // Costo Piattaforma Bus (NRE + TFU) [M€]
    double cost_payload_meur;        // Costo Sviluppo e Qualifica Payload [M€]
    double cost_spacecraft_meur;     // Costo Totale Singolo Satellite C_sc [M€]
    double cost_launch_meur;         // Costo di Lancio C_launch (tariffa per kg di Wet Mass) [M€]
    double cost_ops_5yr_meur;        // Costo Operativo MO&DA su 5 anni (Ground Segment & Flight Ops) [M€]
    double total_mission_cost_meur;  // Costo Totale Missione Singolo Satellite (Lifecycle 5 yr) [M€]
    double constellation_cost_meur;  // Costo Stimato Costellazione Completa Walker (con curva di Wright) [M€]
} SystemEngineeringReport;

// ==========================================
// 5G. SPRINT 18: AUTONOMOUS ONBOARD NAVIGATION & EXTENDED KALMAN FILTER (EKF)
// ==========================================
#define MAX_EKF_SAMPLES 160

typedef struct {
    double gps_sigma_pos_m;      // Rumore di misura GNSS 1-sigma sulla posizione [m] (default 15.0 m)
    double gps_sigma_vel_ms;     // Rumore di misura GNSS 1-sigma sulla velocità [m/s] (default 0.15 m/s)
    double init_sigma_pos_m;     // Incertezza iniziale a priori 1-sigma sulla posizione sigma_(r,0) [m] (default 500.0 m)
    double init_sigma_vel_ms;    // Incertezza iniziale a priori 1-sigma sulla velocità sigma_(v,0) [m/s] (default 2.50 m/s)
    double proc_noise_accel_ms2; // Densità rumore di processo dinamico q_a (drag / J3..J4 non modellati) [m/s^2]
} NavigationEKFSpec;

typedef struct {
    double t_min;              // Tempo di missione cumulativo [min]
    double raw_err_m;          // Errore istantaneo della misura GNSS grezza ||z_r - r_true|| [m]
    double ekf_err_m;          // Errore reale di posizione stimata dall'EKF ||r_hat - r_true|| [m]
    double cov_3sigma_m;       // Inviluppo teorico 3-sigma di posizione 3 * sqrt(P_xx + P_yy + P_zz) [m]
    double vel_err_ms;         // Errore reale di velocità stimata dall'EKF ||v_hat - v_true|| [m/s]
    double cov_vel_3sigma_ms;  // Inviluppo teorico 3-sigma di velocità 3 * sqrt(P_vx + P_vy + P_vz) [m/s]
    double kalman_gain_norm;   // Norma di Frobenius del Guadagno di Kalman ||K_k||
} EKFSample;

typedef struct {
    NavigationEKFSpec spec;
    double raw_pos_rmse_m;          // RMSE di posizione delle misure GNSS non filtrate [m]
    double steady_pos_rmse_m;       // RMSE di posizione EKF a regime (Steady-State) [m]
    double steady_vel_rmse_ms;      // RMSE di velocità EKF a regime (Steady-State) [m/s]
    double steady_vel_rmse_cms;     // RMSE di velocità EKF a regime convertito in [cm/s]
    double final_cov_3sigma_m;      // Limite teorico 3-sigma di posizione a fine missione [m]
    double final_cov_vel_3sigma_ms; // Limite teorico 3-sigma di velocità a fine missione [m/s]
    double convergence_time_min;    // Tempo di convergenza del filtro sotto la soglia 1-sigma GNSS [min]
    double noise_reduction_pct;     // Fattore di attenuazione del rumore GNSS (1 - RMSE_ekf / RMSE_raw) [%]
    double consistency_pct;         // Percentuale di campioni contenuti nell'inviluppo teorico 3-sigma [%]
    double mean_kalman_gain;        // Norma media del Guadagno di Kalman a regime ||K_k||
    int    num_samples;             // Numero di campioni temporali esportati nel profilo EKF
    EKFSample samples[MAX_EKF_SAMPLES];
} NavigationReport;

// ==========================================
// 5H. SPRINT 19: SPACE SITUATIONAL AWARENESS (SSA), B-PLANE CONJUNCTION & CAM
// ==========================================
#define MAX_CAM_SAMPLES 60

typedef struct {
    char   debris_id[64];       // Identificativo oggetto catalogato NORAD/DISCOS (es. "COSMOS-2251 DEB [NORAD #34454]")
    double nominal_miss_m;      // Nominal Miss Distance pre-manovra al TCA [m] (default 180.0 m)
    double rel_vel_kms;         // Velocità relativa d'incontro al TCA v_rel [km/s]
    double hbr_m;               // Combined Hard-Body Radius R_HBR = R_sat + R_debris [m] (default 5.0 m)
    double sigma_xi_m;          // Semiasse 1-sigma covarianza combinata B-Plane asse xi [m]
    double sigma_zeta_m;        // Semiasse 1-sigma covarianza combinata B-Plane asse zeta [m]
    double lead_orbits;         // Anticipo esecuzione manovra CAM prima del TCA [orbite] (default 1.5 orbite)
} ConjunctionEvent;

typedef struct {
    double dv_cam_cms;          // Impulso tangenziale CAM Delta-V [cm/s]
    double dv_cam_ms;           // Impulso tangenziale CAM Delta-V [m/s]
    double xi_m;                // Coordinata B-Plane xi al TCA [m]
    double zeta_m;              // Coordinata B-Plane zeta al TCA [m]
    double miss_dist_m;         // Distanza di separazione risultante nel B-Plane [m]
    double pc;                  // Probabilità di collisione 2D di Foster/Alfriend P_c [-]
    double log10_pc;            // Logaritmo in base 10 della probabilità di collisione log10(P_c)
} BPlaneSample;

typedef struct {
    ConjunctionEvent event;
    double tca_lead_time_min;   // Tempo di anticipo manovra prima del TCA [min]
    double pre_xi_m;            // Coordinata nominale pre-CAM xi_0 nel B-Plane [m]
    double pre_zeta_m;          // Coordinata nominale pre-CAM zeta_0 nel B-Plane [m]
    double pre_miss_m;          // Miss Distance pre-CAM d_miss [m]
    double pre_pc;              // Probabilità di collisione pre-CAM P_c [-]
    double pre_log10_pc;        // log10(P_c) pre-CAM
    double pre_mahalanobis_d;   // Distanza di Mahalanobis 2D pre-CAM [sigma]
    double req_dv_cam_ms;       // Impulso tangenziale richiesto Delta-V_CAM [m/s]
    double req_dv_cam_cms;      // Impulso tangenziale richiesto Delta-V_CAM [cm/s]
    double cam_fuel_g;          // Propellente consumato per la CAM [grammi]
    double delta_a_m;           // Variazione semiasse maggiore indotta dal phasing burn [m]
    double post_xi_m;           // Coordinata post-CAM xi nel B-Plane [m]
    double post_zeta_m;         // Coordinata post-CAM zeta nel B-Plane [m]
    double post_miss_m;         // Distanza di sicurezza post-CAM nel B-Plane [m] (>= 2500 m)
    double post_pc;             // Probabilità di collisione post-CAM P_c (< 1e-7)
    double post_log10_pc;       // log10(P_c) post-CAM
    double post_mahalanobis_d;  // Distanza di Mahalanobis 2D post-CAM [sigma] (>= 5 sigma)
    char   alert_status[48];    // Stato allarme ESA/NASA ("RED ALERT", "YELLOW ADVISORY", "GREEN SAFE")
    int    num_samples;         // Numero di campioni della curva di trade-off CAM
    BPlaneSample samples[MAX_CAM_SAMPLES];
} CAMReport;

// ==========================================
// 5I. SPRINT 20 & SPRINT 27: END-OF-LIFE (EoL) DEORBIT, 3-PHASE EDL, PARACHUTE & POWERED DESCENT
// ==========================================
#define MAX_REENTRY_SAMPLES 160

typedef struct {
    double t_sec;               // Tempo trascorso dall'Entry Interface (h_EI = 120 km) [s]
    double alt_km;              // Quota sopra l'ellissoide terrestre h(t) [km]
    double vel_kms;             // Velocità relativa atmosferica V(t) [km/s]
    double fpa_deg;             // Angolo di rampa di volo (Flight-Path Angle) gamma(t) [deg]
    double decel_g;             // Fattore di carico assiale totale (aerodinamico + paracadute + retrorazzi) [g]
    double heat_flux_wcm2;      // Flusso termico convettivo di ristagno di Sutton-Graves q_dot(t) [W/cm^2]
    double qdyn_kpa;            // Pressione dinamica q_dyn(t) = 0.5 * rho(h) * V^2 [kPa]
    double downrange_km;        // Distanza percorsa lungo traccia a terra (Downrange) s(t) [km]
    double mach;                // Numero di Mach locale M = V / a(h) [-]
    double cd_eff;              // Coefficiente di resistenza aerodinamica effettivo C_D(M) [-]
    double throttle_pct;        // Manetta retrorazzi Powered Descent [%] (0 in AERO/CHUTE, 10-100% in POWERED)
    char   edl_phase[12];       // Fase EDL corrente: "AERO" | "CHUTE" | "POWERED"
} ReentrySample;

typedef struct {
    double entry_fpa_deg;           // Flight-Path Angle all'Entry Interface gamma_EI [deg] (default -1.5 deg)
    double nose_radius_m;           // Raggio di curvatura dello scudo termico R_n [m] (default 0.50 m)
    double chute_deploy_mach;       // Numero di Mach di apertura paracadute Disk-Gap-Band [-] (default 2.1)
    double chute_max_alt_km;        // Quota massima ammissibile apertura paracadute [km] (default 12.0 km)
    double chute_area_m2;           // Area nominale paracadute supersonico dispiegato A_chute [m^2] (default 185.0 m^2)
    double powered_descent_alt_m;   // Quota di sgancio paracadute e accensione retrorazzi [m] (default 2000.0 m)
    double target_touchdown_vel_ms; // Velocità target di soft-touchdown a h = 0 [m/s] (default 1.80 m/s < 2.5 m/s)
} EDLSequenceSpec;

typedef struct {
    EDLSequenceSpec spec;
    // Phase 2: Supersonic Disk-Gap-Band Parachute Deployment
    double chute_deploy_time_sec;   // Istante di apertura paracadute [s]
    double chute_deploy_alt_km;     // Quota effettiva di apertura paracadute [km]
    double chute_deploy_vel_ms;     // Velocità all'apertura del paracadute [m/s]
    double chute_deploy_mach;       // Numero di Mach effettivo all'apertura del paracadute [-]
    double chute_snatch_g;          // Picco impulsivo di decelerazione (Snatch Force Shock) all'apertura [g]
    double chute_snatch_kn;         // Forza di picco sul cavo/riser del paracadute all'apertura [kN]

    // Phase 3: Powered Descent (Gravity Turn / Skycrane Retro-Propulsion)
    double pd_start_time_sec;       // Istante di sgancio paracadute e accensione Powered Descent [s]
    double pd_start_alt_m;          // Quota effettiva di inizio Powered Descent [m]
    double pd_start_vel_ms;         // Velocità all'inizio del Powered Descent [m/s]
    double pd_peak_thrust_kn;       // Spinta massima erogata dai retrorazzi [kN]
    double pd_mean_throttle_pct;    // Manetta media durante il Powered Descent [%]
    double pd_dv_ms;                // Delta-V propulsivo speso nella fase di Powered Descent [m/s]
    double pd_fuel_used_kg;         // Propellente consumato nel Powered Descent [kg]
    double touchdown_vel_ms;        // Velocità finale di Touchdown a h = 0 m [m/s] (< 2.5 m/s)
    double landed_mass_kg;          // Massa finale atterrata al suolo [kg]
    bool   soft_landing_achieved;   // true se V_td < 2.5 m/s
    char   edl_status[64];          // Stato qualifica sequenza EDL NASA JPL
} EDLReport;

typedef struct {
    double entry_alt_km;            // Quota di Entry Interface h_EI [km] (120.0 km)
    double entry_fpa_deg;           // Flight-Path Angle all'Entry Interface gamma_EI [deg] (default -1.5 deg)
    double nose_radius_m;           // Raggio di curvatura effettivo del muso/bus R_n [m] (default 0.50 m)
    double ballistic_coeff_kgm2;    // Coefficiente balistico beta = m / (C_D * A) [kg/m^2]
    double deorbit_perigee_km;      // Quota del perigeo di cattura atmosferica h_p [km] (80.0 km)
    double deorbit_dv_ms;           // Impulso retrogrado di deorbiting attivo Delta-V_deorbit [m/s]
    double deorbit_fuel_kg;         // Propellente richiesto per la manovra di fine vita m_fuel,EoL [kg]
    double passive_decay_years;     // Tempo stimato di decadimento orbitale passivo naturale [anni]
    char   compliance_status[64];   // Conformità ESA Zero-Debris / FCC 5-Year / IADC 25-Year Rule
    double peak_heat_flux_wcm2;     // Picco di flusso termico di ristagno q_dot_max [W/cm^2]
    double peak_heat_alt_km;        // Quota del picco termico aerotermodinamico [km]
    double total_heat_load_jcm2;    // Carico termico integrale totale Q = int(q_dot dt) [J/cm^2]
    double peak_decel_g;            // Decelerazione assiale massima ipersonica n_g,max [g]
    double peak_decel_alt_km;       // Quota del picco di decelerazione meccanica [km]
    double peak_qdyn_kpa;           // Pressione dinamica massima q_dyn,max [kPa]
    double breakup_alt_km;          // Quota di Demise / Peak TPS Ablation [km]
    double breakup_vel_kms;         // Velocità al momento del picco ablativo/breakup [km/s]
    double total_entry_time_sec;    // Durata complessiva della discesa atmosferica da 120 km a 0 km [s]
    double total_downrange_km;      // Downrange complessivo dall'Entry Interface al Touchdown [km]
    EDLReport edl;                  // Sprint 27: Report completo sequenza EDL a 3 fasi (Aero + Chute + Powered)
    int    num_samples;             // Numero di campioni temporali della traiettoria 3-DOF
    ReentrySample samples[MAX_REENTRY_SAMPLES];
} EndOfLifeReport;

// ==========================================
// 5J. SPRINT 21: CISLUNAR ASTRODYNAMICS, CR3BP, LAGRANGE POINTS & LUNAR HALO
// ==========================================
#define MOON_MU           4902.800066  // Costante gravitazionale lunare (km^3/s^2)
#define MOON_RADIUS       1737.4       // Raggio medio lunare (km)
#define CR3BP_MU_EM       0.0121505856 // Parametro di massa adimensionale Terra-Luna mu = m_M / (m_E + m_M)
#define CR3BP_L_STAR_KM   384400.0     // Distanza caratteristica Terra-Luna L* (km)
#define CR3BP_T_STAR_SEC  375190.258   // Tempo caratteristico T* = sqrt(L*^3 / (mu_E + mu_M)) (s)
#define MAX_CR3BP_SAMPLES 200

typedef enum {
    CISLUNAR_L2_HALO     = 0, // L2 Far-Side Southern Halo Orbit (Lunar Gateway / NRHO Regime)
    CISLUNAR_L1_GATEWAY  = 1, // L1 Cislunar Neck Gateway / Lyapunov Orbit
    CISLUNAR_FREE_RETURN = 2  // Lunar Free-Return Circum-Lunar Figure-8 Flyby
} CislunarRegimeType;

typedef struct {
    int    target_regime;      // 0 = L2 Far-Side Halo, 1 = L1 Gateway, 2 = Free-Return Flyby
    double halo_az_km;         // Ampiezza fuori dal piano Halo A_z [km] (default 15000.0 km)
    double leo_parking_alt_km; // Quota orbita terrestre di partenza pre-TLI [km] (default 400.0 km)
} CR3BPSpec;

typedef struct {
    double t_days;        // Tempo di missione nel frame sinodico [giorni]
    double x_syn_km;      // Coordinata X sinodica baricentrica co-rotante [km]
    double y_syn_km;      // Coordinata Y sinodica baricentrica co-rotante [km]
    double z_syn_km;      // Coordinata Z sinodica fuori dal piano [km]
    double jacobi_c;      // Costante di Jacobi istantanea C_J [-]
    double dist_moon_km;  // Distanza dal centro della Luna r_2 [km]
} CR3BPSample;

typedef struct {
    CR3BPSpec spec;
    char   regime_name[64];          // Nome descrittivo del regime cislunare selezionato
    double mu_em;                    // Rapporto di massa adimensionale Terra-Luna (~0.0121505856)
    double l_star_km;                // Distanza caratteristica Terra-Luna L* [km] (384400.0 km)
    double t_star_sec;               // Tempo caratteristico T* [s] (~375190.26 s)
    double v_star_kms;               // Velocità caratteristica V* = L* / T* [km/s] (~1.024548 km/s)
    double earth_x_km;               // Posizione Terra sull'asse X sinodico (-mu * L*) [km]
    double moon_x_km;                // Posizione Luna sull'asse X sinodico ((1 - mu) * L*) [km]

    // Coordinate esatte dei 5 Punti di Librazione di Lagrange (Newton-Raphson sul Quintico di Eulero)
    double l1_x_km;                  // Punto collineare cislunare L1 X [km] (~321710 km dal baricentro)
    double l2_x_km;                  // Punto collineare translunare L2 X [km] (~444244 km dal baricentro)
    double l3_x_km;                  // Punto collineare opposto L3 X [km] (~-386346 km dal baricentro)
    double l4_x_km;                  // Punto triangolare equilatero L4 X [km] ((0.5 - mu) * L*)
    double l4_y_km;                  // Punto triangolare equilatero L4 Y [km] (+sqrt(3)/2 * L*)
    double l5_x_km;                  // Punto triangolare equilatero L5 X [km] ((0.5 - mu) * L*)
    double l5_y_km;                  // Punto triangolare equilatero L5 Y [km] (-sqrt(3)/2 * L*)
    double cj_l1;                    // Energia di soglia di Jacobi al collo di L1 (~3.188341)
    double cj_l2;                    // Energia di soglia di Jacobi al collo di L2 (~3.172160)
    double cj_l3;                    // Energia di soglia di Jacobi a L3 (~3.012147)
    double cj_l4_l5;                 // Energia minima di Jacobi a L4/L5 (~2.987997)

    // Budget Propulsivo Cislunare (Trans-Lunar Injection + Lunar / Halo Orbit Insertion)
    double dv_tli_ms;                // Impulso Trans-Lunar Injection da LEO [m/s] (~3130 m/s)
    double dv_loi_ms;                // Impulso Halo / Lunar Orbit Insertion [m/s] (~420 - 650 m/s)
    double total_cislunar_dv_ms;     // Costo totale Delta-V cislunare (TLI + LOI) [m/s]
    double cislunar_fuel_kg;         // Consumo propellente equivalente per iniezione cislunare [kg]
    double tof_days;                 // Tempo di volo del trasferimento cislunare [giorni] (~4.6 - 5.2 d)
    double halo_period_days;         // Periodo dell'orbita periodica Halo / Lyapunov [giorni] (~14.6 d)
    double halo_ax_km;               // Semiasse in-plane X dell'orbita Halo [km]
    double halo_ay_km;               // Semiasse in-plane Y dell'orbita Halo [km]
    double halo_az_km;               // Ampiezza out-of-plane Z dell'orbita Halo [km]
    double perilune_alt_km;          // Quota minima sopra la superficie lunare [km]

    // Costante di Jacobi e verifica di conservazione numerica RK4 (|dC_J / C_J0| < 1e-8)
    double jacobi_c;                 // Costante di Jacobi del trasferimento cislunare C_J [-]
    double jacobi_halo_c;            // Costante di Jacobi dell'orbita Halo periodica [-]
    double delta_cj_rel;             // Errore relativo massimo di conservazione |dC_J / C_J0| (< 1e-8)

    int    num_transfer_samples;
    CR3BPSample traj_cislunar_transfer[MAX_CR3BP_SAMPLES];
    int    num_halo_samples;
    CR3BPSample traj_halo_orbit[MAX_CR3BP_SAMPLES];
} CislunarReport;

// ==========================================
// 5K. SPRINT 22: FLEXIBLE STRUCTURAL DYNAMICS, PROPELLANT SLOSHING & CSI
// ==========================================
#define MAX_CSI_SAMPLES 140

typedef struct {
    double flex_freq_hz;       // Frequenza 1° modo flessionale pannelli solari f_flex,1 [Hz] (default 0.85 Hz)
    double flex_damping_zeta;  // Smorzamento strutturale modale zeta_flex [-] (default 0.015)
    double slosh_damping_zeta; // Smorzamento setti antisciacquio serbatoio (Baffle Damping) zeta_slosh [-] (default 0.030)
    double tank_radius_m;      // Raggio interno serbatoio propellente R_tank [m] (default 0.38 m)
    double tank_fill_ratio;    // Rapporto riempimento liquido h_f / R_tank [-] (default 1.15)
    double array_span_m;       // Lunghezza singola ala pannello solare dispiegato L_wing [m] (default 3.80 m)
    bool   enable_notch_filter;// Stato filtro digitale Notch GNC anti-risonanza CSI (default true)
} StructuralCSISpec;

typedef struct {
    double t_sec;              // Tempo dal transitorio di accensione motore / Slew ADCS [s]
    double tip_defl_mm;        // Deflessione trasversale punta pannello solare delta_tip(t) [mm]
    double slosh_force_n;      // Forza laterale d'onda di sloshing sul serbatoio F_slosh(t) [N]
    double slosh_torque_mnm;   // Coppia di disturbo di sloshing sull'assetto tau_slosh(t) [mN*m]
    double jitter_arcsec;      // Micro-vibrazione / Pointing Jitter residuo sull'asse ottico [arcsec]
    double bode_freq_hz;       // Frequenza nel diagramma di Bode della risposta CSI [Hz]
    double bode_mag_db;        // Ampiezza Bode ad anello aperto/compensato con Notch Filter [dB]
    double bode_unfilt_db;     // Ampiezza Bode pianta flessibile non filtrata (senza Notch) [dB]
} SloshFlexSample;

typedef struct {
    StructuralCSISpec spec;
    // Parametri Modello Equivalente di Abramson (Sloshing) & Trave a Sbalzo (Pannelli Solari)
    double total_fuel_mass_kg;       // Massa totale di propellente liquido a bordo m_prop [kg]
    double slosh_mass_kg;            // Massa liquida oscillante (1° modo di Abramson ~0.45 * m_prop) [kg]
    double fixed_fuel_mass_kg;       // Massa liquida solidale al serbatoio (~0.55 * m_prop) [kg]
    double eff_accel_ms2;            // Accelerazione efficace di eccitazione a_eff (spinta/slew) [m/s^2]
    double slosh_freq_hz;            // Frequenza naturale di sloshing f_slosh = omega_slosh / (2*pi) [Hz]
    double slosh_omega_rads;         // Pulsazione naturale di sloshing di Abramson omega_slosh [rad/s]
    double slosh_pendulum_len_m;     // Lunghezza del pendolo equivalente di Abramson L_pend [m]
    double flex_freq_hz;             // Frequenza fondamentale 1° modo flessionale f_flex,1 [Hz]
    double flex_freq_mode2_hz;       // Frequenza 2° modo flessionale (torsione/flessione superiore) [Hz]
    double modal_coupling_gamma;     // Coefficiente di partecipazione/accoppiamento modale rotazionale Gamma_flex [kg^0.5*m]

    // Risposta Temporale Smorzata post-Burn/Slew e Pointing Jitter
    double peak_tip_defl_mm;         // Deflessione massima punta pannelli solari max|delta_tip(t)| [mm]
    double peak_slosh_disp_mm;       // Spostamento massimo centro di massa onda sloshing max|x_slosh(t)| [mm]
    double peak_slosh_force_n;       // Forza massima laterale di sloshing sulle pareti/baffle [N]
    double peak_slosh_torque_mnm;    // Coppia massima di disturbo di sloshing tau_slosh,max [mN*m]
    double rms_jitter_arcsec;        // Pointing Jitter RMS durante il transitorio [arcsec]
    double peak_jitter_arcsec;       // Pointing Jitter di picco [arcsec]
    double settling_time_sec;        // Tempo di assestamento al di sotto della soglia di fine-pointing [s]

    // Control-Structure Interaction (CSI) & GNC Digital Notch Filter
    double adcs_bandwidth_hz;        // Banda passante anello di controllo d'assetto ADCS f_bw [Hz]
    double notch_center_hz;          // Frequenza centrale del filtro Notch digitale GNC [Hz]
    double notch_attenuation_db;     // Attenuazione del picco di risonanza introdotta dal Notch Filter [dB]
    double peak_res_unfilt_db;       // Picco di risonanza strutturale Bode senza filtro Notch [dB]
    double peak_res_filt_db;         // Picco residuo Bode dopo compensazione con filtro Notch [dB]
    double csi_gain_margin_db;       // Margine di guadagno di stabilità CSI [dB]
    char   notch_status[64];         // Stato del filtro Notch ADCS e stabilità CSI

    // Verifica di Rigidezza Strutturale al Lancio (Vega-C / Falcon 9 / Ariane 6) & QSL Margin of Safety
    double launch_axial_freq_hz;     // Frequenza fondamentale assiale satellite ripiegato f_ax [Hz] (> 35 Hz)
    double launch_lateral_freq_hz;   // Frequenza fondamentale laterale satellite ripiegato f_lat [Hz] (> 15 Hz)
    double qsl_axial_g;              // Carico Quasi-Statico (QSL) assiale di progetto al lancio [g]
    double qsl_lateral_g;            // Carico Quasi-Statico (QSL) laterale di progetto al lancio [g]
    double max_von_mises_mpa;        // Tensione equivalente massima di von Mises nel cilindro primario [MPa]
    double allowable_stress_mpa;     // Tensione ammissibile lega aerospaziale Al-7075-T7351 [MPa]
    double launch_mos_pct;           // Margin of Safety strutturale ai carichi di lancio MoS [%] (> 0%)
    char   launch_compliance[64];    // Conformità inviluppo lanciatori (Vega-C / Falcon 9 / Ariane 6)

    int    num_samples;
    SloshFlexSample samples[MAX_CSI_SAMPLES];
} StructuralReport;

// ==========================================
// 5L. SPRINT 23: AUTONOMOUS FDIR, NORAD TLE PARSER & CCSDS OEM 2.0 EXPORT
// ==========================================
#define MAX_FDIR_SAMPLES 150

typedef enum {
    FDIR_FAULT_RW2_ACTUATOR = 0, // Guasto tachimetro/driver Ruota di Reazione #2 -> Switchover su Skewed RW #4 B-Side
    FDIR_FAULT_SOLAR_SHUNT  = 1, // Cortocircuito sezione SADA Pannello Solare -> Undervoltage Bus 28V -> Regolatore B-Side
    FDIR_FAULT_STAR_TRACKER = 2  // Accecamento/SEU Star Tracker Primario STR-A -> Commutazione su Hot-Standby STR-B
} FDIRFaultScenario;

typedef enum {
    FDIR_STATE_NOMINAL          = 0,
    FDIR_STATE_FAULT_DETECTED   = 1,
    FDIR_STATE_ISOLATION        = 2,
    FDIR_STATE_REDUNDANT_SWITCH = 3,
    FDIR_STATE_SAFE_SUN_ACQ     = 4,
    FDIR_STATE_NOMINAL_RECOVERY = 5
} FDIRStateId;

typedef struct {
    int    sat_catalog_num;      // Numero di catalogo NORAD (es. 25544 per ISS, 40697 per Sentinel-2A)
    int    epoch_year;           // Anno epoca TLE (4 cifre, es. 2026)
    double epoch_day;            // Giorno frazionario dell'anno (1.0 - 366.99999999)
    double mean_motion_rev_day;  // Moto medio n [rev/day]
    double mean_anomaly_deg;     // Anomalia media M [deg]
    double eccentric_anom_deg;   // Anomalia eccentrica risolta E [deg]
    bool   checksum_valid;       // Validità checksum modulo-10 sulle 2 righe TLE
    KeplerianElements kep;       // Elementi Kepleriani convertiti (a [km], e, i [rad], RAAN [rad], omega [rad], theta [rad])
} NoradTLEData;

typedef struct {
    int    fault_scenario;       // 0 = RW-2 Actuator Failure, 1 = Solar Array Shunt Short, 2 = Star Tracker STR-A Blinded
    double fault_inject_time_s;  // Istante di iniezione dell'anomalia t_fault [s] (default 12.0 s)
    double detection_threshold;  // Soglia normalizzata residuo FDIR (default 3.0 sigma)
    double nominal_bus_v;        // Tensione nominale bus primario non regolato [V] (default 28.0 V)
} FDIRSpec;

typedef struct {
    double t_sec;                // Tempo nel transitorio FDIR [s] (0 .. 120 s)
    double sys_health_pct;       // Indice di salute globale del sistema di bordo [%]
    double att_err_deg;          // Errore di puntamento d'assetto 3-assi [deg]
    double bus_voltage_v;        // Tensione del bus elettrico primario 28V [V]
    int    fdir_state_id;        // Stato macchina a stati OBDH FDIR (0..5)
    double primary_unit_pct;     // Livello di attività/coppia unità primaria A-Side [%]
    double redundant_unit_pct;   // Livello di subentro unità ridondante B-Side [%]
} FDIRSample;

typedef struct {
    FDIRSpec spec;
    char   fault_name[80];           // Nome descrittivo del guasto iniettato
    char   fault_subsystem[48];      // Sottosistema colpito (ADCS / EPS / GNC)
    char   primary_unit_name[48];    // Unità primaria guasta (A-Side)
    char   redundant_unit_name[48];  // Unità ridondante attivata (B-Side)
    double fault_inject_time_s;      // Istante di iniezione del guasto [s]
    double detection_latency_s;      // Latenza di rilevamento FDIR (Fault Detection Latency) [s]
    double isolation_switchover_s;   // Tempo di isolamento LCL + commutazione ramo B-Side [s]
    double safe_mode_duration_s;     // Durata acquisizione Sun-Pointing Safe Mode [s]
    double total_recovery_time_s;    // Tempo totale da iniezione a recupero nominale [s]
    double peak_att_err_deg;         // Escursione massima transitoria d'assetto durante il guasto [deg]
    double steady_post_err_deg;      // Errore di puntamento residuo dopo il recupero [deg]
    double nominal_bus_voltage_v;    // Tensione nominale bus elettrico [V] (28.0 V)
    double min_bus_voltage_v;        // Tensione minima raggiunta durante il transitorio [V]
    double min_sys_health_pct;       // Salute minima del sistema durante l'anomalia [%]
    char   redundancy_status[64];    // Stato della ridondanza hardware (B-Side Active / Hot-Cold Standby)
    char   recovery_outcome[80];     // Esito del recupero missione autonomo
    char   ccsds_oem_filename[48];   // Nome file CCSDS OEM generato ("enterprise_mission.oem")
    char   ccsds_oem_status[64];     // Stato conformità CCSDS 502.0-B-2 (OEM v2.0)
    int    oem_state_vectors;        // Numero di vettori di stato esportati nel file .oem
    int    tle_norad_id;             // ID catalogo NORAD di riferimento testato dal parser TLE
    char   tle_sat_name[48];         // Nome satellite TLE di riferimento

    int    num_samples;
    FDIRSample samples[MAX_FDIR_SAMPLES];
} FDIRReport;

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

// ==========================================
// 5M. SPRINT 24: INTERPLANETARY PATCHED CONICS, PORKCHOP PLOT & B-PLANE TARGETING
// ==========================================
#define SUN_MU            1.3271244e11          // Parametro gravitazionale eliocentrico mu_sun [km^3/s^2]
#define AU_KM             149597870.7           // Unita' Astronomica [km]
#define EARTH_ECLIPTIC_OBLIQUITY_DEG 23.4392911 // Obliquita' dell'eclittica (J2000) [deg]
#define MARS_MU           42828.375             // mu Marte [km^3/s^2]
#define MARS_RADIUS_KM    3389.5                // Raggio medio Marte [km]
#define MARS_SMA_AU       1.523679              // Semiasse maggiore Marte [AU]
#define VENUS_MU          324858.592            // mu Venere [km^3/s^2]
#define VENUS_RADIUS_KM   6051.8                // Raggio Venere [km]
#define VENUS_SMA_AU      0.723332              // Semiasse maggiore Venere [AU]
#define PORK_N_DEP        37                    // Punti griglia partenza: 0..180 d, passo 5 d
#define PORK_N_TOF        35                    // Punti griglia TOF: 150..320 d, passo 5 d
#define PORK_DEP_MAX_DAYS 180.0
#define PORK_TOF_MIN_DAYS 150.0
#define PORK_TOF_MAX_DAYS 320.0
#define MAX_INTERPLANETARY_SAMPLES 181

typedef enum {
    PLANET_MARS  = 0,
    PLANET_VENUS = 1
} InterplanetaryTarget;

typedef struct {
    InterplanetaryTarget target;      // Pianeta di destinazione
    double capture_alt_km;            // Altitudine del periasse di cattura (default 400 km)
    double tof_days;                  // TOF di progetto del trasferimento eliocentrico [days] (default 210)
} InterplanetarySpec;

typedef struct {
    double dep_day;                   // Offset di partenza [days]
    double tof_day;                   // Tempo di volo [days]
    double c3;                        // C3 = v_inf_dep^2 [km^2/s^2]
    double vinf_dep;                  // |v_inf| partenza [km/s]
    double vinf_arr;                  // |v_inf| arrivo [km/s]
    double dv_tmi_ms;                 // Delta-V TMI [m/s]
    double dv_moi_ms;                 // Delta-V MOI/VOI [m/s]
    double dv_tot_ms;                 // TMI + MOI [m/s]
    bool   valid;                     // Lambert convergente e fisicamente sensato
} PorkchopGridPoint;

typedef struct {
    double t_days;                    // Giorni dalla partenza
    double x_au, y_au, z_au;          // Posizione eliocentrica eclittica [AU]
} HelioSample;

typedef struct {
    InterplanetarySpec spec;
    char   target_name[16];           // "Mars" / "Venus"
    double mu_target;                 // [km^3/s^2]
    double r_target_km;               // raggio planetario [km]
    double target_sma_au;             // [AU]
    double soi_earth_km;              // Sfera d'influenza di Laplace Terra [km]
    double soi_target_km;             // Sfera d'influenza di Laplace pianeta target [km]
    double r_park_km;                 // Raggio del perigeo dell'orbita di parcheggio [km]
    double r_capture_km;              // Raggio del periasse di cattura [km]

    // Punto di progetto (TOF selezionato, miglior giorno di partenza = min DeltaV_tot)
    double opt_dep_day;               // [days]
    double opt_tof_days;              // [days]
    double c3_opt;                    // [km^2/s^2]
    double vinf_dep_kms;              // [km/s]
    double dla_deg, rla_deg;          // Declinazione / Ascensione retta dell'asintoto di fuga [deg]
    double dv_tmi_ms;                 // [m/s]
    double vinf_arr_kms;              // [km/s]
    double dv_moi_ms;                 // [m/s]
    double dv_total_ms;               // [m/s]
    double b_plane_km;                // Parametro d'impatto B [km]
    double hyp_e_dep, hyp_e_arr;      // Eccentricita' delle iperboli di partenza / arrivo
    double turn_angle_arr_deg;        // Angolo di deflessione all'arrivo [deg]
    double transfer_a_au, transfer_e; // Elementi dell'ellisse di trasferimento
    double transfer_angle_deg;        // Angolo di trasferimento eliocentrico [deg]
    double lambert_residual_km;       // Residuo posizione di arrivo (verifica propagazione Kepler)

    // Minimo globale della griglia
    double grid_min_c3, grid_min_c3_dep, grid_min_c3_tof;
    double grid_min_dv_ms, grid_min_dv_dep, grid_min_dv_tof;

    PorkchopGridPoint grid[PORK_N_TOF][PORK_N_DEP];
    int num_traj_samples;
    HelioSample traj[MAX_INTERPLANETARY_SAMPLES];
} InterplanetaryReport;

// ==========================================
// 5N. SPRINT 25: MULTI-SATELLITE FORMATION FLYING, J2-INVARIANT ROEs & FLEET RECONFIGURATION
// ==========================================
#define MAX_FORMATION_SAMPLES 145
#define NUM_FORMATION_DEPUTIES 3

typedef enum {
    FORM_PROJECTED_CIRCULAR_PCO = 0, // Projected Circular Orbit (PCO) nel piano T-N con oscillazione R
    FORM_ALONG_TRACK_STRING     = 1, // Stringa Along-Track (Leader-Follower con elica E/I passiva)
    FORM_PENDULUM_CROSS_TRACK   = 2, // Pendolo Cross-Track (oscillazione fuori piano per interferometria)
    FORM_CARTWHEEL_ELLIPSE      = 3  // Cartwheel 2:1 In-Plane Ellipse con separazione E/I sicura
} FormationGeometryType;

// Quasi-Nonsingular Relative Orbital Elements (ROEs di D'Amico) adimensionali e scalati per a_c [m]
typedef struct {
    double da;       // (a_d - a_c) / a_c [-]
    double dlambda;  // (u_d - u_c) + (Omega_d - Omega_c) * cos(i_c) [rad]
    double dex;      // e_d * cos(omega_d) - e_c * cos(omega_c) [-]
    double dey;      // e_d * sin(omega_d) - e_c * sin(omega_c) [-]
    double dix;      // i_d - i_c [rad]
    double diy;      // (Omega_d - Omega_c) * sin(i_c) [rad]
    // Componenti metriche scalate per il semiasse maggiore del Chief a_c [m]
    double da_m;
    double dlambda_m;
    double dex_m;
    double dey_m;
    double dix_m;
    double diy_m;
} RelativeOrbitalElements;

typedef struct {
    FormationGeometryType geometry;   // Topologia geometrica della formazione
    double baseline_radius_m;         // Raggio caratteristico della formazione rho_form [m] (default 500.0 m)
    double tolerance_window_m;        // Finestra di tolleranza della baseline [m] (default 5.0 m)
} FormationSpec;

typedef struct {
    double theta_deg;      // Anomalia orbitale lungo l'orbita del Chief [deg] (0 .. 360)
    double d1_r, d1_t, d1_n; // Coordinate LVLH [x_R, y_T, z_N] di Deputy 1 [m]
    double d2_r, d2_t, d2_n; // Coordinate LVLH [x_R, y_T, z_N] di Deputy 2 [m]
    double d3_r, d3_t, d3_n; // Coordinate LVLH [x_R, y_T, z_N] di Deputy 3 [m]
    double baseline_err_m; // Errore di deriva della baseline indotto da J2 [m] (entro +/-5 m)
} FormationSample;

typedef struct {
    FormationSpec spec;
    char   geometry_name[64];         // Nome descrittivo della topologia (es. "Projected Circular Orbit (PCO)")
    int    num_deputies;              // Numero di satelliti Deputy nel cluster (3)
    double chief_sma_km;              // Semiasse maggiore del Chief a_c [km]
    double chief_inc_deg;             // Inclinazione del Chief i_c [deg]
    double mean_motion_rad_s;         // Moto medio orbitale n_c [rad/s]

    // ROEs (D'Amico) per ciascuno dei 3 satelliti Deputy
    RelativeOrbitalElements roe_deputy[NUM_FORMATION_DEPUTIES];

    // Sicurezza passiva E/I-Vector Separation e metriche di separazione
    double min_safe_rn_sep_m;         // Separazione minima sicura nel piano R-N al passaggio per y_T = 0 [m]
    double min_3d_inter_sat_m;        // Separazione 3D minima inter-satellite lungo l'orbita [m]
    double ei_phase_angle_deg;        // Angolo di fase tra i vettori relativi delta_e e delta_i [deg] (0 = paralleli)
    bool   passive_safe;              // true se la separazione E/I garantisce assenza di collisione passiva

    // Deriva secolare differenziale da J2, Invarianza di Schaub-Alfriend & Costi Delta-V
    double j2_raan_drift_deg_day;     // Deriva secolare differenziale del nodo Delta_Omega_dot_J2 [deg/day]
    double j2_argp_drift_deg_day;     // Deriva secolare differenziale del perigeo Delta_omega_dot_J2 [deg/day]
    double schaub_alfriend_eta;       // Fattore cinematico eta = sqrt(1 - e_c^2) di Schaub-Alfriend [-]
    double reconfig_dv_ms;            // Costo di riconfigurazione / iniezione del cluster Delta-V_reconfig [m/s]
    double fk_dv_yr_ms;               // Costo annuale di Formation-Keeping contro J2 Delta-V_FK [m/s/yr]
    double fk_fuel_yr_kg;             // Consumo annuale di propellente per Formation-Keeping [kg/yr]
    double max_baseline_err_m;        // Massimo scostamento della baseline nella finestra di controllo [m]

    int    num_samples;
    FormationSample samples[MAX_FORMATION_SAMPLES];
} FormationFlyingReport;

// ==========================================
// 5O. SPRINT 26: AUTONOMOUS ONBOARD OPERATIONS SCHEDULER, SSR MEMORY BUFFER & ECSS/NASA PDR COMPLIANCE
// ==========================================
#define MAX_OPS_SAMPLES 160

typedef enum {
    OP_SCIENCE_IMAGING     = 0, // Acquisizione Payload Ottico/SAR (fuori eclissi, fuori manovra/downlink, SSR < 95%)
    OP_GROUND_DOWNLINK     = 1, // Scaricamento radio ad alta velocità sulla Ground Station (Milano-Bovisa)
    OP_ORBIT_MANEUVER_SLEW = 2, // Slew d'assetto ADCS e accensione propulsiva / riconfigurazione orbitale
    OP_SUN_BATTERY_CHARGE  = 3, // Puntamento solare per ricarica batterie Li-Ion o hold anti-overflow SSR (>= 95%)
    OP_ECLIPSE_STANDBY     = 4  // Attraversamento cono d'ombra terrestre (Payload standby, Survival Heaters / Batteria)
} OperationalMode;

typedef struct {
    double ssr_capacity_gb;        // Capacità massima memoria flash NAND SSR con EDAC [GB] (default 128.0 GB)
    double science_data_rate_mbps; // Rateo di generazione dati del Payload scientifico R_gen [Mbps] (default 120.0 Mbps)
    double ssr_overflow_limit_pct; // Soglia di protezione anti-overflow SSR [%] (default 95.0%)
} MissionOpsSpec;

typedef struct {
    double t_min;                  // Tempo cumulativo di missione [min]
    int    op_mode_id;             // Modalità operativa autonoma attiva (0..4, enum OperationalMode)
    double ssr_used_gb;            // Volume di dati istantaneo immagazzinato nel buffer SSR [GB]
    double ssr_fill_pct;           // Percentuale di riempimento del buffer SSR [%]
    double net_data_rate_mbps;     // Flusso netto di memoria dM_SSR/dt = R_gen(t) - R_downlink(t) [Mbps]
} OpsTimelineSample;

typedef struct {
    MissionOpsSpec spec;

    // Bilancio Memoria di Massa a Stato Solido (SSR - Solid-State Recorder)
    double peak_ssr_used_gb;          // Picco di occupazione memoria SSR registrato [GB]
    double peak_ssr_fill_pct;         // Picco percentuale di occupazione memoria SSR [%]
    double total_acquired_gb;         // Volume totale di dati scientifici acquisiti dal Payload [GB]
    double total_downlinked_gb;       // Volume totale di dati scaricati a terra via RF Downlink [GB]
    double net_backlog_gb;            // Backlog residuo in memoria SSR a fine finestra simulata [GB]
    int    overflow_prevented_events; // Numero di interventi automatici della protezione SSR >= 95%
    int    downlink_dump_events;      // Numero di sessioni di scaricamento dati a terra

    // Duty Cycle Operativo su 24 Ore / Finestra di Missione [%]
    double duty_science_pct;          // Percentuale tempo in Science Observation [%]
    double duty_downlink_pct;         // Percentuale tempo in RF Ground Downlink [%]
    double duty_sun_charge_pct;       // Percentuale tempo in Battery Charging / Sun Pointing [%]
    double duty_maneuver_pct;         // Percentuale tempo in GNC Maneuvers & Slew [%]
    double duty_eclipse_pct;          // Percentuale tempo in Eclipse Standby [%]

    // ECSS / NASA Preliminary Design Review (PDR) Compliance Matrix & Readiness Score
    double readiness_score_pct;       // Global Mission System Readiness Score [0-100%]
    char   status_propulsion[16];     // "PASS" / "WARN" (Delta-V, TCM 3-Sigma & Propellant Margin)
    double score_propulsion;
    char   status_eps_battery[16];    // "PASS" / "WARN" (Battery SoC > 20% & 28V Bus Stability)
    double score_eps_battery;
    char   status_thermal_rad[16];    // "PASS" / "WARN" (Bus Temp in [-10, +45] C & Radiation TID)
    double score_thermal_rad;
    char   status_rf_telecom[16];     // "PASS" / "WARN" (Eb/N0 Link Margin > 3.0 dB)
    double score_rf_telecom;
    char   status_ssr_buffer[16];     // "PASS" / "WARN" (SSR Fill < 95% & Downlink Throughput)
    double score_ssr_buffer;
    char   status_gnc_adcs_fdir[16];  // "PASS" / "WARN" (ADCS Pointing, EKF 3-Sigma, CSI Notch, Formation & FDIR)
    double score_gnc_adcs_fdir;
    char   status_ssa_eol[16];        // "PASS" / "WARN" (Post-CAM Pc < 1e-7 & ESA Zero-Debris EoL)
    double score_ssa_eol;
    char   pdr_compliance_verdict[80];// Esito globale revisione ECSS / NASA PDR

    int    num_samples;
    OpsTimelineSample samples[MAX_OPS_SAMPLES];
} MissionOpsReport;

// ==========================================
// 5P. SPRINT 28: SPACE WEATHER, SOLAR CYCLE (F10.7/Kp) & HIGH-FIDELITY ATMOSPHERIC DRAG (NRLMSISE-00 / JACCHIA)
// ==========================================
#define MAX_DECAY_SAMPLES   151
#define MAX_DENSITY_SAMPLES 91

typedef enum {
    SOLAR_CYCLE_MINIMUM = 0, // Solar Minimum (F10.7 = 70 SFU, Quiet Sun)
    SOLAR_CYCLE_AVERAGE = 1, // Moderate / Nominal Solar Activity (F10.7 = 140 SFU)
    SOLAR_CYCLE_MAXIMUM = 2  // Solar Maximum (F10.7 = 220 SFU, Active Sun)
} SolarCyclePhase;

typedef struct {
    SolarCyclePhase cycle_phase;  // Fase del ciclo solare selezionata (0=Min 70, 1=Avg 140, 2=Max 220)
    double f107_sfu;              // Flusso radio solare giornaliero a 10.7 cm (2800 MHz) [SFU = 10^-22 W/(m^2*Hz)]
    double f107_avg_sfu;          // Media centrata a 81 giorni (3 rotazioni di Carrington) di F10.7 [SFU]
    double kp_index;              // Indice planetario di tempesta geomagnetica Kp [0.0 .. 9.0] (default 3.0)
    double diurnal_bulge_lag_deg; // Angolo di ritardo del rigonfiamento diurno (Diurnal Bulge) ad Est del punto subsolare [deg] (~30 deg)
} SpaceWeatherSpec;

typedef struct {
    double t_years;               // Tempo trascorso dall'inizio del decadimento passivo [anni]
    double alt_km;                // Quota media orbitale sotto il profilo Space Weather selezionato [km]
    double solar_flux;            // Flusso solare F10.7(t) modulato dal ciclo undecennale di Schwabe [SFU]
    double alt_solar_min_km;      // Curva di riferimento in condizioni di Solar Minimum (F10.7 = 70 SFU) [km]
    double alt_solar_avg_km;      // Curva di riferimento in condizioni di Solar Average (F10.7 = 140 SFU) [km]
    double alt_solar_max_km;      // Curva di riferimento in condizioni di Solar Maximum (F10.7 = 220 SFU) [km]
} DecaySample;

typedef struct {
    double alt_km;                // Quota geometrica sul livello del mare [km] (120 .. 1000 km)
    double rho_min_kgm3;          // Densità termosferica in Solar Minimum (F10.7=70, Kp=1) [kg/m^3]
    double rho_nom_kgm3;          // Densità termosferica in Solar Nominal (F10.7=140, Kp=3) [kg/m^3]
    double rho_max_kgm3;          // Densità termosferica in Solar Maximum (F10.7=220, Kp=6) [kg/m^3]
    double rho_active_kgm3;       // Densità termosferica per la configurazione attiva selezionata [kg/m^3]
} DensityProfileSample;

typedef struct {
    SpaceWeatherSpec spec;
    char   phase_name[48];            // Nome regime solare (es. "Solar Maximum (F10.7 = 220 SFU)")
    double initial_alt_km;            // Quota iniziale di riferimento dell'orbita target [km]
    double ballistic_coeff_kgm2;      // Coefficiente balistico B = m / (C_D * A) [kg/m^2]
    double ap_index;                  // Ampiezza geomagnetica planetaria equivalente Ap [-]

    // Termosfera Jacchia-Bowman / Harris-Priester & Diurnal Bulge
    double t_night_min_k;             // Temperatura esosferica minima notturna globale T_c [K]
    double t_day_apex_k;              // Temperatura di picco sul Diurnal Bulge diurno (14:00 LST) T_max [K]
    double t_exospheric_k;            // Temperatura esosferica media orbitale T_inf(F10.7, Kp) [K]
    double diurnal_bulge_ratio;       // Rapporto di densità Diurnal Bulge (Day Apex / Night Min) [-]
    double rho_at_target_kgm3;        // Densità atmosferica locale alla quota target [kg/m^3]
    double scale_height_km;           // Altezza di scala termosferica effettiva H(h, T_inf) [km]
    double initial_decay_m_day;       // Rateo iniziale di decadimento del semiasse maggiore |da/dt| [m/day]
    double rev_decay_m;               // Decadimento per singola rivoluzione |Delta a_rev| [m/rev]

    // Vita Orbitale a Lungo Termine (King-Hele / Cowell Secular Decay) & Conformità Detriti Spaziali
    double estimated_lifetime_yr;     // Vita orbitale stimata sotto il profilo attivo (fino a h = 120 km) [anni]
    double lifetime_solar_min_yr;     // Vita orbitale in Solar Minimum (F10.7 = 70 SFU) [anni]
    double lifetime_solar_avg_yr;     // Vita orbitale in Solar Average (F10.7 = 140 SFU) [anni]
    double lifetime_solar_max_yr;     // Vita orbitale in Solar Maximum (F10.7 = 220 SFU) [anni]
    bool   rule_25yr_compliant;       // true se estimated_lifetime_yr <= 25.0 anni (IADC / ISO-24113)
    bool   fcc_5yr_compliant;         // true se estimated_lifetime_yr <= 5.0 anni (FCC / ESA Zero-Debris)
    char   compliance_verdict[80];    // Esito conformità IADC 25-Year / ESA Zero-Debris Rule

    int    num_decay_samples;
    DecaySample decay_samples[MAX_DECAY_SAMPLES];

    int    num_density_samples;
    DensityProfileSample density_samples[MAX_DENSITY_SAMPLES];
} OrbitLifetimeReport;

// ==========================================
// 5Q. SPRINT 29: OPTICAL INTER-SATELLITE LINKS (ISL), LINE-OF-SIGHT & CONSTELLATION MESH NETWORK
// ==========================================
#define MAX_ISL_EDGES 512

typedef struct {
    double max_range_km;          // Portata ottica massima del terminale laser LCT [km] (default 5000.0 km)
    double occ_alt_km;            // Quota atmosferica opaca di occultazione del fascio laser [km] (default 80.0 km)
    double laser_wavelength_nm;   // Lunghezza d'onda laser ottico C-Band / Nd:YAG [nm] (default 1550.0 nm)
    double link_rate_gbps;        // Capacità per singolo link laser coerente BPSK/Homodyne [Gbps] (default 10.0 Gbps)
} ISLNetworkSpec;

typedef struct {
    int    sat_id;                // Indice globale del satellite Walker [0 .. T-1]
    int    plane_id;              // Indice del piano orbitale [0 .. P-1]
    double x_km, y_km, z_km;      // Posizione cartesiana ECI 3D a t = 0 [km]
    double lat_deg, lon_deg;      // Punto sub-satellitare geodetico [deg]
    int    active_degree;         // Grado del nodo (numero di link laser attivi connessi)
} ISLNode3D;

typedef struct {
    int    sat_a;                 // Indice nodo sorgente i
    int    sat_b;                 // Indice nodo destinazione j
    int    is_intra_plane;        // 1 = link intra-piano (anello orbitale), 0 = link inter-piano (cross-plane)
    double distance_km;           // Distanza geometrica euclidea d_ij = ||r_j - r_i|| [km]
    double latency_ms;            // Latenza di propagazione ottica nel vuoto tau = d_ij / c [ms]
    double min_ray_alt_km;        // Quota minima del raggio laser sopra la superficie terrestre r_min - R_E [km]
    double fspl_db;               // Attenuazione ottica in spazio libero FSPL (lambda = 1550 nm) [dB]
    double link_margin_db;        // Margine di potenza del collegamento laser a 10 Gbps [dB]
} ISLEdge;

typedef struct {
    ISLNetworkSpec spec;
    int    num_nodes;                     // Numero totale di satelliti nella costellazione Walker T
    ISLNode3D nodes[MAX_WALKER_SATS];

    int    num_active_links;              // Numero totale di link laser attivi (in LoS e entro portata)
    int    intra_plane_links;             // Link attivi lungo lo stesso piano orbitale (Forward/Aft)
    int    inter_plane_links;             // Link attivi tra piani orbitali adiacenti (Cross-Plane)
    double mean_distance_km;              // Distanza media dei link laser attivi [km]
    double min_distance_km;               // Distanza minima tra satelliti connessi [km]
    double max_distance_km;               // Distanza massima tra satelliti connessi [km]
    double horizon_limit_km;              // Distanza geometrica massima di orizzonte non occultata a h_target [km]
    double mean_latency_ms;               // Latenza ottica media di singolo salto (single-hop) [ms]
    double max_latency_ms;                // Latenza ottica massima di singolo salto [ms]
    double mean_margin_db;                // Margine medio del link budget ottico a 10 Gbps [dB]
    double min_margin_db;                 // Margine minimo del link budget ottico a 10 Gbps [dB]
    double total_capacity_gbps;           // Capacità aggregata della rete Mesh spaziale [Gbps]
    double mean_node_degree;              // Grado medio di connessione per satellite (LCT terminals/sat)
    char   mesh_connectivity_status[64];  // Stato topologico della rete ("FULL SPHERICAL MESH" / "PARTIAL MESH")

    ISLEdge edges[MAX_ISL_EDGES];
} ISLNetworkReport;

// ==========================================
// 5R. SPRINT 30: MULTI-GROUND STATION TRACKING, DEEP SPACE NETWORK (DSN) GEOMETRY & LINK HANDOVERS
// ==========================================
#define NUM_NET_STATIONS         4
#define MAX_GROUND_PASSES        64
#define MAX_GROUND_TRACK_SAMPLES 181

typedef enum {
    NET_NASA_DSN    = 0, // NASA Deep Space Network (Goldstone, Madrid, Canberra + Milano-Bovisa MOC)
    NET_ESA_ESTRACK = 1, // ESA ESTRACK Global Network (Cebreros, Malargue, New Norcia, Kourou)
    NET_POLAR_LEO   = 2  // High-Latitude Polar / LEO Network (Svalbard SvalSat, TrollSat, Fairbanks, Punta Arenas)
} GroundNetworkType;

typedef struct {
    GroundNetworkType network_type; // Rete globale di stazioni di terra selezionata (0=DSN, 1=ESTRACK, 2=Polar LEO)
    double elevation_mask_deg;      // Maschera di elevazione minima sopra l'orizzonte locale [deg] (default 5.0 deg)
} GroundNetworkSpec;

typedef struct {
    int    station_id;              // Indice stazione nella rete [0 .. 3]
    char   name[36];                // Nome della stazione (es. "Goldstone (DSN-14, USA)")
    double lat_deg;                 // Latitudine geodetica WGS84 [deg]
    double lon_deg;                 // Longitudine geodetica WGS84 [deg]
    double alt_km;                  // Quota sul livello del mare [km]
    double min_el_deg;              // Maschera di elevazione minima [deg]
    double contact_time_min_24h;    // Tempo totale di contatto acquisito dalla stazione su 24h [min]
    int    pass_count_24h;          // Numero di passaggi (AOS -> LOS) registrati su 24h
    double max_el_deg;              // Massima elevazione raggiunta sulla stazione [deg]
    double peak_slew_rate_deg_s;    // Velocità angolare massima richiesta ai motori del piedistallo antenna [deg/s]
    double footprint_radius_deg;    // Raggio angolare del cono di visibilità terrestre Lambda [deg]
} GroundStationNode;

typedef struct {
    int    station_id;              // Indice stazione [0 .. 3]
    char   station_name[36];        // Nome della stazione
    double aos_min;                 // Istante di Acquisition of Signal (AOS) [min]
    double los_min;                 // Istante di Loss of Signal (LOS) [min]
    double duration_min;            // Durata della finestra di visibilità [min]
    double max_el_deg;              // Elevazione di culmine durante il passaggio [deg]
    double az_at_aos_deg;           // Azimuth topocentrico all'AOS [deg]
    double az_at_los_deg;           // Azimuth topocentrico al LOS [deg]
    double peak_slew_deg_s;         // Picco di Slew Rate meccanico durante il passaggio [deg/s]
} GroundStationPass;

typedef struct {
    double t_min;                   // Tempo trascorso nella finestra di tracking [min]
    double sat_lat_deg;             // Latitudine sub-satellitare WGS84 [deg]
    double sat_lon_deg;             // Longitudine sub-satellitare WGS84 [deg]
    double sat_alt_km;              // Quota geodetica del satellite [km]
    int    active_station_id;       // ID stazione attiva con massima elevazione (-1 = Comm Gap)
    double azimuth_deg;             // Azimuth topocentrico SEZ della stazione attiva [0 .. 360 deg]
    double elevation_deg;           // Elevazione topocentrica SEZ della stazione attiva [-90 .. +90 deg]
    double range_km;                // Slant Range topocentrico rho [km]
    double range_rate_kms;          // Range-rate radiale d(rho)/dt [km/s]
    double slew_rate_deg_s;         // Velocità angolare istantanea dei motori dell'antenna [deg/s]
} GroundTrackSample;

typedef struct {
    GroundNetworkSpec spec;
    char   network_name[48];              // Nome descrittivo della rete attiva
    int    num_stations;                  // Numero di stazioni nella rete (4)
    GroundStationNode stations[NUM_NET_STATIONS];

    double total_contact_hours_24h;       // Tempo totale di contatto combinato (senza sovrapposizioni) su 24h [ore]
    double coverage_duty_pct;             // Percentuale del giorno con link di terra attivo [%]
    double max_comm_gap_min;              // Massimo intervallo di blackout senza comunicazioni (Max Comm Gap) [min]
    int    total_handovers;               // Numero totale di Handover automatici tra stazioni
    int    best_station_id;               // Indice della stazione con il maggior tempo di contatto
    char   best_station_name[36];         // Nome della stazione con il maggior tempo di contatto
    double best_station_contact_min;      // Tempo di contatto della migliore stazione [min]
    double peak_antenna_slew_deg_s;       // Massima velocità angolare richiesta ai motori d'antenna [deg/s]

    int    num_passes;
    GroundStationPass passes[MAX_GROUND_PASSES];

    int    num_track_samples;
    GroundTrackSample samples[MAX_GROUND_TRACK_SAMPLES];
} GroundSegmentReport;

// ==========================================
// 5S. SPRINT 31: HIGH-FIDELITY ELECTRICAL POWER SYSTEM (EPS), PENUMBRA/UMBRA ECLIPSES & LI-ION DEGRADATION
// ==========================================
#define MAX_EPS_SAMPLES    181
#define SUN_RADIUS_KM      695700.0 // Raggio fotosferico solare [km]

typedef enum {
    BATT_CONFIG_8S4P  = 0, // 8S4P Li-Ion (28V Nominal Bus, 33.6V Max / 24.0V Cutoff, 48 Ah / 1382 Wh)
    BATT_CONFIG_8S6P  = 1, // 8S6P Li-Ion High-Capacity (28V Nominal Bus, 33.6V Max / 24.0V Cutoff, 72 Ah / 2074 Wh)
    BATT_CONFIG_7S4P  = 2, // 7S4P Li-Ion Compact LEO (25.2V Nominal Bus, 29.4V Max / 21.0V Cutoff, 48 Ah / 1210 Wh)
    BATT_CONFIG_14S4P = 3  // 14S4P Li-Ion High-Power (50V Nominal Bus, 58.8V Max / 42.0V Cutoff, 48 Ah / 2419 Wh)
} EPSBatteryConfigType;

typedef enum {
    CELL_GAAS_TRIPLE_JUNCTION = 0, // GaAs Triple-Junction InGaP/GaAs/Ge (eta_0 = 30.2%, beta = -0.20%/K, radiation hard)
    CELL_SILICON_HIGH_EFF     = 1  // High-Efficiency Space Silicon BSFR (eta_0 = 19.8%, beta = -0.45%/K)
} EPSSolarCellType;

typedef struct {
    EPSBatteryConfigType batt_config; // Configurazione pacco batteria Li-Ion (default BATT_CONFIG_8S4P)
    EPSSolarCellType     cell_type;   // Tecnologia celle fotovoltaiche (default CELL_GAAS_TRIPLE_JUNCTION)
    double solar_array_area_m2;       // Area efficace pannelli solari A_sa [m^2] (default 2.85 m^2)
    double mission_design_life_yr;    // Vita operativa di progetto per calcolo SoH EoL [anni] (default 5.0 anni)
} EPSHiFiSpec;

typedef struct {
    double t_min;             // Tempo trascorso nella finestra di missione [min]
    double eclipse_fraction;  // Frazione d'ombra conica terrestre: 0.0 = Luce piena, (0,1) = Penombra, 1.0 = Ombra (Umbra)
    double solar_pwr_w;       // Potenza solare generata dopo MPPT, efficienza termica eta(T) e degradazione TID [W]
    double load_pwr_w;        // Potenza elettrica richiesta dai sottosistemi di bordo (Bus + ADCS + TCS Heaters + RF/Payload) [W]
    double bus_voltage_v;     // Tensione non lineare del bus batteria V_bus(SoC, I_batt) [V]
    double soc_pct;           // State of Charge istantaneo della batteria Li-Ion [%]
    double soh_pct;           // State of Health (capacità residua C_max / C_nom) lungo la vita operativa [%]
} EPSSample;

typedef struct {
    EPSHiFiSpec spec;
    char   batt_config_name[48];      // Nome configurazione batteria (es. "8S4P Li-Ion (28V Unregulated Bus)")
    char   cell_type_name[48];        // Nome tecnologia solare (es. "GaAs Triple-Junction (InGaP/GaAs/Ge)")
    int    series_cells;              // Numero di celle in serie N_s (es. 8 per 28V bus)
    int    parallel_strings;          // Numero di stringhe in parallelo N_p (es. 4)
    double nominal_bus_v;             // Tensione nominale di targa del bus [V]
    double max_bus_voltage_v;         // Tensione massima del bus a piena carica (N_s * 4.2V) [V]
    double min_bus_voltage_v;         // Tensione minima del bus registrata durante la scarica in ombra [V]
    double nominal_cap_wh;            // Capacità energetica nominale a inizio vita (BoL) [Wh]
    double eol_cap_wh;                // Capacità energetica residua a fine vita (EoL) [Wh]
    double mppt_eff_pct;              // Efficienza media del convertitore MPPT (Maximum Power Point Tracking) [%]
    double solar_eff_mean_pct;        // Efficienza fotovoltaica media corretta per temperatura eta_T [%]
    double tid_deg_pct;               // Perdita di potenza dei pannelli solari per dose radiativa TID [%]
    double peak_solar_pwr_w;          // Picco di potenza solare generata (Cold-Array Post-Eclipse) [W]
    double mean_solar_pwr_w;          // Potenza solare media orbitale [W]
    double mean_load_pwr_w;           // Carico elettrico medio richiesto dai sottosistemi [W]
    double peak_dod_pct;              // Massima Depth of Discharge (DoD) registrata per ciclo [%]
    double min_soc_pct;               // Minimo State of Charge (SoC = 100 - Peak DoD) registrato [%]
    double eol_soh_pct;               // State of Health (SoH) della batteria a fine vita (EoL) [%]
    double penumbra_time_min;         // Durata complessiva dei transiti in Penombra nella finestra [min]
    double umbra_time_min;            // Durata complessiva dei transiti in Ombra totale (Umbra) nella finestra [min]
    char   eps_compliance_status[64]; // Verdetto conformità ECSS-E-ST-20C (DoD < 35%, SoH_EoL > 75%, V_bus > V_cut)

    int    num_samples;
    EPSSample samples[MAX_EPS_SAMPLES];
} EPSHiFiReport;

// ==========================================
// 5T. SPRINT 32: EARTH OBSERVATION (EO) RAYCASTING, SYNTHETIC PAYLOAD & COVERAGE HEATMAP
// ==========================================
#define EO_HEATMAP_LAT_BINS      60   // Griglia globale latitudine: 60 celle da 3° (-90° .. +90°)
#define EO_HEATMAP_LON_BINS      120  // Griglia globale longitudine: 120 celle da 3° (-180° .. +180°)
#define MAX_EO_FOOTPRINT_SAMPLES 121  // Campioni temporali del poligono trapezoidale di scansione

typedef struct {
    double fov_x_deg;          // Field of View Along-Track (asse X sensore) [deg] (default 10.0 deg)
    double fov_y_deg;          // Field of View Cross-Track (asse Y sensore) [deg] (default 15.0 deg)
    double off_nadir_deg;      // Angolo di puntamento laterale Boresight Off-Nadir [deg] (default 0.0 deg)
    double opt_aperture_m;     // Diametro apertura ottica per calcolo GSD [m] (default 0.35 m)
    double wavelength_nm;      // Lunghezza d'onda spettrale [nm] (default 550.0 nm)
} EOSensorSpec;

typedef struct {
    double t_min;              // Tempo trascorso nella finestra di osservazione [min]
    double sat_lat_deg;        // Latitudine geodetica sub-satellitare WGS84 [deg]
    double sat_lon_deg;        // Longitudine geodetica sub-satellitare WGS84 [deg]
    double sat_alt_km;         // Quota geodetica del satellite sopra l'ellissoide WGS84 [km]
    double bore_lat_deg;       // Latitudine del punto d'intersezione Boresight sull'ellissoide WGS84 [deg]
    double bore_lon_deg;       // Longitudine del punto d'intersezione Boresight sull'ellissoide WGS84 [deg]
    double slant_range_km;     // Distanza obliqua (Slant Range) lungo il raggio di Boresight [km]
    double incidence_deg;      // Angolo d'incidenza locale sul terreno (inclusa curvatura WGS84) [deg]
    double gsd_m;              // Risoluzione spaziale effettiva a terra (Off-Nadir Smearing GSD) [m]
    double smearing_factor;    // Fattore di allungamento geometrico del pixel (Pixel Smearing) [-]
    double corner_lat_deg[4];  // Latitudini WGS84 dei 4 vertici del poligono del FOV (Swath Polygon) [deg]
    double corner_lon_deg[4];  // Longitudini WGS84 dei 4 vertici del poligono del FOV (Swath Polygon) [deg]
} EOSensorFootprintSample;

typedef struct {
    EOSensorSpec spec;
    char   sensor_mode_name[64];         // Modalità operativa ("Nadir Pushbroom Optical" / "Off-Nadir Side-Looking SAR/EO")
    double global_coverage_pct;          // Percentuale della superficie terrestre scansionata (pesata per cos(lat)) [%]
    double avg_slant_range_km;           // Slant Range medio lungo l'asse ottico di Boresight [km]
    double avg_incidence_deg;            // Angolo d'incidenza locale medio sul terreno [deg]
    double nadir_gsd_m;                  // GSD di riferimento al nadir [m]
    double mean_gsd_m;                   // GSD medio effettivo sull'intero swath (con Pixel Smearing) [m]
    double off_nadir_smearing_factor;    // Fattore medio di distorsione/smearing geometrico fuori dal nadir [-]
    double swath_area_km2;               // Area istantanea del footprint trapezoidale proiettato a terra [km^2]
    int    unique_cells_hit;             // Numero di celle della griglia globale osservate almeno una volta
    int    total_grid_cells;             // Numero totale di celle della griglia globale (60 x 120 = 7200)
    int    max_revisit_hits;             // Massimo conteggio di rivisita registrato su una cella

    int    heatmap_counts[EO_HEATMAP_LAT_BINS][EO_HEATMAP_LON_BINS]; // Matrice 2D conteggi di osservazione lat/lon

    int    num_samples;
    EOSensorFootprintSample samples[MAX_EO_FOOTPRINT_SAMPLES];
} EORaycastReport;

// ==========================================
// 5U. SPRINT 33: INTERPLANETARY GRAVITY ASSIST (PLANETARY FLYBY), V-INFINITY ROTATION & B-PLANE TARGETING
// ==========================================
#define JUPITER_MU          126686534.0 // Parametro gravitazionale di Giove [km^3/s^2]
#define JUPITER_RADIUS_KM   71492.0     // Raggio equatoriale di Giove [km]
#define JUPITER_SMA_AU      5.2044      // Semiasse maggiore eliocentrico di Giove [AU]
#define MAX_FLYBY_SAMPLES   101         // Campioni 3D dell'iperbole planetocentrica di flyby

typedef enum {
    FLYBY_JUPITER = 0, // Slingshot gravitazionale su Giove (Voyager / Cassini / New Horizons class)
    FLYBY_VENUS   = 1, // Gravity assist su Venere (Inner Solar System / Juice / BepiColombo class)
    FLYBY_MARS    = 2  // Gravity assist su Marte (Rosetta / Dawn Deep-Space class)
} FlybyPlanetType;

typedef struct {
    bool            enable_flyby;   // Abilita/disabilita la manovra di Gravity Assist (Flyby)
    FlybyPlanetType flyby_planet;   // Pianeta intermedio per l'effetto fionda (0=Jupiter, 1=Venus, 2=Mars)
    double          flyby_alt_km;   // Altitudine di pericentro del flyby sopra la superficie/nubi h_p [km]
} GravityAssistSpec;

typedef struct {
    double t_hr;     // Tempo relativo al passaggio al pericentro t - t_p [ore]
    double nu_deg;   // Anomalia vera lungo l'iperbole di flyby [deg]
    double x_km;     // Coordinata X planetocentrica 3D [km]
    double y_km;     // Coordinata Y planetocentrica 3D [km]
    double z_km;     // Coordinata Z planetocentrica 3D [km]
    double r_km;     // Distanza radiale dal centro del pianeta [km]
    double v_kms;    // Velocità iperbolica planetocentrica istantanea [km/s]
} FlybyHyperbolaSample;

typedef struct {
    GravityAssistSpec spec;
    char     flyby_planet_name[32];    // Nome del pianeta di flyby ("Jupiter", "Venus", "Mars")
    double   mu_planet;                // Parametro gravitazionale del pianeta mu_p [km^3/s^2]
    double   r_planet_km;              // Raggio medio/equatoriale del pianeta R_p [km]
    double   planet_sma_au;            // Distanza eliocentrica media del pianeta [AU]
    double   soi_planet_km;            // Raggio della Sphere of Influence (SOI) di Laplace [km]
    double   periapsis_radius_km;      // Raggio di pericentro del flyby r_p = R_p + h_flyby [km]

    Vector3D v_planet_vec;             // Velocità orbitale eliocentrica del pianeta V_planet [km/s]
    double   v_planet_kms;             // Modulo velocità eliocentrica del pianeta ||V_planet|| [km/s]
    Vector3D v_sc_in_vec;              // Velocità eliocentrica in ingresso della sonda V_sc,in [km/s]
    double   v_sc_in_kms;              // Modulo velocità eliocentrica in ingresso ||V_sc,in|| [km/s]
    Vector3D vinf_in_vec;              // Eccesso iperbolico in ingresso v_inf,in = V_sc,in - V_planet [km/s]
    double   vinf_kms;                 // Modulo dell'eccesso iperbolico ||v_inf|| [km/s]

    double   deflection_angle_deg;     // Angolo di deflessione (Turning Angle) delta = 2*arcsin(1/e_hyp) [deg]
    double   hyp_ecc;                  // Eccentricità dell'iperbole di flyby e_hyp = 1 + r_p*v_inf^2/mu_p [-]
    double   hyp_sma_km;               // Semiasse maggiore iperbolico a_hyp = -mu_p/v_inf^2 [km]
    double   b_plane_km;               // Parametro d'impatto nel B-Plane B = r_p*sqrt(1 + 2*mu_p/(r_p*v_inf^2)) [km]
    double   b_theta_deg;              // Angolo di clock nel B-Plane theta_B [deg]
    Vector3D b_vec_km;                 // Vettore 3D del parametro d'impatto B nel piano fondamentale [km]
    Vector3D periapsis_vec_km;         // Vettore 3D di pericentro r_p nel riferimento planetocentrico [km]
    double   v_periapsis_kms;          // Velocità massima al pericentro del flyby v_p [km/s]

    Vector3D vinf_out_vec;             // Eccesso iperbolico in uscita ruotato di delta v_inf,out [km/s]
    Vector3D v_sc_out_vec;             // Velocità eliocentrica in uscita V_sc,out = V_planet + v_inf,out [km/s]
    double   v_sc_out_kms;             // Modulo velocità eliocentrica in uscita ||V_sc,out|| [km/s]
    double   free_dv_kms;              // Delta-V vettoriale "gratuito" fornito dal pianeta ||V_sc,out - V_sc,in|| [km/s]
    double   free_dv_ms;               // Delta-V vettoriale "gratuito" fornito dal pianeta [m/s]
    double   helio_speed_gain_kms;     // Guadagno netto di velocità scalare eliocentrica ||V_sc,out|| - ||V_sc,in|| [km/s]
    double   post_flyby_sma_au;        // Semiasse maggiore eliocentrico post-flyby [AU] (<0 se fuga dal Sistema Solare)
    double   post_flyby_ecc;           // Eccentricità eliocentrica post-flyby [-]
    double   prop_fuel_saved_kg;       // Massa equivalente di propellente risparmiata a costo zero [kg]

    int      num_hyperbola_samples;
    FlybyHyperbolaSample hyperbola[MAX_FLYBY_SAMPLES];
} GravityAssistReport;

// ==========================================
// 5V. SPRINT 34: AEROCAPTURE, AEROBRAKING & MARTIAN ATMOSPHERIC PASS
// ==========================================
#define MAX_AEROCAPTURE_SAMPLES 141 // Campioni temporali del passaggio atmosferico 3-DOF di aerocattura

typedef enum {
    INSERTION_PROPULSIVE_BURN = 0, // Cattura orbitale propulsiva classica (Retro-Burn impulsivo dV_MOI)
    INSERTION_AEROCAPTURE     = 1  // Cattura atmosferica (Aerocapture + Periapsis Raise Burn + Aerobraking)
} OrbitInsertionMethod;

typedef struct {
    OrbitInsertionMethod method;          // Metodo di inserzione in orbita target (0=Propulsive, 1=Aerocapture)
    InterplanetaryTarget planet;          // Pianeta dotato di atmosfera per l'aerocattura (0=Mars, 1=Venus)
    double lift_to_drag;                  // Efficienza aerodinamica ipersonica L/D = C_L / C_D (es. 0.25)
    double entry_fpa_deg;                 // Flight-Path Angle all'interfaccia di ingresso h_EI = 125 km [deg] (es. -11.4 deg)
    double ballistic_coeff_kgm2;          // Coefficiente balistico beta = m / (C_D * A) [kg/m^2] (es. 145 kg/m^2)
    double nose_radius_m;                 // Raggio di curvatura dello scudo termico TPS R_n [m] (es. 1.25 m)
    double target_apoapsis_km;            // Quota di apocentro target dell'orbita ellittica post-cattura [km] (es. 1200.0 km)
    double final_periapsis_km;            // Quota di pericentro finale dopo il Periapsis Raise Burn [km] (es. 400.0 km)
} AerocaptureSpec;

typedef struct {
    double t_sec;                  // Tempo trascorso dall'Entry Interface (h_EI = 125 km) [s]
    double alt_km;                 // Quota sopra la superficie planetaria h(t) [km]
    double vel_kms;                // Velocità planetocentrica relativa V(t) [km/s]
    double esc_vel_kms;            // Velocità di fuga locale sqrt(2*mu/r) [km/s]
    double fpa_deg;                // Angolo di rampa (Flight-Path Angle) gamma(t) [deg]
    double bank_angle_deg;         // Angolo di Bank (Rollio) comandato dal guida GNC sigma(t) [deg]
    double rho_kgm3;               // Densità atmosferica locale rho(h) [kg/m^3]
    double heat_flux_wcm2;         // Flusso termico aerotermodinamico totale (Sutton-Graves CO2 + rad) [W/cm^2]
    double decel_g;                // Decelerazione aerodinamica totale in multipli di g0 [g]
    double spec_energy_km2s2;      // Energia orbitale specifica kepleriana E = V^2/2 - mu/r [km^2/s^2]
    double inst_sma_km;            // Semiasse maggiore istantaneo a(t) = -mu / (2*E) [km]
    double inst_ecc;               // Eccentricità orbitale istantanea e(t) [-]
    double inst_apoapsis_km;       // Quota di apocentro istantanea post-cattura (se E < 0) [km]
    double downrange_km;           // Distanza curvilinea percorsa lungo l'arco atmosferico [km]
} AerocaptureSample;

typedef struct {
    AerocaptureSpec spec;
    char   method_name[48];               // "Hypersonic Aerocapture + Periapsis Raise" / "Propulsive MOI Retro-Burn"
    char   planet_name[24];               // "Mars" o "Venus"
    double entry_interface_alt_km;        // Quota dell'Entry Interface atmosferica [km] (125.0 km su Marte)
    double vinf_arr_kms;                  // Eccesso iperbolico di arrivo v_inf,arr [km/s]
    double entry_vel_kms;                 // Velocità all'ingresso atmosferico V_EI [km/s]
    double exit_vel_kms;                  // Velocità all'uscita dall'atmosfera V_exit [km/s]
    double min_altitude_km;               // Altitudine minima raggiunta durante il passaggio atmosferico [km]
    double time_at_min_alt_s;             // Istante del passaggio al pericentro atmosferico [s]
    double pass_duration_s;               // Durata totale dell'attraversamento atmosferico [s]
    double capture_time_s;                // Istante in cui l'energia diventa negativa (E < 0, e < 1) [s]

    double corridor_skip_fpa_deg;         // Limite superiore del corridoio (Skip-Out boundary) [deg]
    double corridor_burn_fpa_deg;         // Limite inferiore del corridoio (Thermal/Impact undershoot) [deg]
    double corridor_width_deg;            // Ampiezza del Corridoio di Aerocattura [deg]
    double corridor_margin_pct;           // Margine di sicurezza dentro il corridoio di volo [%]

    double peak_heat_flux_wcm2;           // Picco di flusso termico sullo scudo TPS [W/cm^2]
    double total_heat_load_jcm2;          // Carico termico integrato sul pass [J/cm^2]
    double peak_decel_g;                  // Massima decelerazione aerodinamica avvertita [g]
    double peak_stag_pressure_kpa;        // Pressione dinamica massima di ristagno q_bar [kPa]

    double initial_hyp_ecc;               // Eccentricità iperbolica di arrivo e_in (> 1)
    double exit_eccentricity;             // Eccentricità dell'orbita ellittica subito all'uscita atmosferica
    double final_eccentricity;            // Eccentricità finale raggiunta dopo innalzamento pericentro e aerobraking
    double exit_sma_km;                   // Semiasse maggiore dell'orbita catturata [km]
    double captured_apoapsis_km;          // Quota di apocentro raggiunta post-passaggio [km]
    double final_periapsis_km;            // Quota di pericentro finale sicura fuori dall'atmosfera [km]

    double drag_dv_bled_kms;              // Velocità cinetica dissipata per attrito atmosferico int(D/m dt) [km/s]
    double propulsive_moi_dv_kms;         // Costo equivalente del burn propulsivo classico dV_MOI [km/s]
    double periapsis_raise_dv_ms;         // Piccolo burn all'apocentro per alzare il pericentro fuori dall'atmosfera [m/s]
    double eq_dv_saved_kms;               // Delta-V propulsivo equivalente risparmiato grazie all'aerocattura [km/s]
    double eq_dv_saved_ms;                // Delta-V propulsivo equivalente risparmiato [m/s]
    double fuel_saved_spacecraft_kg;      // Carburante risparmiato sul veicolo (500 kg dry/wet class) [kg]
    double fuel_saved_heavy_stack_kg;     // Carburante risparmiato su uno stadio pesante marziano (3500 kg class) [kg]
    char   corridor_status[64];           // Verdetto ("NOMINAL AEROCAPTURE — CLOSED ELLIPTICAL ORBIT LOCKED")

    int    num_samples;
    AerocaptureSample samples[MAX_AEROCAPTURE_SAMPLES];
} AerocaptureReport;

// ==========================================
// 5W. SPRINT 35: SOLAR SAIL CONTINUOUS THRUST & SOLAR RADIATION PRESSURE (SRP) DYNAMICS
// ==========================================
#define MAX_SOLAR_SAIL_SAMPLES 161      // Campioni temporali della spirale eliocentrica 3D a vela solare
#define SRP_1AU_N_M2           4.5398e-6 // Pressione di radiazione solare nominale P_sun/c a 1 AU [N/m^2]

typedef struct {
    bool                 enable_sail;      // Abilita/disabilita la propulsione continua a Vela Solare
    InterplanetaryTarget target_planet;    // Pianeta bersaglio della spirale (0=Mars Outward, 1=Venus Inward)
    double               sail_area_m2;     // Area riflettente della vela solare A_sail [m^2] (default 1200.0 m^2)
    double               sailcraft_mass_kg;// Massa totale del veliero solare (bus + boom + membrana CP1) [kg] (default 160.0 kg)
    double               reflectivity;     // Coefficiente di riflettanza speculare della vela eta_refl [-] (default 0.91)
} SolarSailSpec;

typedef struct {
    double t_days;         // Tempo di volo trascorso lungo la spirale eliocentrica [giorni]
    double x_au;           // Coordinata X eliocentrica eclittica [AU]
    double y_au;           // Coordinata Y eliocentrica eclittica [AU]
    double z_au;           // Coordinata Z eliocentrica fuori piano [AU]
    double r_au;           // Distanza radiale dal Sole r(t) [AU]
    double vel_kms;        // Velocità orbitale eliocentrica istantanea [km/s]
    double alpha_deg;      // Angolo di beccheggio (Pitch / Cone Angle) della normale vela rispetto al Sole [deg]
    double psi_deg;        // Angolo di Clock fuori piano comandato dall'ADCS [deg]
    double a_srp_mms2;     // Accelerazione istantanea da Pressione di Radiazione Solare ||a_SRP|| [mm/s^2]
    double thrust_mn;      // Forza di spinta fotonica istantanea [mN]
    double cum_dv_kms;     // Delta-V cumulativo accumulato gratuitamente dalla luce solare [km/s]
    double nx, ny, nz;     // Componenti del versore normale alla vela n_hat nel frame eclittico [-]
} SolarSailSample;

typedef struct {
    SolarSailSpec spec;
    char   target_planet_name[24];         // "Mars" (Outward Spiral) o "Venus" (Inward Spiral)
    char   steering_law_name[64];          // "McInnes Optimal Tangential Gain (alpha* = +35.26 deg)"
    double sail_loading_gm2;               // Carico areale sigma = m / A_sail [g/m^2]
    double lightness_number_beta;          // Numero di leggerezza adimensionale beta = a_c / g_sun(1 AU) [-]
    double char_accel_mms2;                // Accelerazione caratteristica a 1 AU a_c (per alpha = 0) [mm/s^2]
    double eff_accel_1au_mms2;             // Accelerazione efficace a 1 AU all'angolo ottimo cos^2(alpha*) [mm/s^2]
    double peak_srp_thrust_mn;             // Spinta fotonica massima a 1 AU [mN]
    double opt_pitch_angle_deg;            // Angolo di beccheggio ottimo medio alpha* = +/-35.26 deg [deg]
    double max_pitch_excursion_deg;        // Massima escursione di beccheggio ADCS durante la modulazione [deg]
    double initial_r_au;                   // Raggio eliocentrico iniziale (1.000 AU) [AU]
    double final_r_au;                     // Raggio eliocentrico finale raggiunto [AU]
    double target_r_au;                    // Raggio eliocentrico del pianeta target [AU]
    double spiral_tof_days;                // Durata totale del trasferimento a spirale [giorni]
    double spiral_revolutions;             // Numero di rivoluzioni eliocentriche completate [-]
    double total_dv_kms;                   // Delta-V totale accumulato dalla pressione di radiazione solare [km/s]
    double total_dv_ms;                    // Delta-V totale accumulato dalla luce solare [m/s]
    double eq_chemical_fuel_saved_kg;      // Massa di propellente chimico (MMH/NTO) equivalente risparmiata [kg]

    int    num_samples;
    SolarSailSample samples[MAX_SOLAR_SAIL_SAMPLES];
} SolarSailReport;

#endif // ORBIT_TYPES_H
