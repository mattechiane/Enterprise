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
// 5I. SPRINT 20: END-OF-LIFE (EoL) DEORBIT, ATMOSPHERIC RE-ENTRY & SUTTON-GRAVES HEATING
// ==========================================
#define MAX_REENTRY_SAMPLES 160

typedef struct {
    double t_sec;               // Tempo trascorso dall'Entry Interface (h_EI = 120 km) [s]
    double alt_km;              // Quota sopra l'ellissoide terrestre h(t) [km]
    double vel_kms;             // Velocità relativa atmosferica V(t) [km/s]
    double fpa_deg;             // Angolo di rampa di volo (Flight-Path Angle) gamma(t) [deg]
    double decel_g;             // Fattore di carico assiale aerodinamico n_g(t) = D / (m * g0) [g]
    double heat_flux_wcm2;      // Flusso termico convettivo di ristagno di Sutton-Graves q_dot(t) [W/cm^2]
    double qdyn_kpa;            // Pressione dinamica q_dyn(t) = 0.5 * rho(h) * V^2 [kPa]
    double downrange_km;        // Distanza percorsa lungo traccia a terra (Downrange) s(t) [km]
} ReentrySample;

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
    double peak_decel_g;            // Decelerazione assiale massima n_g,max [g]
    double peak_decel_alt_km;       // Quota del picco di decelerazione meccanica [km]
    double peak_qdyn_kpa;           // Pressione dinamica massima q_dyn,max [kPa]
    double breakup_alt_km;          // Quota di Demise / Breakup strutturale e ablazione catastale [km]
    double breakup_vel_kms;         // Velocità al momento del Breakup strutturale [km/s]
    double total_entry_time_sec;    // Durata complessiva della discesa atmosferica da 120 km a 0 km [s]
    double total_downrange_km;      // Downrange complessivo dall'Entry Interface all'impatto/demise [km]
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
