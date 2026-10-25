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