#ifndef ORBIT_TYPES_H
#define ORBIT_TYPES_H

#include <stdbool.h>

// ==========================================
// 1. COSTANTI FISICHE E PLANETARIE (TERRA - WGS84)
// ==========================================
#define EARTH_MU 398600.4418          // Costante gravitazionale terrestre (km^3/s^2)
#define EARTH_RADIUS 6378.137         // Raggio equatoriale terrestre WGS84 (km)
#define EARTH_FLATTENING (1.0 / 298.257223563) // Schiacciamento terrestre WGS84 (f)
#define EARTH_J2 0.00108262668        // Coefficiente armonico zonale J2
#define EARTH_OMEGA 7.2921158553e-5   // Velocità angolare di rotazione terrestre (rad/s)
#define G0 9.80665                    // Accelerazione di gravità standard (m/s^2)

// ==========================================
// 2. MATEMATICA E GEOMETRIA
// ==========================================
typedef struct {
    double x;
    double y;
    double z;
} Vector3D;

// ==========================================
// 3. STATO DEL VEICOLO SPAZIALE E COORDINATE
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

// Coordinate terrestri per Ground Track e stazioni di terra
typedef struct {
    double latitude;  // Latitudine geodetica WGS84 (radianti, -PI/2 a +PI/2)
    double longitude; // Longitudine terrestre (radianti, -PI a +PI)
    double altitude;  // Quota reale sopra l'ellissoide WGS84 (km)
} GeodeticCoord;

// ==========================================
// 4. SOTTOSISTEMI, PROPRIETA' FISICHE E AMBIENTE
// ==========================================
typedef struct {
    bool is_active;      // true se il motore è acceso
    double thrust;       // Spinta in Newton (N)
    double isp;          // Impulso Specifico in secondi
    Vector3D thrust_dir; // Direzione della spinta ECI (versore)
} Thruster;

typedef struct {
    double drag_coeff;   // Coefficiente di resistenza aerodinamica Cd
    double drag_area;    // Area frontale esposta al flusso in m^2
    Thruster thruster;   // Propulsore di bordo
} SpacecraftConfig;

typedef struct {
    bool enable_j2;      // Attiva la Terra non sferica (J2)
    bool enable_drag;    // Attiva l'attrito atmosferico
    Vector3D sun_dir;    // Versore direzione Terra -> Sole in ECI
} ForceModelConfig;

typedef struct {
    bool in_eclipse;     // true se il satellite è nel cono d'ombra terrestre
    double air_density;  // Densità atmosferica locale (kg/m^3)
} SpaceEnvironment;

#endif // ORBIT_TYPES_H