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

#endif // MISSION_ANALYSIS_H