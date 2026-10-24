#ifndef MISSION_ANALYSIS_H
#define MISSION_ANALYSIS_H

#include "orbit_types.h"

#define MIN_SAFE_PERIGEE_ALT_KM 150.0  // Limite minimo atmosferico per orbita chiusa stabile (SLATS/GOCE limit)
#define MAX_HILL_APOGEE_ALT_KM 900000.0 // Limite massimo Sfera di Hill terrestre (orbita chiusa)

Vector3D geodetic_to_ecef(GeodeticCoord geo);
AccessMetrics compute_ground_access(Vector3D sat_ecef, const GroundStation *station);
HohmannTransfer compute_hohmann_transfer(double r1, double r2, double mu);
CartesianState apply_impulsive_maneuver(CartesianState state, Vector3D delta_v_eci, double isp);

PropellantSpec get_propellant_specs(PropellantType type);
KeplerianElements resolve_user_orbit_input(OrbitInputSpec input, double mu);

// Valida e corregge i parametri di un'orbita chiusa (impedisce collisione con la Terra o fuga iperbolica)
bool validate_and_clamp_closed_orbit(KeplerianElements *kep, char *warning_msg);

// Ottimizzatore Globale 3D: include eccentricità (e1, e2), cambio RAAN (Omega1 -> Omega2) e inclinazione (i1 -> i2)
OptimalTransferPlan optimize_orbital_transfer(KeplerianElements orb_init,
                                              KeplerianElements orb_target,
                                              double init_mass,
                                              PropellantType prop_type,
                                              double max_transfer_time_sec,
                                              double mu);

OrbitReport evaluate_orbit_regime(KeplerianElements kep, double mass, const SpacecraftConfig *sc, const GroundStation *station);
void print_orbit_report(const char *title, OrbitReport rep, const char *station_name);

#endif // MISSION_ANALYSIS_H