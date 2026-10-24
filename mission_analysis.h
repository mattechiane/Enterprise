#ifndef MISSION_ANALYSIS_H
#define MISSION_ANALYSIS_H

#include "orbit_types.h"

Vector3D geodetic_to_ecef(GeodeticCoord geo);
AccessMetrics compute_ground_access(Vector3D sat_ecef, const GroundStation *station);
HohmannTransfer compute_hohmann_transfer(double r1, double r2, double mu);
CartesianState apply_impulsive_maneuver(CartesianState state, Vector3D delta_v_eci, double isp);

// Restituisce i dati fisici (Isp, categoria, nome) del propellente scelto dall'utente
PropellantSpec get_propellant_specs(PropellantType type);

// Converte l'input utente (Cartesiano o Kepleriano) in elementi kepleriani standardizzati
KeplerianElements resolve_user_orbit_input(OrbitInputSpec input, double mu);

// Ottimizzatore Globale di Manovra: valuta Hohmann, Optimal Plane-Split, Bi-Ellittica e Fast Transfer
// rispettando eventuali limiti sul tempo massimo (max_transfer_time_sec <= 0 significa nessun limite)
OptimalTransferPlan optimize_orbital_transfer(KeplerianElements orb_init,
                                              KeplerianElements orb_target,
                                              double init_mass,
                                              PropellantType prop_type,
                                              double max_transfer_time_sec,
                                              double mu);

OrbitReport evaluate_orbit_regime(KeplerianElements kep, double mass, const SpacecraftConfig *sc, const GroundStation *station);
void print_orbit_report(const char *title, OrbitReport rep, const char *station_name);

#endif // MISSION_ANALYSIS_H