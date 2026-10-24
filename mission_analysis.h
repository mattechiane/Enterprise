#ifndef MISSION_ANALYSIS_H
#define MISSION_ANALYSIS_H

#include "orbit_types.h"

Vector3D geodetic_to_ecef(GeodeticCoord geo);
AccessMetrics compute_ground_access(Vector3D sat_ecef, const GroundStation *station);
HohmannTransfer compute_hohmann_transfer(double r1, double r2, double mu);
CartesianState apply_impulsive_maneuver(CartesianState state, Vector3D delta_v_eci, double isp);

// Calcola la manovra ottimale 3D (cambio quota + cambio inclinazione combinato) e direzione di spinta
OptimalTransferPlan plan_optimal_transfer(KeplerianElements orb_init, KeplerianElements orb_target, double init_mass, double isp, double mu);

// Analizza un'orbita per 24 ore calcolando Luce/Ombra, Visibilità Stazione e Carburante di Mantenimento (Station-Keeping)
OrbitReport evaluate_orbit_regime(KeplerianElements kep, double mass, const SpacecraftConfig *sc, const GroundStation *station);

// Stampa a schermo la scheda tecnica completa di un'orbita
void print_orbit_report(const char *title, OrbitReport rep, const char *station_name);

#endif // MISSION_ANALYSIS_H