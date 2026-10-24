#ifndef ASTRODYNAMICS_H
#define ASTRODYNAMICS_H

#include "orbit_types.h"

// Costanti matematiche di conversione
#define PI 3.14159265358979323846
#define RAD2DEG (180.0 / PI)
#define DEG2RAD (PI / 180.0)

// ==========================================
// 1. OPERAZIONI VETTORIALI 3D
// ==========================================
Vector3D vec_add(Vector3D a, Vector3D b);
Vector3D vec_sub(Vector3D a, Vector3D b);
Vector3D vec_scale(Vector3D v, double scalar);
double   vec_dot(Vector3D a, Vector3D b);
Vector3D vec_cross(Vector3D a, Vector3D b);
double   vec_norm(Vector3D v);
Vector3D vec_normalize(Vector3D v);

// ==========================================
// 2. CONVERSIONI DI STATO ORBITALE
// ==========================================
// Converte Posizione e Velocità (ECI) nei 6 Elementi Kepleriani
KeplerianElements cartesian_to_keplerian(CartesianState state, double mu);

// Converte i 6 Elementi Kepleriani in Posizione e Velocità (ECI)
CartesianState keplerian_to_cartesian(KeplerianElements kep, double mu, double mass, double epoch);

// Calcola il periodo orbitale in secondi (per orbite chiuse, e < 1)
double orbital_period(double a, double mu);

#endif // ASTRODYNAMICS_H