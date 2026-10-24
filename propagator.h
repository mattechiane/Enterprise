#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "orbit_types.h"

// Struttura interna per le derivate temporali dello stato: d/dt [r, v, m]
typedef struct {
    Vector3D d_position; // Velocità (km/s)
    Vector3D d_velocity; // Accelerazione totale (km/s^2)
    double   d_mass;     // Portata massica dm/dt (kg/s)
} StateDerivative;

// Calcola la densità atmosferica (kg/m^3) e verifica se il satellite è in eclissi
SpaceEnvironment evaluate_environment(Vector3D pos, Vector3D sun_dir);

// Calcola l'accelerazione perturbante dovuta allo schiacciamento terrestre J2 (km/s^2)
Vector3D accel_j2(Vector3D pos);

// Calcola l'accelerazione di frenamento aerodinamico - Drag (km/s^2)
Vector3D accel_drag(CartesianState state, const SpacecraftConfig *sc, double density);

// Equazione differenziale completa del moto: calcola derivate di posizione, velocità e massa
StateDerivative compute_derivatives(CartesianState state, const SpacecraftConfig *sc, const ForceModelConfig *forces);

// Integratore numerico Runge-Kutta 4° Ordine (avanza lo stato di un passo dt in secondi)
CartesianState propagate_step_rk4(CartesianState state, double dt, const SpacecraftConfig *sc, const ForceModelConfig *forces);

#endif // PROPAGATOR_H