#include "propagator.h"
#include "astrodynamics.h"
#include <math.h>

SpaceEnvironment evaluate_environment(Vector3D pos, Vector3D sun_dir) {
    SpaceEnvironment env;
    double r = vec_norm(pos);
    double altitude = r - EARTH_RADIUS; // Quota in km

    // 1. Modello atmosferico esponenziale a fasce (Standard Vallado semplificato)
    if (altitude < 0.0) {
        env.air_density = 1.225;
    } else if (altitude < 1000.0) {
        // Riferimento a h0 = 700 km, rho0 = 3.614e-14 kg/m^3, scala H = 88.667 km
        // Per quote più basse (es. 300-600 km) usiamo una parametrizzazione media realistica
        double rho0 = 1.0e-11; // kg/m^3 a 400 km
        double h0   = 400.0;   // km
        double H    = 60.0;    // Altezza di scala (km)
        env.air_density = rho0 * exp(-(altitude - h0) / H);
    } else {
        env.air_density = 0.0; // Oltre 1000 km il drag è trascurabile
    }

    // 2. Modello geometrico cilindrico di eclissi solare
    Vector3D s_hat = vec_normalize(sun_dir);
    double proj = vec_dot(pos, s_hat); // Proiezione lungo la linea Terra-Sole

    if (proj < 0.0) {
        // Il satellite è nell'emisfero opposto al Sole: calcoliamo la distanza dall'asse d'ombra
        Vector3D perp = vec_sub(pos, vec_scale(s_hat, proj));
        double dist_from_axis = vec_norm(perp);
        env.in_eclipse = (dist_from_axis < EARTH_RADIUS);
    } else {
        env.in_eclipse = false;
    }

    return env;
}

Vector3D accel_j2(Vector3D pos) {
    double r2 = vec_dot(pos, pos);
    double r  = sqrt(r2);
    double r5 = r2 * r2 * r;

    double z2_r2 = (pos.z * pos.z) / r2;
    double factor = -1.5 * EARTH_J2 * EARTH_MU * (EARTH_RADIUS * EARTH_RADIUS) / r5;

    Vector3D a_j2 = {
        factor * pos.x * (1.0 - 5.0 * z2_r2),
        factor * pos.y * (1.0 - 5.0 * z2_r2),
        factor * pos.z * (3.0 - 5.0 * z2_r2)
    };
    return a_j2;
}

Vector3D accel_drag(CartesianState state, const SpacecraftConfig *sc, double density) {
    if (density <= 0.0 || state.mass <= 0.0) {
        Vector3D zero = {0.0, 0.0, 0.0};
        return zero;
    }

    double v_km_s = vec_norm(state.velocity);
    // Conversione coerente: Area in km^2 (m^2 * 1e-6), densità in kg/km^3 (kg/m^3 * 1e9)
    // Il prodotto (1e9 * 1e-6) = 1000.0
    double coeff = -0.5 * density * 1000.0 * (sc->drag_coeff * sc->drag_area / state.mass) * v_km_s;
    return vec_scale(state.velocity, coeff);
}

StateDerivative compute_derivatives(CartesianState state, const SpacecraftConfig *sc, const ForceModelConfig *forces) {
    StateDerivative deriv;
    deriv.d_position = state.velocity;

    double r = vec_norm(state.position);
    double r3 = r * r * r;

    // 1. Gravità kepleriana centrale (Due Corpi)
    Vector3D a_total = vec_scale(state.position, -EARTH_MU / r3);

    // 2. Perturbazione Terra oblata (J2)
    if (forces->enable_j2) {
        a_total = vec_add(a_total, accel_j2(state.position));
    }

    // 3. Perturbazione Atmosferica (Drag)
    if (forces->enable_drag) {
        SpaceEnvironment env = evaluate_environment(state.position, forces->sun_dir);
        a_total = vec_add(a_total, accel_drag(state, sc, env.air_density));
    }

    // 4. Propulsione (Spinta finita e consumo massa)
    deriv.d_mass = 0.0;
    if (sc->thruster.is_active && state.mass > 0.0 && sc->thruster.thrust > 0.0) {
        // Accelerazione spinta: Thrust (N) / Mass (kg) = m/s^2 -> diviso 1000 per avere km/s^2
        double a_thrust_mag = (sc->thruster.thrust / state.mass) / 1000.0;
        Vector3D u_hat = vec_normalize(sc->thruster.thrust_dir);
        a_total = vec_add(a_total, vec_scale(u_hat, a_thrust_mag));

        // Equazione di Tsiolkovsky differenziale: dm/dt = - T / (Isp * g0)
        deriv.d_mass = -sc->thruster.thrust / (sc->thruster.isp * G0);
    }

    deriv.d_velocity = a_total;
    return deriv;
}

static CartesianState add_scaled_derivative(CartesianState s, StateDerivative d, double dt) {
    CartesianState next;
    next.position = vec_add(s.position, vec_scale(d.d_position, dt));
    next.velocity = vec_add(s.velocity, vec_scale(d.d_velocity, dt));
    next.mass     = s.mass + (d.d_mass * dt);
    next.epoch    = s.epoch + dt;
    return next;
}

CartesianState propagate_step_rk4(CartesianState state, double dt, const SpacecraftConfig *sc, const ForceModelConfig *forces) {
    StateDerivative k1 = compute_derivatives(state, sc, forces);
    StateDerivative k2 = compute_derivatives(add_scaled_derivative(state, k1, 0.5 * dt), sc, forces);
    StateDerivative k3 = compute_derivatives(add_scaled_derivative(state, k2, 0.5 * dt), sc, forces);
    StateDerivative k4 = compute_derivatives(add_scaled_derivative(state, k3, dt), sc, forces);

    CartesianState next_state;
    double dt6 = dt / 6.0;

    next_state.position.x = state.position.x + dt6 * (k1.d_position.x + 2.0*k2.d_position.x + 2.0*k3.d_position.x + k4.d_position.x);
    next_state.position.y = state.position.y + dt6 * (k1.d_position.y + 2.0*k2.d_position.y + 2.0*k3.d_position.y + k4.d_position.y);
    next_state.position.z = state.position.z + dt6 * (k1.d_position.z + 2.0*k2.d_position.z + 2.0*k3.d_position.z + k4.d_position.z);

    next_state.velocity.x = state.velocity.x + dt6 * (k1.d_velocity.x + 2.0*k2.d_velocity.x + 2.0*k3.d_velocity.x + k4.d_velocity.x);
    next_state.velocity.y = state.velocity.y + dt6 * (k1.d_velocity.y + 2.0*k2.d_velocity.y + 2.0*k3.d_velocity.y + k4.d_velocity.y);
    next_state.velocity.z = state.velocity.z + dt6 * (k1.d_velocity.z + 2.0*k2.d_velocity.z + 2.0*k3.d_velocity.z + k4.d_velocity.z);

    next_state.mass  = state.mass + dt6 * (k1.d_mass + 2.0*k2.d_mass + 2.0*k3.d_mass + k4.d_mass);
    next_state.epoch = state.epoch + dt;

    return next_state;
}