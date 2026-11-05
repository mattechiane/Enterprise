#include "astrodynamics.h"
#include <math.h>

// Funzione di sicurezza interna per evitare NaN in acos() dovuti a imprecisioni numeriche
static double clamp(double val, double min_val, double max_val) {
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

// ==========================================
// 1. IMPLEMENTAZIONE ALGEBRA VETTORIALE 3D
// ==========================================
Vector3D vec_add(Vector3D a, Vector3D b) {
    Vector3D res = {a.x + b.x, a.y + b.y, a.z + b.z};
    return res;
}

Vector3D vec_sub(Vector3D a, Vector3D b) {
    Vector3D res = {a.x - b.x, a.y - b.y, a.z - b.z};
    return res;
}

Vector3D vec_scale(Vector3D v, double scalar) {
    Vector3D res = {v.x * scalar, v.y * scalar, v.z * scalar};
    return res;
}

double vec_dot(Vector3D a, Vector3D b) {
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

Vector3D vec_cross(Vector3D a, Vector3D b) {
    Vector3D res = {
        (a.y * b.z) - (a.z * b.y),
        (a.z * b.x) - (a.x * b.z),
        (a.x * b.y) - (a.y * b.x)
    };
    return res;
}

double vec_norm(Vector3D v) {
    return sqrt(vec_dot(v, v));
}

Vector3D vec_normalize(Vector3D v) {
    double norm = vec_norm(v);
    if (norm < 1e-12) {
        Vector3D zero = {0.0, 0.0, 0.0};
        return zero;
    }
    return vec_scale(v, 1.0 / norm);
}

// ==========================================
// 2. IMPLEMENTAZIONE CONVERSIONI ORBITALI
// ==========================================
KeplerianElements cartesian_to_keplerian(CartesianState state, double mu) {
    KeplerianElements kep = {0};
    if (mu <= 0.0) return kep;

    Vector3D r_vec = state.position;
    Vector3D v_vec = state.velocity;

    double r = vec_norm(r_vec);
    double v = vec_norm(v_vec);
    if (r < 1e-12) {
        kep.a = EARTH_RADIUS + 400.0;
        return kep;
    }
    double v_radial = vec_dot(r_vec, v_vec) / r;

    // 1. Momento angolare specifico: h = r x v
    Vector3D h_vec = vec_cross(r_vec, v_vec);
    double h = vec_norm(h_vec);
    if (h < 1e-12) {
        kep.a = r;
        return kep;
    }

    // 2. Inclinazione (i)
    kep.i = acos(clamp(h_vec.z / h, -1.0, 1.0));

    // 3. Linea dei nodi: N = K x h = (-hy, hx, 0)
    Vector3D n_vec = {-h_vec.y, h_vec.x, 0.0};
    double n = vec_norm(n_vec);

    // 4. Ascensione Retta del Nodo Ascendente (RAAN)
    if (n > 1e-10) {
        kep.RAAN = acos(clamp(n_vec.x / n, -1.0, 1.0));
        if (n_vec.y < 0.0) {
            kep.RAAN = (2.0 * PI) - kep.RAAN;
        }
    } else {
        kep.RAAN = 0.0; // Orbita equatoriale
    }

    // 5. Vettore eccentricità: e = (1/mu) * [ (v^2 - mu/r)*r - (r.v)*v ]
    double term1 = (v * v) - (mu / r);
    double term2 = vec_dot(r_vec, v_vec);
    Vector3D e_vec = vec_scale(vec_sub(vec_scale(r_vec, term1), vec_scale(v_vec, term2)), 1.0 / mu);
    kep.e = vec_norm(e_vec);

    // 6. Energia meccanica specifica e Semiasse Maggiore (a)
    double energy = (0.5 * v * v) - (mu / r);
    if (fabs(energy) > 1e-12) {
        kep.a = -mu / (2.0 * energy);
    } else {
        kep.a = r / (1.0 - clamp(kep.e, 0.0, 0.99));
    }

    // 7. Argomento del pericentro (omega) - gestisce anche orbite equatoriali (i -> 0)
    if (kep.e > 1e-10) {
        if (n > 1e-10) {
            kep.omega = acos(clamp(vec_dot(n_vec, e_vec) / (n * kep.e), -1.0, 1.0));
            if (e_vec.z < 0.0) {
                kep.omega = (2.0 * PI) - kep.omega;
            }
        } else {
            // Orbita equatoriale eccentrica: longitudine del pericentro
            kep.omega = acos(clamp(e_vec.x / kep.e, -1.0, 1.0));
            if (e_vec.y < 0.0) {
                kep.omega = (2.0 * PI) - kep.omega;
            }
            if (h_vec.z < 0.0) {
                kep.omega = fmod((2.0 * PI) - kep.omega, 2.0 * PI);
            }
        }
    } else {
        kep.omega = 0.0;
    }

    // 8. Anomalia Vera (theta) - gestisce anche orbite circolari (e -> 0) ed equatoriali (i -> 0)
    if (kep.e > 1e-10) {
        kep.theta = acos(clamp(vec_dot(e_vec, r_vec) / (kep.e * r), -1.0, 1.0));
        if (v_radial < 0.0) {
            kep.theta = (2.0 * PI) - kep.theta;
        }
    } else if (n > 1e-10) {
        // Orbita circolare inclinata: argomento di latitudine u
        kep.theta = acos(clamp(vec_dot(n_vec, r_vec) / (n * r), -1.0, 1.0));
        if (r_vec.z < 0.0) {
            kep.theta = (2.0 * PI) - kep.theta;
        }
    } else {
        // Orbita circolare equatoriale: longitudine vera l
        kep.theta = acos(clamp(r_vec.x / r, -1.0, 1.0));
        if (r_vec.y < 0.0) {
            kep.theta = (2.0 * PI) - kep.theta;
        }
        if (h_vec.z < 0.0) {
            kep.theta = fmod((2.0 * PI) - kep.theta, 2.0 * PI);
        }
    }

    return kep;
}

CartesianState keplerian_to_cartesian(KeplerianElements kep, double mu, double mass, double epoch) {
    CartesianState state = {0};
    state.mass = mass;
    state.epoch = epoch;
    if (mu <= 0.0 || kep.a <= 0.0) return state;

    double e_safe = clamp(kep.e, 0.0, 0.9999);
    // Semilato retto: p = a * (1 - e^2)
    double p = kep.a * (1.0 - (e_safe * e_safe));
    if (p < 1e-9) p = 1e-9;
    double denom = 1.0 + e_safe * cos(kep.theta);
    if (fabs(denom) < 1e-9) denom = 1e-9;
    double r = p / denom;

    // Posizione e velocità nel sistema perifocale (PQW)
    double r_p = r * cos(kep.theta);
    double r_q = r * sin(kep.theta);

    double sqrt_mu_p = sqrt(mu / p);
    double v_p = -sqrt_mu_p * sin(kep.theta);
    double v_q = sqrt_mu_p * (e_safe + cos(kep.theta));

    // Pre-calcolo funzioni trigonometriche per la matrice di rotazione PQW -> ECI
    double cos_O = cos(kep.RAAN);
    double sin_O = sin(kep.RAAN);
    double cos_w = cos(kep.omega);
    double sin_w = sin(kep.omega);
    double cos_i = cos(kep.i);
    double sin_i = sin(kep.i);

    // Elementi della matrice di rotazione
    double R11 = (cos_O * cos_w) - (sin_O * sin_w * cos_i);
    double R12 = (-cos_O * sin_w) - (sin_O * cos_w * cos_i);
    double R21 = (sin_O * cos_w) + (cos_O * sin_w * cos_i);
    double R22 = (-sin_O * sin_w) + (cos_O * cos_w * cos_i);
    double R31 = sin_w * sin_i;
    double R32 = cos_w * sin_i;

    // Trasformazione nel sistema inerziale geocentrico (ECI)
    state.position.x = (R11 * r_p) + (R12 * r_q);
    state.position.y = (R21 * r_p) + (R22 * r_q);
    state.position.z = (R31 * r_p) + (R32 * r_q);

    state.velocity.x = (R11 * v_p) + (R12 * v_q);
    state.velocity.y = (R21 * v_p) + (R22 * v_q);
    state.velocity.z = (R31 * v_p) + (R32 * v_q);

    return state;
}

double orbital_period(double a, double mu) {
    if (a <= 0.0) return 0.0;
    return 2.0 * PI * sqrt((a * a * a) / mu);
}