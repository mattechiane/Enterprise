#include "astrodynamics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

// ==========================================
// 3. PARSER STANDARD NORAD TWO-LINE ELEMENTS (TLE)
// ==========================================
static void slice_tle_field(const char *line, int start_col_1based, int end_col_1based, char *buf, size_t buf_size) {
    if (!line || !buf || buf_size == 0) return;
    int len = (int)strlen(line);
    int idx_start = start_col_1based - 1;
    int idx_end   = end_col_1based - 1;
    if (idx_start < 0) idx_start = 0;
    if (idx_start >= len || idx_end < idx_start) {
        buf[0] = '\0';
        return;
    }
    if (idx_end >= len) idx_end = len - 1;
    size_t copy_len = (size_t)(idx_end - idx_start + 1);
    if (copy_len >= buf_size) copy_len = buf_size - 1;
    memcpy(buf, line + idx_start, copy_len);
    buf[copy_len] = '\0';
}

static bool verify_tle_line_checksum(const char *line) {
    if (!line || strlen(line) < 69) return false;
    int sum = 0;
    for (int i = 0; i < 68; i++) {
        char c = line[i];
        if (c >= '0' && c <= '9') {
            sum += (c - '0');
        } else if (c == '-') {
            sum += 1;
        }
    }
    char expected = (char)('0' + (sum % 10));
    return (line[68] == expected);
}

bool parse_norad_tle_extended(const char *line1, const char *line2, NoradTLEData *tle_out) {
    if (!line1 || !line2 || !tle_out) return false;
    memset(tle_out, 0, sizeof(NoradTLEData));

    // Trim eventual leading spaces
    while (*line1 == ' ' || *line1 == '\t' || *line1 == '\r' || *line1 == '\n') line1++;
    while (*line2 == ' ' || *line2 == '\t' || *line2 == '\r' || *line2 == '\n') line2++;

    if (strlen(line1) < 32 || strlen(line2) < 63) return false;
    if (line1[0] != '1' || line2[0] != '2') return false;

    char buf[32];

    // Line 1: Catalog Number (cols 3-7), Epoch Year (cols 19-20), Epoch Day (cols 21-32)
    slice_tle_field(line1, 3, 7, buf, sizeof(buf));
    tle_out->sat_catalog_num = atoi(buf);

    slice_tle_field(line1, 19, 20, buf, sizeof(buf));
    int yy = atoi(buf);
    tle_out->epoch_year = (yy >= 57) ? (1900 + yy) : (2000 + yy);

    slice_tle_field(line1, 21, 32, buf, sizeof(buf));
    tle_out->epoch_day = atof(buf);

    // Line 2:
    // Cols 9-16: Inclination [deg]
    slice_tle_field(line2, 9, 16, buf, sizeof(buf));
    double inc_deg = atof(buf);

    // Cols 18-25: RAAN [deg]
    slice_tle_field(line2, 18, 25, buf, sizeof(buf));
    double raan_deg = atof(buf);

    // Cols 27-33: Eccentricity (7 digits with implied decimal point "0.xxxxxxx")
    slice_tle_field(line2, 27, 33, buf, sizeof(buf));
    char ecc_buf[36];
    ecc_buf[0] = '0';
    ecc_buf[1] = '.';
    size_t k = 2;
    for (size_t i = 0; buf[i] != '\0' && k < sizeof(ecc_buf) - 1; i++) {
        if (isdigit((unsigned char)buf[i])) {
            ecc_buf[k++] = buf[i];
        }
    }
    ecc_buf[k] = '\0';
    double ecc = clamp(atof(ecc_buf), 0.0, 0.99);

    // Cols 35-42: Argument of Perigee [deg]
    slice_tle_field(line2, 35, 42, buf, sizeof(buf));
    double omega_deg = atof(buf);

    // Cols 44-51: Mean Anomaly M [deg]
    slice_tle_field(line2, 44, 51, buf, sizeof(buf));
    double mean_anom_deg = atof(buf);

    // Cols 53-63: Mean Motion n [rev/day]
    slice_tle_field(line2, 53, 63, buf, sizeof(buf));
    double n_rev_day = atof(buf);
    if (n_rev_day <= 0.01) return false;

    // Conversione Mean Motion [rev/day] -> [rad/s] -> Semiasse Maggiore a = (mu / n^2)^(1/3)
    double n_rad_s = n_rev_day * (2.0 * PI) / 86400.0;
    double sma_km  = cbrt(EARTH_MU / (n_rad_s * n_rad_s));

    // Risoluzione Equazione di Keplero M = E - e * sin(E) tramite Newton-Raphson
    double M_rad = fmod(mean_anom_deg * DEG2RAD, 2.0 * PI);
    if (M_rad < 0.0) M_rad += 2.0 * PI;
    double E_rad = (ecc > 0.8) ? PI : M_rad;
    for (int iter = 0; iter < 30; iter++) {
        double f_val = E_rad - ecc * sin(E_rad) - M_rad;
        double f_der = 1.0 - ecc * cos(E_rad);
        if (fabs(f_der) < 1e-12) break;
        double dE = f_val / f_der;
        E_rad -= dE;
        if (fabs(dE) < 1e-12) break;
    }

    // Conversione Anomalia Eccentrica E -> Anomalia Vera theta
    double theta_rad = 2.0 * atan2(sqrt(1.0 + ecc) * sin(0.5 * E_rad),
                                   sqrt(1.0 - ecc) * cos(0.5 * E_rad));
    if (theta_rad < 0.0) theta_rad += 2.0 * PI;

    tle_out->mean_motion_rev_day = n_rev_day;
    tle_out->mean_anomaly_deg    = mean_anom_deg;
    tle_out->eccentric_anom_deg  = E_rad * RAD2DEG;
    tle_out->checksum_valid      = verify_tle_line_checksum(line1) && verify_tle_line_checksum(line2);

    tle_out->kep.a     = sma_km;
    tle_out->kep.e     = ecc;
    tle_out->kep.i     = clamp(inc_deg, 0.0, 180.0) * DEG2RAD;
    tle_out->kep.RAAN  = fmod(raan_deg + 360.0, 360.0) * DEG2RAD;
    tle_out->kep.omega = fmod(omega_deg + 360.0, 360.0) * DEG2RAD;
    tle_out->kep.theta = theta_rad;

    return true;
}

bool parse_norad_tle(const char *line1, const char *line2, KeplerianElements *kep_out) {
    if (!kep_out) return false;
    NoradTLEData tle_data;
    if (!parse_norad_tle_extended(line1, line2, &tle_data)) {
        return false;
    }
    *kep_out = tle_data.kep;
    return true;
}