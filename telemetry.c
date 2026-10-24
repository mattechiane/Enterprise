#include "telemetry.h"
#include "astrodynamics.h"
#include <math.h>

Vector3D eci_to_ecef(Vector3D pos_eci, double epoch, double gmst0) {
    // Angolo di rotazione terrestre (Greenwich Mean Sidereal Time) all'epoca attuale
    double theta_gmst = gmst0 + (EARTH_OMEGA * epoch);
    double cos_t = cos(theta_gmst);
    double sin_t = sin(theta_gmst);

    // Rotazione attorno all'asse Z di +theta_gmst
    Vector3D pos_ecef = {
        ( cos_t * pos_eci.x) + (sin_t * pos_eci.y),
        (-sin_t * pos_eci.x) + (cos_t * pos_eci.y),
        pos_eci.z
    };
    return pos_ecef;
}

GeodeticCoord ecef_to_geodetic(Vector3D pos_ecef) {
    GeodeticCoord geo;
    double x = pos_ecef.x;
    double y = pos_ecef.y;
    double z = pos_ecef.z;

    // 1. Longitudine esatta (-PI a +PI)
    geo.longitude = atan2(y, x);

    // 2. Metodo iterativo di Bowring/Heiskanen per Latitudine e Quota su Ellissoide WGS84
    double p = sqrt((x * x) + (y * y));
    double e2 = 2.0 * EARTH_FLATTENING - (EARTH_FLATTENING * EARTH_FLATTENING); // Quadrato prima eccentricità terrestre

    double lat = atan2(z, p * (1.0 - e2));
    double N = EARTH_RADIUS;

    for (int iter = 0; iter < 5; iter++) {
        double sin_lat = sin(lat);
        N = EARTH_RADIUS / sqrt(1.0 - (e2 * sin_lat * sin_lat));
        lat = atan2(z + (e2 * N * sin_lat), p);
    }

    geo.latitude = lat;
    if (fabs(cos(lat)) > 1e-6) {
        geo.altitude = (p / cos(lat)) - N;
    } else {
        geo.altitude = (fabs(z) / fabs(sin(lat))) - (N * (1.0 - e2));
    }

    return geo;
}

FILE* telemetry_open(const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (fp != NULL) {
        fprintf(fp, "time_s,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,mass_kg,sma_km,ecc,inc_deg,raan_deg,argp_deg,nu_deg,lat_deg,lon_deg,alt_km,in_eclipse,thruster_on\n");
    }
    return fp;
}

void telemetry_write_step(FILE *fp, CartesianState state, KeplerianElements kep, GeodeticCoord geo, SpaceEnvironment env, bool thruster_on) {
    if (fp == NULL) return;
    fprintf(fp, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f,%.4f,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d\n",
            state.epoch,
            state.position.x, state.position.y, state.position.z,
            state.velocity.x, state.velocity.y, state.velocity.z,
            state.mass,
            kep.a, kep.e,
            kep.i * RAD2DEG, kep.RAAN * RAD2DEG, kep.omega * RAD2DEG, kep.theta * RAD2DEG,
            geo.latitude * RAD2DEG, geo.longitude * RAD2DEG, geo.altitude,
            env.in_eclipse ? 1 : 0,
            thruster_on ? 1 : 0);
}

void telemetry_close(FILE *fp) {
    if (fp != NULL) {
        fclose(fp);
    }
}