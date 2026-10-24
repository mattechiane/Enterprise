#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "propagator.h"
#include "telemetry.h"

int main(void) {
    printf("==========================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 3: TELEMETRY & GEODETICS   \n");
    printf("==========================================================\n\n");

    KeplerianElements init_kep = {
        .a     = EARTH_RADIUS + 450.0,
        .e     = 0.001,
        .i     = 51.6 * DEG2RAD,
        .RAAN  = 30.0 * DEG2RAD,
        .omega = 0.0  * DEG2RAD,
        .theta = 0.0  * DEG2RAD
    };

    CartesianState state = keplerian_to_cartesian(init_kep, EARTH_MU, 500.0, 0.0);

    SpacecraftConfig sc = {
        .drag_coeff = 2.2,
        .drag_area  = 4.0,
        .thruster   = {
            .is_active  = false,
            .thrust     = 50.0,
            .isp        = 300.0,
            .thrust_dir = {0.0, 1.0, 0.0}
        }
    };

    ForceModelConfig forces = {
        .enable_j2   = true,
        .enable_drag = true,
        .sun_dir     = {1.0, 0.0, 0.0}
    };

    double dt = 10.0;            // Passo di integrazione: 10 s
    double total_time = 86400.0; // Durata missione: 24 ore
    double eclipse_seconds = 0.0;

    // Apriamo il file di telemetria nella cartella principale del progetto
    const char *csv_path = "../mission_telemetry.csv";
    FILE *csv_file = telemetry_open(csv_path);
    if (csv_file == NULL) {
        printf("[ERRORE] Impossibile creare il file di telemetria!\n");
        return 1;
    }

    int step_count = 0;
    while (state.epoch <= total_time) {
        // Manovra di innalzamento orbitale nei primi 5 minuti (300 s)
        if (state.epoch < 300.0) {
            sc.thruster.is_active  = true;
            sc.thruster.thrust_dir = vec_normalize(state.velocity);
        } else {
            sc.thruster.is_active  = false;
        }

        SpaceEnvironment env = evaluate_environment(state.position, forces.sun_dir);
        if (env.in_eclipse && state.epoch < total_time) {
            eclipse_seconds += dt;
        }

        // Salviamo un punto di telemetria ogni 30 secondi (per avere file leggeri e veloci)
        if (step_count % 3 == 0) {
            KeplerianElements kep_now = cartesian_to_keplerian(state, EARTH_MU);
            Vector3D pos_ecef = eci_to_ecef(state.position, state.epoch, 0.0);
            GeodeticCoord geo = ecef_to_geodetic(pos_ecef);
            telemetry_write_step(csv_file, state, kep_now, geo, env, sc.thruster.is_active);
        }

        if (state.epoch >= total_time) break;

        state = propagate_step_rk4(state, dt, &sc, &forces);
        step_count++;
    }

    telemetry_close(csv_file);

    KeplerianElements final_kep = cartesian_to_keplerian(state, EARTH_MU);
    Vector3D final_ecef = eci_to_ecef(state.position, state.epoch, 0.0);
    GeodeticCoord final_geo = ecef_to_geodetic(final_ecef);

    printf("[TELEMETRIA SALVATA CON SUCCESSO]\n");
    printf("  File generato         : mission_telemetry.csv (%d record registrati)\n\n", (step_count / 3) + 1);

    printf("[STATO FINALE GEODETICO (T = 24.0 ore)]\n");
    printf("  Latitudine WGS84      : %+10.4f deg\n", final_geo.latitude * RAD2DEG);
    printf("  Longitudine           : %+10.4f deg\n", final_geo.longitude * RAD2DEG);
    printf("  Quota Geodetica reale : %10.2f km (Semiasse a: %.2f km)\n", final_geo.altitude, final_kep.a);
    printf("  Precessione Nodo J2   : %+10.4f deg/giorno\n", (final_kep.RAAN - init_kep.RAAN) * RAD2DEG);
    printf("==========================================================\n");

    return 0;
}