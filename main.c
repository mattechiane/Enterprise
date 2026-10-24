#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "propagator.h"
#include "telemetry.h"
#include "mission_analysis.h"

int main(void) {
    printf("==================================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 4: MISSION ANALYSIS & GROUND LOS   \n");
    printf("==================================================================\n\n");

    // 1. PIANIFICAZIONE MANOVRA DI HOHMANN (Da 400 km a 550 km di quota)
    double r_park   = EARTH_RADIUS + 400.0;
    double r_target = EARTH_RADIUS + 550.0;
    HohmannTransfer hohmann = compute_hohmann_transfer(r_park, r_target, EARTH_MU);

    printf("[1] PIANIFICATORE MANOVRA DI HOHMANN (400 km -> 550 km):\n");
    printf("    Impulso 1 (Perigeo)   : %+8.2f m/s\n", hohmann.dv1 * 1000.0);
    printf("    Impulso 2 (Apogeo)    : %+8.2f m/s\n", hohmann.dv2 * 1000.0);
    printf("    Delta-V Totale        : %8.2f m/s\n", hohmann.total_dv * 1000.0);
    printf("    Tempo di Trasferimento: %8.2f min\n\n", hohmann.transfer_time / 60.0);

    // 2. CONFIGURAZIONE ORBITA OPERATIVA (550 km, inclinazione 51.6 deg)
    KeplerianElements init_kep = {
        .a     = r_target,
        .e     = 0.001,
        .i     = 51.6 * DEG2RAD,
        .RAAN  = 30.0 * DEG2RAD,
        .omega = 0.0  * DEG2RAD,
        .theta = 0.0  * DEG2RAD
    };

    CartesianState state = keplerian_to_cartesian(init_kep, EARTH_MU, 500.0, 0.0);

    // Applichiamo a T=0 il consumo dei due impulsi di Hohmann (Isp = 300 s)
    Vector3D dv_vec = vec_scale(vec_normalize(state.velocity), hohmann.total_dv);
    CartesianState post_burn = apply_impulsive_maneuver(state, dv_vec, 300.0);
    state.mass = post_burn.mass; // Aggiorniamo la massa post-manovra di inserimento
    printf("    Massa post-inserimento: %8.2f kg (Carburante Hohmann: %.2f kg)\n\n",
           state.mass, 500.0 - state.mass);

    // 3. DEFINIZIONE STAZIONE DI TERRA: MILANO (Bovisa)
    GroundStation milano = {
        .name = "Milano Ground Station",
        .coord = {
            .latitude  = 45.502 * DEG2RAD,
            .longitude = 9.156  * DEG2RAD,
            .altitude  = 0.140 // 140 metri s.l.m.
        },
        .min_elevation = 10.0 * DEG2RAD // Contatto valido sopra 10 gradi di elevazione
    };

    SpacecraftConfig sc = {
        .drag_coeff = 2.2,
        .drag_area  = 4.0,
        .thruster   = {.is_active = false, .thrust = 0.0, .isp = 300.0}
    };

    ForceModelConfig forces = {
        .enable_j2   = true,
        .enable_drag = true,
        .sun_dir     = {1.0, 0.0, 0.0}
    };

    printf("[2] ANALISI FINESTRE DI VISIBILITA' (Stazione: %s, El_min = 10 deg):\n", milano.name);
    printf("    ----------------------------------------------------------------------\n");
    printf("    Pass # |   AOS (s)   |   LOS (s)   | Durata (min) | Max El (deg) | Min Range (km)\n");
    printf("    ----------------------------------------------------------------------\n");

    double dt = 5.0;             // Passo ad alta risoluzione (5 secondi) per gli accessi
    double total_time = 86400.0; // 24 ore
    bool in_pass = false;
    int pass_count = 0;
    double aos_time = 0.0;
    double max_el = 0.0;
    double min_range = 1e9;

    while (state.epoch <= total_time) {
        Vector3D pos_ecef = eci_to_ecef(state.position, state.epoch, 0.0);
        AccessMetrics access = compute_ground_access(pos_ecef, &milano);

        if (access.is_visible) {
            if (!in_pass) {
                // Inizio di un nuovo passaggio sopra la stazione (AOS)
                in_pass = true;
                pass_count++;
                aos_time = state.epoch;
                max_el = access.elevation;
                min_range = access.range;
            } else {
                if (access.elevation > max_el) max_el = access.elevation;
                if (access.range < min_range)  min_range = access.range;
            }
        } else if (in_pass) {
            // Fine del passaggio sopra la stazione (LOS)
            in_pass = false;
            double los_time = state.epoch;
            double duration_min = (los_time - aos_time) / 60.0;
            printf("      %2d   |  %8.1f   |  %8.1f   |    %6.2f    |    %6.2f    |    %7.1f\n",
                   pass_count, aos_time, los_time, duration_min, max_el * RAD2DEG, min_range);
        }

        if (state.epoch >= total_time) break;
        state = propagate_step_rk4(state, dt, &sc, &forces);
    }

    printf("    ----------------------------------------------------------------------\n");
    printf("    Totale passaggi utili in 24 ore: %d\n", pass_count);
    printf("==================================================================\n");

    return 0;
}