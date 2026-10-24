#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "mission_analysis.h"

int main(void) {
    printf("==================================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 5: END-TO-END MISSION ARCHITECT    \n");
    printf("==================================================================\n\n");

    // SOTTOSISTEMI SATELLITE E STAZIONE DI TERRA
    double initial_mass = 500.0; // kg
    SpacecraftConfig sc = {
        .drag_coeff = 2.2,
        .drag_area  = 4.0, // 4 m^2
        .thruster   = {.is_active = false, .thrust = 50.0, .isp = 300.0}
    };

    GroundStation milano = {
        .name = "Milano-Bovisa",
        .coord = {.latitude = 45.502 * DEG2RAD, .longitude = 9.156 * DEG2RAD, .altitude = 0.140},
        .min_elevation = 10.0 * DEG2RAD
    };

    // 1. INPUT UTENTE: ORBITA DI PARCHEGGIO INIZIALE (400 km, inclinazione 51.6 deg)
    KeplerianElements orb_initial = {
        .a = EARTH_RADIUS + 400.0, .e = 0.001, .i = 51.6 * DEG2RAD,
        .RAAN = 30.0 * DEG2RAD, .omega = 0.0, .theta = 0.0
    };

    // 2. INPUT UTENTE: ORBITA OPERATIVA TARGET (650 km, inclinazione 53.0 deg)
    KeplerianElements orb_target = {
        .a = EARTH_RADIUS + 650.0, .e = 0.001, .i = 53.0 * DEG2RAD,
        .RAAN = 30.0 * DEG2RAD, .omega = 0.0, .theta = 0.0
    };

    // FASE A: ANALISI ORBITA INIZIALE
    OrbitReport rep_init = evaluate_orbit_regime(orb_initial, initial_mass, &sc, &milano);
    print_orbit_report("[FASE 1] REPORT ORBITA INIZIALE (PARKING ORBIT)", rep_init, milano.name);

    // FASE B: CALCOLO MANOVRA OTTIMALE 3D (QUOTA + CAMBIO PIANO)
    OptimalTransferPlan plan = plan_optimal_transfer(orb_initial, orb_target, initial_mass, sc.thruster.isp, EARTH_MU);

    printf("------------------------------------------------------------------\n");
    printf("  [FASE 2] PIANO DI MANOVRA OTTIMALE (400 km, 51.6 deg -> 650 km, 53.0 deg)\n");
    printf("------------------------------------------------------------------\n");
    printf("  * Impulso 1 (Perigeo)   : %7.2f m/s | Direzione RTN: [ %+5.2f, %+5.2f, %+5.2f ]\n",
           plan.dv1_mag * 1000.0, plan.burn1_dir_rtn.x, plan.burn1_dir_rtn.y, plan.burn1_dir_rtn.z);
    printf("  * Impulso 2 (Apogeo 3D) : %7.2f m/s | Direzione RTN: [ %+5.2f, %+5.2f, %+5.2f ]\n",
           plan.dv2_mag * 1000.0, plan.burn2_dir_rtn.x, plan.burn2_dir_rtn.y, plan.burn2_dir_rtn.z);
    printf("  * Delta-V Totale Manovra: %7.2f m/s\n", plan.total_dv * 1000.0);
    printf("  * Durata Trasferimento  : %7.2f minuti\n", plan.transfer_time_sec / 60.0);
    printf("  * Carburante Manovra    : %7.2f kg  (Massa residua all'arrivo: %.2f kg)\n\n",
           plan.propellant_used, plan.final_mass);

    // FASE C: ANALISI NUOVA ORBITA TARGET (con la nuova massa ridotta post-manovra!)
    OrbitReport rep_target = evaluate_orbit_regime(orb_target, plan.final_mass, &sc, &milano);
    print_orbit_report("[FASE 3] REPORT NUOVA ORBITA OPERATIVA (TARGET ORBIT)", rep_target, milano.name);

    printf("==================================================================\n");
    return 0;
}