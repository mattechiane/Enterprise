#include <stdio.h>
#include "orbit_types.h"
#include "astrodynamics.h"
#include "mission_analysis.h"
#include "telemetry.h"

int main(void) {
    printf("======================================================================\n");
    printf("   PROJECT ENTERPRISE - SPRINT 12: TERMINAL RENDEZVOUS & PROXIMITY OPS\n");
    printf("======================================================================\n\n");

    // 1. INPUT UTENTE: Massa, Propellente e Vincolo di Tempo Massimo
    double initial_mass = 500.0; // kg
    PropellantType selected_prop = PROP_LIQUID_MMH_NTO; // Liquido Bipropellente (Isp = 316 s)
    double max_time_limit_sec = 0.0; // 0 = Nessun limite (cerca il minimo Delta-V assoluto)

    // 2. INPUT ORBITA INIZIALE (Scelta tra INPUT_KEPLERIAN o INPUT_CARTESIAN)
    OrbitInputSpec input_init = {
        .mode = INPUT_KEPLERIAN,
        .kep  = {.a = EARTH_RADIUS + 400.0, .e = 0.001, .i = 51.6 * DEG2RAD, .RAAN = 30.0 * DEG2RAD, .omega = 0.0, .theta = 0.0}
    };

    // 3. INPUT ORBITA FINALE TARGET
    OrbitInputSpec input_target = {
        .mode = INPUT_KEPLERIAN,
        .kep  = {.a = EARTH_RADIUS + 650.0, .e = 0.001, .i = 53.0 * DEG2RAD, .RAAN = 30.0 * DEG2RAD, .omega = 0.0, .theta = 0.0}
    };

    // Risoluzione automatica coordinate
    KeplerianElements orb_initial = resolve_user_orbit_input(input_init, EARTH_MU);
    KeplerianElements orb_target  = resolve_user_orbit_input(input_target, EARTH_MU);

    PropellantSpec prop = get_propellant_specs(selected_prop);
    SpacecraftConfig sc = {
        .drag_coeff = 2.2,
        .drag_area  = 4.0,
        .thruster   = {.is_active = false, .thrust = 50.0, .isp = prop.isp}
    };

    GroundStation milano = {
        .name = "Milano-Bovisa",
        .coord = {.latitude = 45.502 * DEG2RAD, .longitude = 9.156 * DEG2RAD, .altitude = 0.140},
        .min_elevation = 10.0 * DEG2RAD
    };

    // Calcolo Report ed Esecuzione Ottimizzatore Globale di Trasferimento
    OrbitReport rep_init = evaluate_orbit_regime(orb_initial, initial_mass, &sc, &milano);
    OptimalTransferPlan plan = optimize_orbital_transfer(orb_initial, orb_target, initial_mass, selected_prop, max_time_limit_sec, EARTH_MU);
    OrbitReport rep_target = evaluate_orbit_regime(orb_target, plan.final_mass, &sc, &milano);

    // 4. SPRINT 12: TERMINAL RENDEZVOUS & DOCKING (CLOHESSY-WILTSHIRE / HILL LVLH)
    // Stato iniziale relativo Chaser vs Target: dx0 = -250 m (R-bar), dy0 = -2000 m (V-bar), dz0 = +100 m (H-bar)
    Vector3D r0_lvlh_m  = {-250.0, -2000.0, 100.0};
    Vector3D v0_lvlh_ms = {0.0, 0.0, 0.0};
    double t_rdv_sec    = 45.0 * 60.0; // 45 minuti di approccio finale

    RendezvousPlan rdv = compute_cw_rendezvous(orb_target.a, r0_lvlh_m, v0_lvlh_ms, t_rdv_sec,
                                               plan.final_mass, plan.prop_info.isp, EARTH_MU);

    printf("[OPTIMAL MANEUVER SELECTED]: %s\n", plan.maneuver_name);
    printf("  * Propellant            : %s [%s] (Isp = %.1f s)\n", plan.prop_info.name, plan.prop_info.category, plan.prop_info.isp);
    printf("  * Burn 1 Point (ECI)    : [ %+9.2f, %+9.2f, %+9.2f ] km (Lat %.1f deg, Lon %.1f deg)\n",
           plan.burn1_pos_eci.x, plan.burn1_pos_eci.y, plan.burn1_pos_eci.z,
           plan.burn1_geo.latitude * RAD2DEG, plan.burn1_geo.longitude * RAD2DEG);
    printf("  * Burn 1 Cost & Dir RTN : %7.2f m/s | Dir: [ %+5.3f, %+5.3f, %+5.3f ]\n",
           plan.dv1_mag * 1000.0, plan.burn1_dir_rtn.x, plan.burn1_dir_rtn.y, plan.burn1_dir_rtn.z);
    printf("  * Burn 2 Point (ECI)    : [ %+9.2f, %+9.2f, %+9.2f ] km (Lat %.1f deg, Lon %.1f deg)\n",
           plan.burn2_pos_eci.x, plan.burn2_pos_eci.y, plan.burn2_pos_eci.z,
           plan.burn2_geo.latitude * RAD2DEG, plan.burn2_geo.longitude * RAD2DEG);
    printf("  * Burn 2 Cost & Dir RTN : %7.2f m/s | Dir: [ %+5.3f, %+5.3f, %+5.3f ]\n",
           plan.dv2_mag * 1000.0, plan.burn2_dir_rtn.x, plan.burn2_dir_rtn.y, plan.burn2_dir_rtn.z);
    printf("  * Total Delta-V / Fuel  : %7.2f m/s | Fuel = %.2f kg | TOF = %.2f min\n\n",
           plan.total_dv * 1000.0, plan.propellant_used, plan.transfer_time_sec / 60.0);

    printf("[TERMINAL RENDEZVOUS & PROXIMITY OPS - CLOHESSY-WILTSHIRE LVLH]\n");
    printf("  * Initial LVLH Offset   : dX(R-bar) = %+.1f m, dY(V-bar) = %+.1f m, dZ(H-bar) = %+.1f m (Range = %.1f m)\n",
           rdv.r0_lvlh_m.x, rdv.r0_lvlh_m.y, rdv.r0_lvlh_m.z, rdv.initial_range_m);
    printf("  * CW Targeting Burn 1   : %6.3f m/s | LVLH: [ %+6.3f, %+6.3f, %+6.3f ] m/s\n",
           rdv.dv_inject_ms, rdv.dv_inject_lvlh_ms.x, rdv.dv_inject_lvlh_ms.y, rdv.dv_inject_lvlh_ms.z);
    printf("  * Soft-Docking Burn 2   : %6.3f m/s | LVLH: [ %+6.3f, %+6.3f, %+6.3f ] m/s\n",
           rdv.dv_brake_ms, rdv.dv_brake_lvlh_ms.x, rdv.dv_brake_lvlh_ms.y, rdv.dv_brake_lvlh_ms.z);
    printf("  * Total Rendezvous Cost : %6.3f m/s | Fuel = %.3f kg | Approach TOF = %.1f min\n\n",
           rdv.total_dv_ms, rdv.propellant_used_kg, rdv.transfer_time_sec / 60.0);

    if (export_web_mission_json("../mission_report.json", rep_init, plan, rep_target, rdv, &milano, &sc)) {
        printf("[SUCCESS] 'mission_report.json' updated with Burn Points, Propellant & CW Rendezvous Data!\n");
    }
    // Esporta copia anche nella cartella corrente nel caso l'eseguibile venga lanciato dalla root
    export_web_mission_json("mission_report.json", rep_init, plan, rep_target, rdv, &milano, &sc);

    printf("======================================================================\n");
    return 0;
}
