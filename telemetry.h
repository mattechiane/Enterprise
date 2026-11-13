#ifndef TELEMETRY_H
#define TELEMETRY_H

#include <stdio.h>
#include <stdbool.h>
#include "orbit_types.h"

// Conversioni di coordinate spaziali e terrestri
Vector3D eci_to_ecef(Vector3D pos_eci, double epoch, double gmst0);
GeodeticCoord ecef_to_geodetic(Vector3D pos_ecef);

// Funzioni per l'esportazione della telemetria in formato CSV (Opzionale/Debug)
FILE* telemetry_open(const char *filename);
void  telemetry_write_step(FILE *fp, CartesianState state, KeplerianElements kep, GeodeticCoord geo, SpaceEnvironment env, bool thruster_on);
void  telemetry_close(FILE *fp);

// Esporta le effemeridi orbitali in formato standard internazionale CCSDS 502.0-B-2 Orbit Ephemeris Message (CCSDS_OEM_VERS = 2.0)
// direttamente compatibile con NASA GMAT, Ansys STK ed ESOC Flight Dynamics
bool export_ccsds_oem_file(const char *filename,
                           OrbitReport rep_init,
                           OptimalTransferPlan plan,
                           OrbitReport rep_target,
                           const SpacecraftConfig *sc);

// Esporta l'intero scenario di missione (Orbita 1 -> Trasferimento -> Orbita 2 -> Rendezvous CW -> Monte Carlo & TCM -> ADCS 6-DOF -> Thermal TCS -> RF Telecom -> Systems Eng -> Onboard EKF -> SSA & CAM -> EoL Re-Entry -> Cislunar CR3BP -> Structural Slosh & CSI -> Autonomous FDIR)
// in formato JSON interattivo per la Web App 3D
bool export_web_mission_json(const char *filename,
                             OrbitReport rep_init,
                             OptimalTransferPlan plan,
                             OrbitReport rep_target,
                             RendezvousPlan rdv,
                             MonteCarloReport mc,
                             ADCSReport adcs,
                             ThermalReport tcs,
                             TelecomReport telecom,
                             SystemEngineeringReport sys,
                             NavigationReport nav,
                             CAMReport cam,
                             EndOfLifeReport eol,
                             CislunarReport cislunar,
                             StructuralReport csi,
                             FDIRReport fdir,
                             const InterplanetaryReport *ipl,
                             const GroundStation *station,
                             const SpacecraftConfig *sc);

#endif // TELEMETRY_H
