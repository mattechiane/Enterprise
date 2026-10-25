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

// Esporta l'intero scenario di missione (Orbita 1 -> Trasferimento -> Orbita 2)
// in formato JSON interattivo per la Web App 3D e i grafici di trade-off
bool export_web_mission_json(const char *filename,
                             OrbitReport rep_init,
                             OptimalTransferPlan plan,
                             OrbitReport rep_target,
                             const GroundStation *station,
                             const SpacecraftConfig *sc);

#endif // TELEMETRY_H