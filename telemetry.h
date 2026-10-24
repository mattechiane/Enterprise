#ifndef TELEMETRY_H
#define TELEMETRY_H

#include <stdio.h>
#include "orbit_types.h"

Vector3D eci_to_ecef(Vector3D pos_eci, double epoch, double gmst0);
GeodeticCoord ecef_to_geodetic(Vector3D pos_ecef);

FILE* telemetry_open(const char *filename);
void  telemetry_write_step(FILE *fp, CartesianState state, KeplerianElements kep, GeodeticCoord geo, SpaceEnvironment env, bool thruster_on);
void  telemetry_close(FILE *fp);

// Esporta l'intero scenario (Orbita 1 -> Trasferimento -> Orbita 2) in formato JSON per la Web App 3D
bool export_web_mission_json(const char *filename,
                             OrbitReport rep_init,
                             OptimalTransferPlan plan,
                             OrbitReport rep_target,
                             const GroundStation *station,
                             const SpacecraftConfig *sc);

#endif // TELEMETRY_H