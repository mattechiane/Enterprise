#ifndef TELEMETRY_H
#define TELEMETRY_H

#include <stdio.h>
#include "orbit_types.h"

// Converte la posizione da Inerziale (ECI) a Solidale con la Terra (ECEF)
// gmst0 è l'angolo di Greenwich a t = 0 (radianti)
Vector3D eci_to_ecef(Vector3D pos_eci, double epoch, double gmst0);

// Calcola Latitudine Geodetica (sull'ellissoide schiacciato WGS84), Longitudine e Quota
GeodeticCoord ecef_to_geodetic(Vector3D pos_ecef);

// Funzioni per la gestione del file CSV di telemetria
FILE* telemetry_open(const char *filename);
void  telemetry_write_step(FILE *fp, CartesianState state, KeplerianElements kep, GeodeticCoord geo, SpaceEnvironment env, bool thruster_on);
void  telemetry_close(FILE *fp);

#endif // TELEMETRY_H