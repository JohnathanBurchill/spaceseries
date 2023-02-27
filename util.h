/*

    spaceseries: util.h

    Copyright (C) 2023  Johnathan K Burchill

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef _UTIL_H
#define _UTIL_H

#include "data.h"

void latLonToXYZ(float lat, float lon, float *x, float *y, float *z);
void xyzToLatLon(float x, float y, float z, float *lat, float *lon);
void interSatelliteDistanceKm(SwarmData *sat1, size_t n1, SwarmData *sat2, size_t n2, float *distance, float *alongTrackDisplacement, float *horizontalCrossTrackDisplacement, float *verticalCrossTrackDisplacement);

int getOptionValue(char *number, float *value);

int setOutputFilename(ProgramState *state);
int openOutputFile(ProgramState *state);
int closeOutputFile(ProgramState *state);

#endif // _UTIL_H
