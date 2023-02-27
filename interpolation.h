/*

    spaceseries: interpolation.h

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

#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

#include "data.h"

#include <stdbool.h>
#include <stdlib.h>

ssize_t closestTimeIndex(SwarmData *swarmData, double time);
void neighboringDistanceIndices(float distance, SwarmData *sat);
void interpolate(Data *data, double tMilliseconds, float requestedDistance, bool constantInterpolation, float *parameter, float *latitude, float *longitude, float *radius, float *qdlat, float *mlt, float *qdargoforbit);
float interpolateAtDistance(Data *data, int satelliteIndex, float requestedDistance, double *time, float *parameter, float * latitude, float *longitude, float *radius, float *qdlat, float *mlt, float *qdargoforbit);
float calculateGradient(float *distances, float *measurements, size_t nMeasurements, int i, int minDistanceIndex);


#endif // _INTERPOLATION_H
