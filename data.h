/*

    spaceseries: data.h

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

#ifndef _DATA_H
#define _DATA_H

#include "programstate.h"

#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#include <cdf.h>

#define OUTPUT_BUFFER_START_SIZE 1024
#define HEADER_NUMBER_OF_FLOATS 13
#define WINDOW_NUMBER_OF_FLOATS 15

int loadData(ProgramState *state);
int loadParameterData(ProgramState *state, char *measurementParameter);
int loadSatelliteData(ProgramState *state, char *measurementParameter, int satellite);

size_t downSample(double newPeriod, double *oldTime, size_t n, double *time, size_t nMeasurements, float **parameter);
void getRecordRange(char *timeVariableName, CDFid cdf, double firstTime, double lastTime, long *firstRec, long *lastRec);
void loadVariable(char *varName, CDFid cdf, long firstRec, long lastRec, void **data, size_t *nMeasurements);
size_t bytesPerRecord(char *varName, CDFid cdf);

void freeData(ProgramState *state);
void freeSatelliteData(SwarmData *data);

int satelliteOrder(const void *first, const void *second);

int loadEphemeres(SwarmData *swarmData, char *modDir, double firstTime, double lastTime, Ephemeres **ephem, size_t *nEphemeres);

void dateTimeStrings(double time, char *dateString, char *timeString);

void calculateAlongTrackDistances(ProgramState *state);
void calculateMagneticLatitudeArgOfOrbit(ProgramState *state);

#endif // _DATA_H
