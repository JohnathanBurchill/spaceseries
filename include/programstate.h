/*

    spaceseries: programstate.h

    Copyright (C) 2024  Johnathan K Burchill

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

#ifndef _PROGRAMSTATE_H
#define _PROGRAMSTATE_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct SwarmData
{
    char *satellite;
    bool initialized;
    size_t nMeasurements;

    double *times; // Start as 2 Hz
    float *latitude;
    float *longitude;
    float *radius;

    float *qdlatitude;
    float *mlt;
    float *qdargoforbit;

    float *vsatn;
    float *vsate;
    float *vsatc;

    float *parameter;

    float staticVal;

    // derived stuff
    float *vsat;
    float *distanceAlongTrack;

    // Helper stuff
    ssize_t closestTimeIndex;
    // Current interpolation indices
    long i1;
    long i2;

} SwarmData;

typedef struct Data 
{

    bool initialized;
    int nSatellites;
    SwarmData *satelliteData[3];
    float initialDistanceLag[3];
    double initialTimeLag[3];
    // Space series data
    float *spaceSeries;
    float *previousSpaceSeries;
    float *latitude;
    float *longitude;
    float *radius;
    float *qdlatitude;
    float *mlt;
    float *qdargoforbit;
    unsigned char *outputBuffer;
    size_t outputBufferSize;

} Data;

typedef struct Ephemeres
{
    double time;
    float latitude;
    float longitude;
    float radius;
    float vsatN;
    float vsatE;
    float vsatC;
} Ephemeres;


typedef struct Directories
{
    char *tctDir;
    char *c7hDir;
    char *modDir;
    char *lpDir;
    char *outDir;
} Directories;

typedef struct programState
{
    Data data;

    char *measurementParameter;
    char *startString;
    char *stopString;
    double firstTime;
    double lastTime;

    int nOptions;

    double samplePeriod;

    Directories directories;
    char *qdDataFile;

    bool binaryOutput;
    int nOutputFiles;
    char *suffix;

    float parameterOffsets[3];
    float scaleFactor;

    bool alongTrackGradient;

    bool constantInterpolation;

    bool splitOutputFilesByOrbit;
    bool labelOutputFileOrbits;
    char outFile[FILENAME_MAX];
    FILE *output;
    int outputfd;

    float minDistance;
    float maxDistance;

    double scanTime;
    double scanDistance;

    char *tctDataset;

} ProgramState;

void initProgramState(ProgramState *state);


#endif // _PROGRAMSTATE_H
