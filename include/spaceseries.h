/*

    spaceseries: spaceseries.h

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

#ifndef _SPACESERIES_H
#define _SPACESERIES_H

#include "programstate.h"

#include <stdint.h>

#define SPACESERIES_MAX_SATELLITE_NAME_LEN 255
#define SPACESERIES_N_TIMESERIES_VALUES 11

enum SPACESERIES_ERRORS
{
    SPACESERIES_OK = 0,
    SPACESERIES_MEMORY,
    SPACESERIES_CDF_READ,
    SPACESERIES_NOT_ENOUGH_MEASUREMENTS,
    SPACESERIES_NOT_ENOUGH_MAG_MEASUREMENTS,
    SPACESERIES_TIME_RANGE,
    SPACESERIES_EPHEM_READ,
    SPACESERIES_PARSE_NUMBER,
    SPACESERIES_POINTERS,
    SPACESERIES_FILE_READ,
    SPACESERIES_FILE_WRITE,
    SPACESERIES_NO_MORE_DATA,
    SPACESERIES_ARGUMENTS
};

typedef struct timeSeriesPoint
{
    float distanceAlongTrack;
    float latitude;
    float longitude;
    float radius;
    float qdLatitude;
    float mlt;
    float qdArgOfOrbit;
    float alongTrackDisplacement;
    float horizontalCrossTrackDisplacement;
    float verticalCrossTrackDisplacement;
    float staticParam;

} TimeSeriesPoint;

typedef struct spaceSeriesHeader
{
    double timeStamp;
    int nSpaceSeriesPoints;
    int nTimeSeriesPoints;
    TimeSeriesPoint *satelliteTimeSeriesPoint;
} SpaceSeriesHeader;

typedef struct spaceSeriesDataPoint
{
    float distanceAlongTrack;
    float latitude;
    float longitude;
    float radius;
    float qdLatitude;
    float mlt;
    float qdArgOfOrbit;
    float param;
    float dparamdt;
    int nErrors;
    float *staticAssumptionErrors;
} SpaceSeriesDataPoint;

typedef struct SpaceSeries
{
    int nDimensions;
    SpaceSeriesHeader header;
    SpaceSeriesDataPoint *points;
} SpaceSeries;

typedef struct satelliteInfo
{
    char *name;
    double timeLag;
    float distanceLag;
} SatelliteInfo;

typedef struct spaceSeriesScanHeader
{
    int nSatellites;
    SatelliteInfo *satellites;
} SpaceSeriesScanHeader;

int generateSpaceSeries(ProgramState *state);

int readSpaceSeriesScanHeader(int binaryfiledescriptor, SpaceSeriesScanHeader *header);
char *readSatelliteName(int binaryfiledescriptor);

int readSpaceSeriesHeader(int binaryfiledescriptor, SpaceSeriesHeader *header);

int readSpaceSeries(int binaryfiledescriptor, SpaceSeries *series);


#endif // _SPACESERIES_H
