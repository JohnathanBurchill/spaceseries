/*

    spaceseries: spaceseries.c

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

#include "spaceseries.h"

#include "interpolation.h"

#include "data.h"
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int generateSpaceSeries(ProgramState *state)
{
    int status = SPACESERIES_OK;

    float alongTrackDisplacement = 0.0;
    float horizontalCrossTrackDisplacement = 0.0;
    float verticalCrossTrackDisplacement = 0.0;
    float d0 = 0.0;
    float d1 = 0.0;
    float d2 = 0.0;
    float param = 0.0;
    float param1 = 0.0;
    float param2 = 0.0;
    float latitude = 0.0;
    float longitude = 0.0;
    float radius = 0.0;
    float qdlat = 0.0;
    float mlt = 0.0;
    float qdargoforbit = 0.0;
    char epoch[200] = {0};
    int minDistanceIndex = 0;
    double timeIntervalToDisplay = 60.0 * 1000; // 60 s
    double timePassed = 0.0;
    void *mem = NULL;

    ssize_t closestTimeInd = 0;

    SwarmData *sat = NULL;

    sat = state->data.satelliteData[0];
    double tStart = sat->times[0];
    sat = state->data.satelliteData[state->data.nSatellites-1];
    double tStop = sat->times[sat->nMeasurements-1];

    sat = state->data.satelliteData[0];

    for (int i = 0; i < sat->nMeasurements; i++)
        state->data.previousSpaceSeries[i] = nan("");

    float gradient = 0.0;
    float dparamdt = 0.0;
    float staticVal = 0.0;
    float staticVal1 = 0.0;
    float staticVal2 = 0.0;
    float diffVal = 0.0;
    float floatBuf = 0.0;
    float interpFrac = 0.0;
    int windowPoints = 0;
    size_t requiredBufferSize = 0;

    size_t nToWrite = 0;
    ssize_t nWritten = 0;

    SwarmData *tmpSat = NULL;

    unsigned char *obi = state->data.outputBuffer; // output buffer pointer

    // Estimate the space series at each time
    for (double time = tStart; time < tStop; time += 1000 * state->samplePeriod)
    {
        timePassed += state->samplePeriod;
        if ((timePassed >= 3600.0))
        {
            toEncodeEPOCH(time, 4, epoch);
            fprintf(stdout, "\r%s", epoch);
            fflush(stdout);
            timePassed = 0.0;
        }
        // Init some indexes into the data arrays
        for (int s = 0; s < state->data.nSatellites; s++)
        {
            closestTimeIndex(state->data.satelliteData[s], time);
            state->data.satelliteData[s]->i1 = (long)minDistanceIndex;
            state->data.satelliteData[s]->i2 = (long)minDistanceIndex;
        }

        // Option to split output files by orbit
        if (state->splitOutputFilesByOrbit && (sat->closestTimeIndex > 0 && (fmodf(sat->qdlatitude[sat->closestTimeIndex], 360.0) >= 0.0 && fmodf(sat->qdlatitude[sat->closestTimeIndex - 1], 360.0) < 0.)))
        {
            state->lastTime = time - 1000.0 * state->samplePeriod;
            status = closeOutputFile(state);
            if (status != SPACESERIES_OK)
            {
                fprintf(stderr, "Unable to close output file.\n");
                return SPACESERIES_FILE_WRITE;
            }
            state->firstTime = time;
            state->lastTime = tStop;
            status = openOutputFile(state);
            if (status != SPACESERIES_OK)
            {
                fprintf(stderr, "Unable to open output file for writing.\n");
                return SPACESERIES_FILE_WRITE;
            }
        }

        // First pass through the space series
        // Calculate basic space series for the requested parameter, or its anti-derivative
        // in the case of an along-track gradient

        // Keep track of the number of data points in the space series
        // Here the term "window" is a synonym for "space series"
        windowPoints = 0;
        for (int i = minDistanceIndex; i < sat->nMeasurements && sat->distanceAlongTrack[i] < state->maxDistance; i++)
        {
            // Getting measurements and ephemeres at positions interpolated between the three satellites
            // sat->distanceAlongTrack is that of the first satellite.
            d0 = sat->distanceAlongTrack[i];
            interpolate(&state->data, time, d0, state->constantInterpolation, &param, &latitude, &longitude, &radius, &qdlat, &mlt, &qdargoforbit);
            state->data.spaceSeries[i] = param;
            state->data.latitude[i] = latitude;
            state->data.longitude[i] = longitude;
            state->data.radius[i] = radius;
            state->data.qdlatitude[i] = qdlat;
            state->data.mlt[i] = mlt;
            state->data.qdargoforbit[i] = qdargoforbit;

            // Result will not be finite if the distance along track is not bounded by the space series
            if (isfinite(param))
                windowPoints++;
            if (isinff(param) && signbit(param) != 0)
            {
                // Leading satellite has not reached any of the following points. Can go to the next time
                break;
            }
        }

        // Write space series header information
        if (!state->binaryOutput)
            fprintf(state->output, "%lf %d", time, windowPoints);
        else
        {
            obi = state->data.outputBuffer;
            *((double*)obi) = time;
            obi += sizeof(double);
            *((int*)obi) = windowPoints;
            obi += sizeof(int);
        }
        d0 = state->data.satelliteData[0]->distanceAlongTrack[closestTimeInd];
        for (int s = 0; s < state->data.nSatellites; s++)
        {
            tmpSat = state->data.satelliteData[s];
            closestTimeInd = tmpSat->closestTimeIndex;
            // Satellite position at the current time, and the value obtained
            // under the static assumption at the current position (visited by
            // the satellite at another time), and the difference
            // between the static and dynamic values
            // Redundant displacement data for lead satellite
            interSatelliteDistanceKm(state->data.satelliteData[0], state->data.satelliteData[0]->closestTimeIndex, tmpSat, closestTimeInd, NULL, &alongTrackDisplacement, &horizontalCrossTrackDisplacement, &verticalCrossTrackDisplacement);

            // Static assumption value
            if (state->alongTrackGradient)
                staticVal = calculateGradient(tmpSat->distanceAlongTrack, tmpSat->parameter, tmpSat->nMeasurements, closestTimeInd, minDistanceIndex);
            else
                staticVal = tmpSat->parameter[closestTimeInd];

            if (!state->binaryOutput)
                fprintf(state->output, " %.3f %.3f %.1f %.1f %.2f %.2f %.2f %.2f %.2f %.2f %lg", tmpSat->latitude[closestTimeInd], tmpSat->longitude[closestTimeInd], tmpSat->radius[closestTimeInd], tmpSat->distanceAlongTrack[closestTimeInd], tmpSat->qdlatitude[closestTimeInd], tmpSat->mlt[closestTimeInd], fmodf(tmpSat->qdargoforbit[closestTimeInd], 360.0), alongTrackDisplacement, horizontalCrossTrackDisplacement, verticalCrossTrackDisplacement, staticVal);
            else
            {
                // For efficient binary output, store points in a buffer
                // and write out a full space series in one go
                *((float*)obi) = tmpSat->latitude[closestTimeInd];
                obi += sizeof(float);
                *((float*)obi) = tmpSat->longitude[closestTimeInd];
                obi += sizeof(float);
                *((float*)obi) = tmpSat->radius[closestTimeInd];
                obi += sizeof(float);
                *((float*)obi) = tmpSat->distanceAlongTrack[closestTimeInd];
                obi += sizeof(float);
                *((float*)obi) = tmpSat->qdlatitude[closestTimeInd];
                obi += sizeof(float);
                *((float*)obi) = tmpSat->mlt[closestTimeInd];
                obi += sizeof(float);
                *((float*)obi) = fmodf(tmpSat->qdargoforbit[closestTimeInd], 360.0);
                obi += sizeof(float);
                *((float*)obi) = alongTrackDisplacement;
                obi += sizeof(float);
                *((float*)obi) = horizontalCrossTrackDisplacement;
                obi += sizeof(float);
                *((float*)obi) = verticalCrossTrackDisplacement;
                obi += sizeof(float);
                *((float*)obi) = staticVal;
                obi += sizeof(float);
            }
            // Reset distance indices for pass 2
            tmpSat->i1 = minDistanceIndex;
            tmpSat->i2 = minDistanceIndex;
        }
        if (!state->binaryOutput)
            fprintf(state->output, "\n");
        else
        {
            nToWrite = obi - state->data.outputBuffer;
            nWritten = write(state->outputfd, state->data.outputBuffer, nToWrite);
            if (nWritten != (ssize_t)nToWrite)
            {
                fprintf(stderr, "Could not write binary file.\n");
                exit(EXIT_FAILURE);
            }
            // Update output Buffer size if needed
            requiredBufferSize = ((size_t)HEADER_NUMBER_OF_FLOATS + (size_t)WINDOW_NUMBER_OF_FLOATS) * (size_t)windowPoints * sizeof(float);
            if (state->data.outputBufferSize < requiredBufferSize)
            {
                mem = realloc(state->data.outputBuffer, requiredBufferSize);
                if (mem == NULL)
                {
                    printf("Could not allocate memory.\n");
                    exit(EXIT_FAILURE);
                }
                state->data.outputBuffer = mem;
                state->data.outputBufferSize = requiredBufferSize;
            }
            obi = state->data.outputBuffer;
        }
        
        // 2nd pass through the space series: estimate gradient (if requested)
        // estimate the total time derivative (in a frame fixed to the Earth)
        // Write results (if ASCII output) or store for binary write later
        for (int i = minDistanceIndex; i < sat->nMeasurements && sat->distanceAlongTrack[i] < state->maxDistance; i++)
        {
            d0 = sat->distanceAlongTrack[i];
            if (isfinite(state->data.spaceSeries[i]))
            {
                if (state->alongTrackGradient)
                    param = calculateGradient(sat->distanceAlongTrack, state->data.spaceSeries, sat->nMeasurements, i, minDistanceIndex);
                else
                    param = state->data.spaceSeries[i];

                dparamdt = (param - state->data.previousSpaceSeries[i]) / state->samplePeriod;
                state->data.previousSpaceSeries[i] = param;

                // Space series data point
                if (!state->binaryOutput)
                    fprintf(state->output, "%.2f %.3f %.3f %.1f %.2f %.2f %.2f %g %g", sat->distanceAlongTrack[i], state->data.latitude[i], state->data.longitude[i], state->data.radius[i], state->data.qdlatitude[i], state->data.mlt[i], fmodf(state->data.qdargoforbit[i], 360.0), param, dparamdt);
                else
                {
                    *((float*)obi) = sat->distanceAlongTrack[i];
                    obi += sizeof(float);
                    *((float*)obi) = state->data.latitude[i];
                    obi += sizeof(float);
                    *((float*)obi) = state->data.longitude[i];
                    obi += sizeof(float);
                    *((float*)obi) = state->data.radius[i];
                    obi += sizeof(float);
                    *((float*)obi) = state->data.qdlatitude[i];
                    obi += sizeof(float);
                    *((float*)obi) = state->data.mlt[i];
                    obi += sizeof(float);
                    *((float*)obi) = fmod(state->data.qdargoforbit[i], 360.0);
                    obi += sizeof(float);
                    *((float*)obi) = param;
                    obi += sizeof(float);
                    *((float*)obi) = dparamdt;
                    obi += sizeof(float);
                }

                // estimate difference from static assumption
                for (int s = 0; s < state->data.nSatellites; s++)
                {
                    tmpSat = state->data.satelliteData[s];
                    neighboringDistanceIndices(d0, tmpSat);
                    if (state->alongTrackGradient)
                    {
                        staticVal1 = calculateGradient(tmpSat->distanceAlongTrack, tmpSat->parameter, tmpSat->nMeasurements, tmpSat->i1, minDistanceIndex);
                        staticVal2 = calculateGradient(tmpSat->distanceAlongTrack, tmpSat->parameter, tmpSat->nMeasurements, tmpSat->i2, minDistanceIndex);
                    }
                    else
                    {
                        staticVal1 = tmpSat->parameter[tmpSat->i1];
                        staticVal2 = tmpSat->parameter[tmpSat->i2];
                    }
                    d1 = tmpSat->distanceAlongTrack[tmpSat->i1];
                    d2 = tmpSat->distanceAlongTrack[tmpSat->i2];
                    interpFrac = (d0 - d1) / (d2 - d1);
                    if (interpFrac <= 0.0)
                        interpFrac = 1.0;
                    staticVal = staticVal1 + (staticVal2 - staticVal1) * interpFrac;
                    diffVal = staticVal - param;

                    if (!state->binaryOutput)
                        fprintf(state->output, " %g", diffVal);
                    else
                    {
                        *((float*)obi) = diffVal;
                        obi += sizeof(float);
                    }
                }
                if (!state->binaryOutput)
                    fprintf(state->output, "\n");
            }
            else if (isinff(state->data.spaceSeries[i]) && signbit(state->data.spaceSeries[i]) != 0)
            {
                // Leading satellite has not reached any of the following points. Go to the next time
                break;
            }
            else if (isinff(state->data.spaceSeries[i]))
            {
                // Trailing satellite already visited this point, none of the points behind have values
                minDistanceIndex = i + 1;
            }
        }
        // Write the space series data points to binary file if requested
        if (state->binaryOutput)
        {
            nToWrite = obi - state->data.outputBuffer;
            nWritten = write(state->outputfd, state->data.outputBuffer, nToWrite);
            if (nWritten != nToWrite)
            {
                fprintf(stderr, "Could not write binary file.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    return status;
}


// Reads a binary space series file header 
int readSpaceSeriesScanHeader(int binaryfiledescriptor, SpaceSeriesScanHeader *header)
{
    if (binaryfiledescriptor == -1)
        return SPACESERIES_FILE_READ;

    ssize_t bytesRead = -1;

    bytesRead = read(binaryfiledescriptor, &header->nSatellites, sizeof header->nSatellites);
    if (bytesRead != sizeof header->nSatellites)
        return SPACESERIES_CDF_READ;

    header->satellites = calloc(sizeof *header->satellites, (size_t)header->nSatellites);
    if (header->satellites == NULL)
        return SPACESERIES_MEMORY;

    for (int i = 0; i < header->nSatellites; i++)
    {
        header->satellites[i].name = readSatelliteName(binaryfiledescriptor);
        if (header->satellites[i].name == NULL)
            return SPACESERIES_FILE_READ;
    }
    for (int i = 0; i < header->nSatellites; i++)
    {
        bytesRead = read(binaryfiledescriptor, &header->satellites[i].timeLag, sizeof header->satellites[i].timeLag);
        if (bytesRead != sizeof header->satellites[i].timeLag)
            return SPACESERIES_FILE_READ;

        bytesRead = read(binaryfiledescriptor, &header->satellites[i].distanceLag, sizeof header->satellites[i].distanceLag);
        if (bytesRead != sizeof header->satellites[i].distanceLag)
            return SPACESERIES_FILE_READ;
    }

    return SPACESERIES_OK;
}

// Get the 0-terminated satellite name from a binary file
char *readSatelliteName(int binaryfiledescriptor)
{
    int stringLength = 0;
    ssize_t bytesRead = -1;

    char name[SPACESERIES_MAX_SATELLITE_NAME_LEN+1] = {0};
    do
    {
        bytesRead = read(binaryfiledescriptor, name + stringLength, sizeof *name);
        if (bytesRead != sizeof *name)
            return NULL;
        stringLength += bytesRead;
    }
    while (name[stringLength - 1] != 0);

    if (stringLength == 0)
        return NULL;

    return strdup(name);

}

// Read header info for a space series
int readSpaceSeriesHeader(int binaryfiledescriptor, SpaceSeriesHeader *header)
{
    ssize_t bytesRead = 0;
    float floatBuf[11] = {0};
    int offset = 0;

    bytesRead = read(binaryfiledescriptor, &header->timeStamp, sizeof header->timeStamp);
    if (bytesRead != sizeof header->timeStamp)
        return bytesRead == 0 ? SPACESERIES_NO_MORE_DATA : SPACESERIES_FILE_READ;

    bytesRead = read(binaryfiledescriptor, &header->nSpaceSeriesPoints, sizeof header->nSpaceSeriesPoints);
    if (bytesRead != sizeof header->nSpaceSeriesPoints)
        return bytesRead == 0 ? SPACESERIES_NO_MORE_DATA : SPACESERIES_FILE_READ;

    for (int i = 0; i < header->nTimeSeriesPoints; i++)
    {

        bytesRead = read(binaryfiledescriptor, floatBuf, SPACESERIES_N_TIMESERIES_VALUES * sizeof *floatBuf);
        if (bytesRead != SPACESERIES_N_TIMESERIES_VALUES * sizeof *floatBuf)
            return bytesRead == 0 ? SPACESERIES_NO_MORE_DATA : SPACESERIES_FILE_READ;

        offset = 0;
        header->satelliteTimeSeriesPoint[i].latitude = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].longitude = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].radius = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].distanceAlongTrack = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].qdLatitude = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].mlt = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].qdArgOfOrbit = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].alongTrackDisplacement = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].horizontalCrossTrackDisplacement = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].verticalCrossTrackDisplacement = floatBuf[offset++];
        header->satelliteTimeSeriesPoint[i].staticParam = floatBuf[offset++];
    }

    return SPACESERIES_OK;
}

// Read space series data points
int readSpaceSeries(int binaryfiledescriptor, SpaceSeries *series)
{
    int status = SPACESERIES_OK;

    void *mem = NULL;

    int currentSpaceSeriesNPoints = series->header.nSpaceSeriesPoints;
    status = readSpaceSeriesHeader(binaryfiledescriptor, &series->header);
    if (status != SPACESERIES_OK)
        return status;

    if (series->header.nSpaceSeriesPoints > currentSpaceSeriesNPoints)
    {
        mem = realloc(series->points, series->header.nSpaceSeriesPoints * sizeof *series->points);
        if (mem == NULL)
            return SPACESERIES_MEMORY;
        series->points = mem;
        for (int i = currentSpaceSeriesNPoints; i < series->header.nSpaceSeriesPoints; i++)
            series->points[i].staticAssumptionErrors = NULL;
        for (int i = 0; i < series->header.nSpaceSeriesPoints; i++)
        {
            mem = realloc(series->points[i].staticAssumptionErrors, series->header.nTimeSeriesPoints * sizeof *series->points[i].staticAssumptionErrors);
            if (mem == NULL)
                return SPACESERIES_MEMORY;

            series->points[i].staticAssumptionErrors = mem;
        }
        
    }

    ssize_t bytesRead = 0;

    float floatBuf[9] = {0};
    int offset = 0;
    int nBytesToRead = 0;

    nBytesToRead = series->header.nTimeSeriesPoints * sizeof *series->points[0].staticAssumptionErrors;
    for (int i = 0; i < series->header.nSpaceSeriesPoints; i++)
    {
        bytesRead = read(binaryfiledescriptor, floatBuf, 9 * sizeof *floatBuf);
        if (bytesRead == 0)
            return SPACESERIES_NO_MORE_DATA;
        if (bytesRead != 9 * sizeof *floatBuf)
            return SPACESERIES_FILE_READ;

        offset = 0;
        series->points[i].distanceAlongTrack = floatBuf[offset++];
        series->points[i].latitude = floatBuf[offset++];
        series->points[i].longitude = floatBuf[offset++];
        series->points[i].radius = floatBuf[offset++];
        series->points[i].qdLatitude = floatBuf[offset++];
        series->points[i].mlt = floatBuf[offset++];
        series->points[i].qdArgOfOrbit = floatBuf[offset++];
        series->points[i].param = floatBuf[offset++];
        series->points[i].dparamdt = floatBuf[offset++];
        bytesRead = read(binaryfiledescriptor, series->points[i].staticAssumptionErrors, nBytesToRead);
        if (bytesRead == 0)
            return SPACESERIES_NO_MORE_DATA;
        if (bytesRead != nBytesToRead)
            return SPACESERIES_FILE_READ;
    }

    return status;
}