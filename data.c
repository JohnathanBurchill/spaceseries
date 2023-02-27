/*

    spaceseries: data.c

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

#include "data.h"
#include "util.h"
#include "spaceseries.h"
#include "settings.h"

#include <stdio.h>
#include <stdlib.h>
#include <fts.h>
#include <string.h>
#include <math.h>

#include <cdf.h>
#include <quasidipole.h>

int loadData(ProgramState *state)
{
    int status = SPACESERIES_OK;

    if (strcmp("jz", state->measurementParameter) == 0)
    {
        state->alongTrackGradient = true;
        state->scaleFactor = -0.7957747; // µA/m^2
        status = loadParameterData(state, "dby");
    }
    else
        status = loadParameterData(state, state->measurementParameter);

    return status;

}

int loadParameterData(ProgramState *state, char *measurementParameter)
{

    if (state == NULL || state->parameterOffsets == NULL || measurementParameter == NULL || state->directories.c7hDir == NULL || state->directories.modDir == NULL || state->directories.tctDir == NULL || state->directories.lpDir == NULL)
        return SPACESERIES_POINTERS;

    // TII cross-track ion drift and magnetic field residual files should be
    // in the working directory.

    // Load data from as many of the satellites as possible

    int status = SPACESERIES_OK;
    char *satelliteNames[3] = {"Swarm A", "Swarm B", "Swarm C"};
    float alongTrackDisplacement = 0.0;

    for (int s = 0; s < 3; s++)
    {
        state->data.satelliteData[s] = (SwarmData *)calloc(1, sizeof(SwarmData));
        if (state->data.satelliteData[s] == NULL)
        {
            fprintf(stderr, "Unable to allocate satellite data memory.\n");
            return SPACESERIES_MEMORY;
        }
        state->data.satelliteData[s]->satellite = strdup(satelliteNames[s]);
        status = loadSatelliteData(state, measurementParameter, s);
        if (status != SPACESERIES_OK)
            return status;

    }

    for (int i = 0; i < 3; i++)
    {
        state->data.nSatellites += (state->data.satelliteData[i]->nMeasurements > 1);
    }
    state->data.initialized = state->data.nSatellites > 1;

    if (state->data.initialized)
    {
        // Sort all satellite data pointers, even if some data are missing
        qsort((void*)&(state->data.satelliteData[0]), 3, sizeof(SwarmData*), &satelliteOrder);
        for (int s = 0; s < state->data.nSatellites; s++)
        {
            interSatelliteDistanceKm(state->data.satelliteData[s], 0, state->data.satelliteData[0], 0, NULL, &alongTrackDisplacement, NULL, NULL);
            state->data.initialDistanceLag[s] = alongTrackDisplacement;
            state->data.initialTimeLag[s] = (1000.0 * alongTrackDisplacement) / (float)state->data.satelliteData[s]->vsat[0];

            // Remove requested parameter offsets
            if (state->parameterOffsets[s] != 0.0)
                for (size_t i = 0; i < state->data.satelliteData[s]->nMeasurements; i++)
                    state->data.satelliteData[s]->parameter[i] -= state->parameterOffsets[s];

            // Apply scale factor
            if (state->scaleFactor != 1.0)
                for (size_t i = 0; i < state->data.satelliteData[s]->nMeasurements; i++)
                    state->data.satelliteData[s]->parameter[i] *= state->scaleFactor;

        }

    }

    // Space series is calculated along the track of the leading satellite
    // Buffer for estimating time derivatives
    state->data.spaceSeries = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.spaceSeries);
    state->data.previousSpaceSeries = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.previousSpaceSeries);
    state->data.outputBuffer = calloc(OUTPUT_BUFFER_START_SIZE, sizeof *state->data.outputBuffer);
    state->data.outputBufferSize = OUTPUT_BUFFER_START_SIZE * sizeof *state->data.outputBuffer;
    state->data.latitude = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.latitude);
    state->data.longitude = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.longitude);
    state->data.radius = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.radius);
    state->data.qdlatitude = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.qdlatitude);
    state->data.mlt = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.mlt);
    state->data.qdargoforbit = calloc(state->data.satelliteData[0]->nMeasurements, sizeof *state->data.qdargoforbit);

    if (state->data.spaceSeries == NULL || state->data.previousSpaceSeries == NULL || state->data.latitude == NULL || state->data.longitude == NULL || state->data.radius == NULL || state->data.qdlatitude == NULL || state->data.mlt == NULL)
        status = SPACESERIES_MEMORY;

    return status;

}

int loadSatelliteData(ProgramState *state, char *measurementParameter, int satellite)
{
    if (state == NULL || measurementParameter == NULL)
        return SPACESERIES_POINTERS;

    if (satellite < 0)
        return SPACESERIES_ARGUMENTS;

    int status = SPACESERIES_OK;

    char firstDateStr[9] = {0};
    char lastDateStr[9] = {0};

    char firstTimeStr[7] = {0};
    char lastTimeStr[7] = {0};

    long year = 0;
    long month = 0;
    long day = 0;
    long hour = 0;
    long minute = 0;
    long second = 0;
    long msec = 0;

    double unixTime = 0;

    EPOCHbreakdown(state->firstTime, &year, &month, &day, &hour, &minute, &second, &msec);
    sprintf(firstDateStr, "%4ld%02ld%02ld", year, month, day);
    sprintf(firstTimeStr, "%02ld%02ld%02ld", hour, minute, second);
    EPOCHbreakdown(state->lastTime, &year, &month, &day, &hour, &minute, &second, &msec);
    sprintf(lastDateStr, "%4ld%02ld%02ld", year, month, day);
    sprintf(lastTimeStr, "%02ld%02ld%02ld", hour, minute, second);

    // Directory information
    char* dir[2] = {state->directories.tctDir, NULL};
    // Load TII CT data
    FTSENT *e = NULL;

    FTS *fts = fts_open(dir, FTS_LOGICAL | FTS_NOCHDIR | FTS_NOSTAT, NULL);

    CDFid cdf = NULL;
    CDFstatus cdfstatus = CDF_OK;
    bool gotEphemeres = false;
    bool gotParameter = false;
    size_t nValues = 0;

    long firstRec = 0;
    long lastRec = 0;

    bool magParameter = strcmp("dbx", measurementParameter) == 0 || strcmp("dby", measurementParameter) == 0 || strcmp("dbz", measurementParameter) == 0;

    SwarmData *swarmData = state->data.satelliteData[satellite];

    e = fts_read(fts);
    while (e)
    {
        if ((e->fts_namelen == 59) && (strcmp(".cdf", e->fts_name + e->fts_namelen - 4) == 0) && (*(e->fts_name + 11) == swarmData->satellite[strlen(swarmData->satellite)-1]) && (strncmp("_TCT02", e->fts_name + 12, 6) == 0) && (strncmp(firstDateStr, e->fts_name + 19, 8) == 0) && (strncmp(firstTimeStr, e->fts_name+28, 6)>=0) && (strncmp(lastTimeStr, e->fts_name+44, 6)<=0))
        {
            cdfstatus = CDFopen(e->fts_path, &cdf);
            if (cdfstatus != CDF_OK)
            {
                fprintf(stderr, "Could not open CDF %s\n", e->fts_path);
                fts_close(fts);
                return SPACESERIES_CDF_READ;
            }
            // Errors reading variables are indicated by setting nMeasurements to -1
            getRecordRange("Timestamp", cdf, state->firstTime, state->lastTime, &firstRec, &lastRec);
            loadVariable("Timestamp", cdf, firstRec, lastRec, (void**)&(swarmData->times), &nValues);
            loadVariable("Latitude", cdf, firstRec, lastRec, (void**)&(swarmData->latitude), &nValues);
            loadVariable("Longitude", cdf, firstRec, lastRec, (void**)&(swarmData->longitude), &nValues);
            loadVariable("Radius", cdf, firstRec, lastRec, (void**)&(swarmData->radius), &nValues);
            loadVariable("QDLatitude", cdf, firstRec, lastRec, (void**)&(swarmData->qdlatitude), &nValues);
            loadVariable("MLT", cdf, firstRec, lastRec, (void**)&(swarmData->mlt), &nValues);
            loadVariable("VsatN", cdf, firstRec, lastRec, (void**)&(swarmData->vsatn), &nValues);
            loadVariable("VsatE", cdf, firstRec, lastRec, (void**)&(swarmData->vsate), &nValues);
            loadVariable("VsatC", cdf, firstRec, lastRec, (void**)&(swarmData->vsatc), &nValues);

            // If the parameter is a TCT param, read it too
            if (CDFconfirmzVarExistence(cdf, measurementParameter) == CDF_OK)
            {
                loadVariable(measurementParameter, cdf, firstRec, lastRec, (void**)&(swarmData->parameter), &(swarmData->nMeasurements));
                gotParameter = true;
            }
            CDFclose(cdf);

            break;
        }

        e = fts_read(fts);
    }
    fts_close(fts);

    if (swarmData->nMeasurements == -1)
    {
        fprintf(stderr, "Error loading TCT data.\n");
        return SPACESERIES_CDF_READ;
    }

    if (gotParameter && swarmData->nMeasurements < 2)
        return SPACESERIES_OK;
    
    // If there are no TII data and we are reading a mag residual parameter,
    // read in ephemeres from MOD file and calculate QD Lat and MLT  
    // and interpolate them to 2 Hz at the UT second and UT second + 0.5 s.
    if (!gotEphemeres)
    {
        Ephemeres *ephem = NULL;
        size_t nEphem = 0;
        // MOD files do not start at UT 00:00:00. Load two days' worth of measurements
        // Load previous day
        status = loadEphemeres(swarmData, state->directories.modDir, state->firstTime-86400000, state->lastTime-86400000, &ephem, &nEphem);
        // Try anyway for the requested day even on an error.
        status = loadEphemeres(swarmData, state->directories.modDir, state->firstTime, state->lastTime, &ephem, &nEphem);

        if (status != SPACESERIES_OK || nEphem == 0)
        {
            if (ephem != NULL)
                free(ephem);
            return status;
        }
        // Interpolate to times and store data
        double sampleIntervalMs = 500; // 2 Hz like TCT dataset
        ssize_t nSamples = (ssize_t)floor((state->lastTime - state->firstTime) / sampleIntervalMs);
        if (nSamples <= 1)
        {
            if (ephem != NULL)
                free(ephem);
            return SPACESERIES_TIME_RANGE;
        }
        swarmData->times = calloc((size_t)nSamples, sizeof(double));
        swarmData->latitude = calloc((size_t)nSamples, sizeof(float));
        swarmData->longitude = calloc((size_t)nSamples, sizeof(float));
        swarmData->radius = calloc((size_t)nSamples, sizeof(float));
        swarmData->qdlatitude = calloc((size_t)nSamples, sizeof(float));
        swarmData->mlt = calloc((size_t)nSamples, sizeof(float));
        swarmData->vsatn = calloc((size_t)nSamples, sizeof(float));
        swarmData->vsate = calloc((size_t)nSamples, sizeof(float));
        swarmData->vsatc = calloc((size_t)nSamples, sizeof(float));
        if (swarmData->times == NULL || swarmData->latitude == NULL || swarmData->longitude == NULL || swarmData->radius == NULL || swarmData->qdlatitude == NULL || swarmData->mlt == NULL || swarmData->vsatn == NULL || swarmData->vsate == NULL || swarmData->vsatc == NULL)
        {
            if (ephem != NULL)
                free(ephem);
            return SPACESERIES_MEMORY;
        }

        for (ssize_t i = 0; i < nSamples; i++)
            swarmData->times[i] = state->firstTime + (double)i * sampleIntervalMs;

        // Interpolate other ephemeres
        double lastEphemTime = ephem[0].time;
        double nextEphemTime = ephem[0].time;
        ssize_t lastEphemIndex = 0;
        ssize_t ephemIndex = 0;
        ssize_t timeIndex = 0;
        double dt = 0.0;
        float interpFraction = 0.0;
        float x0 = 0.0;
        float y0 = 0.0;
        float z0 = 0.0;
        float x1 = 0.0;
        float y1 = 0.0;
        float z1 = 0.0;
        float x2 = 0.0;
        float y2 = 0.0;
        float z2 = 0.0;
        float lat2 = 0.0;
        float lon2 = 0.0;        
        double lat = 0.0;
        double lon = 0.0;
        double alt = 0.0;
        double qdlat = 0.0;
        double qdlon = 0.0;
        double mlt = 0.0;
        for (ssize_t timeIndex = 0; timeIndex < nSamples; timeIndex++)
        {

            while (ephemIndex < nEphem && timeIndex < nSamples && swarmData->times[timeIndex] > ephem[ephemIndex].time)
            {
                lastEphemIndex = ephemIndex++;
                lastEphemTime = ephem[lastEphemIndex].time;
            }

            if (ephemIndex == nEphem)
                ephemIndex--;
            nextEphemTime = ephem[ephemIndex].time;
            dt = nextEphemTime - lastEphemTime;
            if (dt < 500) // ephems are from 1000 ms samples
                dt = 1.0; // arbitrary
            interpFraction = (float)(swarmData->times[timeIndex] - lastEphemTime) / dt;

            // Calculate ephem interpolations
            latLonToXYZ(ephem[lastEphemIndex].latitude, ephem[lastEphemIndex].longitude, &x0, &y0, &z0);
            latLonToXYZ(ephem[ephemIndex].latitude, ephem[ephemIndex].longitude, &x1, &y1, &z1);
            x2 = x0 + (x1 - x0) * interpFraction;
            y2 = y0 + (y1 - y0) * interpFraction;
            z2 = z0 + (z1 - z0) * interpFraction;
            xyzToLatLon(x2, y2, z2, &lat2, &lon2);
            swarmData->latitude[timeIndex] = lat2;
            swarmData->longitude[timeIndex] = lon2;
            swarmData->radius[timeIndex] = ephem[lastEphemIndex].radius + (ephem[ephemIndex].radius - ephem[lastEphemIndex].radius) * interpFraction;

            swarmData->vsatn[timeIndex] = ephem[lastEphemIndex].vsatN + (ephem[ephemIndex].vsatN - ephem[lastEphemIndex].vsatN) * interpFraction;
            swarmData->vsate[timeIndex] = ephem[lastEphemIndex].vsatE + (ephem[ephemIndex].vsatE - ephem[lastEphemIndex].vsatE) * interpFraction;
            swarmData->vsatc[timeIndex] = ephem[lastEphemIndex].vsatC + (ephem[ephemIndex].vsatC - ephem[lastEphemIndex].vsatC) * interpFraction;
            
            // QD Latitude, MLT
            EPOCHtoUnixTime(&swarmData->times[timeIndex], &unixTime, 1);
            lat = (double)swarmData->latitude[timeIndex];
            lon = (double)swarmData->longitude[timeIndex];
            alt = ((double)swarmData->radius[timeIndex] - 6371200.0) / 1000.0;
            geographicToQuasiDipole(state->qdDataFile, unixTime, lat, lon, alt, &qdlat, &qdlon);
            quasiDipoleMagneticLocalTime(state->qdDataFile, unixTime, qdlat, qdlon, &mlt);
            swarmData->qdlatitude[timeIndex] = (float)qdlat;
            swarmData->mlt[timeIndex] = (float)mlt;
        }
        
        nValues = nSamples;

        if (ephem != NULL)
            free(ephem);

    }

    // Allocate memory for derived parameters
    swarmData->vsat = (float*) calloc(nValues, sizeof(float));
    swarmData->distanceAlongTrack = (float*) calloc(nValues, sizeof(float));
    swarmData->qdargoforbit = (float *)calloc(nValues, sizeof(float));
    if (swarmData->distanceAlongTrack == NULL || swarmData->distanceAlongTrack == NULL || swarmData->qdargoforbit == NULL)
    {
        fprintf(stderr, "Could not allocate memory for derived parameters.\n");
        return SPACESERIES_MEMORY;
    }
    for (size_t i = 0; i < nValues; i++)
    {
        swarmData->vsat[i] = sqrtf(swarmData->vsatn[i]*swarmData->vsatn[i] + swarmData->vsate[i]*swarmData->vsate[i] + swarmData->vsatc[i]*swarmData->vsatc[i]);
    }

    // if the parameter to read is a mag residual parameter, read it now
    if (!gotParameter && magParameter)
    {
        dir[0] = state->directories.c7hDir;
        fts = fts_open(dir, FTS_LOGICAL | FTS_NOCHDIR | FTS_NOSTAT, NULL);
        e = fts_read(fts);
        while (e)
        {
            // Mag residuals may be processed in short intervals, giving multiple files per day. Check that the time range covers the requested interval
            if ((e->fts_namelen == 59) && (strcmp(".cdf", e->fts_name + e->fts_namelen - 4) == 0) && (*(e->fts_name + 11) == swarmData->satellite[strlen(swarmData->satellite)-1]) && (strncmp("C7H_2_", e->fts_name + 12, 6) == 0) && (strncmp(firstDateStr, e->fts_name + 19, 8) == 0) && (strncmp(firstTimeStr, e->fts_name+28, 6)>=0) && (strncmp(lastTimeStr, e->fts_name+44, 6)<=0))
            {
                cdfstatus = CDFopen(e->fts_path, &cdf);
                if (cdfstatus != CDF_OK)
                {
                    fprintf(stderr, "Could not open CDF %s\n", e->fts_path);
                    fts_close(fts);
                    return SPACESERIES_CDF_READ;
                }
                // Data are high res: load then resample
                double *btime = NULL;
                double *db = NULL; // 3-vectors
                size_t n = 0;
                // Get 0.5 s extra at each end of interval.
                getRecordRange("Timestamp", cdf, state->firstTime-500, state->lastTime+500, &firstRec, &lastRec);
                loadVariable("Timestamp", cdf, firstRec, lastRec, (void**)&btime, &n);
                loadVariable("dB_nec", cdf, firstRec, lastRec, (void**)&db, &n);
                CDFclose(cdf);

                if (btime == NULL || db == NULL)
                {
                    fts_close(fts);
                    return SPACESERIES_MEMORY;
                }

                // Might be able to continue without measurements from this satellite
                if (n == 0 || btime[0] > state->firstTime || btime[n-1] < state->lastTime)
                {
                    fts_close(fts);
                    free(btime);
                    free(db);
                    return SPACESERIES_OK;
                }


                // Transform dB from NEC to Satellite-track-aligned frame

                float *dbn = calloc(n, sizeof(float));
                float *dbe = calloc(n, sizeof(float));
                float *dbc = calloc(n, sizeof(float));
                swarmData->parameter = calloc(n, sizeof(float));
                if (dbn == NULL || dbe == NULL || dbc == NULL || swarmData->parameter == NULL)
                {
                    fprintf(stderr, "Could not allocate memory.\n");
                    fts_close(fts);
                    return SPACESERIES_MEMORY;
                }
                for (size_t i = 0; i < n; i++)
                {
                    dbn[i] = (float) (db[i*3]);
                    dbe[i] = (float) (db[i*3 + 1]);
                    dbc[i] = (float) (db[i*3 + 2]);
                }

                // Resample by averaging 50 Hz to 2 Hz centered on each TII time
                int ndbn = downSample(0.5, btime, n, swarmData->times, nValues, &dbn);
                int ndbe = downSample(0.5, btime, n, swarmData->times, nValues, &dbe);
                int ndbc = downSample(0.5, btime, n, swarmData->times, nValues, &dbc);

                if (ndbn != nValues || ndbe != nValues || ndbc != nValues)
                {
                    fprintf(stderr, "Resampled residual field does not have same number of measuremets as TCT data.\n");
                    free(btime);
                    free(db);
                    free(dbn);
                    free(dbe);
                    free(dbc);
                    fts_close(fts);
                    return SPACESERIES_OK;
                }

                // transform
                float v = 0.0;
                float xn = 0.0;
                float xe = 0.0;
                float xc = 0.0;
                float yn = 0.0;
                float ye = 0.0;
                float yc = 0.0;
                float zn = 0.0;
                float ze = 0.0;
                float zc = 0.0;
                float ymag = 0.0;
                float zmag = 0.0;
                for (size_t i = 0; i < nValues; i++)
                {
                    v = swarmData->vsat[i];
                    xn = swarmData->vsatn[i] / v;
                    xe = swarmData->vsate[i] / v;
                    xc = swarmData->vsatc[i] / v;
                    yn = xe * -1.0;
                    ye = -xn * -1.0;
                    yc = 0.0;
                    ymag = sqrtf(yn * yn + ye * ye + yc * yc);
                    yn /= ymag;
                    ye /= ymag;
                    yc /= ymag;
                    zn = xe * yc - ye * xc;
                    ze = -xn * yc + xc * yn;
                    zc = xn * ye - xe * yn;
                    zmag = sqrtf(zn * zn + ze * ze + zc * zc);
                    zn /= zmag;
                    ze /= zmag;
                    zc /= zmag;
                    switch(measurementParameter[2])
                    {
                        case 'x':
                            swarmData->parameter[i] = xn * dbn[i] + xe * dbe[i] + xc * dbc[i];
                            break;
                        case 'y':
                            swarmData->parameter[i] = yn * dbn[i] + ye * dbe[i] + yc * dbc[i];
                            break;
                        case 'z':
                            swarmData->parameter[i] = zn * dbn[i] + ze * dbe[i] + zc * dbc[i];
                            break;
                        default:
                            swarmData->parameter[i] = NAN;
                    }
                }

                gotParameter = true;
                swarmData->nMeasurements = nValues;
                free(btime);
                free(db);
                free(dbn);
                free(dbe);
                free(dbc);

                break;
            }

            e = fts_read(fts);
        }

        fts_close(fts);

        if (swarmData->nMeasurements == -1)
        {
            fprintf(stderr, "Error loading TCT data.\n");
            return SPACESERIES_CDF_READ;
        }
        else if (swarmData->nMeasurements == 0)
        {
            return SPACESERIES_OK;
        }
    }

    // if the parameter to read is an LP EXTD parameter, read it now
    if (!gotParameter)
    {
        dir[0] = state->directories.lpDir;
        fts = fts_open(dir, FTS_LOGICAL | FTS_NOCHDIR | FTS_NOSTAT, NULL);
        e = fts_read(fts);
        while (e)
        {
            //  residuals may be processed in short intervals, giving multiple files per day. Check that the time range covers the requested interval
            if ((e->fts_namelen == 59) && (strcmp(".cdf", e->fts_name + e->fts_namelen - 4) == 0) && (*(e->fts_name + 11) == swarmData->satellite[strlen(swarmData->satellite)-1]) && (strncmp("LP_HM", e->fts_name + 13, 5) == 0) && (strncmp(firstDateStr, e->fts_name + 19, 8) == 0) && (strncmp(firstTimeStr, e->fts_name+28, 6)>=0) && (strncmp(lastTimeStr, e->fts_name+44, 6)<=0))
            {
                cdfstatus = CDFopen(e->fts_path, &cdf);
                if (cdfstatus != CDF_OK)
                {
                    fprintf(stderr, "Could not open CDF %s\n", e->fts_path);
                    fts_close(fts);
                    return SPACESERIES_CDF_READ;
                }

                // Check that file has this parameter
                cdfstatus = CDFgetVarNum(cdf, measurementParameter);
                if (cdfstatus < CDF_OK)
                {
                    e = fts_read(fts);
                    continue;
                }

                // Data are 2 Hz: load then resample
                double *lptime = NULL;
                double *lpVal = NULL;
                size_t n = 0;
                // Get 0.5 s extra at each end of interval.
                getRecordRange("Timestamp", cdf, state->firstTime-500, state->lastTime+500, &firstRec, &lastRec);
                loadVariable("Timestamp", cdf, firstRec, lastRec, (void**)&lptime, &n);
                loadVariable(measurementParameter, cdf, firstRec, lastRec, (void**)&lpVal, &n);
                CDFclose(cdf);
                swarmData->parameter = calloc(nValues, sizeof(float));

                if (lptime == NULL || lpVal == NULL || swarmData->parameter == NULL)
                {
                    fts_close(fts);
                    return SPACESERIES_MEMORY;
                }

                // Might be able to continue without measurements from this satellite
                if (n == 0 || lptime[0] > state->firstTime || lptime[n-1] < state->lastTime)
                {
                    fts_close(fts);
                    free(lptime);
                    free(lpVal);
                    return SPACESERIES_OK;
                }

                double lastLpTime = lptime[0];
                double nextLpTime = lptime[0];
                ssize_t lastLpIndex = 0;
                ssize_t lpIndex = 0;
                double dt = 0.0;
                float interpFraction = 0.0;
                float p0 = 0.0;
                float p1 = 0.0;

                for (ssize_t timeIndex = 0; timeIndex < nValues; timeIndex++)
                {

                    while (lpIndex < n && timeIndex < nValues && swarmData->times[timeIndex] > lptime[lpIndex])
                    {
                        lastLpIndex = lpIndex++;
                        lastLpTime = lptime[lastLpIndex];
                    }

                    if (lpIndex == n)
                        lpIndex--;
                    nextLpTime = lptime[lpIndex];
                    dt = nextLpTime - lastLpTime;
                    if (dt < 100) // ephems are from 1000 ms samples
                        dt = 1.0; // arbitrary
                    interpFraction = (float)(swarmData->times[timeIndex] - lastLpTime) / dt;
                    swarmData->parameter[timeIndex] = (float)(lpVal[lastLpIndex] + (lpVal[lpIndex] - lpVal[lastLpIndex]) * interpFraction);
                }
                
                gotParameter = true;
                swarmData->nMeasurements = nValues;
                free(lptime);
                free(lpVal);
                break;
            }

            e = fts_read(fts);
        }

        fts_close(fts);

        if (swarmData->nMeasurements == -1)
        {
            fprintf(stderr, "Error loading TCT data.\n");
            return SPACESERIES_CDF_READ;
        }
        else if (swarmData->nMeasurements == 0)
        {
            return SPACESERIES_OK;
        }

    }



    swarmData->initialized = true;

    return status;
}

size_t downSample(double newPeriod, double *oldTime, size_t n, double *time, size_t nMeasurements, float **parameter)
{
    size_t newSampleIndex = 0;
    float sampleTotal = 0.0;
    size_t nSampled = 0;
    double dT = newPeriod * 1000; // milliseconds
    size_t oldSampleIndex = 0;
    while (oldSampleIndex < n && newSampleIndex < nMeasurements)
    {
        while (oldSampleIndex < n && oldTime[oldSampleIndex] < time[newSampleIndex] - dT / 2.0)
            oldSampleIndex++;
        while (oldSampleIndex < n && oldTime[oldSampleIndex] < time[newSampleIndex] + dT / 2.0)
        {
            sampleTotal += (*parameter)[oldSampleIndex];
            oldSampleIndex++;
            nSampled++;
        }
        if (nSampled > 0)
            (*parameter)[newSampleIndex] = sampleTotal / (float) nSampled;                        
        newSampleIndex++;
        nSampled = 0;
        sampleTotal = 0.0;
    }

    // Resize array
    float * newPointer = NULL;
    newPointer = (float *) realloc(*parameter, sizeof(float) * nMeasurements);
    if (newPointer == NULL)
    {
        fprintf(stderr, "Unable to reallocate parameter memory.\n");
        return 0;
    }

    *parameter = newPointer;


    return newSampleIndex;

}

void getRecordRange(char *timeVariableName, CDFid cdf, double firstTime, double lastTime, long *firstRec, long *lastRec)
{
    if (firstRec == NULL || lastRec == NULL)
        return;

    CDFstatus status = CDF_OK;

    CDFdata data = NULL;
    long numRecs = 0;
    long dataType = 0;
    long numElems = 0;
    long numDims = 0;
    long dimSizes[CDF_MAX_DIMS] = {0};
    long recVary = 0;
    long dimVarys[CDF_MAX_DIMS] = {0};

    status = CDFreadzVarAllByVarName(cdf, timeVariableName, &numRecs, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys, &data);
    if (status != CDF_OK)
        return;

    double *time = (double *)data;

    *firstRec = 0;
    *lastRec = 0;

    while (*firstRec < numRecs && time[*firstRec] < firstTime)
        (*firstRec)++;
    if (*firstRec == numRecs)
    {
        *firstRec = 0;
        return;
    }

    *lastRec = *firstRec;
    while (*lastRec < numRecs && time[*lastRec] <= lastTime)
        (*lastRec)++;

    if (*lastRec == numRecs)
    {
        *lastRec = numRecs -1;
    }

    CDFdataFree(data);

    return;

}

void loadVariable(char *varName, CDFid cdf, long firstRec, long lastRec, void **data, size_t *nMeasurements)
{

    CDFstatus status = CDF_OK;

    if (data == NULL)
        return;
    if (nMeasurements == NULL)
        return;

    *nMeasurements = 0;

    if (firstRec < 0 || lastRec < 0)
        return;

    size_t totalBytes = bytesPerRecord(varName, cdf) * (lastRec - firstRec + 1);
    if (totalBytes == 0)
        return;

    *data = malloc(totalBytes);
    if (*data == NULL)
        return;

    status = CDFgetVarRangeRecordsByVarName(cdf, varName, firstRec, lastRec, *data);
    if (status == CDF_OK)
        *nMeasurements = lastRec - firstRec + 1;
    else
        *nMeasurements = 0;

    return;

}


size_t bytesPerRecord(char *varName, CDFid cdf)
{
    CDFstatus  status = CDF_OK;

    // From CDF 3.8.1 C reference manual
    // Get number of bytes per record
    long varNumber = CDFgetVarNum(cdf, varName);
    if (varNumber < 0)
        return 0;

    long nDims = 0;
    status = CDFgetzVarNumDims(cdf, varNumber, &nDims);
    if (status != CDF_OK)
        return 0;

    long dimSizes[CDF_MAX_DIMS];
    status = CDFgetzVarDimSizes(cdf, varNumber, dimSizes);
    if (status != CDF_OK)
        return 0;

    long dataType = 0;
    status = CDFgetzVarDataType(cdf, varNumber, &dataType);
    if (status != CDF_OK)
        return 0;

    long bytesPerRecord = 1;
    status = CDFgetDataTypeSize(dataType, &bytesPerRecord);
    if (status != CDF_OK)
        return 0;

    for (int i = 0; i < nDims; i++)
    {
        bytesPerRecord *= dimSizes[i];
    }

    return bytesPerRecord;

}

void freeData(ProgramState *state)
{
    if (state == NULL)
        return;

    for (int s = 0; s < 3; s++)
    {
        freeSatelliteData(state->data.satelliteData[s]);
        free(state->data.satelliteData[s]);
    }
    free(state->data.spaceSeries);
    free(state->data.previousSpaceSeries);
    free(state->data.latitude);
    free(state->data.longitude);
    free(state->data.radius);
    free(state->data.qdlatitude);
    free(state->data.mlt);
    free(state->data.qdargoforbit);
    free(state->data.outputBuffer);

    return;
}

void freeSatelliteData(SwarmData *data)
{

    if (data == NULL)
        return;

    free(data->times);
    free(data->latitude);
    free(data->longitude);
    free(data->radius);
    free(data->qdlatitude);
    free(data->mlt);
    free(data->vsatn);
    free(data->vsate);
    free(data->vsatc);
    free(data->parameter);
    free(data->vsat);
    free(data->distanceAlongTrack);
    free(data->qdargoforbit);

    return;
}

int satelliteOrder(const void *first, const void *second)
{
    // Push missing data to the end
    if (first == NULL || second == NULL)
        return 1;

    SwarmData *sat1 = (SwarmData*)(*(SwarmData**)first);
    SwarmData *sat2 = (SwarmData*)(*(SwarmData**)second);

    if (sat1->nMeasurements < 2)
        return 1;
    if (sat2->nMeasurements < 2)
        return -1;

    // qsort arranges items in ascending order
    // Order here means distance behind the first satellite
    // Return -1 if first is ahead of second
    // Return +1 if first is behind second
    // Return 0 otherwise

    // Assess positions using first two measurements
    // First satellite is reference
    // Using along-track component of displacement from sat 1 to sat 2
    // Only works when satellites are not to far apart. Good, say, up to ~5000 km separation.
    float d1 = 0.0;
    float d2 = 0.0;
    interSatelliteDistanceKm(sat1, 0, sat2, 0, NULL, &d1, NULL, NULL);
    interSatelliteDistanceKm(sat1, 0, sat2, 1, NULL, &d2, NULL, NULL);
    // If the second satellite's second position is closer to the first satellite's first position
    // The second satellite is behind the first
    if (fabs(d2) < fabs(d1))
        return -1;
    else if (fabs(d2) > fabs(d1))
        return 1;
    else
        return 0;

}

int highestVersionFirst(const FTSENT **first, const FTSENT **second)
{
    if (first == NULL)
        return 1;
    if (second == NULL)
        return -1;

    FTSENT *e1 = (FTSENT *)first;
    FTSENT *e2 = (FTSENT *)second;
    if (e1 == NULL)
        return 1;
    if (e2 == NULL)
        return -1;

    if (e1->fts_namelen != 59)
        return 1;
    if (e2->fts_namelen != 59)
        return -1;

    int comparison = strncmp(e1->fts_name + e1->fts_namelen - 8, e2->fts_name + e2->fts_namelen - 8, 4);
    return -comparison;
}

int loadEphemeres(SwarmData *swarmData, char *modDir, double firstTime, double lastTime, Ephemeres **ephem, size_t *nEphemeres)
{

    if (ephem == NULL || nEphemeres == NULL)
        return SPACESERIES_EPHEM_READ;

    char* dir[2] = {modDir, NULL};
    FTSENT *e = NULL;
    FTS *fts = NULL;

    char firstDateStr[9] = {0};
    char firstTimeStr[7] = {0};
    dateTimeStrings(firstTime, firstDateStr, firstTimeStr);
    char lastDateStr[9] = {0};
    char lastTimeStr[7] = {0};
    dateTimeStrings(lastTime, lastDateStr, lastTimeStr);

    // For dates with multiple versions, reads the highest first.
    fts = fts_open(dir, FTS_LOGICAL | FTS_NOCHDIR | FTS_NOSTAT, &highestVersionFirst);
    e = fts_read(fts);

    int lastFileVersion = 0;
    int fileVersion = 0;
    char fileVersionString[5] = {0};

    while (e)
    {
        // TODO read in previous day file as well to get full coverage, if needed.
        if ((e->fts_namelen == 59) && ((strcmp(".sp3", e->fts_name + e->fts_namelen - 4) == 0) || (strcmp(".DBL", e->fts_name + e->fts_namelen - 4) == 0)) && (*(e->fts_name + 11) == swarmData->satellite[strlen(swarmData->satellite)-1]) && (strncmp("_SC_1B", e->fts_name + 12, 6) == 0) && (strncmp(firstDateStr, e->fts_name + 19, 8) == 0) && (strncmp(firstTimeStr, e->fts_name+28, 6)>=0) && (strncmp(lastTimeStr, e->fts_name+44, 6)<=0))
        {
            strncpy(fileVersionString, e->fts_name + e->fts_namelen - 8, 4);
            fileVersion = atoi(fileVersionString);
            if (fileVersion > lastFileVersion)
            {
                lastFileVersion = fileVersion;
                // Read and store ephemeres
                FILE *file = fopen(e->fts_path, "r");
                if (file == NULL)
                {
                    e = fts_read(fts);
                    continue;
                }
                // read number of points and GPS time offset
                char buffer[256] = {0};
                char *res = fgets(buffer, 255, file);
                if (res != buffer)
                {
                    fclose(file);
                    fts_close(fts);
                    return SPACESERIES_EPHEM_READ;
                }
                char nEpochsString[8];
                strncpy(nEpochsString, buffer+32, 8);
                int nEpochs = atoi(nEpochsString);
                if (nEpochs < 1 || nEpochs > 86400)
                {
                    fclose(file);
                    fts_close(fts);
                    return SPACESERIES_EPHEM_READ;
                }

                // ra-allocate memory for ephemeres
                void *tmp = realloc(*ephem,(size_t) (*nEphemeres + nEpochs) * sizeof(Ephemeres));
                if (tmp == NULL)
                {
                    fclose(file);
                    fts_close(fts);
                    return SPACESERIES_MEMORY;
                }
                *ephem = tmp;
                size_t count = 0;
                int year = 0;
                int month = 0;
                int day = 0;
                int hour = 0;
                int minute = 0;
                double second = 0.0;
                int sec = 0;
                int msec = 0;
                bool gotTime = false;
                bool gotPosition = false;
                float x = 0.0;
                float y = 0.0;
                float z = 0.0;
                float ufo = 0.0; // unidentified formatted object
                float vx = 0.0;
                float vy = 0.0;
                float vz = 0.0;

                float chatx = 0.0;
                float chaty = 0.0;
                float chatz = 0.0;
                float ehatx = 0.0;
                float ehaty = 0.0;
                float ehatz = 0.0;
                float nhatx = 0.0;
                float nhaty = 0.0;
                float nhatz = 0.0;
                float cmag = 0.0;
                float emag = 0.0;
                float nmag = 0.0;

                float r = 0.0;
                float degrees = M_PI / 180.0;
                char d1 = 0;
                char d2 = 0;

                int nRead = 0;
                nRead = sscanf(buffer, "#cV%4d%3d%3d%3d%3d%12lf", &year, &month, &day, &hour, &minute, &second);
                sec = (int)floor(second);
                msec = (int)floor(1000.0 * (second - (double)sec));
                double gpsEpoch = computeEPOCH(year, month, day, hour, minute, sec, msec);
                double utEpoch = computeEPOCH(year, month, day, 0, 0, 0, 0);
                double gpsTimeOffset = gpsEpoch - utEpoch;
                double gpsTime = 0.0;
                while ((count < nEpochs) && (fgets(buffer, 255, file) != NULL))
                {
                    if ((buffer[0] == '*') && (gotPosition == false))
                    {
                        nRead = sscanf(buffer, "*  %4d%3d%3d%3d%3d%12lf", &year, &month, &day, &hour, &minute, &second);
                        if (nRead == 6)
                        {
                            sec = (int)floor(second);
                            msec = (int)floor(1000.0 * (second - (double)sec));
                            gpsTime = computeEPOCH(year, month, day, hour, minute, second, msec);
                            (*ephem)[count].time = gpsTime - gpsTimeOffset;
                            gotTime = true;
                        }
                        else
                            gotTime = false;
                    }
                    else if ((strncmp(buffer, "PL", 2) == 0) && (gotTime == true))
                    {
                        nRead = sscanf(buffer, "PL%c%c%14f%14f%14f%14f", &d1, &d2, &x, &y, &z, &ufo);
                        if (nRead == 6)
                        {
                            gotPosition = true;
                        }
                        else
                        {
                            gotPosition = false;
                        }
                    }
                    else if ((strncmp(buffer, "VL", 2) == 0) && (gotTime == true && gotPosition == true))
                    {
                        nRead = sscanf(buffer, "VL%c%c%14f%14f%14f%14f", &d1, &d2, &vx, &vy, &vz, &ufo);
                        if (nRead == 6)
                        {
                            // Calculate ephemeres
                            r = sqrtf(x*x + y*y + z*z);
                            (*ephem)[count].latitude = 90.0 - acosf(z / r) / degrees;
                            (*ephem)[count].longitude = atan2f(y, x) / degrees;
                            (*ephem)[count].radius = 1000.0 * r;
                            chatx = -x / r;
                            chaty = -y / r;
                            chatz = -z / r;
                            ehatx = chaty;
                            ehaty = -chatx;
                            ehatz = 0.0;
                            emag = sqrtf(ehatx*ehatx + ehaty*ehaty);
                            ehatx /= emag;
                            ehaty /= emag;
                            nhatx = ehaty * chatz - ehatz * chaty;
                            nhaty = -ehatx * chatz + ehatz * chatx;
                            nhatz = ehatx * chaty - ehaty * chatx;
                            nmag = sqrtf(nhatx*nhatx + nhaty*nhaty + nhatz*nhatz);
                            nhatx /= nmag;
                            nhaty /= nmag;
                            nhatz /= nmag;
                            (*ephem)[count].vsatN = (vx * nhatx + vy * nhaty + vz * nhatz) / 10.0;
                            (*ephem)[count].vsatE = (vx * ehatx + vy * ehaty + vz * ehatz) / 10.0;
                            (*ephem)[count].vsatC = (vx * chatx + vy * chaty + vz * chatz) / 10.0;

                            // Reset for next epoch
                            gotTime = false;
                            gotPosition = false;
                            count++;
                        }
                    }
                    else
                    {
                        continue;
                    }
                    
                }
                *nEphemeres += (size_t) nEpochs;
                
                fclose(file);
                break;
            }            
        }
        e = fts_read(fts);
    }
    fts_close(fts);

    return SPACESERIES_OK;

}

void dateTimeStrings(double time, char *dateString, char *timeString)
{
    if (dateString == NULL || timeString == NULL)
        return;

    long year = 0;
    long month = 0;
    long day = 0;
    long hour = 0;
    long minute = 0;
    long second = 0;
    long msec = 0;

    EPOCHbreakdown(time, &year, &month, &day, &hour, &minute, &second, &msec);
    sprintf(dateString, "%4ld%02ld%02ld", year, month, day);
    sprintf(timeString, "%02ld%02ld%02ld", hour, minute, second);

    return;
}

void calculateAlongTrackDistances(ProgramState *state)
{
    // Get distance in km along the leading satellite's track
    SwarmData *sat = state->data.satelliteData[0];
    // Estimate distance lags of the second and third satellites
    // Use latitude and longitude, and calculate starting distances
    float alongTrackDisplacement = 0.0;

    state->data.satelliteData[0]->distanceAlongTrack[0] = alongTrackDisplacement;
    double totalDistance = 0.0;
    for (int n = 1; n < sat->nMeasurements; n++)
    {
        interSatelliteDistanceKm(state->data.satelliteData[0], n-1, sat, n, NULL, &alongTrackDisplacement, NULL, NULL);
        totalDistance += alongTrackDisplacement;
        sat->distanceAlongTrack[n] = (float)totalDistance;
    }
    // Along-track distances for trailing satellites with respect to Swarm B
    // This method takes into account differences and variations in satellite velocities
    // so that distance along track has a common meaning for all satellites.
    for (int i = 1; i < state->data.nSatellites; i++)
    {
        sat = state->data.satelliteData[i];
        alongTrackDisplacement = 0.0;
        interSatelliteDistanceKm(state->data.satelliteData[0], 0, sat, 0, NULL, &alongTrackDisplacement, NULL, NULL);
        sat->distanceAlongTrack[0] = alongTrackDisplacement;
        totalDistance = (double)sat->distanceAlongTrack[0];
        for (int n = 1; n < sat->nMeasurements; n++)
        {
            interSatelliteDistanceKm(state->data.satelliteData[0], n, sat, n, NULL, &alongTrackDisplacement, NULL, NULL);
            totalDistance = (double)state->data.satelliteData[0]->distanceAlongTrack[n] + alongTrackDisplacement;
            sat->distanceAlongTrack[n] = (float)totalDistance;
        }
    }

    state->minDistance = state->data.satelliteData[0]->distanceAlongTrack[0];
    state->maxDistance = state->data.satelliteData[state->data.nSatellites-1]->distanceAlongTrack[state->data.satelliteData[state->data.nSatellites-1]->nMeasurements-1];

    return;
}

void calculateMagneticLatitudeArgOfOrbit(ProgramState *state)
{
    // Calculate Quasi-dipole magnetic latitude argument of orbit
    // Let arg of orbit increase indefinitely to make it
    // easy to interpolate later
    float nOrbits = 0;

    SwarmData *sat = NULL;

    for (int i = 0; i < state->data.nSatellites; i++)
    {
        sat = state->data.satelliteData[i];
        sat->qdargoforbit[0] = sat->qdlatitude[0];
        nOrbits = 0;
        if (sat->nMeasurements > 1)
        {
            sat->qdargoforbit[0] = sat->qdlatitude[0];
            if (sat->qdlatitude[1] - sat->qdlatitude[0] < 0)
                sat->qdargoforbit[0] = 180.0 - sat->qdargoforbit[0];
            else if (sat->qdlatitude[0] < 0)
                sat->qdargoforbit[0] = 360.0 + sat->qdargoforbit[0];
        }
        for (int n = 1; n < sat->nMeasurements; n++)
        {
            if (sat->qdlatitude[n] >= 0 && sat->qdlatitude[n-1] < 0)
                nOrbits++;
            if (sat->qdlatitude[n] - sat->qdlatitude[n-1] < 0)
                sat->qdargoforbit[n] = 360.0 * nOrbits + 180.0 - sat->qdlatitude[n];
            else if (sat->qdlatitude[n] < 0)
                sat->qdargoforbit[n] = 360.0 * (nOrbits + 1) + sat->qdlatitude[n];
            else
                sat->qdargoforbit[n] = 360.0 * nOrbits + sat->qdlatitude[n];
        }
    }

    return;
}