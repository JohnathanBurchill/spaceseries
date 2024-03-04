/*

    spaceseries: util.c

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

#include "util.h"

#include "data.h"
#include "spaceseries.h"

#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cdf.h>

void latLonToXYZ(float lat, float lon, float *x, float *y, float *z)
{
    if (x == NULL || y == NULL || z == NULL)
        return;

    float degrees = M_PI / 180.0;
    *x = sinf((90.0 - lat) * degrees) * cosf(lon * degrees);
    *y = sinf((90.0 - lat) * degrees) * sinf(lon * degrees);
    *z = cosf((90.0 - lat) * degrees);

    return;
}

void xyzToLatLon(float x, float y, float z, float *lat, float *lon)
{
    if (lat == NULL || lon == NULL)
        return;

    float degrees = M_PI / 180.0;
    float r = sqrt(x * x + y * y + z * z);
    float rho = sqrt(x * x + y * y);
    *lat = 90.0 - acosf(z / r) / degrees;
    *lon = atan2f(y, x) / degrees;

    return;
}

void interSatelliteDistanceKm(SwarmData *sat1, size_t n1, SwarmData *sat2, size_t n2, float *distance, float *alongTrackDisplacement, float *horizontalCrossTrackDisplacement, float *verticalCrossTrackDisplacement)
{
    if (n1 >= sat1->nMeasurements || n2 >= sat2->nMeasurements)
    {
        if (distance != NULL)
            *distance = NAN;
        if (alongTrackDisplacement != NULL)
            *alongTrackDisplacement = NAN;
        if (horizontalCrossTrackDisplacement != NULL)
            *horizontalCrossTrackDisplacement = NAN;
        if (verticalCrossTrackDisplacement != NULL)
            *verticalCrossTrackDisplacement = NAN;
        return;
    }

    float dtor = M_PI / 180.0;

    float theta1 = (90.0 - sat1->latitude[n1]) * dtor;
    float phi1 = sat1->longitude[n1] * dtor;
    float r1 = sat1->radius[n1];

    float theta2 = (90.0 - sat2->latitude[n2]) * dtor;
    float phi2 = sat2->longitude[n2] * dtor;
    float r2 = sat2->radius[n2];

    float x1 = r1 * sinf(theta1) * cosf(phi1);
    float y1 = r1 * sinf(theta1) * sinf(phi1);
    float z1 = r1 * cosf(theta1);

    float x2 = r2 * sinf(theta2) * cosf(phi2);
    float y2 = r2 * sinf(theta2) * sinf(phi2);
    float z2 = r2 * cosf(theta2);

    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    float mag = 0.0;

    if (distance != NULL)
        *distance = sqrtf(dx*dx + dy*dy + dz*dz) / 1000.0;

    if (alongTrackDisplacement != NULL || horizontalCrossTrackDisplacement != NULL || verticalCrossTrackDisplacement != NULL)
    {
        // Use first satellite's velocity for track-aligned coordinates
        // Probably better to use actual trajectories and arc lengths
        float vn = sat1->vsatn[n1];
        float ve = sat1->vsate[n1];
        float vc = sat1->vsatc[n1];

        // Transform VSatNEC to XYZ frame at first satellite
        float r = sqrtf(x1*x1 + y1*y1 + z1*z1);

        // C is minus the normalized position vector
        float cx = - x1 / r1;
        float cy = - y1 / r1;
        float cz = - z1 / r1;

        // E is C-cross-z, normalized
        float ex = cy;
        float ey = -cx;
        float ez = 0.0;
        mag = sqrtf(ex*ex + ey*ey + ez*ez);
        ex /= mag;
        ey /= mag;
        ez /= mag;

        // N is E cross C
        float nx = ey * cz - ez *cy;
        float ny = -ex * cz + ez * cx;
        float nz = ex * cy - ey * cx;
        mag = sqrtf(nx*nx + ny*ny + nz*nz);
        nx /= mag;
        ny /= mag;
        nz /= mag;

        // sat velocity in XYZ
        float vx = vn * nx + ve * ex + vc * cx;
        float vy = vn * ny + ve * ey + vc * cy;
        float vz = vn * nz + ve * ez + vc * cz;
        float vsat = sqrtf(vx*vx + vy*vy + vz*vz);

        // along-track x: VSat
        float sxx = vx / vsat;
        float sxy = vy / vsat;
        float sxz = vz / vsat;
        // horizontal cross-track sy: C cross sx
        float syx = cy * sxz - cz * sxy;
        float syy = -cx * sxz + cz * sxx;
        float syz = cx * sxy - cy * sxx;
        mag = sqrtf(syx*syx + syy*syy + syz*syz);
        syx /= mag;
        syy /= mag;
        syz /= mag;
        // vertical (down) cross-track z: sx cross sy
        float szx = sxy * syz - sxz * syy;
        float szy = -sxx * syz + sxz * syx;
        float szz = sxx * syy - sxy * syx;
        // Probably unnecessary to renormalize z
        mag = sqrtf(szx*szx + szy*szy + szz*szz);
        szx /= mag;
        szy /= mag;
        szz /= mag;

        // Approximate, using mean track-aligned coordinate system at the midpoint between the sats
        if (alongTrackDisplacement != NULL)
        {
            // In the direction of motion
            *alongTrackDisplacement = (dx * sxx + dy * sxy + dz * sxz) / 1000.0;
        }
        if (horizontalCrossTrackDisplacement != NULL)
        {
            // To the right (facing forward)
            *horizontalCrossTrackDisplacement = (dx * syx + dy * syy + dz * syz) / 1000.0;
        }
        if (verticalCrossTrackDisplacement != NULL)
        {
            // Downward, use difference in altitudes, since the displacement
            // does not take into account curvature of the track
            *verticalCrossTrackDisplacement = (r2 - r1) / 1000.0;
        }
    }

    return;

}

int getOptionValue(char *number, float *value)
{
    if (value == NULL)
        return SPACESERIES_PARSE_NUMBER;

    char *lastPos = number;
    float val = (float) strtod(number, &lastPos);
    int status = SPACESERIES_OK;
    if (lastPos != number)
        *value = val;
    else
        status = SPACESERIES_PARSE_NUMBER;

    return status;
}

int setOutputFilename(ProgramState *state)
{
    if (state == NULL)
        return SPACESERIES_POINTERS;

    char format[EPOCHx_FORMAT_MAX] = "<year><mm.02><dom.02>T<hour><min><sec>";
    char t1[EPOCHx_STRING_MAX] = {0};
    char t2[EPOCHx_STRING_MAX] = {0};
    encodeEPOCHx(floor(state->firstTime/1000.0)*1000.0, format, t1);
    encodeEPOCHx(ceil(state->lastTime/1000.0)*1000.0, format, t2);

    char outputFileNumber[15] = {0};
    if (state->splitOutputFilesByOrbit && state->labelOutputFileOrbits)
        sprintf(outputFileNumber, "_orbit_%02d", state->nOutputFiles);

    if (state->binaryOutput)
        (void)snprintf(state->outFile, FILENAME_MAX, "%s/%s_%s_%s%s%s.bin", state->directories.outDir, state->measurementParameter, t1, t2, state->suffix, outputFileNumber);
    else
        (void)snprintf(state->outFile, FILENAME_MAX, "%s/%s_%s_%s%s%s.txt", state->directories.outDir, state->measurementParameter, t1, t2, state->suffix, outputFileNumber);

    return SPACESERIES_OK;
}

int openOutputFile(ProgramState *state)
{
    if (state == NULL)
        return SPACESERIES_POINTERS;

    int status = SPACESERIES_OK;

    char nul = '\0';

    Data *data = &state->data;

    state->nOutputFiles++;

    status = setOutputFilename(state);
    if (status != SPACESERIES_OK)
        return status;

    if (state->binaryOutput)
    {
        ssize_t nWritten = 0;
        state->outputfd = open(state->outFile, O_WRONLY | O_CREAT | O_TRUNC, 0755);
        if (state->outputfd == -1)
            return SPACESERIES_FILE_WRITE;
        nWritten = write(state->outputfd, &data->nSatellites, sizeof data->nSatellites);
        if (nWritten != sizeof data->nSatellites)
            return SPACESERIES_FILE_WRITE;
        for (int i = 0; i < data->nSatellites; i++)
        {
            nWritten = write(state->outputfd, data->satelliteData[i]->satellite, sizeof *data->satelliteData[i]->satellite * strlen(data->satelliteData[i]->satellite));
            if (nWritten != sizeof *data->satelliteData[i]->satellite * strlen(data->satelliteData[i]->satellite))
                return SPACESERIES_FILE_WRITE;
            nWritten = write(state->outputfd, &nul, sizeof nul);
            if (nWritten != sizeof nul)
                return SPACESERIES_FILE_WRITE;
        }
        for (int i = 0; i < data->nSatellites; i++)
        {
            nWritten = write(state->outputfd, &data->initialTimeLag[i], sizeof data->initialTimeLag[i]);
            if (nWritten != sizeof data->initialTimeLag[i])
                return SPACESERIES_FILE_WRITE;
            nWritten = write(state->outputfd, &data->initialDistanceLag[i], sizeof data->initialDistanceLag[i]);
            if (nWritten != sizeof data->initialDistanceLag[i])
                return SPACESERIES_FILE_WRITE;
        }
    }
    else
    {
        state->output = fopen(state->outFile, "w");
        if (state->output == NULL)
            return SPACESERIES_FILE_WRITE;
        fprintf(state->output, "Satellite order:");
        for (int i = 0; i < data->nSatellites; i++)
            fprintf(state->output, " %s", data->satelliteData[i]->satellite);
        fprintf(state->output, "\n");
        fprintf(state->output, "Satellite lags:");
        for (int i = 0; i < data->nSatellites; i++)
            fprintf(state->output, " %s -%.1f s -%.1f km;", data->satelliteData[i]->satellite, data->initialTimeLag[i], data->initialDistanceLag[i]);
        fprintf(state->output, "\n");
        fprintf(state->output, "Space series consist of a header for each time, and one or more space series points.\n");
        fprintf(state->output, "Headers: <time (ms from 0000-01-01T00:00:00.000)> <number of points in series>, and the following on the same line:\n");
        for (int s = 0; s < data->nSatellites; s++)
        {
            fprintf(state->output, " %s:", data->satelliteData[s]->satellite);
            fprintf(state->output, " <Lat (deg)>");
            fprintf(state->output, " <Lon (deg)>");
            fprintf(state->output, " <Rad (m)>");
            fprintf(state->output, " <Along track distance (km)>");
            fprintf(state->output, " <QDLat (deg)>");
            fprintf(state->output, " <MLT (hours)>");
            fprintf(state->output, " <QDArgOfOrbit (deg)>");
            fprintf(state->output, " <Along-track displacement (km)>");
            fprintf(state->output, " <Horizontal cross-track displacement (km)>");
            fprintf(state->output, " <Vertical cross-track displacement (km)>");
            fprintf(state->output, " <static %s (km)>", state->measurementParameter);
        }
        fprintf(state->output, "Series point, one per line: <time (ms from 0000-01-01T00:00:00.000)> <distance along track (km)> <latitude (deg)> <longitude (deg)> <radius (m)> <QDLatitude (deg)> <MLT (hours)> <QDArgOfORbit (deg)> <%s (units from source CDF)> <d(%s)/dt (units per second)>", state->measurementParameter, state->measurementParameter);
        for (int s = 0; s < data->nSatellites; s++)
        {
            fprintf(state->output, " <delta %s (dynamic - static) %s (km)>", state->measurementParameter, state->data.satelliteData[s]->satellite);
        }
        fprintf(state->output, "\n");
    }

    return status;

}

int closeOutputFile(ProgramState *state)
{
    int status = SPACESERIES_OK;

    char *currentFilename = strdup(state->outFile);
    if (currentFilename == NULL)
        return SPACESERIES_MEMORY;

    if (state->binaryOutput)
        close(state->outputfd);
    else
        fclose(state->output);

    // Update filename to actual start and stop times
    setOutputFilename(state);
    if (strcmp(currentFilename, state->outFile) != 0)
        status = rename(currentFilename, state->outFile);

    free(currentFilename);

    return status;
}
