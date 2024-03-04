/*

    spaceseries: interpolation.c

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

#include "interpolation.h"

#include "data.h"
#include "util.h"

#include <math.h>
#include <float.h>

ssize_t closestTimeIndex(SwarmData *swarmData, double time)
{
    if (swarmData == NULL)
        return 0;

    double dt = 0.0;
    double lastdt = DBL_MAX;

    if (swarmData->closestTimeIndex >= swarmData->nMeasurements)
        swarmData->closestTimeIndex = swarmData->nMeasurements-1;
    if (swarmData->closestTimeIndex < 0)
        swarmData->closestTimeIndex = 0;

    // Probably a faster way involves cutting the invertal by two on each check?

    if (swarmData->times[swarmData->closestTimeIndex]<time)
    {
        while(swarmData->closestTimeIndex < swarmData->nMeasurements && (dt = fabs(swarmData->times[swarmData->closestTimeIndex] - time)) < lastdt)
        {
            lastdt = dt;
            swarmData->closestTimeIndex++;
        }
        swarmData->closestTimeIndex--;
    }
    else if (swarmData->times[swarmData->closestTimeIndex] > time)
    {
        while(swarmData->closestTimeIndex > -1 && (dt = fabs(swarmData->times[swarmData->closestTimeIndex] - time)) < lastdt)
        {
            lastdt = dt;
            swarmData->closestTimeIndex--;
        }
    }
    swarmData->closestTimeIndex++;

    if (swarmData->closestTimeIndex < 0)
        swarmData->closestTimeIndex = 0;
    else if (swarmData->closestTimeIndex > swarmData->nMeasurements - 1)
        swarmData->closestTimeIndex = swarmData->nMeasurements - 1;

    // Return for convenience. Can discard return value just to update the closest index for later use.
    return swarmData->closestTimeIndex;
    
}

void neighboringDistanceIndices(float distance, SwarmData *sat)
{
    if (sat == NULL)
        return;

    long i1 = (long)sat->i1 - 1;
    long i2 = (long)sat->i2 - 1;
    long n = (long)sat->nMeasurements;
    if (i1 < 0)
        i1 = 0;
    if (i2 < 0)
        i2 = 0;
    while(i2 < n && sat->distanceAlongTrack[i2] <= distance)
        i2++;

    if (i2 == 0 || i2 == n)
    {
        sat->i1 = -1;
        sat->i2 = -1;
        return;
    }

    // The satellite's along-track distance at i2 is greater than or equal to the reference distance.
    // This means the along-track distrance at i2 - 1 is less than the reference distance
    i1 = i2 - 1;

    sat->i1 = i1;
    sat->i2 = i2;

    return;
}

void interpolate(Data *data, double tMilliseconds, float requestedDistance, bool constantInterpolation, float *parameter, float *latitude, float *longitude, float *radius, float *qdlat, float *mlt, float *qdargoforbit)
{

    if (parameter == NULL)
        return;

    double t0 = -1.0;
    double t1 = -1.0;
    double t2 = -1.0;
    double tmid = -1.0;
    double tmax = -1.0;
    float interpFraction = 0.0;
    float param0 = 0.0;
    float param1 = 0.0;
    float param2 = 0.0;
    float param = 0.0;
    float latitude0 = 0.0;
    float latitude1 = 0.0;
    float latitude2 = 0.0;
    float lat = 0.0;
    float longitude0 = 0.0;
    float longitude1 = 0.0;
    float longitude2 = 0.0;
    float lon = 0.0;
    float radius0 = 0.0;
    float radius1 = 0.0;
    float radius2 = 0.0;
    float rad = 0.0;
    float qdlat0 = 0.0;
    float qdlat1 = 0.0;
    float qdlat2 = 0.0;
    float qdl = 0.0;
    float mlt0 = 0.0;
    float mlt1 = 0.0;
    float mlt2 = 0.0;
    float ml = 0.0;
    float qdargoforbit0 = 0.0;
    float qdargoforbit1 = 0.0;
    float qdargoforbit2 = 0.0;
    float qdarg = 0.0;

    SwarmData *sat = NULL;

    // Interpolate the parameter and time of measurement at the requested along-track distance 
    // for each satellite.
    interpolateAtDistance(data, 0, requestedDistance, &t0, &param0, &latitude0, &longitude0, &radius0, &qdlat0, &mlt0, &qdargoforbit0);
    interpolateAtDistance(data, 1, requestedDistance, &t1, &param1, &latitude1, &longitude1, &radius1, &qdlat1, &mlt1, &qdargoforbit1);
    if (data->nSatellites == 3)
        interpolateAtDistance(data, 2, requestedDistance, &t2, &param2, &latitude2, &longitude2, &radius2, &qdlat2, &mlt2, &qdargoforbit2);

    if (t0 < 0.0 || t1 < 0.0 || (data->nSatellites == 3 && t2 < 0.0))
    {
        *parameter = NAN;
        if (latitude != NULL)
            *latitude = NAN;
        if (longitude != NULL)
            *longitude = NAN;
        if (radius != NULL)
            *radius = NAN;
        if (qdlat != NULL)
            *qdlat = NAN;
        if (mlt != NULL)
            *mlt = NAN;
        if (qdargoforbit != NULL)
            *qdargoforbit = NAN;
        return;
    }

    if (data->nSatellites == 3)
        tmax = t2;
    else
        tmax = t1;

    // No extrapolation in time...
    // -INFINITY means Swarm B has not visited the point yet
    // INFINITY means Swarm C already passed the point
    if (tMilliseconds < t0)
    {
        *parameter = -INFINITY;
        if (latitude != NULL)
            *latitude = -INFINITY;
        if (longitude != NULL)
            *longitude = -INFINITY;
        if (radius != NULL)
            *radius = -INFINITY;
        if (qdlat != NULL)
            *qdlat = -INFINITY;
        if (mlt != NULL)
            *mlt = -INFINITY;
        if (qdargoforbit != NULL)
            *qdargoforbit = -INFINITY;
        return;
    }
    if (tMilliseconds > tmax)
    {
        *parameter = INFINITY;
        if (latitude != NULL)
            *latitude = INFINITY;
        if (longitude != NULL)
            *longitude = INFINITY;
        if (radius != NULL)
            *radius = INFINITY;
        if (qdlat != NULL)
            *qdlat = INFINITY;
        if (mlt != NULL)
            *mlt = INFINITY;
        if (qdargoforbit != NULL)
            *qdargoforbit = INFINITY;
        return;
    }

    if (t0 <= tMilliseconds && tMilliseconds < t1)
    {
        interpFraction = (tMilliseconds - t0) / (t1 - t0);
        if (constantInterpolation)
        {
            tmid = (t0 + t1) / 2.0;
            if (tMilliseconds <= tmid)
            {
                param = param0;
                lat = latitude0;
                lon = longitude0;
                rad = radius0;
                qdl = qdlat0;
                ml = mlt0;
                qdarg = qdargoforbit0;
            }
            else
            {
                param = param1;
                lat = latitude1;
                lon = longitude1;
                rad = radius1;
                qdl = qdlat1;
                ml = mlt1;
                qdarg = qdargoforbit1;
            }
        }
        else
        {
            param = param0 + (param1 - param0) * interpFraction;

            float x0 = 0.0, y0 = 0.0, z0 = 0.0;
            float x1 = 0.0, y1 = 0.0, z1 = 0.0;
            float x = 0.0, y = 0.0, z = 0.0;
            latLonToXYZ(latitude0, longitude0, &x0, &y0, &z0);
            latLonToXYZ(latitude1, longitude1, &x1, &y1, &z1);
            x = x0 + (x1 - x0) * interpFraction;
            y = y0 + (y1 - y0) * interpFraction;
            z = z0 + (z1 - z0) * interpFraction;
            xyzToLatLon(x, y, z, &lat, &lon);

            rad = radius0 + (radius1 - radius0) * interpFraction; 

            latLonToXYZ(qdlat0, mlt0 / 24.0 * 360.0, &x0, &y0, &z0);
            latLonToXYZ(qdlat1, mlt1 / 24.0 * 360.0, &x1, &y1, &z1);
            x = x0 + (x1 - x0) * interpFraction;
            y = y0 + (y1 - y0) * interpFraction;
            z = z0 + (z1 - z0) * interpFraction;
            xyzToLatLon(x, y, z, &qdl, &ml);
            ml = ml /360.0 * 24.0;
            if (ml < 0.0)
                ml = 24.0 + ml;
            qdarg = qdargoforbit0 + (qdargoforbit1 - qdargoforbit0) * interpFraction;
        }
    }
    else if (data->nSatellites == 3)
    {
        interpFraction = (tMilliseconds - t1) / (t2 - t1);
        if (constantInterpolation)
        {
            tmid = (t1 + t2) / 2.0;
            if (tMilliseconds <= tmid)
            {
                param = param1;
                lat = latitude1;
                lon = longitude1;
                rad = radius1;
                qdl = qdlat1;
                ml = mlt1;
                qdarg = qdargoforbit1;
            }
            else
            {
                param = param2;
                lat = latitude2;
                lon = longitude2;
                rad = radius2;
                qdl = qdlat2;
                ml = mlt2;
                qdarg = qdargoforbit2;
            }
        }
        else
        {
            param = param1 + (param2 - param1) * interpFraction;

            float x0 = 0.0, y0 = 0.0, z0 = 0.0;
            float x1 = 0.0, y1 = 0.0, z1 = 0.0;
            float x = 0.0, y = 0.0, z = 0.0;
            latLonToXYZ(latitude1, longitude1, &x0, &y0, &z0);
            latLonToXYZ(latitude2, longitude2, &x1, &y1, &z1);
            x = x0 + (x1 - x0) * interpFraction;
            y = y0 + (y1 - y0) * interpFraction;
            z = z0 + (z1 - z0) * interpFraction;
            xyzToLatLon(x, y, z, &lat, &lon);

            rad = radius1 + (radius2 - radius1) * interpFraction; 

            latLonToXYZ(qdlat1, mlt1 / 24.0 * 360.0, &x0, &y0, &z0);
            latLonToXYZ(qdlat2, mlt2 / 24.0 * 360.0, &x1, &y1, &z1);
            x = x0 + (x1 - x0) * interpFraction;
            y = y0 + (y1 - y0) * interpFraction;
            z = z0 + (z1 - z0) * interpFraction;
            xyzToLatLon(x, y, z, &qdl, &ml);
            ml = ml /360.0 * 24.0;
            if (ml < 0.0)
                ml = 24.0 + ml;
            qdarg = qdargoforbit1 + (qdargoforbit2 - qdargoforbit1) * interpFraction;
        }
    }

    *parameter = param;
    if (latitude != NULL)
        *latitude = lat;
    if (longitude != NULL)
        *longitude = lon;
    if (radius != NULL)
        *radius = rad;
    if (qdlat != NULL)
        *qdlat = qdl;
    if (mlt != NULL)
        *mlt = ml;
    if (qdargoforbit != NULL)
        *qdargoforbit = qdarg;

    return;

}

float interpolateAtDistance(Data *data, int satelliteIndex, float requestedDistance, double *time, float *parameter, float *latitude, float *longitude, float *radius, float *qdlat, float *mlt, float *qdargoforbit)
{
    SwarmData *sat = NULL;
    sat = data->satelliteData[satelliteIndex];
    neighboringDistanceIndices(requestedDistance, sat);
    ssize_t i1 = sat->i1;
    ssize_t i2 = sat->i2;
    if (i1 < 0 || i2 < 0)
    {
        *time = -1.0;
        *parameter = 0.0;
        return 0.0;
    }
    float d1 = sat->distanceAlongTrack[i1];
    float d2 = sat->distanceAlongTrack[i2];
    float interpFraction = ((requestedDistance - d1) / (d2 - d1));
    if (time != NULL)
        *time = sat->times[i1] + (sat->times[i2] - sat->times[i1]) * interpFraction;
    if (parameter != NULL)
        *parameter = sat->parameter[i1] + (sat->parameter[i2] - sat->parameter[i1]) * interpFraction;


    float degrees = M_PI / 180.0;
    if (latitude != NULL || longitude != NULL)
    {
        float x1 = 0.0, y1 = 0.0, z1 = 0.0;
        float x2 = 0.0, y2 = 0.0, z2 = 0.0;
        float x3 = 0.0, y3 = 0.0, z3 = 0.0;
        float lat = 0.0;
        float lon = 0.0;
        latLonToXYZ(sat->latitude[i1], sat->longitude[i1], &x1, &y1, &z1);
        latLonToXYZ(sat->latitude[i2], sat->longitude[i2], &x2, &y2, &z2);
        x3 = x1 + (x2 - x1) * interpFraction;
        y3 = y1 + (y2 - y1) * interpFraction;
        z3 = z1 + (z2 - z1) * interpFraction;
        xyzToLatLon(x3, y3, z3, &lat, &lon);
        if (latitude != NULL)
             *latitude = lat;
        if (longitude != NULL)
            *longitude = lon;
    }

    if (radius != NULL)
        *radius = sat->radius[i1] + (sat->radius[i2] - sat->radius[i1]) * interpFraction;

    if (qdlat != NULL || mlt != NULL)
    {
        float x1 = 0.0, y1 = 0.0, z1 = 0.0;
        float x2 = 0.0, y2 = 0.0, z2 = 0.0;
        float x3 = 0.0, y3 = 0.0, z3 = 0.0;
        float qdlat1 = 0.0;
        float lon1 = 0.0;
        float mlt1 = 0.0;
        latLonToXYZ(sat->qdlatitude[i1], sat->mlt[i1]/24.0 * 360.0, &x1, &y1, &z1);
        latLonToXYZ(sat->qdlatitude[i2], sat->mlt[i2]/24.0 * 360.0, &x2, &y2, &z2);
        x3 = x1 + (x2 - x1) * interpFraction;
        y3 = y1 + (y2 - y1) * interpFraction;
        z3 = z1 + (z2 - z1) * interpFraction;
        xyzToLatLon(x3, y3, z3, &qdlat1, &lon1);

        if (qdlat != NULL)
            *qdlat = qdlat1;
        if (mlt != NULL)
        {
            mlt1 = lon1 / 360.0 * 24.0;
            if (mlt1 < 0.0)
                mlt1 = 24.0 + mlt1;
            *mlt = mlt1;
        }
    }
    if (qdargoforbit != NULL)
    {
        float ao1 = sat->qdargoforbit[i1];
        float ao2 = sat->qdargoforbit[i2];
        *qdargoforbit = ao1 + (ao2 - ao1) * interpFraction;
    }

    return interpFraction;

}

float calculateGradient(float *distances, float *measurements, size_t nMeasurements, int i, int minDistanceIndex)
{
    // Average derivatives using measurements ahead of and behind current position
    // Except at edge of window, where a single derivative is used
    // Derivative is finite difference derivative 

    float gradient = 0.0;
    float gradientBehind = 0.0;
    float gradientAhead = 0.0;

    if (i <= minDistanceIndex)
    {
        if (i >= nMeasurements - 1 || i < minDistanceIndex)
            gradient = nan("");
        else
            gradient = (measurements[i + 1] - measurements[i]) / (distances[i + 1] - distances[i]);
    }
    else if (i >= nMeasurements)
    {
        if (i == minDistanceIndex || i > nMeasurements - 1)
            gradient = nan("");
        else
            gradient = (measurements[i] - measurements[i - 1]) / (distances[i] - distances[i - 1]);
    }
    else
    {
        gradientBehind = (measurements[i] - measurements[i - 1]) / (distances[i] - distances[i - 1]);
        gradientAhead = (measurements[i + 1] - measurements[i]) / (distances[i + 1] - distances[i]);
        gradient = (gradientBehind + gradientAhead) / 2.0;
    }

    return gradient;

}
