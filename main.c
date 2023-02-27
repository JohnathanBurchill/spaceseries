/*

    spaceseries: main.c

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

#include "main.h"

#include "programstate.h"
#include "data.h"
#include "parseargs.h"
#include "spaceseries.h"
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{

    int status = SPACESERIES_OK;

    ProgramState state = {0};
    initProgramState(&state);
    // Can exit program
    parseArgs(&state, argc, argv);

    status = loadData(&state);
    if (status != SPACESERIES_OK || !state.data.initialized)
    {
        freeData(&state);
        fprintf(stderr, "%s: failed to load data.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    calculateAlongTrackDistances(&state);

    calculateMagneticLatitudeArgOfOrbit(&state);

    status = openOutputFile(&state);
    if (status != SPACESERIES_OK)
    {
        fprintf(stderr, "Unable to open output file for writing.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stdout, "Calculating %d-satellite space series of %s from %.2f Mm to %.1f Mm along track, with respect to position of %s at the first measurement time.\n", state.data.nSatellites, state.measurementParameter, state.minDistance / 1000.0, state.maxDistance / 1000.0, state.data.satelliteData[0]->satellite);

    // Generate space series from first valid time to last valid time at the requested sample period
    generateSpaceSeries(&state);
    fprintf(stdout, "\n");

    closeOutputFile(&state);
    
    freeData(&state);

    return EXIT_SUCCESS;

}
