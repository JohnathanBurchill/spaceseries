/*

    spaceseries: parseargs.c

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

#include "parseargs.h"
#include "spaceseries.h"
#include "util.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cdf.h>

void parseArgs(ProgramState *state, int argc, char *argv[])
{
    int status = SPACESERIES_OK;
    
    for (int i = 0; i < argc; i++)
    {
        if (strcmp("--about", argv[i])==0)
        {
            // call about function then exit
            about();
            exit(EXIT_SUCCESS);
        }
        else if (strcmp("--help", argv[i])==0)
        {
            usage(argv[0]);
            exit(EXIT_SUCCESS);
        }
        else if (strncmp("--qd-data-file=", argv[i], 15) == 0)
        {
            state->nOptions++;
            state->qdDataFile = argv[i] + 15;
        }
        else if (strcmp("--ascii-output", argv[i])==0)
        {
            state->nOptions++;
            state->binaryOutput = false;
        }
        else if (strcmp("--single-output-file", argv[i])==0)
        {
            state->nOptions++;
            state->splitOutputFilesByOrbit = false;
        }
        else if (strcmp("--no-orbit-numbers", argv[i])==0)
        {
            state->nOptions++;
            state->labelOutputFileOrbits = false;
        }
        else if (strncmp("--output-directory=", argv[i], 19) == 0)
        {
            state->nOptions++;
            state->directories.outDir = argv[i] + 19;
        }
        else if (strncmp("--suffix=", argv[i], 9) == 0)
        {
            state->nOptions++;
            state->suffix = argv[i] + 9;
        }
        else if (strcmp("--along-track-gradient", argv[i])==0)
        {
            state->nOptions++;
            state->alongTrackGradient = true;
        }
        else if (strcmp("--constant-interpolation", argv[i])==0)
        {
            state->nOptions++;
            state->constantInterpolation = true;
        }
        else if (strncmp("--sample-period=", argv[i], 16) == 0)
        {
            state->nOptions++;
            char *lastPos = argv[i] + 16;
            double period = strtod(argv[i] + 16, &lastPos);
            if (lastPos != argv[i] + 16)
            {
                state->samplePeriod = period;
            }
            else
            {
                fprintf(stderr, "Unable to parse requested %s.\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp("--tct-directory=", argv[i], 16) == 0)
        {
            state->nOptions++;
            state->directories.tctDir = argv[i] + 16;
        }
        else if (strncmp("--c7h-directory=", argv[i], 16) == 0)
        {
            state->nOptions++;
            state->directories.c7hDir = argv[i] + 16;
        }
        else if (strncmp("--mod-directory=", argv[i], 16) == 0)
        {
            state->nOptions++;
            state->directories.modDir = argv[i] + 16;
        }
        else if (strncmp("--lp-directory=", argv[i], 15) == 0)
        {
            state->nOptions++;
            state->directories.lpDir = argv[i] + 15;
        }
        else if (strncmp("--scale-factor=", argv[i], 15) == 0)
        {
            state->nOptions++;
            status = getOptionValue(argv[i] + 15, &state->scaleFactor);
            if (status != SPACESERIES_OK)
            {
                fprintf(stderr, "Could not parse option %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp("--sat1-offset=", argv[i], 14) == 0)
        {
            state->nOptions++;
            status = getOptionValue(argv[i] + 14, &(state->parameterOffsets[0]));
            if (status != SPACESERIES_OK)
            {
                fprintf(stderr, "Could not parse option %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp("--sat2-offset=", argv[i], 14) == 0)
        {
            state->nOptions++;
            status = getOptionValue(argv[i] + 14, &(state->parameterOffsets[1]));
            if (status != SPACESERIES_OK)
            {
                fprintf(stderr, "Could not parse option %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp("--sat3-offset=", argv[i], 14) == 0)
        {
            state->nOptions++;
            status = getOptionValue(argv[i] + 14, &(state->parameterOffsets[2]));
            if (status != SPACESERIES_OK)
            {
                fprintf(stderr, "Could not parse option %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp("--", argv[i], 2) == 0)
        {
            fprintf(stderr, "Unknown option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if ((argc - state->nOptions) == 2)
    {
        if (strncmp(argv[1]+strlen(argv[1]) - 47, "C7H", 3) == 0)
        {
            fprintf(stdout, "Available parameters:\n");
            fprintf(stdout, " dbx\n");
            fprintf(stdout, " dby\n");
            fprintf(stdout, " dbz\n");
            fprintf(stdout, " jz\n");
            exit(EXIT_SUCCESS);
        }
        CDFid cdf = NULL;
        CDFstatus s = CDFopen(argv[1], &cdf);
        if (s != CDF_OK)
        {
            fprintf(stderr, "Unable to read %s\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        long nVars = 0;
        char varName[CDF_VAR_NAME_LEN+1] = {0};
        s = CDFgetNumzVars(cdf, &nVars);
        if (s != CDF_OK)
        {
            fprintf(stderr, "Unable to get number of variables in %s\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        fprintf(stdout, "Available parameters:\n");
        for (int i = 0; i < nVars; i++)
        {
            s = CDFgetzVarName(cdf, i, varName);
            if (s == CDF_OK)
            {
                if (strcmp(varName, "Timestamp") != 0 && strcmp(varName, "Latitude") && strcmp(varName, "Longitude") != 0 && strcmp(varName, "Height") != 0 && strcmp(varName, "Radius") != 0 && strcmp(varName, "SZA") != 0 && strcmp(varName, "SAz") != 0 && strcmp(varName, "ST") != 0 && strcmp(varName, "Diplat") != 0 && strcmp(varName, "Diplon") != 0 && strcmp(varName, "MLat") != 0 && strcmp(varName, "QDLatitude") != 0 && strcmp(varName, "MLT") != 0 && strcmp(varName, "AACGMLat") != 0 && strcmp(varName, "AACGMLon") != 0 && strcmp(varName, "Flagbits") != 0 && strcmp(varName, "Quality_flags") != 0 && strcmp(varName, "Calibration_flags") != 0 && strcmp(varName, "VsatN") != 0 && strcmp(varName, "VsatE") != 0 && strcmp(varName, "VsatC") != 0)
                    fprintf(stdout, " %s\n", varName);
            }
        }
        exit(EXIT_SUCCESS);
    }

    if (argc - state->nOptions != 4)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (access(state->qdDataFile, F_OK) != 0)
    {
        fprintf(stderr, "QuasiDipole spherical harmonic data file %s not found.\n", state->qdDataFile);
        fprintf(stderr, "Place something like \"apexsh.dat\" in the working directory, or specify a path with \"--qd-data-file=<filename>\".\n");
        exit(EXIT_FAILURE);
    }


    state->measurementParameter = argv[1];
    state->startString = argv[2];
    state->stopString = argv[3];
    state->firstTime = parseEPOCH4(state->startString);
    state->lastTime = parseEPOCH4(state->stopString);

    if (state->firstTime == ILLEGAL_EPOCH_VALUE || state->lastTime == ILLEGAL_EPOCH_VALUE)
    {
        fprintf(stderr, "%s: could not parse time format as yyyy-mm-ddThh:mm:ss.ccc\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    return;
}

void usage(char *commandName)
{
    fprintf(stdout, "Usage:\n");
    fprintf(stdout, "First form: %s <parameter> <firstTime> <lastTime> [--ascii-output] [--single-output-file] [--no-orbit-numbers] [--output-directory=<dir>] [--sufix=<suffix>] [--tct-directory=<dir>] [--qd-data-file=<filename>] [--c7h-directory=<dir>] [--mod-directory=<dir>] [--lp-directory=<dir>] [--scale-factor=<factor>] [--sat1-offset=<value>] [--sat2-offset=<value>] [--sat3-offset=<value>] [--sample-period=seconds] [--along-track-gradient] [--about] [--help]\n", commandName);
    fprintf(stdout, "Reduces data read from Swarm CDF files in working directory for estimating the requested parameter's space series.\n");    
    fprintf(stdout, " <parameter> One of the EXPT TCT file variables or one of dbx, dby, dbz.\n");
    fprintf(stdout, " time format: yyyy-mo-ddThh:mm:ss.ccc\n");
    fprintf(stdout, " --ascii-output\tFormatted ascii output file. Default is \"unformatted binary output\"\n");
    fprintf(stdout, " --single-output-file\tGenerate a single output file. Default is to split output file by orbits starting at QD Lat = 0 deg ascending.\n");
    fprintf(stdout, " --no-orbit-numbers\tSuppress orbit numbers. Default is to append an orbit number to the filename.\n");
    fprintf(stdout, " --output-directory=<dir>\tSet location of output file. Default \".\"\n");
    fprintf(stdout, " --suffix=<name>\tAppend <name> to filename. Default no suffix\n");
    fprintf(stdout, " --qd-data-file=<filename>\tSet path to apexsh.dat file used for QuasiDipole magnetic coordinate calculations. Default path is \"./apexsh.dat\".\n");
    fprintf(stdout, " --tct-directory=<dir>\tSet location of TCT files. Default \".\"\n");
    fprintf(stdout, " --c7h-directory=<dir>\tSet location of C7H CHAOS 7 magnetic residual files. Default \".\"\n");
    fprintf(stdout, " --mod-directory=<dir>\tSet location of MOD files. Default \".\"\n");
    fprintf(stdout, " --lp-directory=<dir>\tSet location of EXTD LP_HM files. Default \".\"\n");
    fprintf(stdout, " --scale-factor=<factor>\tMultiply parameter by <factor> after removing offsets.\n");
    fprintf(stdout, " --sat1-offset=<value>\tRemove <value> from each raw measurement parameter for the leading satellite.\n");
    fprintf(stdout, " --sat2-offset=<value>\tRemove <value> from each raw measurement parameter for the second satellite.\n");
    fprintf(stdout, " --sat3-offset=<value>\tRemove <value> from each raw measurement parameter for the third satellite.\n");
    fprintf(stdout, " --sample-period=seconds\tgenerates the interpolated space series at \"seconds\" intervals.");
    fprintf(stdout, " --along-track-gradient\tEstimate along-track gradient of requested measurement parameter for individual satellites assuming Doppler-shifted stationary structure.\n");

    fprintf(stdout, "Second form: %s <cdfFile>\n", commandName);
    fprintf(stdout, "Prints parameter names for <cdfFile>\n");
    fprintf(stdout, "General options:\n");

    fprintf(stdout, " --about\tprint program summary and license info.\n");
    fprintf(stdout, " --help\tprint this message.\n");

    return;
}

void about(void)
{
    fprintf(stdout, "swarm_scan\n");
    fprintf(stdout, "Copyright (C) 2023 Johnathan K. Burchill.\n");
    fprintf(stdout, " Estimates dynamic space series from Swarm electric and magnetic\n field measurements from all three satellites during the\n pearls-on-a-string phase of the mission.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
    fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
    fprintf(stdout, "under the terms of the GNU General Public License.\n");
    fprintf(stdout, "See the file LICENSE in the source repository for details.\n");
    return;
}
