/*

    spaceseries: analyze_space_series.c

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

// Equal-area bin statistics adapted from TII-Ion-Drift-Processor
// https://github.com/JohnathanBurchill/TII-Ion-Drift-Processor/tree/tracis_flagging

#include "spaceseries.h"
#include "statistics.h"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fts.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <libgen.h>
#include <errno.h>

#include <cdf.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

typedef struct processingParameters
{
    int nOptions;
    int referenceSatellite;
    bool timeDerivative;
    bool staticAssumption;
    bool staticAssumptionError;
    bool absoluteValue;

    bool verbose;
    bool showFileProgress;

    char *binaryDirectory;
    char *inputFile;

    char *parameter;
    char *statistic;
    
    BinningState binningState;

    char *firstTimeString;
    char *lastTimeString;
    double firstTime;
    double lastTime;
    bool processAllSpaceSeries;

    SpaceSeriesScanHeader scanHeader;
    SpaceSeries spaceSeries;
    long nSpaceSeries;
    long nFileSpaceSeries;

} ProcessingParameters;


void usage(char*);
void about(void);
void parseCommandLine(ProcessingParameters *params, int argc, char *argv[]);
bool fileMatch(FTSENT *e, ProcessingParameters *params);
int processFile(char *filename, ProcessingParameters *params);

int main(int argc, char *argv[])
{

    int status = 0;

    ProcessingParameters params = {0};

    params.referenceSatellite = -1;
    params.timeDerivative = false;
    params.staticAssumption = false;
    params.staticAssumptionError = false;
    params.staticAssumptionError = false;
    params.verbose = false;
    params.showFileProgress = true;
    params.binaryDirectory = ".";
    params.inputFile = NULL;

    // Defaults to equal-area binning
    params.binningState.equalArea = true;
    params.binningState.qdlatmin = 50.0;
    params.binningState.qdlatmax = 90.0;
    params.binningState.deltaqdlat = 5.0;
    params.binningState.mltmin = 0.0;
    params.binningState.mltmax = 24.0;
    params.binningState.deltamlt = 8.0;
    params.binningState.flipParamWhenDescending = false;

    // check options and arguments
    parseCommandLine(&params, argc, argv);

    if (argc - params.nOptions == 4)
    {
        // Process a single file
        params.inputFile = argv[3];
        if (strlen(params.inputFile) == 0 || access(params.inputFile, F_OK) != 0)
        {
            fprintf(stderr, "Input file %s not found.\n", params.inputFile);
            exit(EXIT_FAILURE);
        }
        params.processAllSpaceSeries = true;
    }
    else
    {
        // Process files in directory based on time range
        params.firstTimeString = argv[3];
        params.lastTimeString = argv[4];
        params.firstTime = parseEPOCH4(params.firstTimeString);
        params.lastTime = parseEPOCH4(params.lastTimeString);
    }

    status = initBinningState(&params.binningState);
    if (status != BIN_OK)
        exit(EXIT_FAILURE);

    // Turn off GSL failsafe error handler. Check the GSL return codes.
    gsl_set_error_handler_off();

    long processedFiles = 0;
    float percentDone = 0.0;

    void *mem = NULL;

    if (params.inputFile != NULL)
    {
        processFile(params.inputFile, &params);
    }
    else
    {
        // Analyze files in requested directory
        char *dir[2] = {params.binaryDirectory, NULL};
        long nFiles = 0;
        char *filename = NULL;
        char fullpath[FILENAME_MAX + 1] = {0};
        char *cwdresult = getcwd(fullpath, FILENAME_MAX);
        size_t pathlen = strlen(fullpath);
        if (cwdresult == NULL)
        {
            fprintf(stderr, "Unable to get current directory.\n");
            exit(EXIT_FAILURE);
        }

        // Count files
        FTS *f = fts_open(dir, FTS_PHYSICAL | FTS_NOSTAT, NULL);
        FTSENT *e = fts_read(f);
        while (e != NULL)
        {
            if (fileMatch(e, &params) == true)
                nFiles++;            
            e = fts_read(f);
        }
        fts_close(f);

        int percentCheck = (int) ceil(0.01 * (float)nFiles);

        // Reopen directory listing to do the analysis
        f = fts_open(dir, FTS_PHYSICAL | FTS_NOSTAT, NULL);
        e = fts_read(f);
        while (e != NULL)
        {
            if (fileMatch(e, &params) == true)
            {
                if (params.verbose)
                    fprintf(stderr, "\nAnalyzing %s\n", e->fts_name);

                snprintf(fullpath + pathlen, FILENAME_MAX - pathlen, "/%s", e->fts_path);
                processFile(fullpath, &params);
                processedFiles++;

                if (params.verbose)
                {
                    fprintf(stderr, "\nnSats: %d\n", params.scanHeader.nSatellites);
                    fprintf(stderr, "\n%10s\t%11s\t%11s\n", "Name", "timelag (s)", "distlag (km)\n");
                    for (int i = 0; i < params.scanHeader.nSatellites; i++)
                        fprintf(stderr, "\n%10s\t%11.1lf\t%11.1f\n", params.scanHeader.satellites[i].name, params.scanHeader.satellites[i].timeLag, params.scanHeader.satellites[i].distanceLag);
                }
                if (params.showFileProgress)
                {
                    percentDone = (float)processedFiles / (float)nFiles * 100.0;
                    if (processedFiles % percentCheck == 0)
                        fprintf(stderr, "\r%s: %ld of %ld files processed (%3.0f%%)", argv[0], processedFiles, nFiles, percentDone);
                }
            }

            e = fts_read(f);
        }
        fts_close(f);

        if (params.showFileProgress)
            fprintf(stderr, "\r\n");
    }

    printBinningResults(&params.binningState, params.parameter, params.statistic);

    freeBinStorage(&params.binningState);
    free(params.spaceSeries.header.satelliteTimeSeriesPoint);
    free(params.spaceSeries.points);

    return EXIT_SUCCESS;
}

void usage(char *name)
{
    fprintf(stdout, "usage: %s <parameter> <statistic> <startDate> <stopDate>\n", name);
    fprintf(stdout, "%35s - %s\n", "--help", "print this message");
    fprintf(stdout, "%35s - %s\n", "--about", "print program and license info");
    fprintf(stdout, "%35s - %s\n", "--verbose", "extra processing information");
    fprintf(stdout, "%35s - %s\n", "--time-derivative", "get statistic of the time derivative");
    fprintf(stdout, "%35s - %s\n", "--static-assumption=<satnum>", "get statistic of the result from assuming <satnmu> traversed a static structure");
    fprintf(stdout, "%35s - %s\n", "--static-assumption-error=<satnum>", "get statistic of the change in static estimate for <satnum> with respect to the space series value");
    fprintf(stdout, "%35s - %s\n", "--available-statistics", "print list of supported statistics");
    fprintf(stdout, "%35s - %s\n", "--no-file-progress", "suppress printing file progress");
    fprintf(stdout, "%35s - %s\n", "--equal-length-bins", "use standing binning rather than equal-area");
    fprintf(stdout, "%35s - %s\n", "--qdlatmin=<value>", "minimum quasi-dipole magnetic latitude");
    fprintf(stdout, "%35s - %s\n", "--qdlatmax=<value>", "maximum quasi-dipole magnetic latitude");
    fprintf(stdout, "%35s - %s\n", "--deltaqdlat=<value>", "quasi-dipole magnetic latitude bin width");
    fprintf(stdout, "%35s - %s\n", "--mltmin=<value>", "minimum magnetic local time");
    fprintf(stdout, "%35s - %s\n", "--mltmax=<value>", "maximum magnetic local time");
    fprintf(stdout, "%35s - %s\n", "--deltamlt=<value>", "magnetic local time bin width (at the polar cap if for equal-area binning)");
    fprintf(stdout, "%35s - %s\n", "--flip-when-descending", "change sign of value when on descending part of orbit");
    fprintf(stdout, "%35s - %s\n", "--binary-input-directory=<dir>", "path to directory containing binary input files");

    return;
}

void about(void)
{
    fprintf(stdout, "analyze_space_series\n");
    fprintf(stdout, "Copyright (C) 2023 Johnathan K. Burchill.\n");
    fprintf(stdout, " Analyzes dynamic space series.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
    fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
    fprintf(stdout, "under the terms of the GNU General Public License.\n");
    fprintf(stdout, "See the file LICENSE in the source repository for details.\n");
    return;
}

void parseCommandLine(ProcessingParameters *params, int argc, char *argv[])
{

    for (int i = 0; i < argc; i++)
    {
        if (strcmp("--help", argv[i]) == 0)
        {
            usage(argv[0]);
            exit(EXIT_SUCCESS);
        }
        else if (strcmp("--about", argv[i])==0)
        {
            // call about function then exit
            about();
            exit(EXIT_SUCCESS);
        }
        else if (strcmp("--verbose", argv[i])==0)
        {
            params->nOptions++;
            params->verbose = true;
        }
        else if (strcmp("--time-derivative", argv[i])==0)
        {
            params->nOptions++;
            params->timeDerivative = true;
        }
        else if (strcmp("--absolute-value", argv[i])==0)
        {
            params->nOptions++;
            params->absoluteValue = true;
        }
        else if (strncmp("--static-assumption=", argv[i], 20)==0)
        {
            params->nOptions++;
            params->staticAssumption = true;
            if (strlen(argv[i]) < 21)
            {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            params->referenceSatellite = atoi(argv[i] + 20);
        }
        else if (strncmp("--static-assumption-error=", argv[i], 26)==0)
        {
            params->nOptions++;
            params->staticAssumptionError = true;
            if (strlen(argv[i]) < 27)
            {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            params->referenceSatellite = atoi(argv[i] + 26);
        }
        else if (strcmp(argv[i], "--available-statistics") == 0)
        {
            fprintf(stdout, "Available statistics:\n");
            printAvailableStatistics(stdout);
            exit(EXIT_SUCCESS);
        }
        else if (strcmp(argv[i], "--no-file-progress") == 0)
        {
            params->nOptions++;
            params->showFileProgress = false;
        }
        else if (strcmp(argv[i], "--equal-length-bins") == 0)
        {
            params->nOptions++;
            params->binningState.equalArea = false;
        }
        else if (strncmp(argv[i], "--qdlatmin=", 11) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 11)
                params->binningState.qdlatmin = atof(argv[i]+11);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--qdlatmax=", 11) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 11)
                params->binningState.qdlatmax = atof(argv[i]+11);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--deltaqdlat=", 13) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 13)
                params->binningState.deltaqdlat = atof(argv[i]+13);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--mltmin=", 9) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 9)
                params->binningState.mltmin = atof(argv[i]+9);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--mltmax=", 9) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 9)
                params->binningState.mltmax = atof(argv[i]+9);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "--deltamlt=", 11) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 11)
                params->binningState.deltamlt = atof(argv[i]+11);
            else
            {
                fprintf(stderr, "Could not parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--flip-when-descending") == 0)
        {
            params->nOptions++;
            params->binningState.flipParamWhenDescending = true;
        }
        else if (strncmp(argv[i], "--binary-input-directory=", 25) == 0)
        {
            params->nOptions++;
            if (strlen(argv[i]) > 25)
                params->binaryDirectory = argv[i] + 25;
            else
            {
                fprintf(stderr, "Unable to parse %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strncmp("--", argv[i], 2) == 0)
        {
            fprintf(stderr, "Unrecognized option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if ((argc - params->nOptions != 4) && (argc - params->nOptions != 5))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (params->staticAssumption == true && params->staticAssumptionError == true)
    {
        fprintf(stderr, "Choose --static-assumption=<satIndex> or --static-assumption-error=<satIndex> but not both.\n");
        exit(EXIT_FAILURE);
    }

    params->parameter = argv[1];
    params->statistic = argv[2];

    if (!validStatistic(params->statistic))
    {
        fprintf(stderr, "Invalid statistic '%s'\n", params->statistic);
        fprintf(stderr, "Must be one of:\n");
        printAvailableStatistics(stderr);
        exit(EXIT_FAILURE);
    }

    return;
}

bool fileMatch(FTSENT *e, ProcessingParameters *params)
{
    // Filename pattern is a <parameter> '_' <startDateTime> '_' <stopDateTime> ['_orbit_??'] '.bin'
    // Split filename by "_"
    // Adapted from the strsep manual page example
    char *inputstring = strdup(e->fts_name);
    if (inputstring == NULL)
        return false;

    char **ap = NULL;
    char *argv[10] = {0};

    for (ap = argv; (*ap = strsep(&inputstring, "_")) != NULL;)
        if (**ap != '\0')
                if (++ap >= &argv[10])
                    break;
   
    if (argv[0] == NULL || argv[1] == NULL || argv[2] == NULL)
        return false;

    if (strcmp(argv[0], params->parameter) != 0)
        return false;

    long y1, m1, d1, h1, min1, s1;
    long y2, m2, d2, h2, min2, s2;
    int nAssigned = 0;
    nAssigned = sscanf(argv[1], "%4ld%2ld%2ldT%2ld%2ld%2ld", &y1, &m1, &d1, &h1, &min1, &s1);
    if (nAssigned != 6)
        return false;
    nAssigned = sscanf(argv[2], "%4ld%2ld%2ldT%2ld%2ld%2ld", &y2, &m2, &d2, &h2, &min2, &s2);
    if (nAssigned != 6)
        return false;
    
    double fileFirstTime = computeEPOCH(y1, m1, d1, h1, min1, s1, 0);
    if (fileFirstTime == ILLEGAL_EPOCH_VALUE)
        return false;
    double fileLastTime = computeEPOCH(y2, m2, d2, h2, min2, s2, 0);
    if (fileLastTime == ILLEGAL_EPOCH_VALUE)
        return false;

    double firstTime = params->firstTime;
    double lastTime = params->lastTime;

    bool match = ((strcmp(".bin", e->fts_name + e->fts_namelen - 4) == 0) && ((firstTime >= fileFirstTime && firstTime <= fileLastTime) || (lastTime >= fileFirstTime && lastTime <= fileLastTime) || (firstTime < fileFirstTime && lastTime > fileLastTime)));

    free(inputstring);

    return match;
}

int processFile(char *filename, ProcessingParameters *params)
{
    if (filename == NULL || params == NULL)
        return SPACESERIES_POINTERS;

    void *mem = NULL;

    float mlt = 0.0;
    float qdlat = 0.0;
    float lastQdLat = 0.0;
    float value = 0.0;
    bool includeValue = false;

    ssize_t bytesRead = 0;

    int binfd = open(filename, O_RDONLY);
    int status = readSpaceSeriesScanHeader(binfd, &params->scanHeader);
    if (status != 0)
        return status;
    if (params->referenceSatellite > params->scanHeader.nSatellites - 1)
    {
        fprintf(stderr, "Unknown requested reference satellite %d\n", params->referenceSatellite);
        return SPACESERIES_ARGUMENTS;
    }
    if (params->scanHeader.nSatellites > params->spaceSeries.header.nTimeSeriesPoints)
    {
        mem = realloc(params->spaceSeries.header.satelliteTimeSeriesPoint, params->scanHeader.nSatellites * sizeof *params->spaceSeries.header.satelliteTimeSeriesPoint);
        if (mem == NULL)
        {
            fprintf(stderr, "Could not allocate memory.\n");
            exit(EXIT_FAILURE);
        }
        params->spaceSeries.header.satelliteTimeSeriesPoint = mem;
    }
    params->spaceSeries.header.nTimeSeriesPoints = params->scanHeader.nSatellites;

    status = SPACESERIES_OK;
    float qdDirection = 0.0;

    // TODO allow thresholds to exclude values
    includeValue = true;

    while (true)
    {
        status = readSpaceSeries(binfd, &params->spaceSeries);
        if (status != SPACESERIES_OK)
            break;

        // Check if this space series's time is within the requested time interval
        if (!params->processAllSpaceSeries && params->spaceSeries.header.timeStamp < params->firstTime)
            continue;
        else if (!params->processAllSpaceSeries && params->spaceSeries.header.timeStamp > params->lastTime)
            break;

        params->nSpaceSeries++;

        // Calculate space series statistic
        if (params->staticAssumption && params->referenceSatellite < params->spaceSeries.header.nTimeSeriesPoints)
        {
            value = params->spaceSeries.header.satelliteTimeSeriesPoint[params->referenceSatellite].staticParam;
            qdlat = params->spaceSeries.header.satelliteTimeSeriesPoint[params->referenceSatellite].qdLatitude;
            mlt = params->spaceSeries.header.satelliteTimeSeriesPoint[params->referenceSatellite].mlt;
            if (params->nFileSpaceSeries > 1)
                qdDirection = qdlat - lastQdLat;
            else
                qdDirection = 0.0;
            lastQdLat = qdlat;

            if (params->absoluteValue)
                value = fabs(value);
            else if (params->binningState.flipParamWhenDescending && qdDirection < 0.0)
                value = -value;
            params->binningState.nValsRead++;
            binData(&params->binningState, qdlat, mlt, value, includeValue);
        }
        else
            for (int i = 0; i < params->spaceSeries.header.nSpaceSeriesPoints; i++)
            {
                if (params->timeDerivative)
                {
                    value = params->spaceSeries.points[i].dparamdt;
                }
                else if (params->staticAssumptionError)
                {
                    if (params->referenceSatellite < params->spaceSeries.header.nTimeSeriesPoints)
                        value = params->spaceSeries.points[i].staticAssumptionErrors[params->referenceSatellite];
                    else
                        value = NAN;
                }
                else
                    value = params->spaceSeries.points[i].param;

                qdlat = params->spaceSeries.points[i].qdLatitude;
                mlt = params->spaceSeries.points[i].mlt;
                if (params->nFileSpaceSeries > 1)
                    qdDirection = qdlat - lastQdLat;
                else
                    qdDirection = 0.0;
                lastQdLat = qdlat;

                if (params->binningState.flipParamWhenDescending && qdDirection < 0.0)
                    value = -value;
                params->binningState.nValsRead++;
                binData(&params->binningState, qdlat, mlt, value, includeValue);
            }

    }
next:
    if (binfd != -1)
        close(binfd);

    return SPACESERIES_OK;
}
