/*

    spaceseries: programstate.c

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

#include "programstate.h"

void initProgramState(ProgramState *state)
{
    state->samplePeriod = 0.5;
    state->directories.tctDir = ".";
    state->directories.c7hDir = ".";
    state->directories.modDir = ".";
    state->directories.lpDir = ".";
    state->directories.outDir = ".";
    state->suffix = "";
    state->scaleFactor = 1.0;
    state->binaryOutput = true;
    state->splitOutputFilesByOrbit = true;
    state->labelOutputFileOrbits = true;
    state->qdDataFile = "apexsh.dat";
    state->tctDataset = "_TCT02";
    return;
}

