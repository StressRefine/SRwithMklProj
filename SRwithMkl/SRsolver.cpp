/*
Copyright (c) 2020 Richard King

The stressRefine analysis executable "SRwithMkl" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"SRwithMkl" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING.txt,
also available at <https://www.gnu.org/licenses/>

Note that in it's current form "SRwithMkl" make use of the pardiso sparse solver
in the Intel MKL library, with which it must be linked.
Copyright (c) 2018 Intel Corporation.
You may use and redistribute the Intel MKL library, without modification, provided the conditions of
the Intel Simplified Software license are met:
https://software.intel.com/en-us/license/intel-simplified-software-license

It is perfectly permissable to replace the use of the pardiso software from the MKL library
with an equivalent open-source solver
*/

//////////////////////////////////////////////////////////////////////
//
// SRsolver.cpp: implementation of the SRsolver class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SRlinux.h"
#include "SRanalysis.h"
#include "SRinput.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRanalysis analysis;

#ifdef _DEBUG
static bool solverecho=true;
#endif

void SRsolver::Cleanup()
{
#ifndef NOSOLVER
	parDisoSolver.clear();
#endif
}

void SRsolver::DoSolution()
{
	//bookkeep, assmble, and solve for the global stiffness matrix
	//version using the Intel pardiso solver:
#ifndef NOSOLVER
	int neq = analysis.GetNumEquations();
	LOGPRINT("Solution Setup: bookkeep, assemble, preanalyze stiffness matrix\n");
	parDisoSolver.bookkeep();
	parDisoSolver.assemble();
	parDisoSolver.solve();
#endif
}
