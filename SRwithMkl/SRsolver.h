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
// SRsolver.h: interface for the SRsolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRSOLVER_INCLUDED)
#define SRSOLVER_INCLUDED



#include "SRmklUtil.h"
#include "SRutil.h"

class SRFaceForceGroup;
class SRforce;

class SRpardiso
{
#ifndef NOSOLVER
	friend class SRsolver;
public:
	SRpardiso();
	void fillRowNZForSmoothing(int row, int neq, SRmklIntVector* colindr);
	void smoothfillColind(int neq);
	void smoothBookkeep();
	void clear();
	void assemble();
	void smoothAssemble(int elnum, double *elstiff);
	void solve(bool smoothing = false);
	//simpler routines for debugging:
	void fillRowNZSimple(int row, int neq, SRmklIntVector* colindr);
	void bookkeep(bool smoothing = false);
	void fillColindSimple(int neq);
	//for multifaceforce smoothing:
	void fillColindFFG(SRFaceForceGroup *ffg, int neq);
	void fillRowNZFFG(SRFaceForceGroup *ffg, int row, int neq, SRmklIntVector* colindr);
	void bookkeepFFG(SRFaceForceGroup *ffg);
	void assembleFFG(SRFaceForceGroup* ffg, int elnum);
	void solveFFG(SRFaceForceGroup *ffg);
	int getFFGElemStiffLocation(SRforce* faf, int row, int col);

private:
	MKL_INT64 nzTotal;
	SRmklIntVector rowNumels;
	SRmklIntVector rowNZ;
	SRmklIntVector colIndexNonZeroforOneRow;
	SRmklIntVector rowIndex;
	SRmklIntVector colIndex;
	SRmklDoubleVector globalStiff;
	SRvector <SRmklIntVector> AllElEquationNumbers;
	SRvector <SRmklIntVector> AllElSortedFunNumbers;
	SRvector <SRmklIntVector> rowColinds;
	SRvector <SRmklIntVector> rowEls;
	SRmklIntVector nzStore;
	SRvector <SRmklIntVector> AllElpackedEquationNumbers;
	SRmklIntVector elmaxrow;
	SRmklIntVector elminrow;
	double maxInt32;
#endif
};

class SRsolver
{

public:
	void Cleanup();
	void DoSolution();
	SRpardiso* parDisoPtr(){ return &parDisoSolver; };

private:
	SRpardiso parDisoSolver;
};

#endif //!defined(SRSOLVER_INCLUDED)
