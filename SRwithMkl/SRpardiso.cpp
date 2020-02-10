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
// SRpardiso.cpp: implementation of the SRpardiso class
//for the intel mkl solver
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SRmodel.h"
#include "SRanalysis.h"
#include "SRinput.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern SRmodel model;
extern SRanalysis analysis;

#ifdef _DEBUG
static bool solverecho = true;
#endif

#define MAXINT32 2E9

SRpardiso::SRpardiso()
{
	maxInt32 = MAXINT32;
}


void SRpardiso::smoothBookkeep()
{
	//bookkeeping for the Intel pardiso solver.
	//version for solving for strain smoothing 

	nzTotal = 0;
	int nfunelmax = model.GetmaxNumElementFunctions();
	SRmklIntVector packedtmp;
	packedtmp.Allocate(nfunelmax);
	int nel = model.GetNumElements();
	AllElEquationNumbers.Allocate(nel);//global equation numbers for each local equation in element, includes "-1" for constrained dofs
	AllElSortedFunNumbers.Allocate(nel); //global equation numbers for each local equation in element, excludes constrained dofs

	int neq = analysis.GetNumSmoothEquations();
	rowNumels.Allocate(neq);
	rowEls.Allocate(neq);

	int rowfun, roweq;

	//first pass through elements: count number of elements that own each equation (rowNumels),
	// and fill up ElEquationNumbers and ElSortedNumbers
	for (int e = 0; e < nel; e++)
	{
		if (!analysis.post.getElSmooth(e))
			continue;
		SRelement* elem = model.GetElement(e);
		int nfun = elem->GetNumFunctions();
		int leq = 0;
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		ElEquationNumbers->Allocate(nfun);

		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = elem->GetFunctionNumber(lfun);
			int geq = model.GetSmoothFunctionEquation(gfun);
			ElEquationNumbers->Put(leq, geq);
			packedtmp.Put(lfun, geq);
			rowNumels.PlusPlus(geq);
			leq++;
		}

		SRmklIntVector* ElSortedFunNumbers = AllElSortedFunNumbers.GetPointer(e);
		ElSortedFunNumbers->Allocate(nfun);
		ElSortedFunNumbers->Copy(packedtmp, nfun);
		ElSortedFunNumbers->Sort();
	}

	packedtmp.Free();

	SRmklIntVector* oneRowEls;

	for (int row = 0; row < neq; row++)
	{
		oneRowEls = rowEls.GetPointer(row);
		int rownel = rowNumels.Get(row);
		oneRowEls->Allocate(rownel);
	}
	rowNumels.Zero();

	//second pass through elements:fill list of elements that own each equation (smoothing) or function (!smoothing) and store as rowEls
	for (int e = 0; e < nel; e++)
	{
		if (!analysis.post.getElSmooth(e))
			continue;
		SRelement* elem = model.GetElement(e);

		int nfun = elem->GetNumFunctions();
		for (int lfun = 0; lfun < nfun; lfun++)
		{
			int gfun = elem->GetFunctionNumber(lfun);
			int geq = model.GetSmoothFunctionEquation(gfun);
			oneRowEls = rowEls.GetPointer(geq);
			int c = rowNumels.Get(geq);
			oneRowEls->Put(c, e);
			rowNumels.PlusPlus(geq);
		}
	}

	smoothfillColind(neq);

	AllElSortedFunNumbers.Free();

	rowNumels.Free();
	for (int row = 0; row < neq; row++)
	{
		oneRowEls = rowEls.GetPointer(row);
		oneRowEls->Free();
	}
	rowEls.Free();
}

void SRpardiso::smoothfillColind(int neq)
{
	//support routine. fill up the column index array
	SRmklIntVector* oneRowEls;
	rowColinds.Allocate(neq);
	rowNZ.Allocate(neq);
	SRmklIntVector* colind1 = &colIndex;
	colind1->Allocate(neq);
	nzStore.Allocate(neq);

	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		fillRowNZForSmoothing(row, neq, colindr);
	}

	nzStore.Free();

	colind1->Free();
	for (int row = 0; row < neq; row++)
		nzTotal += rowNZ.Get(row);

	globalStiff.Allocate((int)nzTotal);
	rowIndex.Allocate(neq + 1);
	rowIndex.Put(0, 1);
	for (int row = 1; row <= neq; row++)
	{
		int rit0 = rowIndex.Get(row - 1);
		int nzt = rowNZ.Get(row - 1);
		int rit = rit0 + nzt;
		rowIndex.Put(row, rit);
	}
	colIndex.Allocate(nzTotal);
	int colg = 0;
	int colind;
	for (int row = 0; row < neq; row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		for (int cnz = 0; cnz < colindr->GetNum(); cnz++)
		{
			colind = (int)colindr->Get(cnz);
			colIndex.Put(colg, colind + 1);//colIndex is passed into parDiso solver, which expects 1-based
			colg++;
		}
	}

}

void SRpardiso::clear()
{
	rowColinds.Free();
	rowIndex.Free();
	colIndex.Free();
	rowNZ.Free();
	AllElEquationNumbers.Free();
}

void SRpardiso::fillRowNZForSmoothing(int row, int neq, SRmklIntVector* colindr)
{
	//fill number of nonzero elements for each row (1 row per global equation)
	//version for solving for smoothing 
	//also fill the nonzero column indices for each row
	//input:
		//row = global equation number
		//neq = number global function
	//output:
		//colindr = column indices for this row
	//notes:
		//fills class variable rownz with number of nonzero elements for the rows associated with rowfun

	SRmklIntVector* colind1 = &colIndex;
	colind1->Zero();
	int *colindOneData = colind1->d;
	int nz = 0;
	SRmklIntVector* oneRowEls = rowEls.GetPointer(row);
	int* nzstoreData = nzStore.d;
	for (int et = 0; et < rowNumels.d[row]; et++)
	{
		int e = oneRowEls->Get(et);
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* ElSortedNumbers = AllElSortedFunNumbers.GetPointer(e);
		int elnz = ElSortedNumbers->GetNum();
		int elroweq = ElSortedNumbers->Find(row);
		int coleq;
		for (int elcoleq = elroweq; elcoleq < elnz; elcoleq++)
		{
			coleq = ElSortedNumbers->Get(elcoleq);
			if (colindOneData[coleq] == 0)
			{
				colindOneData[coleq] = 1;
				nzstoreData[nz] = coleq;
				nz++;
			}
		}
	}
	rowNZ.Put(row, nz);
	colindr->Allocate(nz);
	int* colindrData = colindr->d;
	int coleq;
	for (int c = 0; c < nz; c++)
	{
		coleq = nzstoreData[c];
		colindrData[c] = coleq;
	}
	colindr->Sort();
}

void SRpardiso::smoothAssemble(int elnum, double *elstiff)
{
	//bookkeeping for solver assembly of the element stiffness matrices
	//this version is for solving for strain smoothing 

	SRelement* elem = model.GetElement(elnum);
	SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(elnum);

	int neqel = elem->GetNumFunctions();
	elem->FillStiffDiag(neqel);
	int rowColLocal;
	int rowColGlobal;
	for (int row = 0; row < neqel; row++)
	{
		int roweq = ElEquationNumbers->Get(row);
		if (roweq < 0)
			continue;
		int col0 = rowIndex.Get(roweq) - 1;//row index is one-based
		SRmklIntVector* colindr = rowColinds.GetPointer(roweq);
		for (int col = 0; col < neqel; col++)
		{
			int coleq = ElEquationNumbers->Get(col);
			if (coleq < roweq)
				continue;
			rowColLocal = elem->GetStiffnessLocation(row, col);
			double kelrowcol = elstiff[rowColLocal];
			rowColGlobal = colindr->Find(coleq) + col0;
			globalStiff.PlusAssign(rowColGlobal, kelrowcol);
		}
	}
}


void SRpardiso::assemble()
{
	//bookkeeping for solver assembly of the element stiffness matrices

	double* elstiff = NULL;
	int nel = model.GetNumElements();
	for (int e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		SRmklIntVector* ElEquationNumbers = AllElEquationNumbers.GetPointer(e);
		elstiff = analysis.ReadElementStiffness(elem);
		int neqel = ElEquationNumbers->GetNum();
		int rowColLocal;
		int rowColGlobal;
		for (int row = 0; row < neqel; row++)
		{
			int roweq = ElEquationNumbers->Get(row);
			if (roweq < 0)
				continue;
			int col0 = rowIndex.Get(roweq) - 1;//row index is one-based
			SRmklIntVector* colindr = rowColinds.GetPointer(roweq);
			for (int col = 0; col < neqel; col++)
			{
				int coleq = ElEquationNumbers->Get(col);
				if (coleq < roweq)
					continue;
				rowColLocal = elem->GetStiffnessLocation(row, col);
				double kelrowcol = elstiff[rowColLocal];
				rowColGlobal = colindr->Find(coleq) + col0;
				globalStiff.PlusAssign(rowColGlobal, kelrowcol);
			}
		}
	}

	for (int row = 0; row < analysis.GetNumEquations(); row++)
	{
		SRmklIntVector* colindr = rowColinds.GetPointer(row);
		colindr->Free();
	}
	rowColinds.Free();

}

void SRpardiso::solve(bool smoothing)
{
	//perform decomposition and backsolve of the global stiffness matrix using the Intel pardiso sparse solve
	//input:
	//smoothing = true to do solution on matrix for strain smoothing, false for stiffness matrix

	_MKL_DSS_HANDLE_t parSolveHandle;
	MKL_INT iparm[64];//solver control parameters. iparm[0] = 0 means use defaults
	SRmklIntVector perm;
	MKL_INT neq;
	if (smoothing)
		neq = analysis.GetNumSmoothEquations();
	else
		neq = analysis.GetNumEquations();
	SRdoubleVector solutionTmp;
	solutionTmp.Allocate(neq);
	perm.Allocate(neq);
	MKL_INT *permp = perm.GetVector();
	MKL_INT *ja = colIndex.GetVector();
	MKL_INT *ia = rowIndex.GetVector();
	int ptv[64];
	for (int i = 0; i < 64; i++)
		ptv[i] = 0;
	parSolveHandle = (_MKL_DSS_HANDLE_t)ptv;
	MKL_INT maxfct = 1;
	MKL_INT mnum = 1;
	MKL_INT mtype = 2; // sym posdef
	MKL_INT solvphase;
	MKL_INT nrhs = 1;

	iparm[0] = 0;
	MKL_INT msglvl = 0;//message level, 0 for no stats, 1 for stats to screen;
	MKL_INT error = 0;
	double *a = globalStiff.GetVector();
	double *x = solutionTmp.GetVector();
	double *b = NULL;

	bool oocNeeded = false;

	solvphase = 11; //analyze the matrix;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	int memNeeded = iparm[14]; //see mkl documentation
	int mem2 = iparm[15] + iparm[16];
	if (mem2 > memNeeded)
		memNeeded = mem2;

	double MbytesNeeded = ((double)memNeeded) / 1024.0;

	if (!smoothing)
		LOGPRINT("Memory needed for equation solving (Mbytes): %lg\n", MbytesNeeded);

	double availmem = SRmachDep::availMemCheck();

	LOGPRINT("Pardiso sparse equation solving ");

	if (memNeeded > availmem)
		oocNeeded = true;

	if (oocNeeded)
	{
		LOGPRINT(" not enough memory for incore solution, switching to out of core");
		//not enough memeory. try switching to OOC:
		iparm[59] = 2;//flag for ooc
		solvphase = 11; //reanalyze the matrix in OOC mode:
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	}

	solvphase = 22; //factor the matrix;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	if (error == -2 && !oocNeeded)
	{
		//try again with ooc on. redo analysis and factor in one pass with solvphase 12
		oocNeeded = true;
		iparm[59] = 2;//flag for ooc
		solvphase = 12; //analysize and factor the matrix;
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
	}
	if (error != 0)
	{
		LOGPRINT("pardiso factor phase. error = %d", error);
		if (!smoothing)
		{
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			REPPRINT("Poorly constrained model. Please check constraints and material properties.");
			REPPRINT("If loads are in equilibrium, soft spring stabilization may fix the problem.");
		}
		ERROREXIT;
	}

	solvphase = 33; //backsolve;
	if (smoothing)
	{
		for (int comp = 0; comp < 6; comp++)
		{
			b = analysis.post.getSmoothedStrainVec(comp);
			pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
			for (int i = 0; i < neq; i++)
				b[i] = x[i];
			if (error != 0)
			{
				LOGPRINT("pardiso backsolve phase. error = %d", error);
				ERROREXIT;
			}
		}

	}
	else
	{
		b = analysis.GetSolutionVector();
		pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
		if (error != 0)
		{
			LOGPRINT("pardiso backsolve phase. error = %d", error);
			LOGPRINT("Probable poorly constrained model. Please check constraints and material properties", error);
			REPPRINT("Poorly constrained model. Please check constraints and material properties.");
			REPPRINT("If loads are in equilibrium, soft spring stabilization may fix the problem.");
			ERROREXIT;
		}
		analysis.CopyToSolutionVector(solutionTmp);
	}
	colIndex.Free();

#ifndef _DEBUG
	//release pardiso memory:
	solvphase = -1;
	pardiso(parSolveHandle, &maxfct, &mnum, &mtype, &solvphase, &neq, (void*)a, ia, ja, permp, &nrhs, iparm, &msglvl, (void*)b, (void*)x, &error);
#endif
}
