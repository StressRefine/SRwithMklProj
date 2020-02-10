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
// SRmklUtil.cpp: implementation of the SRmklUtil class.
//
// this class supports use of the Intel mlk pardiso solver
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdlib.h>
#include "SRmklUtil.h"
#include "SRutil.h"

void SRmklIntVector::PushBack(int v)
{
	//append integer v to end of this vector
	SRmklIntVector tmp;
	tmp.Copy(*this);
	num++;
	int nt = num;
	Free();
	num = nt;
	Allocate(num);
	for (int i = 0; i < tmp.GetNum(); i++)
		d[i] = tmp.d[i];
	d[num - 1] = v;
}

int SRmklIntVector::Find(int intIn)
{
	//find an integer in this vector using binary search
	//input:
		//intIn = integer value
	//return:
		//location where intIn resides in this vector, -1 if not found
	//note:
		//Sort must be called before using this routine
	int* intOutPtr = (int *)bsearch(&intIn, d, num, sizeof(int), SRmklIntCompareFunc);
	if (intOutPtr == NULL)
		return -1;
	else
		return intOutPtr - d;
}

void SRmklIntVector::Sort()
{
	//binary sort this vector
	qsort((void *)d, num, sizeof(int), SRmklIntCompareFunc);
}


MKL_INT64 SRmklIntVector64::Find(MKL_INT64 intIn)
{
	//find an integer in this vector using binary search
	//input:
		//intIn = integer value
	//return:
		//location where intIn resides in this vector, -1 if not found
	//note:
		//Sort must be called before using this routine

	MKL_INT64* intOutPtr = (MKL_INT64 *)bsearch(&intIn, d, num, sizeof(MKL_INT64), SRmklIntCompareFunc);
	return intOutPtr - d;
}

void SRmklIntVector64::Sort()
{
	//binary sort this vector
	qsort(d, num, sizeof(MKL_INT64), SRmklInt64CompareFunc);
}


void SRmklDoubleVector::PushBack(double v)
{
	//append integer v to end of this vector
	SRmklDoubleVector tmp;
	tmp.Copy(*this);
	num++;
	Allocate(num);
	for (int i = 0; i < tmp.GetNum(); i++)
		d[i] = tmp.d[i];
	d[num - 1] = v;
};

int SRmklInt64CompareFunc(const void *v1, const void *v2)
{
	//comparison function for binary search
	return *((MKL_INT64 *) v1) - *((MKL_INT64*) v2);
}
int SRmklIntCompareFunc(const void *v1, const void *v2)
{
	//comparison function for binary search
	return *((int *)v1) - *((int *)v2);
}
