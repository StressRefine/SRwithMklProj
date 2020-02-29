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
// SRmklUtil.h: interface for the SRmklUtil class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMKLUTIL_INCLUDED)
#define SRMKLUTIL_INCLUDED
#include "SRlinux.h"
#ifndef NOSOLVER
#include "mkl.h"
#include "SRutil.h"

int SRmklIntCompareFunc(const void *v1, const void *v2);
int SRmklInt64CompareFunc(const void *v1, const void *v2);


class SRmklIntVector
{
public:
	bool isEmpty(){ return (d == NULL); };
	int GetNum(){ return num; };
	void SetNum(int n){ num = n; };

	void Free()
	{
		if(num == 0)
			return;
		mkl_free(d);
		d = NULL;
		num = 0;
	};
	void Allocate(int nt)
	{
		if(nt == 0)
			return;
		if (d != NULL)
			Free();
		num = nt;
		d = (int *)mkl_calloc(num, sizeof(int), 32);
		//note: it would be a bit more efficient to use mkl_realloc, but
		//I had run time problems with it.
		Zero();

	};
	void PushBack(int v);
	void Zero()
	{
		for(int i = 0; i < num; i++)
			d[i] = 0;
	};
	void Set(int val)
	{
		for (int i = 0; i < num; i++)
			d[i] = val;
	};
	void Copy(SRmklIntVector& v2, int numt = 0)
	{
		int len;
		if(numt == 0)
			len = v2.GetNum();
		else
			len = numt;
		if(num < len)
			Allocate(len);
		for(int i = 0; i < len; i++)
			d[i] = v2.d[i];
	};
	int* GetVector(){ return d; };
	inline int Get(int i){ return d[i];};
	inline void Put(int i,int di){ d[i] = di; };
	inline void PlusAssign(int i, int di){ d[i] += di; };
	inline void PlusPlus(int i){ d[i] ++; };
	int operator [] (int i) { return Get(i); };
	void Sort();
	int Find(int intIn);

	SRmklIntVector(int nt){ Allocate(nt); };
	SRmklIntVector(){ num = 0; d = NULL; };
	~SRmklIntVector(){	Free();	};

	int* d;
	int num;
};

class SRmklIntVector64
{
public:
	bool isEmpty(){ return (d == NULL); };
	MKL_INT64 GetNum(){ return num; };
	void Free()
	{
		if (num == 0)
			return;
		mkl_free(d);
		d = NULL;
		num = 0;
	};
	void Allocate(int nt)
	{
		MKL_INT64 n64 = (MKL_INT64)nt;
		Allocate(n64);
	}
	void Allocate(MKL_INT64 nt)
	{
		if (nt == 0)
			return;
		if (d != NULL)
			Free();
		num = nt;
		d = (MKL_INT64 *)mkl_calloc(num, sizeof(MKL_INT64), 64);
		//note: it would be a bit more efficient to use mkl_realloc, but
		//I had run time problems with it.
		Zero();
	};
	void Zero()
	{
		for (MKL_INT64 i = 0; i < num; i++)
			d[i] = 0;
	};
	void Set(MKL_INT64 val)
	{
		MKL_INT64 n;
		for (MKL_INT64 i = 0; i < num; i++)
			d[i] = val;
	};
	void Copy(SRmklIntVector64& v2, MKL_INT64 numt = 0)
	{
		MKL_INT64 len;
		if (numt == 0)
			len = v2.GetNum();
		else
			len = numt;
		if (num < len)
			Allocate(len);
		for (MKL_INT64 i = 0; i < len; i++)
			d[i] = v2.d[i];
	};

	MKL_INT64 * GetVector(){ return d; };
	inline MKL_INT64 Get(MKL_INT64 i){ return d[i]; };
	inline void Put(MKL_INT64 i, MKL_INT64 di){ d[i] = di; };
	inline void PlusAssign(MKL_INT64 i, MKL_INT64 di){ d[i] += di; };
	MKL_INT64 operator [] (MKL_INT64 i) { return Get(i); };
	void Sort();
	MKL_INT64 Find(MKL_INT64 intIn);


	SRmklIntVector64(MKL_INT64 nt){ Allocate(nt); };
	SRmklIntVector64(){ num = 0; d = NULL; };
	~SRmklIntVector64(){ Free(); };

	MKL_INT64 * d;
	MKL_INT64 num;
};

class SRmklDoubleVector
{
public:
	bool isEmpty(){ return (d == NULL); };
	void EquateVector(double* v){ d = v; };
	int GetNum(){ return num; };
	void Copy(SRmklDoubleVector& v2){ Copy(v2.d, v2.GetNum()); };
	void PushBack(double v);
	void Free()
	{
		if(num == 0)
			return;
		mkl_free(d);
		d = NULL;
		num = 0;
	};
	void Allocate(int nt)
	{
		if(nt == 0)
			return;
		if (d != NULL)
			Free();
		num = nt;
		d = (double*)mkl_calloc(num, sizeof(double), 64);
		//note: it would be a bit more efficient to use mkl_realloc, but
		//I had run time problems with it.
		Zero();
	};
	void Zero()
	{
		for(int i = 0; i < num; i++)
			d[i] = 0.0;
	};
	void Set(double val)
	{
		for (int i = 0; i < num; i++)
			d[i] = val;
	};
	void Copy(double* v2, int len)
	{
		if(num < len)
			Allocate(len);
		for(int i = 0; i < len; i++)
			d[i] = v2[i];
	};
	double* GetVector(){ return d; };
	inline double Get(int i){ return d[i]; };
	inline void Put(int i, double di){ d[i] = di; };
	inline void PlusAssign(int i, double di){ d[i] += di; };
	double operator [] (int i) { return Get(i); };

	SRmklDoubleVector(int nt){ Allocate(nt); };
	SRmklDoubleVector(){ num = 0; d = NULL; };
	~SRmklDoubleVector(){ Free(); };

private:
	double* d;
	int num;
};

class SRmklDouble64Vector
{
public:
	bool isEmpty(){ return (d == NULL); };
	void EquateVector(double* v){ d = v; };
	MKL_INT64 GetNum(){ return num; };
	void Copy(SRmklDouble64Vector& v2){ Copy(v2.d, v2.GetNum()); };
	void PushBack(double v);
	void Free()
	{
		if (num == 0)
			return;
		mkl_free(d);
		d = NULL;
		num = 0;
	};
	void Allocate(MKL_INT64 nt)
	{
		if (nt == 0)
			return;
		if (d != NULL)
			Free();
		num = nt;
		d = (double*)mkl_calloc(num, sizeof(double), 64);
		//note: it would be a bit more efficient to use mkl_realloc, but
		//I had run time problems with it.
		Zero();
	};
	void Zero()
	{
		for (MKL_INT64 i = 0; i < num; i++)
			d[i] = 0.0;
	};
	void Set(double val)
	{
		for (MKL_INT64 i = 0; i < num; i++)
			d[i] = val;
	};
	void Copy(double* v2, MKL_INT64 len)
	{
		if (num < len)
			Allocate(len);
		for (MKL_INT64 i = 0; i < len; i++)
			d[i] = v2[i];
	};
	double* GetVector(){ return d; };
	inline double Get(MKL_INT64 i){ return d[i]; };
	inline void Put(MKL_INT64 i, double di){ d[i] = di; };
	inline void PlusAssign(MKL_INT64 i, double di){ d[i] += di; };
	double operator [] (MKL_INT64 i) { return Get(i); };

	SRmklDouble64Vector(MKL_INT64 nt){ Allocate(nt); };
	SRmklDouble64Vector(){ num = 0; d = NULL; };
	~SRmklDouble64Vector(){ Free(); };

private:
	double* d;
	MKL_INT64 num;
};

#endif// NOSOLVER

#endif //!defined(SRMKLUTIL_INCLUDED)


