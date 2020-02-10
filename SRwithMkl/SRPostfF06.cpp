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

////////////////////////////////////////////////////////////////////////////////////////////
//
// SRpostF06.cpp: implementation of the f06 output routines for the  SRpostProcess class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"
#include "SRanalysis.h"
#include "SRelement.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern SRmodel model;
extern SRanalysis analysis;

static char bdfbuf[20];
static void bdfEwrite(double r)
{
	sprintf_s(bdfbuf, "%13.6lE", r);
	int l = strlen(bdfbuf);
	char e1 = bdfbuf[l - 2];
	char e2 = bdfbuf[l - 1];
	//strip off extra 0, e.g E+001 -> e+01
	bdfbuf[l - 3] = e1;
	bdfbuf[l - 2] = e2;
	bdfbuf[l - 1] = '\0';
}

static int brickBdftoSR[20] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15 };
static bool bricks = false;

void SRpostProcess::OutputF06()
{
	//output stress results to .f06 file
	//this only works for all tets or all hexas
	//ttd: implement wedges, then test on bricks and wedges and mixed br and wej	if (analysis.GetAnyWedges())
		return;
	if (analysis.GetAnyBricks())
	{
		if (analysis.GetAnyTets())
			return;
		bricks = true;
	}

	//ttd: to generalize to mixed
	//change "bricks" flag to f06Type (same as SRelementType) and promote to class variable
	//if(anybricks), set f06Type to brick and write stressheader, then write out stresses for all brick elements
	//if(anywedges), set f06Type to wedge and write stressheader, then write out stresses for all wedge elements
	//if(anytets), set f06Type to tet and write stressheader, then write out stresses for all tet elements
	analysis.f06outFile.Open(SRoutputMode);

	int e;
	SRelement* elem = model.GetElement(0);
	int UidAdd = 0;
	if (elem->GetType() == brick)
		bricks = true;
	if (elem->GetUserid() == 0)
		UidAdd = 1;

	int nlines = 0;
	analysis.f06outFile.PrintLine("1    NX NASTRAN STATIC ANALYSIS SET                                            JULY  12, 2014  NX NASTRAN  8/25/13   PAGE    1");
	analysis.f06outFile.PrintLine("                                                                                                                                    ");
	analysis.f06outFile.PrintLine("0                                                                                                                                   ");
	analysis.f06outFile.PrintLine("1    NX NASTRAN STATIC ANALYSIS SET                                            JULY  12, 2014  NX NASTRAN  8/25/13   PAGE    2");
	analysis.f06outFile.PrintLine("                                                                                                                                    ");
	analysis.f06outFile.PrintLine("0                                                                                                                                   ");
	analysis.f06outFile.PrintLine(" ");

	double stressConv = 1.0;
	double lengthConv = 1.0;

	double svm;
	double stress[6];

	//stress header:
	writeHeaderF06();

	double r, s, t, rf, sf;

	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		elem->SetBasisData();
		nlines++;
		int euid = elem->GetUserid() + UidAdd;
		if (bricks)
			analysis.f06outFile.PrintLine("0    %5d           0GRID CS 8 GP", euid);
		else
			analysis.f06outFile.PrintLine("0    %5d           0GRID CS 4 GP", euid);
		nlines++;
		int nn = elem->GetNumNodes();
		int n1, n2;

		model.map.ElementCentroid(elem, r, s, t);
		elem->GetStress(r, s, t, stress);
		svm = model.math.GetSvm(stress);
		svm *= stressConv;
		for (int ss = 0; ss < 6; ss++)
			stress[ss] *= stressConv;
		analysis.f06outFile.Print("0                CENTER");
		bdfEwrite(stress[xxComponent]);
		analysis.f06outFile.Print("  X  %13s", bdfbuf);
		bdfEwrite(stress[xyComponent]);
		analysis.f06outFile.Print("  XY  %13s", bdfbuf);
		analysis.f06outFile.Print("   A   0.0           LX 0.00 1.00 0.00            0.0");
		bdfEwrite(svm);
		analysis.f06outFile.PrintLine("   %13s", bdfbuf);
		nlines++;
		bdfEwrite(stress[yyComponent]);
		analysis.f06outFile.Print("                         Y  %13s", bdfbuf);
		bdfEwrite(stress[yzComponent]);
		analysis.f06outFile.Print("  YZ  %13s", bdfbuf);
		analysis.f06outFile.PrintLine("   B   0.0           LY 0.00 1.00 0.00");
		nlines++;
		bdfEwrite(stress[zzComponent]);
		analysis.f06outFile.Print("                         Z  %13s", bdfbuf);
		bdfEwrite(stress[xzComponent]);
		analysis.f06outFile.Print("  ZX  %13s", bdfbuf);
		analysis.f06outFile.PrintLine("   C   0.0           LZ 0.00 1.00 0.00");
		pageCheckF06(nlines);

		for (int i = 0; i < nn; i++)
		{
			int nid = i;
			if (bricks)
				nid = brickBdftoSR[i];
			int id = elem->GetNodeId(nid);
			SRnode* node = model.GetNode(id);
			for (int c = 0; c < 6; c++)
				stress[c] = stressConv*nodalStress.Get(id, c);
			svm = model.math.GetSvm(stress);
			int nuid = node->GetUserid() + UidAdd;
			analysis.f06outFile.Print("0                %6d", nuid);
			bdfEwrite(stress[xxComponent]);
			analysis.f06outFile.Print("  X  %13s", bdfbuf);
			bdfEwrite(stress[xyComponent]);
			analysis.f06outFile.Print("  XY  %13s", bdfbuf);
			analysis.f06outFile.Print("   A   0.0           LX 0.00 1.00 0.00            0.0");
			bdfEwrite(svm);
			analysis.f06outFile.PrintLine("   %13s", bdfbuf);
			nlines++;
			bdfEwrite(stress[yyComponent]);
			analysis.f06outFile.Print("                         Y  %13s", bdfbuf);
			bdfEwrite(stress[yzComponent]);
			analysis.f06outFile.Print("  YZ  %13s", bdfbuf);
			analysis.f06outFile.PrintLine("   B   0.0           LY 0.00 1.00 0.00");
			nlines++;
			bdfEwrite(stress[zzComponent]);
			analysis.f06outFile.Print("                         Z  %13s", bdfbuf);
			bdfEwrite(stress[xzComponent]);
			analysis.f06outFile.Print("  ZX  %13s", bdfbuf);
			analysis.f06outFile.PrintLine("   C   0.0           LZ 0.00 1.00 0.00");
			pageCheckF06(nlines);
		}
	}
	analysis.f06outFile.Close();

}

static int npage = 2;//skip the fixed two pages in the header
void SRpostProcess::pageCheckF06(int& nlines, bool disp)
{
	//see if output file needs a header for a new "page"
	//this dates back to the impact printers for IBM mainframes.
	//Nastran has been around a while!

	nlines++;
	if (nlines > 39)
	{
		npage++;
		if (disp)
			writeDispHeaderF06();
		else
			writeHeaderF06();
		nlines = 0;
	}
}

void SRpostProcess::writeHeaderF06()
{
	//write stress header for a new "page"

	analysis.f06outFile.PrintLine("  ");
	analysis.f06outFile.PrintLine("1    NX NASTRAN STATIC ANALYSIS SET                                            JULY  12, 2014  NX NASTRAN  8/25/13   PAGE%6d", npage);
	analysis.f06outFile.PrintLine("  ");
	analysis.f06outFile.PrintLine("0 ");
	analysis.f06outFile.PrintLine("  ");
	//stress header:
	if (bricks)
		analysis.f06outFile.PrintLine("                   S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )");
	else
		analysis.f06outFile.PrintLine("                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )");

	analysis.f06outFile.PrintLine("0                CORNER------CENTER AND CORNER POINT STRESSES-------- - DIR.COSINES       MEAN");
	analysis.f06outFile.PrintLine("ELEMENT - ID    GRID - ID        NORMAL              SHEAR             PRINCIPAL - A - -B - -C - PRESSURE       VON MISES");

}

void SRpostProcess::writeDispHeaderF06()
{
	//write displacement header for a new "page"
	analysis.f06outFile.PrintLine("  ");
	analysis.f06outFile.PrintLine("1    NX NASTRAN STATIC ANALYSIS SET                                            JULY  12, 2014  NX NASTRAN  8/25/13   PAGE%6d", npage);
	analysis.f06outFile.PrintLine("  ");
	analysis.f06outFile.PrintLine("0 ");
	analysis.f06outFile.PrintLine("  ");
	analysis.f06outFile.PrintLine("                                             D I S P L A C E M E N T   V E C T O R");
	analysis.f06outFile.PrintLine("                                                                                                                                    ");
	analysis.f06outFile.PrintLine("      POINT ID.TYPE          T1             T2             T3                R1             R2             R3");
}
