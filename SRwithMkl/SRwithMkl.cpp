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

//
// SRwithMkl.cpp : Defines the entry point for the console application.
//

#include <vector>
#include "SRanalysis.h"
#include "globalWrappers.h"

#include <stdlib.h>
#ifndef linux
#include <direct.h>
#include <process.h>
#endif

using namespace std;

//One global instance of "analysis""
SRanalysis analysis;

extern SRmodel model;

using namespace std;

#define LOOP_COUNT 10

int main(int argc, char *argv[])
{
	SRstring line;

	char buf[256];
	MDGETCWD(buf, sizeof(buf));
	SCREENPRINT(" _getcwd: %s", buf);
	analysis.wkdir = buf;
	analysis.wkdir += slashStr;

	SCREENPRINT(" wkdir: %s", analysis.wkdir.getStr());
	line = analysis.wkdir;
	line += "engine.log";
	model.logFile.setFileName(line);
	model.logFile.Delete();
	LOGPRINT("wkdir: %s", analysis.wkdir.getStr());

	SRstring tok;

	SRstring tail, filename, foldername;

	line = analysis.wkdir;
	line += "ModelFileName.txt";
	SRfile modelF;
	if (!modelF.Open(line, SRinputMode))
	{
		LOGPRINT("ModelFileName file %s not found\n", line.getStr());
		ERROREXIT;
	}

	modelF.GetLine(foldername);
	modelF.Close();

	foldername.Right(slashChar, tail);
	if (tail.getLength() == 0)
	{
		LOGPRINT("Error in ModelFileName file ");
		ERROREXIT;
	}
	analysis.fileNameTail = tail;
	filename = foldername;
	filename += slashStr;
	filename += tail;
	filename += ".msh";
	if (!analysis.inputFile.Existcheck(filename.getStr()))
	{
		LOGPRINT("input file %s not found\n", filename.getStr());
		exit(0);
	}

	analysis.outdir = foldername;
	analysis.inputFile.setFileName(filename);
	LOGPRINT("\nStressRefine\n model: %s\n", tail.getStr());
	SCREENPRINT("\nStressRefine\n model: %s\n", tail.getStr());

	model.repFile.setFileName(analysis.outdir);
	model.repFile.filename.Cat(slashStr);
	model.repFile.filename.Cat("report.txt");
	if (model.repFile.Existcheck())
		model.repFile.Delete();
	model.repFile.Open(SRoutputMode);
	model.repFile.PrintLine("stressRefine Copyright (C) 2020 Richard B King.");
	model.repFile.PrintLine("This software uses the Intel MKL library Copyright (c) 2018 Intel Corporation.");
	model.repFile.PrintLine("Model: %s\n", analysis.fileNameTail.getStr());
	model.repFile.Close();

	filename.Left('.', analysis.settingsFile.filename);
	analysis.settingsFile.filename.Cat(".srs");

	analysis.CreateModel();

	analysis.Run();

	analysis.CleanUp();

	model.repFile.Open(SRappendMode);
	model.repFile.PrintReturn();
	model.repFile.PrintLine("Adaptive Solution Completed\n");
	model.repFile.Close();
	LOGPRINT("\nStress Refine Run Completed\n");

	LOGPRINT("SUCCESSFUL COMPLETION");
	SCREENPRINT("SUCCESSFUL_COMPLETION");

	return 0;
}

