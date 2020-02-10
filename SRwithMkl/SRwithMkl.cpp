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

#include "stdafx.h"
#include <vector>
#include "SRanalysis.h"
#include <direct.h>
#include <stdlib.h>
#include <process.h>


using namespace std;


#define LOOP_COUNT 10


//One global instance of "analysis" and "model"
SRanalysis analysis;
SRmodel model;

using namespace std;


#define LOOP_COUNT 10


int main(int argc, char *argv[])
{
	_set_output_format(_TWO_DIGIT_EXPONENT);

	SRstring line;
	bool batchMode = false;

	char buf[256];
	_getcwd(buf, sizeof(buf));
	SCREENPRINT(" _getcwd: %s", buf);
	analysis.wkdir = buf;
	analysis.wkdir += "\\";
	SCREENPRINT(" wkdir: %s", analysis.wkdir.str);
	line = analysis.wkdir;
	line += "engine.log";
	model.logFile.setFileName(line);
	model.logFile.Delete();
	SRmachDep::GetTime(line);
	LOGPRINT("%s", line.str);
	LOGPRINT("wkdir: %s", analysis.wkdir.str);

	SRstring tok;

	SRstring tail, filename, foldername;

	line = analysis.wkdir;
	line += "ModelFileName.txt";
	SRfile modelF;
	if (!modelF.Open(line, SRinputMode))
	{
		LOGPRINT("ModelFileName file %s not found\n", line.str);
		ERROREXIT;
	}

	modelF.GetLine(foldername);
	modelF.Close();

	foldername.Right('\\', tail);
	if (tail.len == 0)
	{
		LOGPRINT("Error in ModelFileName file ");
		ERROREXIT;
}
	analysis.fileNameTail = tail;
	filename = foldername;
	filename += "\\";
	filename += tail;
	filename += ".msh";
	if (!analysis.inputFile.Existcheck(filename.str))
	{
		LOGPRINT("input file %s not found\n", filename.str);
		exit(0);
	}

	analysis.outdir = foldername;
	analysis.inputFile.setFileName(filename);
	LOGPRINT("\nStressRefine\n model: %s\n", tail.str);
	SCREENPRINT("\nStressRefine\n model: %s\n", tail.str);

	model.reportFile.setFileName(analysis.outdir);
	model.reportFile.catFileName("\\report.txt");
	if (model.reportFile.Existcheck())
		model.reportFile.Delete();
	model.reportFile.Open(SRoutputMode);
	SRmachDep::GetTime(line);
	model.reportFile.PrintLine("stressRefine Copyright (C) 2020 Richard B King.");
	model.reportFile.PrintLine("This software uses the Intel MKL library Copyright (c) 2018 Intel Corporation.");
	model.reportFile.PrintLine("%s   Model: %s\n", line.str, analysis.fileNameTail.str);
	model.reportFile.Close();

	filename.Left('.', analysis.settingsFile.getFileName());
	analysis.settingsFile.catFileName(".srs");
	filename.Left('.', line);
	line.Cat(".out");
	model.setOutFileName(line.str);
	if (model.outFile.Existcheck())
		model.outFile.Delete();

	analysis.CreateModel();

	analysis.Run();

	analysis.CleanUp();

	model.reportFile.Open(SRappendMode);
	model.reportFile.PrintReturn();
	model.reportFile.PrintLine("Adaptive Solution Completed\n");
	model.reportFile.Close();
	LOGPRINT("\nStress Refine Run Completed\n");
	SRanalysis::TimeStamp();

	LOGPRINT("SUCCESSFUL COMPLETION");

	if (!batchMode)
	{
		SCREENPRINT("hit any char to exit");
		int c = getchar();
	}

	return 0;
}

