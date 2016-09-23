/*
 * programExecution.h
 *
 *  Created on: Sep 22 2016
 *      Author: Kyle Kloberdanz
 */

#ifndef PROGRAM_EXECUTION_H
#define PROGRAM_EXECUTION_H

#include "defs.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <sys/wait.h>
#include <stdlib.h>

#include <thread>
#include <mutex>

using namespace std;


void startExecutables(int);

int executeBwaIndex(string sReferenceFile);

int executeBwaAligner(string sReferenceFile, string sReadsFile, string sOutputFile);

int getReads(string sUnalignedFile, string sBwaOutputFile, string sFlag1, string sFlag2, bool bFlagFirst);

int convertSAMtoFASTA(string sUnalignedFile);

int filterOut (string, string, string, string);
#endif /* PROGRAM_EXECUTION_H */
