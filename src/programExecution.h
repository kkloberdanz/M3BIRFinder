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
#include <cstdint>

#include <thread>
#include <mutex>

using namespace std;


void startExecutables(int64_t);

int64_t executeBwaIndex(string sReferenceFile);

int64_t executeBwaAligner(string sReferenceFile, string sReadsFile, string sOutputFile);

int64_t getReads(string sUnalignedFile, string sBwaOutputFile, string sFlag1, string sFlag2, bool bFlagFirst);

int64_t convertSAMtoFASTA(string sUnalignedFile);

int64_t filterOut (string, string, string, string);
#endif /* PROGRAM_EXECUTION_H */
