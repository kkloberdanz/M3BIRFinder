//============================================================================
// Name        : main.cpp
// Author      : Matt Segar
// Version     :
// Copyright   : 
// Description : Main program
//============================================================================

#define MAIN_CPP_

#include "defs.h"
#include "globals.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <unistd.h>

#include <ctime>
#include <thread>

#include "programExecution.h"
#include "BWA.hpp" 

unsigned int NUM_THREADS = 1;
std::string  PROGRAM_NAME;

using namespace std;

/*
void test_threads(int i) {
    startExecutables(i);
}
*/

double get_end_time_ms(clock_t start_time) {
    return 1000*(double)(clock() - start_time)/CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {
    PROGRAM_NAME = argv[0];
    clock_t start_time = clock();
    // start timing
    time0 = time(NULL);

    // get command line arguments
    start_time = clock();
    int cl = processCommandLine(argc, argv);
    if (cl != 0)
        cout << "Command line exited with: " << cl << endl;

    cout << "COMPLETED: processCommandLine(): " << get_end_time_ms(start_time) << "ms" << endl;

    start_time = clock();
    sBaseFileName = setupBaseFileName();
    cout << "COMPLETED: setupBaseFileName(): " << get_end_time_ms(start_time) << "ms" << endl;

    // Set up general log file and project directory
    start_time = clock();
    prepareLogFile(fLogFileOut,sBaseFileName,confDB.getKey("sLogDirName").stringVal);
    sProjectDirectory = confDB.getKey("projectDirectory").stringVal;
    cout << "COMPLETED: prepareLogFile(): " << get_end_time_ms(start_time) << "ms" << endl;

    // Print the method configuration parameters
    printConfig(fLogFileOut);

    start_time = clock();
    // execute the main alignment programs (BWA and samtools)
    
    /* Changes start here */

    // run this over several threads
    /*
    const int NUM_THREADS = sReadsFile_v.size();
    */
    //bir::BWA bwa_align;
    //bwa_align.startExecutables(0);
    startExecutables(1);

    /*
    cout << "Making threads" << endl;
    thread t[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; ++i) {
        t[i] = thread(test_threads, i);
    }

    //int exec = startExecutables();

    // join threads here
    for (int i = 0; i < NUM_THREADS; ++i) {
        t[i].join();
    } 
    cout << "Threads rejoined" << endl;
    */

    /*
    if (exec != 0)
        cout << "Exec exited with: " << exec << endl;
    cout << "COMPLETED: startExecutables(): " << get_end_time_ms(start_time) << "ms" << endl;
    */

    /* Changes end here */

    start_time = clock();
    // read in reference genome
    readInReferenceGenome(); 
    cout << "COMPLETED: readInReferenceGenome(): " << get_end_time_ms(start_time) << "ms" << endl;

    start_time = clock();
    // find reads with one half anchored and the other unaligned
    int cand = startCandidateReads();
    if (cand != 0)
        cout << "startCandidateReads exited with: " << cand << endl;
    cout << "COMPLETED: startCandidateReads(): " << get_end_time_ms(start_time) << "ms" << endl;

    start_time = clock();
    // consolidate reads to narrow down BIR locations
    int cons = startConsolidate();
    if (cons != 0)
        cout << "startConsolidate exited with: " << cons << endl;
    cout << "COMPLETED: startConsolidate(): " << get_end_time_ms(start_time) << "ms" << endl;

    start_time = clock();
    // perform alignments to get candidate BIR locations
    int bir = startBirFinder();
    if (bir != 0)
        cout << "startBirFinder exited with: " << bir << endl;
    cout << "COMPLETED: startBirFinder(): " << get_end_time_ms(start_time) << "ms" << endl;

    start_time = clock();
    // we have the possible bir strings and their respective locations so let's find the template to confirm
    int temp = startTemplateFinder();
    if (temp != 0)
        cout << "startTemplateFinder exited with: " << temp << endl;
    cout << "COMPLETED: startTemplateFinder(): " << get_end_time_ms(start_time) << "ms" << endl;

    return 0;
}


int processCommandLine(int argc, char* argv[])
{
    char *pch;
    string sName, sVal;
    string::size_type pos;

    if (argc < 2) {
        cout << "\nMissing configuration file name!\n";
        printUsageAndExit(argv[0]);
    }

    confDB.initializeDB(argv[1]);

    // Now, any additional command line options are designed to override options
    // are in the configuration file
    for(int i = 2; i < argc; i++)
    {
        pch = argv[i];
        if (*pch++ != '-')
            printUsageAndExit(argv[0]);
        if (*pch == '\0')
            printUsageAndExit(argv[0]);

        switch(*pch++)
        {
        case '?':
            printUsageAndExit(argv[0]);
            break;

        case 'O':
            sName = pch;
            pos = sName.find('=');
            if (pos == string::npos || pos + 1 == sName.length())
                break;
            sVal = sName.substr(pos+1);
            sName.erase(pos);
            confDB.setKey(sName,sVal);
            cout << "addConfigOption(" << sName << "," << sVal << ")\n";
            break;

        case 'j':
            i++;
            if (i >= argc)
                printUsageAndExit(argv[0]);
            sJobId = argv[i];
            break;

        case 'n':
            if (i >= argc)
                printUsageAndExit(argv[0]);
            NUM_THREADS = atoi(argv[++i]);
            cout << "Spliting over '" << NUM_THREADS << "' files" << endl;
            //TODO split over argv[i+1] files
            break;

        default: 
            cerr << "Unknown Option: " << (pch - 1) << endl;
            printUsageAndExit(argv[0]);
            break;
        }
    }

    return 0;
}

void printUsageAndExit(char *sName)
{
    cout << "\nUsage: " << sName << " configFile [options]" << endl
            << " -?                    Command line options\n"
            << " configFile            Program configuration file name\n"
            << " -j jobIdString        Identifier to be used for this job\n"
            << " -n number of threads  Default = 1, (Reccomended: " << std::thread::hardware_concurrency() << ")\n"
            << " -OconfOption=Value    Override options in the configuration file here\n\n";

    exit(1);
}

/**
 * Program fails out somewhere in this function
 */
void readInReferenceGenome(){
    cout << "\nReading in reference genome..." << endl;
    fLogFileOut << "\nReading in reference genome..." << endl;
    ifstream input;
    t_chromosome chromosome;
    bool first = true;
    int iChr = confDB.getKey("chromosome").intVal;
    int currChr = 1;

    input.open(confDB.getKey("referenceFile").stringVal.c_str());

    // check if file is open
    if (input.is_open() == false){
        cout << "\nThe reference file could not be found: " << confDB.getKey("referenceFile").stringVal << endl;
        exit(1);
    }

    for (string row; getline(input, row, '\n');){
        if (row[0] == '>'){
            if (first){
                chromosome.fastaHeader = row.substr(1);
                first = false;
            } else {
                vReferenceGenome.push_back(chromosome);
                chromosome.fastaHeader = row.substr(1);
                chromosome.sequence = "";
                ++currChr;
            }
        } else if (row[0] == ' ')
            continue;
        else {
            if (currChr == iChr || iChr == 0)
                chromosome.sequence += row;
        }
    }
    vReferenceGenome.push_back(chromosome);

    /*ofstream output;
    output.open("chromosome.fasta");
    for (int i = 0; i < vReferenceGenome.size(); ++i){
        if (vReferenceGenome[i].sequence.length() < 2)
            continue;
        output << ">" << vReferenceGenome[i].fastaHeader << endl;
        int lngth = vReferenceGenome[i].sequence.length();
        for (int j = 0; j < lngth; ++j){
            output << vReferenceGenome[i].sequence.substr(j, 80) << endl;
            j += 79;
        }
    }

    output.close();
    exit(1);*/
    // print reference genome information
    printReferenceGenomeInfo(fLogFileOut);
}


/*
 * void setupBaseFileName()
 *
 * Sets up a base filename to use throughout the program for all log files,
 * etc... The file name will be based on the current time/date stamp.
 */
string setupBaseFileName()
{
    string sBaseFileName;
    char sTime[50];
    time_t t;
    tm *tmptr;

    t = time(NULL);
    tmptr = localtime(&t);
    strftime(sTime,50,"%Y%m%d_%H%M",tmptr);

    sBaseFileName = sTime;
    if (sJobId == "")
        sJobId = "0";

    sBaseFileName += "_" + sJobId + "_";

    cout << "BASE FILE NAME: " << sBaseFileName << endl;
    return sBaseFileName;
}

void prepareLogFile(ofstream& fout, string sBaseFileName, string sLogDirName)
{
#define DIRNAME_LEN 1024

    char curwd[DIRNAME_LEN];
    // Get the current working directory
    if (getcwd(curwd, DIRNAME_LEN) == NULL)
        perror("getcwd error");

    string sFile = sLogDirName + "/" + sBaseFileName + "Log.txt";

    // Turn on exceptions for file output stream
    fout.exceptions(ios::failbit|ios::badbit);
    try {
        fout.open(sFile.c_str());
    }
    catch (ios_base::failure& e) {
        throw ios_base::failure(string("Failed to open ") + sFile);
    }

    time_t t;
    t = time(NULL);
    fout << "Log File Started: " << ctime(&t) << "\r\n";
    fout << "Base File Name: " << sBaseFileName << "\r\n";
}

