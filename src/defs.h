/*
 * defs.h
 *
 *  Created on: Dec 12, 2012
 *      Author: msegar
 */

#ifndef DEFS_H_
#define DEFS_H_

#include "ConfigDB.h"
#include "Fragment.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <cstdint>

using namespace std;

struct t_chromosome {
    string fastaHeader;
    string samHeader;
    string sequence;
};

struct t_consolidated{
    string sReadName;
    string sParentRead;
    int64_t iParentStart;
    int64_t iParentEnd;
    int64_t iBirStart;
    int64_t iBirEnd;
    string sBir;
    int64_t iBirLength;
    int64_t iTemplateStart;
    int64_t iTemplateEnd;
    string sTemplate;
    int64_t iTemplateLength;
    bool bBirCandidateFound;
    int64_t iChromosome;
    bool bAnchorLeft;
    int64_t iFlag;
    bool bBadRead;
};


struct t_alignment_struct{
    int64_t iEndPosI;
    int64_t iStartPosI;
    int64_t iEndPosJ;
    int64_t iStartPosJ;
    string sAlignedRegionI;
    string sAlignedRegionJ;
};

// globals
#ifdef MAIN_CPP_
ConfigDB confDB;
vector<t_consolidated> vCandidateReads; // map the read name to the children fragments
vector<vector <t_consolidated> > vConsolidated; // vector of CONSOLIDATED BIR locations
vector<t_consolidated> vCandidateRegions;
vector<t_consolidated> vFinalBirLocs;
vector<t_chromosome> vReferenceGenome; // the reference genome
int64_t time0;
string sJobId;
string sBaseFileName;
string sProjectDirectory;
ofstream fLogFileOut;
size_t num_candidate_reads;

#else
extern ConfigDB confDB;
extern vector<t_consolidated> vCandidateReads; // map the read name to the children fragments
extern vector<vector <t_consolidated> > vConsolidated; // vector of CONSOLIDATED BIR locations
extern vector<t_consolidated> vCandidateRegions;
extern vector<t_consolidated> vFinalBirLocs;
extern vector<t_chromosome> vReferenceGenome; // the reference genome
extern int64_t time0;
extern string sJobId;
extern string sBaseFileName;
extern string sProjectDirectory;
extern ofstream fLogFileOut;
extern size_t num_candidate_reads;

#endif /* MAIN_CPP_ */


// main
void readInReferenceGenome();
int64_t startExecutables();
int64_t processCommandLine(int argc, char* argv[]);
void printUsageAndExit(char *sName);
string setupBaseFileName();
void prepareLogFile(ofstream& fout, string sBaseFileName, string sLogDirName);

// programExecution
int64_t startExecutables();
int64_t executeBwaIndex(string s);
int64_t executeBwaAligner(string s1, string s2, string s3);
int64_t getReads(string s1, string s2, string s3, string s4, bool b1);
int64_t filterOut(string s1, string s2, string s3, string s4);
int64_t convertSAMtoFASTA(string s);

// birFinder
int64_t startCandidateReads();
int64_t createSplitReadDatabase(string s1);

// birAligner
int64_t startBirFinder();
int64_t getLastPosition(string &sReference, string &sAligned);
void getBirLoc(t_alignment_struct tAligned, int64_t st, int64_t i);

// birConsolidate
int64_t startConsolidate();
void consolidateLocations();
int64_t createParentReadsFromMySQL();
int64_t createParentReads();
void getConsensus();
void consolidateBaseCalls();
bool compareStart(const t_consolidated &a, const t_consolidated &b);

// templateFinder
int64_t startTemplateFinder();
string getReverseComplement(string &str);

// alignment
t_alignment_struct getLocalAlignment(string seq1, string seq2, double d1, double d2);
t_alignment_struct getGlobalAlignment(string &seq1, string &seq2, double d1, double d2);
double getSimilarityScore(char a, char b);
int64_t getMaxArrayValue(double array[], int64_t length);

// printOptions
void printConfig(ostream& out);
void printMatches(ostream& out, t_alignment_struct &a);
void printError(ostream& out, string &bir, string &ref, t_alignment_struct tAligned);
void printAlignment(ostream& out, string &bir, string &ref, int64_t readStart, t_alignment_struct tAligned);
void printFinal();
void printReferenceGenomeInfo(ostream& out);
void printCandidateReadsInfo(ostream& out);

#endif /* DEFS_H_ */
