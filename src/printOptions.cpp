/*
 * printOptions.cpp
 *
 *  Created on: Jan 30, 2013
 *      Author: msegar
 */


#include "defs.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <sys/wait.h>
#include <stdlib.h>

using namespace std;


/*
 * Prints out the options from the config.txt file
 */
void printConfig(ostream& out) {
    out << "jobId:\t" << sJobId << endl;

    out << "projectDirectory:\t" << confDB.getKey("projectDirectory").stringVal << endl;
    out << "referenceGenome:\t" << confDB.getKey("referenceFile").stringVal << endl;
    out << "chromosome:\t" << confDB.getKey("chromosome").stringVal << endl;
    out << "readsFile:\t" << confDB.getKey("readsFile").stringVal << endl;
    out << "alignedFile:\t" << confDB.getKey("alignedFile").stringVal << endl;
    out << "unalignedFile:\t" << confDB.getKey("unalignedFile").stringVal << endl;
    out << "outputFile:\t" << confDB.getKey("outputFile").stringVal << endl;
    out << "locationsFile:\t" << confDB.getKey("locationsFile").stringVal << endl;
    out << "pairedEnd:\t" << confDB.getKey("pairedEnd").stringVal << endl;
    out << endl;
    out << "runAlignment:\t" << confDB.getKey("runAlignment").stringVal << endl;
    out << "onlyAlign:\t" << confDB.getKey("onlyAlign").stringVal << endl;
    out << "indexGenome:\t" << confDB.getKey("indexGenome").stringVal << endl;
    out << "fullAlign:\t" << confDB.getKey("fullAlign").stringVal << endl;
    out << "bamFile:\t" << confDB.getKey("bamFile").stringVal << endl;
    out << "extractUnalignedReads:\t" << confDB.getKey("extractUnalignedReads").stringVal << endl;
    out << "halfAlign:\t" << confDB.getKey("halfAlign").stringVal << endl;
    out << "extractHalfReads:\t" << confDB.getKey("extractHalfReads").stringVal << endl;
    out << "filterOut:\t" << confDB.getKey("filterOut").stringVal << endl;
    out << endl;
    out << "performClustering:\t" << confDB.getKey("performClustering").stringVal << endl;
    out << "minConsolidate:\t" << confDB.getKey("minConsolidate").stringVal << endl;
    out << "mysql:\t" << confDB.getKey("mysql").stringVal << endl;
    out << endl;
    out << "searchLength:\t" << confDB.getKey("searchLength").stringVal << endl;
    out << "minSeqLength:\t" << confDB.getKey("minSeqLength").stringVal << endl;
    out << "minBirLength:\t" << confDB.getKey("minBirLength").stringVal << endl;
    out << "minAlignedLength:\t" << confDB.getKey("minAlignedLength").stringVal << endl;
    out << "missCount:\t" << confDB.getKey("missCount").stringVal << endl;
    out << "hitCount:\t" << confDB.getKey("hitCount").stringVal << endl;
    out << "tolerance:\t" << confDB.getKey("tolerance").stringVal << endl;
    out << endl;
}

void printMatches(ostream& out, t_alignment_struct &tAligned){
    int64_t size = tAligned.sAlignedRegionI.length();
    int64_t counter = 0;
    string astericks = "";
    int64_t iLongestTag = 7; // TODO maybe add a config option for longest tag?
    int64_t iMinTagLength = 3;

    // get the asterick's string
    for (int64_t i = 0; i < size; ++i){
        if (tAligned.sAlignedRegionI[i] == tAligned.sAlignedRegionJ[counter])
            astericks += "*";
        else
            astericks += " ";
        ++counter;
    }

    int64_t astericks_length = astericks.length();
    int64_t tagLength = 0;
    int64_t iTagStart = 0;

    // get the tag distance from end
    for (int64_t i = astericks_length - 1; i >= 0; --i){
        if (astericks[i] == '*')
            ++tagLength;
        else {
            if (tagLength >= iMinTagLength && tagLength <= iLongestTag && i >= astericks_length - iLongestTag - iMinTagLength){
                iTagStart = i;
                break;
            }
            tagLength = 0;
        }
    }
    out << endl;


    // print64_t sAlignedRegionI
    for (int64_t i = 0; i < size; ++i){
        if (i == iTagStart + 1 && iTagStart != 0)
            out << '-';
        out << tAligned.sAlignedRegionI[i];
    }
    out << endl;

    // print64_t sAlignedRegionJ
    for (int64_t i = 0; i < size; ++i){
        if (i == iTagStart + 1 && iTagStart != 0)
            out << '-';
        out << tAligned.sAlignedRegionJ[i];
    }
    out << endl;

    // print64_t astericks
    for (int64_t i = 0; i < size; ++i){
        if (i == iTagStart + 1 && iTagStart != 0)
            out << '-';
        out << astericks[i];
    }
    out << endl;
}

void printError(ostream& out, string &bir, string &ref, t_alignment_struct tAligned){
    out << "\n*** There was an error with the BIR alignment" << endl;
    out << "r: " << ref << endl;
    out << "b: " << bir << endl;
    printMatches(out, tAligned);
    out << "length: " << bir.length() << endl;
    out << "aligned st: " << tAligned.iStartPosJ << endl;
    out << "aligned en: " << tAligned.iEndPosJ << endl;
    out << "st percent: " << (tAligned.iStartPosJ * 100 / bir.length()) << "%" << endl;
    out << "en percent: " << (tAligned.iEndPosJ * 100 / bir.length()) << "%" << endl;
}

void printAlignment(ostream& out, string &bir, string &ref, int64_t readStart, t_alignment_struct tAligned){
    out << "readStart: " << readStart << endl;
    out << "ref: " << ref << endl;
    out << "bir: " << bir << endl;
    out << tAligned.iStartPosJ << ", " << tAligned.iEndPosJ << ", " << bir.length() << endl;
    printMatches(out, tAligned);
}

void printFinal(){
    ofstream output;
    output.open("final_bir_locations.txt");

    for (uint64_t i = 0; i < vFinalBirLocs.size(); ++i){
        output << "Num: " << i+1 << "/" << vFinalBirLocs.size() << endl;
        output << "BIR: " << vFinalBirLocs[i].sBir << endl;
        output << " bs: " << vFinalBirLocs[i].iBirStart << endl;
        output << " be: " << vFinalBirLocs[i].iBirEnd << endl;
        output << " bl: " << vFinalBirLocs[i].sBir.length() << endl;
        output << "TEM: " << vFinalBirLocs[i].sTemplate << endl;
        output << " ts: " << vFinalBirLocs[i].iTemplateStart << endl;
        output << " te: " << vFinalBirLocs[i].iTemplateEnd << endl;
        output << " tl: " << vFinalBirLocs[i].sTemplate.length() << endl;
        output << "CHR: " << vFinalBirLocs[i].iChromosome << endl;
        output << "INS: " << vFinalBirLocs[i].iBirStart - vFinalBirLocs[i].iTemplateEnd << endl;
        output << endl;
    }

    output.close();
}

void printReferenceGenomeInfo(ostream& out){
    out << "\n----- REFERENCE GENOME INFORMATION -----" << endl;
    out << "Number of chromosomes: " << vReferenceGenome.size() << endl;
    for (uint64_t i = 0; i < vReferenceGenome.size(); ++i)
        out << "Chr " << i+1 << ": name = " << vReferenceGenome[i].fastaHeader << ", length = " << vReferenceGenome[i].sequence.length() << endl;
}

void printCandidateReadsInfo(ostream& out){
    out << "\n----- CANDIDATE READS INFORMATION -----" << endl;
    out << "Size: " << vCandidateReads.size() << endl;
    int64_t count = 0;
    for (uint64_t i = 0; i <vCandidateReads.size(); ++i){
        if (vCandidateReads[i].bAnchorLeft)
            ++count;
    }
    out << "Anchored left: " << count << endl;
    out << "        right: " << vCandidateReads.size() - count << endl;
    for (uint64_t i = 0; i < vReferenceGenome.size(); ++i){
        count = 0;
        for (uint64_t j = 0; j < vCandidateReads.size(); ++j){
            if (vCandidateReads[j].iChromosome == i)
                ++count;
        }
        out << "Chr " << i+1 << ": " << count << endl;
    }
}
