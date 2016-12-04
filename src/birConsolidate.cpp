/*
 * birConsolidate.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: msegar
 */


#include "defs.h"
#include "CSVRow.hpp"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/wait.h>
#include <stdlib.h>
#include <algorithm>
#include <limits.h>
#include <ctime>
#include <cmath>
#include <cstdint>
#include <map>

using namespace std;

std::string remove_char(std::string s, char bad_c) {
    std::string ret_s;
    for (const char& c : s) {
        if (c != bad_c) {
            ret_s += c;
        }
    }
    return ret_s;
}

template <typename T>
void print_vector(std::vector<T> v) {
    std::cout << "[";
    for (size_t i = 0; i < v.size() - 1; ++i) { 
        std::cout << v.at(i) << ", ";
    }
    std::cout << v.back() << "]" << std::endl;
}

void print_t_consolidated(t_consolidated f) {
    std::cout << "----------------------" << std::endl;
    std::cout << f.sReadName << std::endl;
    std::cout << f.sParentRead << std::endl;
    std::cout << f.iParentStart << std::endl;
    std::cout << f.iParentEnd << std::endl;
    std::cout << f.iBirStart << std::endl;
    std::cout << f.iBirEnd << std::endl;
    std::cout << f.sBir << std::endl;
    std::cout << f.iBirLength << std::endl;
    std::cout << f.iTemplateStart << std::endl;
    std::cout << f.iTemplateEnd << std::endl;
    std::cout << f.sTemplate << std::endl;
    std::cout << f.iTemplateLength << std::endl;
    std::cout << f.bBirCandidateFound << std::endl;
    std::cout << f.iChromosome << std::endl;
    std::cout << f.bAnchorLeft << std::endl;
    std::cout << f.iFlag << std::endl;
    std::cout << f.bBadRead << std::endl; 
    std::cout << "----------------------" << std::endl;
}

// birConsolidate
int64_t startConsolidate(){
    cout << "\nstartConsolidate start..." << endl;
    fLogFileOut << "\nstartConsolidate start..." << endl;
    ofstream output;

    // only perform clustering if user option is specified
    if (confDB.getKey("performClustering").boolVal == true){
        cout << "start candidate read sorting..." << endl;
        fLogFileOut << "start candidate read sorting..." << endl;
        // sort the parent starting locations for each chromosome from smallest to largest
        //sort(vCandidateReads.begin(), vCandidateReads.end(), compareStart);

        //fLogFileOut << "Before consolidation: " << vCandidateReads.size() << endl;
        fLogFileOut << "Before consolidation: " << num_candidate_reads << endl;
        cout << "Printing unconsolidated_bir_locations.txt..." << endl;

        /*
        output.open((sProjectDirectory + "unconsolidated_bir_locations.txt").c_str());
        for (uint64_t i = 0; i < vCandidateReads.size(); ++i){
            if (vCandidateReads[i].bBadRead)
                continue;
            output << vCandidateReads[i].iChromosome << ", " << vCandidateReads[i].iParentStart << "," << vCandidateReads[i].iParentEnd << endl;
        }
        output.close();
        */


        // Finally, let's consolidated the possible BIR locations to speed up the alignment steps. This is based on overlapping regions
        consolidateLocations();

        fLogFileOut << "\nAfter consolidation: " << vConsolidated.size() << endl;

        if (confDB.getKey("mysql").boolVal){
            vCandidateReads.clear();
            vConsolidated.clear();
            std::cout << "EXITING SUCCESS" << std::endl;
            exit(0);
        }
    }

    // Let's destruct all the elements in the database to save space
    vCandidateReads.clear();
    //vConsolidated.clear();

    /*
     * *************************************************************************************************************************
     * ******************            THIS IS WHERE THE PROGRAM SPLITS IN TWO. PREVIOUS TO THIS POINT,          ******************
     * ******************           THE PROGRAM WORKS AS SPECIFIED. THE PROGRAM NOW REQUIRES ANOTHER          ******************
     * ******************           SCRIPT THAT UTILIZES MYSQL TO SEARCH FOR THE HALF-READS QUICKLY.          ******************
     * ******************          THE PROGRAM NOW WILL IMPORT THE RESULTS FROM THE SCRIPT AND CONTINUE.      ******************
     * *************************************************************************************************************************
     */

    // Find the other size of the half reads
    if (confDB.getKey("mysql").boolVal)
        createParentReadsFromMySQL();
    else
        createParentReads();

    // Make all the reads in a cluster the same length
    getConsensus();

    // Create a consensus read that takes int64_to consideration all the base calls at a specific location
    consolidateBaseCalls();

    // Let's destruct all the elements in the database to save space
    vConsolidated.clear();

    fLogFileOut << "\nConsolidation time = " << (int64_t)time(NULL)-time0 << endl;

    output.open((sProjectDirectory+ "consolidated_bir_locations.txt").c_str());
    for (uint64_t i = 0; i < vCandidateRegions.size(); ++i)
        output << vCandidateRegions[i].iChromosome << ", " << vCandidateRegions[i].iParentStart << "," << vCandidateRegions[i].iParentStart + vCandidateRegions[i].sParentRead.length() << endl;
    output.close();

    return 0;
}


/*
 * ************************************************************
 * consolidateLocations
 *
 * This function takes all the candidate regions and groups them together by the starting location of the unaligned read. If
 * there is any overlap between the reads, then they will be grouped together.
 *
 * A config.txt parameter (minConsolidate) determines how many reads need to be used to consolidate.
 */
void consolidateLocations(){
    cout << "\nConsolidating read locations start..." << endl;
    fLogFileOut << "\nConsolidating read locations start..." << endl;

    // TODO Find better way to consolidate, maybe based on stdev and discard ones that aren't (ie through out outliers)
    int64_t iLastEndLoc = 0;
    uint64_t iMinConsolidate = confDB.getKey("minConsolidate").intVal; // the minimum number of reads to consolidate, otherwise don't add to vector
    vector<t_consolidated> curr;
    int64_t iNumReads = 0;
    int64_t iSkippedCluster = 0;
    bool first = true;
    int64_t iChr = (confDB.getKey("chromosome").intVal == 0 ? 0 : confDB.getKey("chromosome").intVal - 1); // get the user specified chromoome
    bool bOnlyOneChr = (confDB.getKey("chromosome").intVal != 0 ? true :false);
    //int64_t sizeCan = vCandidateReads.size();
    int64_t sizeCan = num_candidate_reads;
    std::cout << "sizeCan: " << sizeCan << std::endl;
    int64_t fivePercent = sizeCan / 20;
    string sHalfReadClusters = confDB.getKey("clusterFile").stringVal;

    cout << "Chromsome: " << iChr << endl;

    std::string CandidateReads_filename;
    CandidateReads_filename = sProjectDirectory + "CandidateReads.csv";
    std::ifstream candidate_reads_file;
	candidate_reads_file.open(CandidateReads_filename); 

    CSVRow row;
    std::vector<std::string> row_v;
    t_consolidated frag, prev;

    try {
        for (size_t i = 0; candidate_reads_file >> row; ++i) {

            //row.print();

            //std::cout << "Itteration: " << i << std::endl;

            for (size_t j = 0; j < row.size(); ++j) {
                row_v.push_back(remove_char(row[j], ' '));
            }

            //print_vector(row_v);

            if (first && (i > 0)) {
                frag = prev;
            } else {
                frag.sReadName            = row_v[3];
                frag.sParentRead          = row_v[4];
                frag.iParentStart         = atoi(row_v[1].c_str());
                frag.iParentEnd           = atoi(row_v[2].c_str());
                frag.iBirStart            = atoi(row_v[5].c_str());
                frag.iBirEnd              = atoi(row_v[6].c_str());
                frag.sBir                 = row_v[7];
                frag.iBirLength           = atoi(row_v[8].c_str());
                frag.iTemplateStart       = atoi(row_v[9].c_str());
                frag.iTemplateEnd         = atoi(row_v[10].c_str());
                frag.sTemplate            = row_v[11];
                frag.iTemplateLength      = atoi(row_v[12].c_str());
                frag.bBirCandidateFound   = atoi(row_v[13].c_str()) != 0;
                frag.iChromosome          = atoi(row_v[0].c_str());
                frag.bAnchorLeft          = atoi(row_v[14].c_str()) != 0;
                frag.iFlag                = atoi(row_v[15].c_str());
                frag.bBadRead             = atoi(row_v[16].c_str()) != 0; 
            }
            row_v.clear();

            //print_t_consolidated(frag);

            if (i % fivePercent == 0)
                cout << "\tcandidateRead: " << i << " of " << sizeCan << " on chromosome " << iChr << " (" << ((i * 100) / sizeCan) << "%)" << endl;
            if (frag.iChromosome != iChr && frag.iChromosome != iChr+1)
                continue;
            else if (frag.iChromosome != iChr && frag.iChromosome == iChr+1){
                if (bOnlyOneChr) // if the user only wants to cluster one chromosome, we break out here
                    break;
                ++iChr;
                cout << "\nChromsome: " << iChr << endl;
                if (curr.size() >= iMinConsolidate && curr.size() <= 200){ // TODO make 200 a config variable
                    vConsolidated.push_back(curr);
                    iNumReads += curr.size();
                } else {
                    ++iSkippedCluster;
                }
                curr.clear();

                // now instead of copying the code below, we set first to true (since it's the start of a new chromosome) and repeat the loop
                first = true;
                --i;
                prev = frag;
                continue;
            }

            if (first){
                curr.push_back(prev);
                first = false;
                iLastEndLoc = frag.iParentStart + frag.sParentRead.length() - 1;
                continue;
            }

            if (frag.iParentStart > iLastEndLoc){
                if (curr.size() >= iMinConsolidate && curr.size() <= 200){ // TODO make 200 a config variable
                    vConsolidated.push_back(curr);
                    iNumReads += curr.size();
                } else {
                    ++iSkippedCluster;
                }
                curr.clear();
                curr.push_back(frag);
            } else {
                curr.push_back(frag);
            }
            iLastEndLoc = frag.iParentStart + frag.sParentRead.length() - 1;
        }
    } catch (exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    candidate_reads_file.close();

    fLogFileOut << "Number clusters skipped: " << iSkippedCluster << endl;
    cout << "\nNumber reads after consolidation: " << iNumReads << endl;
    fLogFileOut << "Number reads after consolidation: " << iNumReads << endl;
    cout << "Number of clusters: " << vConsolidated.size() << endl;
    fLogFileOut << "Number of clusters: " << vConsolidated.size() << endl;

    // Now the program will print64_t out the consolidate read locations
    // A MySQL and Perl script is required to parse the data to further the program
    // The program will now exit
    fLogFileOut << "\nPrinting half_read_clusters.txt" << endl;
    ofstream output;
    output.open((sProjectDirectory + sHalfReadClusters.c_str()).c_str());
    int64_t size2 = 0;
    for (int64_t i = 0; i < vConsolidated.size(); ++i){
        size2 = vConsolidated[i].size();
        for (int64_t j = 0; j < size2; ++j){
            output << vConsolidated[i].at(j).sReadName << "\t" << vConsolidated[i].at(j).sParentRead << "\t" << vConsolidated[i].at(j).iParentStart << "\t" << vConsolidated[i].at(j).iChromosome << "\t"<< i << endl;
        }
    }
    output.close();

        /*
    { // Clean up Candidate Reads File
        std::string command;
        command = "rm -f " + sProjectDirectory + "CandidateReads.csv";
        system(command.c_str());
    }
        */

    // print64_t out the half-read clustered locations
    /*fLogFileOut << "Printing half_read_clustered_locations.txt" << endl;
    output.open("half_read_clustered_locations.txt");
    output << "i\tchr\tsize\tmin\tmax\tlength" << endl;
    for (int64_t i = 0; i < iNumReads; ++i){
        int64_t min = 0;
        int64_t max = 0;
        size2 = vConsolidated[i].size();
        for (int64_t j = 0; j < size2; ++j){
            if (j == 0)
                min = vConsolidated[i].at(j).iParentStart;
            else if (vConsolidated[i].at(j).iParentStart < min)
                min = vConsolidated[i].at(j).iParentStart;
            if (max < vConsolidated[i].at(j).iParentEnd)
                max = vConsolidated[i].at(j).iParentEnd;
        }
        output << i << "\t" << vConsolidated[i].at(0).iChromosome << "\t" << size2 << "\t" << min << "\t" << max << "\t" << (max - min + 1) << endl;
    }
    output.close();*/

    fLogFileOut << "\nCluster time = " << (int64_t)time(NULL)-time0 << endl;

    // clean up memory before exiting
    //cout << "\nCleaning vCandidateReads..." << endl;
    //vCandidateReads.clear();
    //cout << "Cleaning vConsolidated..." << endl;
    //vConsolidated.clear();

    //exit(0);
}

int64_t createParentReadsFromMySQL(){
    cout << "\nImporting clusters start... " << endl;
    fLogFileOut << "\nImporting clusters start... " << endl;

    ifstream input;
    input.open((sProjectDirectory + confDB.getKey("mysqlFile").stringVal).c_str());
    // check if file is open
    if (input.is_open() == false){
        cout << "\nThe mysql_results file could not be found" << endl;
        exit(1);
    }

    t_consolidated cons;
    vector<string> curr;
    vector<t_consolidated> vCluster;
    int64_t iCluster = 0;
    int64_t iCurrCluster = 0;
    int64_t iClustersImported = 0;
    int64_t iReadsImported = 0;
    int64_t iChromosome = confDB.getKey("chromosome").intVal;
    bool first = true;
    //uint64_t iMinConsolidate = confDB.getKey("minConsolidate").intVal; // the minimum number of reads to consolidate, otherwise don't add to vector
    int64_t iSkipped = 0;

    for (string row; getline(input, row, '\n');){
        istringstream ss(row);
        curr.clear();
        for(string word; getline(ss, word, '\t');)
            curr.push_back(word);

        cons.sParentRead = curr[0];
        cons.sReadName = curr[1].substr(0, curr[1].length()-2);
        cons.iParentStart = atoi(curr[3].c_str());
        cons.iChromosome = atoi(curr[4].c_str());
        iCurrCluster = atoi(curr[5].c_str());

        if (iChromosome != (cons.iChromosome+1) && iChromosome != 0)
            continue;

        // catch reads that happen to have a starting location less than 0
        if (cons.iParentStart <= 0){
            ++iSkipped;
            continue;
        }

        if (first){
            iCluster = iCurrCluster;
            first = false;
        }

        if (iCurrCluster == iCluster)
            vCluster.push_back(cons);
        else {
            //if (vCluster.size() >= iMinConsolidate){
                vConsolidated.push_back(vCluster);
                iReadsImported += vCluster.size();
                ++iClustersImported;
            //} else
            //    ++iSkipped;
            vCluster.clear();
            vCluster.push_back(cons);
            iCluster = iCurrCluster;
        }
    }
    if (vCluster.size() > 0){
        vConsolidated.push_back(vCluster);
        iReadsImported += vCluster.size();
        vCluster.clear();
        ++iClustersImported;
    }

    input.close();

    cout << "Reads skipped (start < 0): " << iSkipped << endl;
    fLogFileOut << "\nReads skipped (start < 0): " << iSkipped << endl;
    cout << "Clusters imported: " << iClustersImported << endl;
    fLogFileOut << "Clusters imported: " << iClustersImported << endl;
    cout << "Reads imported: " << iReadsImported << endl;
    fLogFileOut << "Reads imported: " << iReadsImported << endl;
    fLogFileOut << "Cluster import time = " << (int64_t)time(NULL)-time0 << endl;

    ofstream output;
    output.open((sProjectDirectory+ "clustered_locations.txt").c_str());

    output << "i\tj\tstart\tchr\tread" << endl;
    int64_t size = vConsolidated.size();
    int64_t size2 = 0;
    for (int64_t i = 0; i < size; ++i) {
        size2 = vConsolidated[i].size();
        for (int64_t j = 0; j < size2; ++j)
            output << i << ", " << j << ", " << vConsolidated[i].at(j).iParentStart << ", " << vConsolidated[i].at(j).iChromosome << ", " << vConsolidated[i].at(j).sParentRead << endl;
    }

    output.close();

    return 0;
}

int64_t createParentReads(){
    cout << "\nCreating parent reads start... " << endl;
    fLogFileOut << "\nCreating parent reads start... " << endl;

    string sUnalignedFile = confDB.getKey("unalignedFile").stringVal + "_2.sam"; // the name of the unaligned FASTA/SAM file
    map<string, string> mUnaligned;
    map<string, string>::iterator iter;

    ifstream input;
    input.open(sUnalignedFile.c_str());
    vector<string> curr; // temporary vector to hold SAM file line

    // read in unaligned file
    for (string row; getline(input, row, '\n');){
        if (row[0] == '@') // skip lines that start with '@'
            continue;
        curr.clear();
        istringstream ss(row);
        for (string word; getline(ss, word, '\t');)
            curr.push_back(word);

        mUnaligned[curr[0].erase(curr[0].size() - 2)] = curr[9]; // remove last 2 characters from read name
    }
    input.close();

    cout << "Number of unaligned reads in " << sUnalignedFile << ": " << mUnaligned.size() << endl;
    fLogFileOut << "Number of unaligned reads in " << sUnalignedFile << ": " << mUnaligned.size() << endl;

    string sReadName = "";
    string sTemp = "";

    // Now loop through vCanidateReads appending unaligned info to the parent read

    //vector<t_consolidated>& temp;

    int64_t iFound = 0;
    int64_t iSize = vConsolidated.size();
    int64_t iSize2 = 0;
    char last;
    for (int64_t i = 0; i < iSize; ++i){
        vector<t_consolidated>& temp = vConsolidated[i];
        iSize2 = temp.size();
        for (int64_t j = 0; j < iSize2; ++j){
            sReadName = temp[j].sReadName;
            //sReadName = vConsolidated[i].at(j).sReadName;
            last = sReadName[sReadName.size() - 1];
            sReadName.erase(sReadName.size() - 2); // remove last 2 characters

            iter = mUnaligned.find(sReadName);
            if (iter != mUnaligned.end()){
                if (last == '1'){
                    temp[j].sParentRead += iter->second;
                } else {
                    //cout << iter->second << " + " << temp[j].sParentRead << endl;
                    //sTemp = iter->second + temp[j].sParentRead;
                    temp[j].sParentRead = iter->second + temp[j].sParentRead;
                }
                ++iFound;
            }
        }
    }

    fLogFileOut << "Found parent reads: " << iFound << endl;

    return 0;
}

/*
 * ************************************************************
 * getConsensus
 *
 * This function takes the grouped (consolidated) reads and makes them all the same length by adding '-' to the beginning and end.
 *
 * For example, if the reads are:
 *
 *   AACGTCGA                                      --AACGTCGA---
 * CGAACGTC            will be transformed int64_to    CGAACGTC-----
 *     CGTCGACGT                                ----CGTCGACGT
 *
 * This helps for the the next step which is calculate the base calls at each location.
 */
void getConsensus(){
    cout << "\nMaking reads same length..." << endl;
    fLogFileOut << "\nMaking reads same length..." << endl;

    int64_t iLowestStart;
    int64_t iLongest;
    int64_t iCurrStart;
    int64_t diff;
    int64_t size = 0;
    string currString;
    int64_t sizeCons = vConsolidated.size();
    int64_t fivePercent = sizeCons / 20;

    // make all the reads in each consolidation have the same starting point
    for (uint64_t i = 0; i < sizeCons; ++i){
        if (i % fivePercent == 0)
            cout << "consolidated: " << (i+1) << " of " << sizeCons << " (" << ((i * 100) / sizeCons) << "%)" << endl;
        size = vConsolidated[i].size();
        iLowestStart = vConsolidated[i].at(0).iParentStart;
        for (int64_t j = 1; j < size; ++j){
            iCurrStart = vConsolidated[i].at(j).iParentStart;
            if (iLowestStart > iCurrStart)
                iLowestStart = iCurrStart;
        }

        // Now we have the lowest starting point, let's add '-' to the beginning so they all match up
        for (int64_t j = 0; j < size; ++j){
            iCurrStart =  vConsolidated[i].at(j).iParentStart;
            currString = vConsolidated[i].at(j).sParentRead;
            diff = iCurrStart - iLowestStart;
            for (int64_t k = 0; k < diff; ++k){
                currString = "-" + currString;
                ++iCurrStart;
            }
            vConsolidated[i].at(j).sParentRead = currString;
            currString = "";
        }

        // Let's repeat only adding '-' to the end
        iLongest = vConsolidated[i].at(0).sParentRead.length();
        for (int64_t j = 1; j < size; ++j){
            iCurrStart = vConsolidated[i].at(j).sParentRead.length();
            if (iLongest < iCurrStart)
                iLongest = iCurrStart;
        }

        for (int64_t j = 0; j < size; ++j){
            iCurrStart = vConsolidated[i].at(j).sParentRead.length();
            currString = vConsolidated[i].at(j).sParentRead;
            for (int64_t k = iCurrStart; k < iLongest; ++k){
                currString = currString + "-";
                ++iCurrStart;
            }
            vConsolidated[i].at(j).sParentRead = currString;
            currString = "";
        }
    }
    cout << "exiting..." << endl;
    ofstream output;
    output.open((sProjectDirectory + "clustered_locations2.txt").c_str());

    output << "i\tj\tstart\tchr\tread" << endl;
    int64_t size2 = vConsolidated.size();
    int64_t size3 = 0;
    for (int64_t i = 0; i < size2; ++i) {
        size3 = vConsolidated[i].size();
        for (int64_t j = 0; j < size3; ++j)
            output << i << ", " << j << ", " << vConsolidated[i].at(j).iParentStart << ", " << vConsolidated[i].at(j).iChromosome << ", " << vConsolidated[i].at(j).sParentRead << endl;
        output << endl;
    }

    output.close();
}


/*
 * ************************************************************
 * consolidateBaseCalls
 *
 * This is basically a modified de novo assembly method that takes the base call frequency at each location to create a new read.
 */
void consolidateBaseCalls(){
    cout << "\nCalculating consensus base calls..." << endl;
    fLogFileOut << "\nCalculating consensus base calls..." << endl;

    int64_t iConsolidatedSize = vConsolidated.size();
    int64_t iCurrSize = 0;
    int64_t arr[85];
    int64_t numTot = 0;
    int64_t iReadLength = 0;
    string blank = "";
    string sConsensus = ""; // the consensus string
    int64_t iMinStartPos = 0;
    string sCurrChar;
    int64_t iHighestProb;
    int64_t fivePercent = iConsolidatedSize / 20;
    t_consolidated consensusRead;

    /*
     * Memory is starting to become an issue. To compensate for this, we will start deleting the vConsolidated cluster after the vCandidateRegions
     * struct is added. However, erasing in a vector is inefficient since all the memory is copied. There is one way to compensate for this, start from
     * the back. Since it doesn't matter which cluster we start on, and reallocation won't be affected erasing from the back, let's reverse traverse the
     * vector.
     *
     * However, we need to start rearranging previous methods to start from the back as well. Start with the sorting algorithm and make it so the largest
     * chromosome is at the front of the vector.
     */
    // TODO erase starting from the back.

    for (int64_t i = 0; i < iConsolidatedSize; ++i){
        if (i % fivePercent == 0)
            cout << "consolidated: " << i << " of " << iConsolidatedSize << " (" << ((i * 100) / iConsolidatedSize) << "%)" << endl;

        // get the address to make lookups faster
        vector<t_consolidated>& curr = vConsolidated[i];
        iReadLength = curr[0].sParentRead.length();
        sConsensus = "";

        iCurrSize = curr.size();

        iMinStartPos = curr[0].iParentStart;
        for (int64_t j = 1; j < iCurrSize; ++j){
            if (curr[j].iParentStart < iMinStartPos)
                iMinStartPos = curr[j].iParentStart;
        }

        for (int64_t j = 0; j < iReadLength; ++j){
            // initialize array and variables to 0s
            for (int64_t k = 0; k < 85; ++k)
                arr[k] = 0;
            numTot = 0;

            for (int64_t k = 0; k < iCurrSize; ++k){
                ++arr[(int64_t) curr[k].sParentRead.at(j)];
                if (curr[k].sParentRead.at(j) != '-')
                    ++numTot;
            }
            //curr.clear();
            iHighestProb = 0;

            // -s
            if (numTot == 0)
                continue;
            // As
            if ((arr[65] * 10) / numTot > iHighestProb){
                iHighestProb = (arr[65] * 10) / numTot;
                sCurrChar = "A";
            }
            // Cs
            if ((arr[67] * 10) / numTot > iHighestProb){
                iHighestProb = (arr[67] * 10) / numTot;
                sCurrChar = "C";
            }
            // Gs
            if ((arr[71] * 10) / numTot > iHighestProb){
                iHighestProb = (arr[71] * 10) / numTot;
                sCurrChar = "G";
            }
            // Ts
            if ((arr[84] * 10) / numTot > iHighestProb){
                iHighestProb = (arr[84] * 10) / numTot;
                sCurrChar = "T";
            }

            sConsensus += sCurrChar;
            sCurrChar = "-";
        }

        curr.clear();
        if (sConsensus.length() > 900) // TODO THIS NEEDS TO BE FIXED, BUT FOR NOW WE'LL SEE
            continue;
        consensusRead.sParentRead = sConsensus;
        consensusRead.iParentStart = iMinStartPos;
        consensusRead.iChromosome = curr[0].iChromosome;
        consensusRead.bBirCandidateFound = false;
        vCandidateRegions.push_back(consensusRead);
    }

    // sort the consensus reads
    sort(vCandidateRegions.begin(), vCandidateRegions.end(), compareStart);
}


bool compareStart(const t_consolidated &a, const t_consolidated &b){
    if (a.iChromosome < b.iChromosome) return true;
    if (a.iChromosome == b.iChromosome) return a.iParentStart < b.iParentStart;
    return false;
}
