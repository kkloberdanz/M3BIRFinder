/*
 * birFinder.cpp
 *
 *  Created on: Dec 12, 2012
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

#define WRITE_BUFF_SIZE 10000

using namespace std;

struct t_read{
	int iLength;
	int iStart;
	int iChromosome;
	string sSequence;
};

void sort_consolidated_file() {
    std::cout << "Sorting unconsolidated_bir_locations.txt ..." << std::endl;

    // Uses GNU sort with the project directory as the tmp directory
    std::string command = "sort -nk1,1 -nk2,2 -nk3,3 -T ./" + sProjectDirectory + " " + sProjectDirectory + "unconsolidated_bir_locations.txt > " + sProjectDirectory + "unconsolidated_bir_locations.txt.tmp";
    std::cout << command << std::endl;
    if (system(command.c_str()) != 0) {
        std::cout << "Error: sort not successful. Ensure there is enough space on disk to perform sort" << std::endl;
        std::cerr << "Error: sort not successful. Ensure there is enough space on disk to perform sort" << std::endl;
        exit(EXIT_FAILURE);
    }
    command = "mv " + sProjectDirectory + "unconsolidated_bir_locations.txt.tmp " + sProjectDirectory + "unconsolidated_bir_locations.txt";
    system(command.c_str());
    std::cout << command << std::endl;
    std::cout << "Done" << std::endl;

    //command = "sort -n -T ./" + sProjectDirectory + " " + sProjectDirectory + "CandidateReads.csv > " + sProjectDirectory + "CandidateReads.csv.tmp";
    command = "sort  -nk1,1 -nk2,2 -nk3,3 -T ./" + sProjectDirectory + " " + sProjectDirectory + "CandidateReads.csv > " + sProjectDirectory + "CandidateReads.csv.tmp";
    std::cout << command << std::endl;
    if (system(command.c_str()) != 0) {
        std::cout << "Error: sort not successful. Ensure there is enough space on disk to perform sort" << std::endl;
        std::cerr << "Error: sort not successful. Ensure there is enough space on disk to perform sort" << std::endl;
        exit(EXIT_FAILURE);
    }
    command = "mv " + sProjectDirectory + "CandidateReads.csv.tmp " + sProjectDirectory + "CandidateReads.csv";
    system(command.c_str());
    std::cout << command << std::endl;
    std::cout << "Done" << std::endl;
}

// birFinder
int startCandidateReads(){
	string sFinalAlignedFile = confDB.getKey("finalAlignedFile").stringVal;

	if (confDB.getKey("performClustering").boolVal == false)
		return 0;

	cout << "\nDB before: " << vCandidateReads.size() << endl;
	fLogFileOut << "\nDB before: " << vCandidateReads.size() << endl;

	// Let's take the final unaligned reads and find the locations of the matched pair reads
	createSplitReadDatabase(sFinalAlignedFile + ".sam");

	cout << "\nDB after: " << vCandidateReads.size() << endl;
	fLogFileOut << "\nDB after: " << vCandidateReads.size() << endl;

	//createParentReads();

	fLogFileOut << "\nRead database time = " << (int)time(NULL)-time0 << endl;

	return 0;
}

int createSplitReadDatabase(string sAlignedFilename){
	cout << "\nCreating candidate read database starting..." << endl;
	fLogFileOut << "\nCreating candidate read database starting..." << endl;
	ifstream input;
    ofstream unconsolidated_output, candidate_reads_output;
	vector<string> curr;
	char row_delim = '\n';
	char field_delim = '\t';
	string sReadName;
	t_consolidated frag;
	int iChr = 0;
	int chromosome = confDB.getKey("chromosome").intVal;
	int iExludedReads = 0;
	int iBadReads = 0;
	int index = 0;


    t_consolidated f;
    std::vector<t_consolidated> good_write_buffer;
    std::vector<t_consolidated> total_write_buffer;

    size_t good_index = 0;
    size_t total_index = 0;

	fLogFileOut << "Aligned file: " << sProjectDirectory + sAlignedFilename << endl;

	// Read in the aligned SAM file and create a database of anchored half reads and their locations
	// We are treating these as candidate reads
    std::cout << "Opening: " << sProjectDirectory + sAlignedFilename << std::endl;
	input.open((sProjectDirectory + sAlignedFilename).c_str());
    if (!input.is_open()) {
        std::cout << "Could not open file: '" << sProjectDirectory + sAlignedFilename << "'" << std::endl;
    }

    ofstream output;
    candidate_reads_output.open((sProjectDirectory + "CandidateReads.csv").c_str());
    unconsolidated_output.open((sProjectDirectory + "unconsolidated_bir_locations.txt").c_str());
    num_candidate_reads = 0;
	for (string row; getline(input, row, row_delim);){ 

        try {
            ++index;
            if (index % 10000000 == 0) {
                cout << "half-read: " << index << endl;
                /*
                std::cout << "Size of curr: " << curr.size() << std::endl;
                std::cout << "Size of vCandidateReads: " << vCandidateReads.size() << std::endl;
                std::cout << "capacity: " << vCandidateReads.capacity() << std::endl;
                std::cout << "max_size: " << vCandidateReads.max_size() << std::endl;
                */
            }
            // reset vector
            curr.clear();

            istringstream ss(row);
            for(string word; getline(ss, word, field_delim);)
                curr.push_back(word);


            if ((row.size() > 2) && row[0]=='@' && row[1]=='S' && row[2]=='Q'){
                // Get the header information
                vReferenceGenome.at(iChr).samHeader = curr.at(1).substr(3);
                ++iChr;
            } else if ((row.size() > 0) && (row[0]!='@')) {
                // We store the read name as the key and the position, length, and sequence as the elements in a struct
                sReadName = curr.at(0);
                frag.sReadName = curr.at(0);
                frag.iParentStart = atoi(curr.at(3).c_str());
                frag.sParentRead = curr.at(9);
                frag.iParentEnd = frag.iParentStart + frag.sParentRead.length() - 1;
                frag.iFlag = atoi(curr.at(1).c_str());

                if (sReadName[sReadName.length()-1] == '1'){ // the first half is anchored
                    frag.bAnchorLeft = true;
                } else if (sReadName[sReadName.length()-1] == '2'){ // the second half is anchored, first half is unaligned
                    frag.bAnchorLeft = false;
                } else {
                    fLogFileOut << "\n ** ERROR ** There was an error when searching for the other half of an anchored read" << endl;
                    continue;
                }

                for (unsigned int i = 0; i < vReferenceGenome.size(); ++i){
                    //cout << curr[2] << " ?= " << vReferenceGenome[i].fastaHeader << endl;
                    if (curr.at(2).find(vReferenceGenome.at(i).fastaHeader) != string::npos)
                        frag.iChromosome = i;
                    else if (vReferenceGenome.at(i).fastaHeader.find(curr.at(2)) != string::npos) // TODO find better way of searching for fasta header
                        frag.iChromosome = i;
                }

                if (frag.iParentStart < 0 || frag.iFlag == 4){
                    frag.bBadRead = true;
                    ++iBadReads;
                } else
                    frag.bBadRead = false;

                // limit to single chromosome
                if (chromosome != (frag.iChromosome+1) && chromosome > 0){
                    ++iExludedReads;
                    continue;
                }


                /*
                 * If read is good, store in write buffer.
                 * Write buffer to a file after some time
                 */
                //vCandidateReads.push_back(frag);
                num_candidate_reads++;
                if (not frag.bBadRead) {
                    //good_write_buffer.at(good_index) = frag;
                    //good_index++;
                    //good_write_buffer.push_back(frag);
                        unconsolidated_output << frag.iChromosome << ", " << frag.iParentStart << ", " << frag.iParentEnd << std::endl;
                }

                /*
                //if (good_index >= WRITE_BUFF_SIZE) {
                if (good_write_buffer.size() >= WRITE_BUFF_SIZE) {
                    //for (size_t i = 0; i < good_index; ++i) { 
                    for (const auto& f : good_write_buffer) { 
                        //f = good_write_buffer.at(i);
                        unconsolidated_output << f.iChromosome << ", " 
                               << f.iParentStart << ", " 
                               << f.iParentEnd  << ", "
                               << std::endl;
                    } 
                    //good_index = 0;
                    good_write_buffer.clear();
                }
                */


                //total_write_buffer.push_back(frag);
                //total_write_buffer.at(total_index) = frag; 
                //total_index++;
                //if (total_index >= WRITE_BUFF_SIZE) {
                //if (total_write_buffer.size() >= WRITE_BUFF_SIZE) {
                    //for (size_t i = 0; i < total_index; ++i) {
                    //for (const auto& f : total_write_buffer) {
                        //f = total_write_buffer.at(i);

                        /* Sort based on these */
                        candidate_reads_output 
                               << frag.iChromosome << ", "  // 0
                               << frag.iParentStart << ", " // 1
                               << frag.iParentEnd  << ", "  // 2

                        /* Store these */
                               << frag.sReadName << ", "    // 3
                               << frag.sParentRead << ", "  // 4
                               << frag.iBirStart << ", "    // 5

                               << frag.iBirEnd << ", "      // 6
                               << frag.sBir << ", "         // 7
                               << frag.iBirLength << ", "   // 8 

                               << frag.iTemplateStart << ", " // 9
                               << frag.iTemplateEnd << ", "   // 10
                               << frag.sTemplate << ", "      // 11

                               << frag.iTemplateLength << ", "    // 12
                               << frag.bBirCandidateFound << ", " // 13
                               << frag.bAnchorLeft << ", "        // 14

                               << frag.iFlag << ", "              // 15
                               << frag.bBadRead                   // 16

                               << endl;
                    //}
                    //total_index = 0;
                    //total_write_buffer.clear();
                //}

                }
        } catch (const std::bad_alloc& e) {
            std::cout << "Allocation failed: " << e.what() << std::endl;
            std::cout << "Size of vCandidateReads: " << vCandidateReads.size() << std::endl;
            std::cout << "capacity: " << vCandidateReads.capacity() << std::endl;
            std::cout << "max_size: " << vCandidateReads.max_size() << std::endl;
            std::cout << "sizeof t_consolidated: " << sizeof(t_consolidated) << std::endl;
            output.close();
            //std::cout << "Allocation failed: " << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }

	}


    /*
    //for (size_t i = 0; i < good_index; ++i) {
    for (const auto& f : good_write_buffer) {
        //f = good_write_buffer.at(i);
        unconsolidated_output << f.iChromosome << ", " 
               << f.iParentStart << ", " 
               << f.iParentEnd  << ", "
               << std::endl;
    } 
    good_write_buffer.clear();

    //for (size_t i = 0; i < total_index; ++i) {
    for (const auto & f : total_write_buffer) {
        //f = total_write_buffer.at(i);

        // Sort based on these
        candidate_reads_output << f.iChromosome << ", " 
               << f.iParentStart << ", " 
               << f.iParentEnd  << ", " 

        // Store these
               << f.sReadName << ", " 
               << f.sParentRead << ", " 
               << f.iBirStart << ", " 

               << f.iBirEnd << ", " 
               << f.sBir << ", " 
               << f.iBirLength << ", " 

               << f.iTemplateStart << ", " 
               << f.iTemplateEnd << ", " 
               << f.sTemplate << ", " 

               << f.iTemplateLength << ", " 
               << f.bBirCandidateFound << ", " 
               << f.bAnchorLeft << ", " 

               << f.iFlag << ", " 
               << f.bBadRead

               << endl;
    }
    total_write_buffer.clear();
    */

    sort_consolidated_file();

    output.close(); 
	input.close();

	fLogFileOut << "Excluded non-chromosome specific reads: " << iExludedReads << endl;
	fLogFileOut << "Bad reads: " << iBadReads << endl;

	cout << "Number candidate reads: " << vCandidateReads.size() << endl;
	fLogFileOut << "Number candidate reads: " << vCandidateReads.size() << endl;

	return 0;
}
// samtools view -S unaligned_1.sam | awk '{OFS="\t"; print ">"$1"-1\n"substr($10,1,length($10)/2)"\n>"$1"-2\n"substr($10,length($10)/2+1,length($10))}' - > filename.fasta

