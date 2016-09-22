#!/usr/bin/python3

import sys
import os

def get_extention(filename):
    for i in range(len(filename)):
        if filename[i] == '.':
            index = i
    return filename[index+1 :]

def remove_extension(filename):
    for i in range(len(filename)):
        if filename[i] == '.':
            index = i
    filename = filename[ : index]
    return filename 

def split_fastaq(input_filename, output_filename_list):

    print("Calculating size of input file...")
    num_lines = sum(1 for line in open(input_filename))
    print("Input file has " + str(num_lines) + " lines")

    num_files = len(output_filename_list)
    file_index = 0
    num_lines_curr_file = 0

    input_file = open(input_filename, "r")
    output_file = open(output_filename_list[file_index], "w")

    apprx_lines_per_file = num_lines // num_files
    print("Writing approxamately " + str(apprx_lines_per_file) + " lines per each output file")

    for line in input_file: 

        # split text to next file
        if num_lines_curr_file > apprx_lines_per_file and line[0] =='@':
            file_index += 1
            output_file.close()
            output_file = open(output_filename_list[file_index], "w")
            num_lines_curr_file = 0

        output_file.write(line) 
        num_lines_curr_file += 1

    return 0

def split_sam(input_filename, output_filename_list):
    return 0

number_of_output_files = 0
input_filename = ""

i = 1
while i < len(sys.argv):

    if sys.argv[i] == '-n':
        number_of_output_files = int(sys.argv[i+1])
        i += 1

    elif sys.argv[i] == '-i':
        input_filename = sys.argv[i+1]
        i += 1

    else:
        print("Error: invalid argument '" + sys.argv[i] + "'")
        quit()

    i += 1

if not number_of_output_files:
    print("Error: specify a non-zero number of output files")
    quit()

if number_of_output_files < 0:
    print("Error: specify a positive number of ouput files")
    quit()

if not input_filename:
    print("Error: specify an input file")
    quit()

print(input_filename, number_of_output_files)
print(get_extention(input_filename))
print(remove_extension(input_filename))


# fill output_filename_list with names of output files
# The number of output files will be number_of_output_files 
output_filename_list = []
for i in range(number_of_output_files):
    output_filename_list.append(remove_extension(input_filename) + ".part_" + str(i) + "." + get_extention(input_filename))

print(output_filename_list)

if get_extention(input_filename) == "fastq":
    split_fastaq(input_filename, output_filename_list)
else:
    print("Error: filetype '" + get_extention(input_filename) + "' not supported")
