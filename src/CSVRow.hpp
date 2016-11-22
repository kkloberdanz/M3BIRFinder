/**
 * Credit goes to Loki Astari on stackoverflow
 */
#pragma once

#ifndef CSV_ROW_HPP
#define CSV_ROW_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>


class CSVRow {
    private:
        std::vector<std::string> row_v;

    public:
        std::string const& operator[](std::size_t index) const {
            return row_v.at(index);
        }

        std::size_t size() const {
            return row_v.size();
        }

        void readNextRow(std::istream& str) {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            row_v.clear();
            while(std::getline(lineStream, cell, ',')) {
                row_v.push_back(cell);
            }
        }

        void print() {
            std::cout << "[";
            for (size_t i = 0; i < row_v.size() - 1; ++i) {
                std::cout << row_v.at(i) << ", ";
            }
            std::cout << row_v.back() << "]" << std::endl;
        }

        /*
        void remove_spaces() { 
            for (size_t j = 0; j < row_v.size(); ++j) {
                std::remove_if(row_v.at(j).begin(), row_v.at(j).end(), isspace);
            }
        } 
        */
};

std::istream& operator>>(std::istream& str, CSVRow& data) {
    data.readNextRow(str);
    return str;
}

#endif // CSV_ROW_HPP 
