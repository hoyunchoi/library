#pragma once
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

//* Write and Read CSV file
namespace CSV {
//* Write
template <typename T, typename TT>
void write(const std::string& t_writeFileName, const std::map<T, TT>& t_data, const char t_seperate = ',', const int& precision = 0) {
    std::ofstream writeFile(t_writeFileName);
    if (precision) {
        writeFile.precision(precision);
    }
    for (auto& row : t_data) {
        writeFile << row.first << t_seperate << row.second << "\n";
    }
}

template <typename T, typename TT, typename TTT>
void write(const std::string& t_writeFileName, const std::map<std::pair<T, TT>, TTT>& t_data, const char t_seperate = ',', const int& precision = 0) {
    std::ofstream writeFile(t_writeFileName);
    if (precision) {
        writeFile.precision(precision);
    }
    for (auto& row : t_data) {
        writeFile << row.first.first << t_seperate << row.first.second << t_seperate << row.second << "\n";
    }
}

template <typename T>
void write(const std::string& t_writeFileName, const std::vector<T>& t_data, const int& precision = 0) {
    std::ofstream writeFile(t_writeFileName);
    if (precision) {
        writeFile.precision(precision);
    }
    std::copy(t_data.begin(), t_data.end(), std::ostream_iterator<T>(writeFile, "\n"));
}

template <typename T>
void write(const std::string t_writeFileName, const std::vector<std::vector<T>>& t_data, const char t_seperate = ',', const int& precision = 0) {
    std::ofstream writeFile(t_writeFileName);
    if (precision) {
        writeFile.precision(precision);
    }
    for (const std::vector<T>& row : t_data) {
        for (unsigned i = 0; i < row.size() - 1; ++i) {
            writeFile << row[i] << t_seperate;
        }
        writeFile << row.back() << "\n";
    }
}

//* Read
void read(const std::string& readFileName, std::vector<double>& t_data) {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row;
    while (getline(readFile, row)) {
        t_data.emplace_back(std::stod(row));
    }
}
void read(const std::string& readFileName, std::vector<int>& t_data) {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row;
    while (getline(readFile, row)) {
        t_data.emplace_back(std::stoi(row));
    }
}
void read(const std::string& readFileName, std::vector<unsigned>& t_data) {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row;
    while (getline(readFile, row)) {
        t_data.emplace_back(std::stoul(row));
    }
}

void read(const std::string& readFileName, std::vector<std::vector<double>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row, column;
    std::vector<double> rowvector;
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        rowvector.clear();
        while (getline(rowstream, column, t_seperate)) {
            rowvector.emplace_back(std::stod(column));
        }
        t_data.emplace_back(rowvector);
    }
}

void read(const std::string& readFileName, std::vector<std::vector<int>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row, column;
    std::vector<int> rowvector;
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        rowvector.clear();
        while (getline(rowstream, column, t_seperate)) {
            rowvector.emplace_back(std::stoi(column));
        }
        t_data.emplace_back(rowvector);
    }
}

void read(const std::string& readFileName, std::vector<std::vector<unsigned>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row, column;
    std::vector<unsigned> rowvector;
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        rowvector.clear();
        while (getline(rowstream, column, t_seperate)) {
            rowvector.emplace_back(std::stoul(column));
        }
        t_data.emplace_back(rowvector);
    }
}

void read(const std::string& readFileName, std::map<double, double>& t_data, const char t_seperate = ',') {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row, key, value;
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        int i = 0;
        while (getline(rowstream, value, t_seperate)) {
            if (i % 2 == 1) {
                t_data[std::stod(key)] = std::stod(value);
            }
            key = value;
            ++i;
        }
    }
}

void read(const std::string& readFileName, std::map<int, double>& t_data, const char t_seperate = ',') {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file : " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row, key, value;
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        int i = 0;
        while (getline(rowstream, value, t_seperate)) {
            if (i % 2 == 1) {
                t_data[std::stoi(key)] = std::stod(value);
            } else {
                key = value;
            }
            ++i;
        }
    }
}

void read(const std::string& readFileName, std::map<std::pair<int, int>, double>& t_data, const char t_seperate = ',') {
    std::ifstream readFile(readFileName);
    if (readFile.fail()) {
        std::cout << "No such file: " << readFileName << std::endl;
        exit(1);
    }
    t_data.clear();
    std::string row, firstKey, secondKey, value;
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        int i = 0;
        while (getline(rowstream, value, t_seperate)) {
            if (i % 3 == 0) {
                firstKey = value;
            } else if (i % 3 == 1) {
                secondKey = value;
            } else {
                t_data[std::make_pair(std::stoi(firstKey), std::stoi(secondKey))] = std::stod(value);
            }
            ++i;
        }
    }
}

//* Delete File
void deleteFile(const std::string& t_deleteFileName) {
    if (std::remove(t_deleteFileName.c_str())) {
        std::cout << "Error in removing " << t_deleteFileName << "\n";
    }
}
}  // namespace CSV