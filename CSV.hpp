#pragma once
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
//* Write and Read CSV file
namespace CSV {
//* Write
template <typename T, typename TT>
void write(const std::string& t_writeFile, const std::map<T, TT>& t_data, const int& t_precision = -1, const char t_seperate = ',') {
    std::ofstream writeFile(t_writeFile);
    t_precision < 0 ? writeFile.precision(std::numeric_limits<double>::digits10 + 1) : writeFile.precision(t_precision);
    for (const auto& row : t_data) {
        writeFile << row.first << t_seperate << row.second << "\n";
    }
}

template <typename T, typename TT, typename TTT>
void write(const std::string& t_writeFile, const std::map<std::pair<T, TT>, TTT>& t_data, const int& t_precision = -1, const char t_seperate = ',') {
    std::ofstream writeFile(t_writeFile);
    t_precision < 0 ? writeFile.precision(std::numeric_limits<double>::digits10 + 1) : writeFile.precision(t_precision);
    for (const auto& row : t_data) {
        writeFile << row.first.first << t_seperate << row.first.second << t_seperate << row.second << "\n";
    }
}

template <typename T, typename TT>
void write(const std::string& t_writeFile, const std::map<T, std::vector<TT>>& t_data, const int& t_precision = -1, const char t_seperate = ',') {
    std::ofstream writeFile(t_writeFile);
    t_precision < 0 ? writeFile.precision(std::numeric_limits<double>::digits10 + 1) : writeFile.precision(t_precision);
    for (const auto& row : t_data) {
        writeFile << row.first << t_seperate;
        for (unsigned i = 0; i < row.second.size() - 1; ++i) {
            writeFile << row.second[i] << t_seperate;
        }
        writeFile << row.second.back() << "\n";
    }
}

template <typename T>
void write(const std::string& t_writeFile, const std::vector<T>& t_data, const int& t_precision = -1) {
    std::ofstream writeFile(t_writeFile);
    t_precision < 0 ? writeFile.precision(std::numeric_limits<double>::digits10 + 1) : writeFile.precision(t_precision);
    std::copy(t_data.begin(), t_data.end(), std::ostream_iterator<T>(writeFile, "\n"));
}

template <typename T>
void write(const std::string t_writeFile, const std::vector<std::vector<T>>& t_data, const int& t_precision = -1, const char t_seperate = ',') {
    std::ofstream writeFile(t_writeFile);
    t_precision < 0 ? writeFile.precision(std::numeric_limits<double>::digits10 + 1) : writeFile.precision(t_precision);
    for (const std::vector<T>& row : t_data) {
        for (unsigned i = 0; i < row.size() - 1; ++i) {
            writeFile << row[i] << t_seperate;
        }
        writeFile << row.back() << "\n";s
    }
}

//* Read
std::ifstream check(const std::string& t_readFile) {
    std::ifstream readFile(t_readFile);
    if (readFile.fail()) {
        std::cout << "No such file : " << t_readFile << std::endl;
        exit(1);
    }
    return readFile;
}

void read(const std::string& t_readFile, std::vector<double>& t_data) {
    std::ifstream readFile = check(t_readFile);
    std::string row;
    t_data.clear();
    while (getline(readFile, row)) {
        t_data.emplace_back(std::stod(row));
    }
}
void read(const std::string& t_readFile, std::vector<int>& t_data) {
    std::ifstream readFile = check(t_readFile);
    std::string row;
    t_data.clear();
    while (getline(readFile, row)) {
        t_data.emplace_back(std::stoi(row));
    }
}
void read(const std::string& t_readFile, std::vector<unsigned>& t_data) {
    std::ifstream readFile = check(t_readFile);
    std::string row;
    t_data.clear();
    while (getline(readFile, row)) {
        t_data.emplace_back(std::stoul(row));
    }
}

void read(const std::string& t_readFile, std::vector<std::vector<double>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, column;
    t_data.clear();
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        std::vector<double> rowvector;
        while (getline(rowstream, column, t_seperate)) {
            rowvector.emplace_back(std::stod(column));
        }
        t_data.emplace_back(rowvector);
    }
}
void read(const std::string& t_readFile, std::vector<std::vector<int>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, column;
    t_data.clear();
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        std::vector<int> rowvector;
        while (getline(rowstream, column, t_seperate)) {
            rowvector.emplace_back(std::stoi(column));
        }
        t_data.emplace_back(rowvector);
    }
}
void read(const std::string& t_readFile, std::vector<std::vector<unsigned>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, column;
    t_data.clear();
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        std::vector<unsigned> rowvector;
        while (getline(rowstream, column, t_seperate)) {
            rowvector.emplace_back(std::stoul(column));
        }
        t_data.emplace_back(rowvector);
    }
}

void read(const std::string& t_readFile, std::map<double, double>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, key, value;
    t_data.clear();
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

void read(const std::string& t_readFile, std::map<int, double>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, key, value;
    t_data.clear();
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

void read(const std::string& t_readFile, std::map<std::pair<int, int>, double>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, firstKey, secondKey, value;
    t_data.clear();
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
    return;
}

void read(const std::string& t_readFile, std::map<int, std::vector<double>>& t_data, const char t_seperate = ',') {
    std::ifstream readFile = check(t_readFile);
    std::string row, key, value;
    t_data.clear();
    while (getline(readFile, row)) {
        std::stringstream rowstream(row);
        std::vector<double> valueVec;
        int i = 0;
        while (getline(rowstream, value, t_seperate)) {
            if (i == 0) {
                key = value;
            } else {
                valueVec.emplace_back(std::stod(value));
            }
            ++i;
        }
        t_data[std::stoi(key)] = valueVec;
    }
    return;
}

//* Define new directory
void generateDirectory(const std::string& t_directory) {
    if (!std::filesystem::exists(t_directory)) {
        std::filesystem::create_directories(t_directory);
    }
    return;
}

//* Delete File
void deleteFile(const std::string& t_deleteFile) {
    if (std::remove(t_deleteFile.c_str())) {
        std::cout << "Error in removing " << t_deleteFile << "\n";
    }
    return;
}
}  // namespace CSV