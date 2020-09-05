#pragma once
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <iterator>
#include <iostream>

//Write CSV file
template<typename T, typename TT>
void writeCSV(const std::string& outfilename, const std::map<T,TT>& data, const int& precision=0){
    std::ofstream outfile(outfilename);
    if (precision){
        outfile.precision(precision);
    }
    for (auto& row : data){
        outfile<<row.first<<","<<row.second<<"\n";
    }
}

template<typename T>
void writeCSV(const std::string &outfilename, const std::vector<T>&data, const int& precision=0){
    std::ofstream outfile(outfilename);
    if (precision){
        outfile.precision(precision);
    }
    std::copy(data.begin(), data.end(), std::ostream_iterator<T>(outfile,"\n"));
}

template<typename T>
void writeCSV(const std::string outfilename, const std::vector<std::vector<T>>& data, const int& precision=0){
    std::ofstream outfile(outfilename);
    if (precision){
        outfile.precision(precision);
    }
    for (auto &row : data){
        std::copy(row.begin(), row.end(), std::ostream_iterator<T>(outfile,","));
        outfile<<"\n";
    }
}


// Read CSV file
template<typename T>
void readCSV(const std::string& readfilename, std::vector<T>& data){
    std::ifstream readfile(readfilename);
    if (readfile.fail()){
        std::cout<<"No such file : "<<readfilename<<std::endl;
        exit(1);
    }
    data.clear();
    std::string row;
    while(getline(readfile,row)){
        data.emplace_back(std::stod(row));
    }
}

void readCSV(const std::string& readfilename, std::vector<std::vector<double>>& data){
    std::ifstream readfile(readfilename);
    if (readfile.fail()){
        std::cout<<"No such file : "<<readfilename<<std::endl;
        exit(1);
    }
    data.clear();
    std::string row, column;
    std::vector<double> rowvector;
    while(getline(readfile,row)){
        std::stringstream rowstream(row);
        rowvector.clear();
        while(getline(rowstream, column,',')){
            rowvector.emplace_back(std::stod(column));
        }
        data.emplace_back(rowvector);
    }
}

template<typename T, typename TT>
void readCSV(const std::string& readfilename, std::map<T,TT>& data){
    std::ifstream readfile(readfilename);
    if (readfile.fail()){
        std::cout<<"No such file : "<<readfilename<<std::endl;
        exit(1);
    }
    data.clear();
    std::string row,key,value;
    while(getline(readfile,row)){
        std::stringstream rowstream(row);
        int i=0;
        while(getline(rowstream, value,',')){
            if (i%2==1){
                data[std::stod(key)]=std::stod(value);
            }
            key=value;
            ++i;
        }
    }

}