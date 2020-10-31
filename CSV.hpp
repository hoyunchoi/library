#pragma once
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <iterator>
#include <iostream>
#include <cstdio>

//* Write and Read CSV file
namespace CSV{
    //* Write
    template<typename T, typename TT>
    void write(const std::string& writeFileName, const std::map<T,TT>& data, const int& precision=0){
        std::ofstream outFile(writeFileName);
        if (precision){
            outFile.precision(precision);
        }
        for (auto& row : data){
            outFile<<row.first<<","<<row.second<<"\n";
        }
    }

    template<typename T>
    void write(const std::string &writeFileName, const std::vector<T>&data, const int& precision=0){
        std::ofstream outFile(writeFileName);
        if (precision){
            outFile.precision(precision);
        }
        std::copy(data.begin(), data.end(), std::ostream_iterator<T>(outFile,"\n"));
    }

    template<typename T>
    void write(const std::string writeFileName, const std::vector<std::vector<T>>& data, const int& precision=0){
        std::ofstream outFile(writeFileName);
        if (precision){
            outFile.precision(precision);
        }
        for (auto &row : data){
            std::copy(row.begin(), row.end(), std::ostream_iterator<T>(outFile,","));
            outFile<<"\n";
        }
    }


    //* Read
    template<typename T>
    void read(const std::string& readFileName, std::vector<T>& data){
        std::ifstream readFile(readFileName);
        if (readFile.fail()){
            std::cout<<"No such file : "<<readFileName<<std::endl;
            exit(1);
        }
        data.clear();
        std::string row;
        while(getline(readFile,row)){
            data.emplace_back(std::stod(row));
        }
    }

    template<typename T>
    void read(const std::string& readFileName, std::vector<std::vector<T>>& data){
        std::ifstream readFile(readFileName);
        if (readFile.fail()){
            std::cout<<"No such file : "<<readFileName<<std::endl;
            exit(1);
        }
        data.clear();
        std::string row, column;
        std::vector<T> rowvector;
        while(getline(readFile,row)){
            std::stringstream rowstream(row);
            rowvector.clear();
            while(getline(rowstream, column,',')){
                rowvector.emplace_back(std::stod(column));
            }
            data.emplace_back(rowvector);
        }
    }

    template<typename T, typename TT>
    void read(const std::string& readFileName, std::map<T,TT>& data){
        std::ifstream readFile(readFileName);
        if (readFile.fail()){
            std::cout<<"No such file : "<<readFileName<<std::endl;
            exit(1);
        }
        data.clear();
        std::string row,key,value;
        while(getline(readFile,row)){
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

    //* Delete File
    void deleteFile(const std::string& t_deleteFileName){
        if (std::remove(t_deleteFileName.c_str())){
            std::cout << "Error in removing " << t_deleteFileName << "\n";
        }
    }




}//* End of namespace CSV