#pragma once

#include <cmath>
#include <vector>
#include <map>
#include <utility>

#include "../library/linearAlgebra.hpp"

/*
    Data process code for binning
*/


namespace dataProcess{
    //* Log Binning
    std::map<double, double> logBin(const std::vector<double>& t_raw, const double& t_minExponent=0.0, const double& t_maxExponent=10.0, const double& t_deltaExponent=0.1);

    template <typename T>
    std::map<double, double> logBin(const std::map<T, double>& t_raw, const double& t_minExponent=0.0, const double& t_maxExponent=10.0, const double& t_deltaExponent=0.1);

    //* Log-Log Binning
    std::map<std::pair<double, double>, double> logLogBin(const std::map<std::pair<int, int>, double>& t_raw, const double& t_minExponent=0.0, const double& t_maxExponent=10.0, const double& t_deltaExponent=0.1);

    //* Linear Binning
    template <typename T>
    std::map<double, double> linBin(const std::map<T, double>& t_raw, const double& t_minVal=0.0, const double& t_maxVal=1.0, const double& t_delta=5e-4);

}//* End of namespace dataProcess


std::map<double, double> dataProcess::logBin(const std::vector<double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent){
    //* Setup values for integer log binning
    const std::vector<double> minBin = linearAlgebra::elementPow(10.0, linearAlgebra::arange(t_minExponent, t_maxExponent, t_deltaExponent));
    std::vector<double> value(minBin.size()-1, 0.0);
    for (unsigned i=0; i<minBin.size()-1; ++i){
        value[i] = std::sqrt(minBin[i] * minBin[i+1]);
    }

    //* Bin the data
    std::map<double, double> binned;
    std::map<double, unsigned> sampled;
    for (unsigned i=0; i<t_raw.size(); ++i){
        for (unsigned j=0; j<value.size(); ++j){
            if (i < minBin[j+1]){
                binned[value[j]] += t_raw[j];
                ++sampled[value[j]];
                break;
            }
        }
    }
    for (const auto& e : sampled){
        binned.at(e.first) /= e.second;
    }
    return binned;
}//* End of function dataProcess::logBin(vector)

template <typename T>
std::map<double, double> dataProcess::logBin(const std::map<T, double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent){
    //* Setup values for integer log binning
    const std::vector<double> minBin = linearAlgebra::elementPow(10.0, linearAlgebra::arange(t_minExponent, t_maxExponent, t_deltaExponent));
    std::vector<double> value(minBin.size()-1, 0.0);
    for (unsigned i=0; i<minBin.size()-1; ++i){
        value[i] = std::sqrt(minBin[i] * minBin[i+1]);
    }

    //* Bin the data
    std::map<double, double> binned;
    std::map<double, unsigned> sampled;
    for (const auto& raw : t_raw){
        for (unsigned j=0; j<value.size(); ++j){
            if (raw.first < minBin[j+1]){
                binned[value[j]] += raw.second;
                ++sampled[value[j]];
                break;
            }
        }
    }
    for (const auto& e : sampled){
        binned.at(e.first) /= e.second;
    }
    return binned;
}//* End of function dataProcess::logBin(map)


std::map<std::pair<double, double>, double> dataProcess::logLogBin(const std::map<std::pair<int, int>, double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent){
    //* Setup values for integer log binning
    const std::vector<double> minBin = linearAlgebra::elementPow(10.0, linearAlgebra::arange(t_minExponent, t_maxExponent, t_deltaExponent));
    std::vector<double> value(minBin.size()-1, 0.0);
    for (unsigned i=0; i<minBin.size()-1; ++i){
        value[i] = std::sqrt(minBin[i] * minBin[i+1]);
    }

    //* Bin the data
    std::map<std::pair<double, double>, double> binned;
    std::map<std::pair<double, double>, unsigned> sampled;
    for (const auto& raw : t_raw){
        for (unsigned j=0; j<value.size(); ++j){
            if (raw.first.first < minBin[j+1]){
                for (unsigned k=0; k<value.size(); ++k){
                    if (raw.first.second < minBin[k+1]){
                        binned[std::make_pair(value[j], value[k])] += raw.second;
                        ++sampled[std::make_pair(value[j], value[k])];
                        break;
                    }
                }
            }
        }
    }
    for (const auto& e : sampled){
        binned.at(e.first) /= e.second;
    }
    return binned;
}//* End of function dataProcess::logLogBin

template <typename T>
std::map<double, double> dataProcess::linBin(const std::map<T, double>& t_raw, const double& t_minVal, const double& t_maxVal, const double& t_delta){
    //* Setup values for integer log binning
    const std::vector<double> minBin = linearAlgebra::arange(t_minVal, t_maxVal, t_delta);
    std::vector<double> value(minBin.size()-1, 0.0);
    for (unsigned i=0; i<minBin.size(); ++i){
        value[i] = (minBin[i] + minBin[i+1])/2.0;
    }

    //* Bin the data
    std::map<double, double> binned;
    std::map<double, unsigned> sampled;
    for (const auto& raw : t_raw){
        for (unsigned j=0; j<value.size(); ++j){
            if (raw.first < minBin[j+1]){
                binned[value[j]] += raw.second;
                ++sampled[value[j]];
                break;
            }
        }
    }
    for (const auto& e : sampled){
        binned.at(e.first) /= e.second;
    }
    return binned;
}//* End of function dataProcess::linBin



