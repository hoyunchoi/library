#pragma once

#include <cmath>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "../library/linearAlgebra.hpp"

/*
    Data process code for binning
*/

namespace dataProcess {
//! Setup variables for binning
//* Return tuple (minBin, value, binSize)
const std::tuple<std::vector<int>, std::vector<double>, std::vector<unsigned>> setLogBin_(const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent);

//* Return tuple(minBin, value)
const std::tuple<std::vector<double>, std::vector<double>> setLinBin_(const double& t_minVal, const double& t_maxVal, const double& t_delta);

//! Log Binning
//* Log Binning distribution
std::map<double, double> distLogBin(const std::vector<double>& t_raw, const double& t_minExponent = 0.0, const double& t_maxExponent = 10.0, const double& t_deltaExponent = 0.1);
template <typename T>
std::map<double, double> distLogBin(const std::map<T, double>& t_raw, const double& t_minExponent = 0.0, const double& t_maxExponent = 10.0, const double& t_deltaExponent = 0.1);

//* Log Binning average
std::map<double, double> avgLogBin(const std::vector<double>& t_raw, const double& t_minExponent = 0.0, const double& t_maxExponent = 10.0, const double& t_deltaExponent = 0.1);
template <typename T>
std::map<double, double> avgLogBin(const std::map<T, double>& t_raw, const double& t_minExponent = 0.0, const double& t_maxExponent = 10.0, const double& t_deltaExponent = 0.1);

//* Log-Log Binning distribution
std::map<std::pair<double, double>, double> distLogLogBin(const std::map<std::pair<int, int>, double>& t_raw, const double& t_minExponent = 0.0, const double& t_maxExponent = 10.0, const double& t_deltaExponent = 0.1);
std::map<std::pair<double, double>, double> avgLogLogBin(const std::map<std::pair<int, int>, double>& t_raw, const double& t_minExponent = 0.0, const double& t_maxExponent = 10.0, const double& t_deltaExponent = 0.1);

//* Linear Binning
template <typename T>
std::map<double, double> distLinBin(const std::map<T, double>& t_raw, const double& t_minVal = 0.0, const double& t_maxVal = 1.0, const double& t_delta = 5e-4);
template <typename T>
std::map<double, double> avgLinBin(const std::map<T, double>& t_raw, const double& t_minVal = 0.0, const double& t_maxVal = 1.0, const double& t_delta = 5e-4);
}  // namespace dataProcess

const std::tuple<std::vector<int>, std::vector<double>, std::vector<unsigned>> dataProcess::setLogBin_(const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    const std::vector<double> interval = linearAlgebra::elementPow(10.0, linearAlgebra::arange(t_minExponent, t_maxExponent, t_deltaExponent));
    std::vector<int> minBin;
    minBin.reserve(interval.size() - 1);
    std::vector<double> value;
    value.reserve(interval.size() - 1);
    std::vector<unsigned> binSize;
    binSize.reserve(interval.size() - 1);
    for (unsigned i = 0; i < interval.size() - 1; ++i) {
        const unsigned size = std::ceil(interval[i + 1]) - std::ceil(interval[i]);
        if (size) {
            minBin.emplace_back(std::ceil(interval[i]));
            value.emplace_back(std::sqrt((double)std::ceil(interval[i]) * (double)std::ceil(interval[i + 1])));
            binSize.emplace_back(size);
        }
    }
    return std::make_tuple(minBin, value, binSize);
}  //* End of function dataProcess::setLogBin_

const std::tuple<std::vector<double>, std::vector<double>> dataProcess::setLinBin_(const double& t_minVal, const double& t_maxVal, const double& t_delta) {
    const std::vector<double> minBin = linearAlgebra::arange(t_minVal, t_maxVal, t_delta);
    std::vector<double> value;
    value.reserve(minBin.size() - 1);
    for (unsigned i = 0; i < minBin.size(); ++i) {
        value.emplace_back((minBin[i] + minBin[i + 1]) / 2.0);
    }
    return std::make_tuple(minBin, value);
}  //* End of function dataProcess::setLinBin_

std::map<double, double> dataProcess::distLogBin(const std::vector<double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    using namespace linearAlgebra;
    std::map<double, double> binned;
    const auto [minBin, value, binSize] = setLogBin_(t_minExponent, t_maxExponent, t_deltaExponent);
    for (unsigned i = 0; i < t_raw.size(); ++i) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if ((int)i < minBin[j + 1]) {
                binned[value[j]] += t_raw[i] / (double)binSize[j];
                break;
            }
        }
    }
    for (auto it = binned.begin(); it != binned.end();) {
        it->second == 0 ? binned.erase(it++) : ++it;
    }
    binned /= linearAlgebra::accumulate(binned);
    return binned;
}

template <typename T>
std::map<double, double> dataProcess::distLogBin(const std::map<T, double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    using namespace linearAlgebra;
    std::map<double, double> binned;
    const auto [minBin, value, binSize] = setLogBin_(t_minExponent, t_maxExponent, t_deltaExponent);
    for (const auto& raw : t_raw) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if (raw.first < minBin[j + 1]) {
                binned[value[j]] += raw.second / (double)binSize[j];
                break;
            }
        }
    }
    binned /= linearAlgebra::accumulate(binned);
    return binned;
}

std::map<double, double> dataProcess::avgLogBin(const std::vector<double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    std::map<double, double> binned;
    std::map<double, unsigned> sampled;
    const auto [minBin, value, binSize] = setLogBin_(t_minExponent, t_maxExponent, t_deltaExponent);
    (void)binSize;
    for (unsigned i = 0; i < t_raw.size(); ++i) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if ((int)i < minBin[j + 1]) {
                binned[value[j]] += t_raw[j];
                ++sampled[value[j]];
                break;
            }
        }
    }
    for (const auto& e : sampled) {
        binned.at(e.first) /= e.second;
    }
    return binned;
}  //* End of function dataProcess::logBin(vector)

template <typename T>
std::map<double, double> dataProcess::avgLogBin(const std::map<T, double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    std::map<double, double> binned;
    std::map<double, unsigned> sampled;
    const auto [minBin, value, binSize] = setLogBin_(t_minExponent, t_maxExponent, t_deltaExponent);
    (void)binSize;
    for (const auto& raw : t_raw) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if (raw.first < minBin[j + 1]) {
                binned[value[j]] += raw.second;
                ++sampled[value[j]];
                break;
            }
        }
    }
    for (const auto& e : sampled) {
        binned.at(e.first) /= e.second;
    }
    return binned;
}  //* End of function dataProcess::logBin(map)

std::map<std::pair<double, double>, double> dataProcess::distLogLogBin(const std::map<std::pair<int, int>, double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    using namespace linearAlgebra;
    std::map<std::pair<double, double>, double> binned;
    const auto [minBin, value, binSize] = setLogBin_(t_minExponent, t_maxExponent, t_deltaExponent);
    for (const auto& raw : t_raw) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if (raw.first.first < minBin[j + 1]) {
                for (unsigned k = 0; k < value.size(); ++k) {
                    if (raw.first.second < minBin[k + 1]) {
                        binned[std::make_pair(value[j], value[k])] += raw.second / (binSize[j] * binSize[k]);
                        break;
                    }
                }
            }
        }
    }
    binned /= linearAlgebra::accumulate(binned);
    return binned;
}  //* End of function dataProcess::distlogLogBin

std::map<std::pair<double, double>, double> dataProcess::avgLogLogBin(const std::map<std::pair<int, int>, double>& t_raw, const double& t_minExponent, const double& t_maxExponent, const double& t_deltaExponent) {
    std::map<std::pair<double, double>, double> binned;
    std::map<std::pair<double, double>, unsigned> sampled;
    const auto [minBin, value, binSize] = setLogBin_(t_minExponent, t_maxExponent, t_deltaExponent);
    (void)binSize;
    for (const auto& raw : t_raw) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if (raw.first.first < minBin[j + 1]) {
                for (unsigned k = 0; k < value.size(); ++k) {
                    if (raw.first.second < minBin[k + 1]) {
                        binned[std::make_pair(value[j], value[k])] += raw.second;
                        ++sampled[std::make_pair(value[j], value[k])];
                        break;
                    }
                }
            }
        }
    }
    for (const auto& e : sampled) {
        binned.at(e.first) /= e.second;
    }
    return binned;
}  //* End of function dataProcess::avgLogLogBin

template <typename T>
std::map<double, double> dataProcess::distLinBin(const std::map<T, double>& t_raw, const double& t_minVal, const double& t_maxVal, const double& t_delta) {
    using namespace linearAlgebra;
    std::map<double, double> binned;
    const auto [minBin, value] = setLinBin_(t_minVal, t_maxVal, t_delta);
    for (const auto& raw : t_raw) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if (raw.first < minBin[j + 1]) {
                binned[value[j]] += raw.second;
                break;
            }
        }
    }
    binned /= linearAlgebra::accumulate(binned);
    return binned;
}  //* End of function dataProcess::distLinBin

template <typename T>
std::map<double, double> dataProcess::avgLinBin(const std::map<T, double>& t_raw, const double& t_minVal, const double& t_maxVal, const double& t_delta) {
    std::map<double, double> binned;
    std::map<double, unsigned> sampled;
    const auto [minBin, value] = setLinBin_(t_minVal, t_maxVal, t_delta);
    for (const auto& raw : t_raw) {
        for (unsigned j = 0; j < value.size(); ++j) {
            if (raw.first < minBin[j + 1]) {
                binned[value[j]] += raw.second;
                ++sampled[value[j]];
                break;
            }
        }
    }
    for (const auto& e : sampled) {
        binned.at(e.first) /= e.second;
    }
    return binned;
}  //* End of function dataProcess::avgLinBin
