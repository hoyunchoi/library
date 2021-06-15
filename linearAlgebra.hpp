#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

namespace linearAlgebra {
//* transpose matrix
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& t_mat) {
    std::vector<std::vector<T>> result;
    result.resize(t_mat[0].size(), std::vector<T>(t_mat.size()));
    for (unsigned i = 0; i < t_mat.size(); ++i) {
        for (unsigned j = 0; j < t_mat[0].size(); ++j) {
            result[j][i] = t_mat[i][j];
        }
    }
    return result;
}

//* vector + vector -> vector
template <typename T>
std::vector<T> operator+(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    if (t_vec1.size() != t_vec2.size()) {
        std::cout << "In plus operation, two vectors have different length" << std::endl;
        exit(1);
    }

    std::vector<T> result(t_vec1.size());
    for (unsigned i = 0; i < t_vec1.size(); ++i) {
        result[i] = t_vec1[i] + t_vec2[i];
    }
    return result;
}
template <typename T>
std::vector<T>& operator+=(std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    if (t_vec1.empty()) {
        t_vec1.resize(t_vec2.size());
    } else if (t_vec1.size() != t_vec2.size()) {
        std::cout << "In plus operation, two vectors have different length" << std::endl;
        exit(1);
    }
    for (unsigned i = 0; i < t_vec1.size(); ++i) {
        t_vec1[i] += t_vec2[i];
    }
    return t_vec1;
}

//* vector - vector -> vector
template <typename T>
std::vector<T> operator-(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    if (t_vec1.size() != t_vec2.size()) {
        std::cout << "In plus operation, two vectors have different length" << std::endl;
        exit(1);
    }
    std::vector<T> result(t_vec1.size());
    for (unsigned i = 0; i < t_vec1.size(); ++i) {
        result[i] = t_vec1[i] - t_vec2[i];
    }
    return result;
}
template <typename T>
std::vector<T>& operator-=(std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    if (t_vec1.size() != t_vec2.size()) {
        std::cout << "In plus operation, two vectors have different length" << std::endl;
        exit(1);
    }
    for (unsigned i = 0; i < t_vec1.size(); ++i) {
        t_vec1[i] -= t_vec2[i];
    }
    return t_vec1;
}

//* matrix + matrix -> matrix
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    if (t_mat1.size() != t_mat2.size() || t_mat1[0].size() != t_mat2[0].size()) {
        std::cout << "In plus operation, two matrices have different size" << std::endl;
        exit(1);
    }
    std::vector<std::vector<T>> result;
    result.resize(t_mat1.size(), std::vector<T>(t_mat1[0].size()));
    for (unsigned i = 0; i < t_mat1.size(); ++i) {
        for (unsigned j = 0; j < t_mat1[0].size(); ++j) {
            result[i][j] = t_mat1[i][j] + t_mat2[i][j];
        }
    }
    return result;
}
template <typename T>
std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    if (t_mat1.size() != t_mat2.size() || t_mat1[0].size() != t_mat2[0].size()) {
        std::cout << "In plus operation, two matrices have different size" << std::endl;
        exit(1);
    }
    for (unsigned i = 0; i < t_mat1.size(); ++i) {
        for (unsigned j = 0; j < t_mat1[0].size(); ++j) {
            t_mat1[i][j] += t_mat2[i][j];
        }
    }
    return t_mat1;
}

//* matrix - matrix -> matrix
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    if (t_mat1.size() != t_mat2.size() || t_mat1[0].size() != t_mat2[0].size()) {
        std::cout << "In plus operation, two matrices have different size" << std::endl;
        exit(1);
    }

    std::vector<std::vector<T>> result;
    result.resize(t_mat1.size(), std::vector<T>(t_mat1[0].size()));
    for (unsigned i = 0; i < t_mat1.size(); ++i) {
        for (unsigned j = 0; j < t_mat1[0].size(); ++j) {
            result[i][j] = t_mat1[i][j] - t_mat2[i][j];
        }
    }
    return result;
}
template <typename T>
std::vector<std::vector<T>>& operator-=(std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    if (t_mat1.size() != t_mat2.size() || t_mat1[0].size() != t_mat2[0].size()) {
        std::cout << "In plus operation, two matrices have different size" << std::endl;
        exit(1);
    }
    for (unsigned i = 0; i < t_mat1.size(); ++i) {
        t_mat1[i] -= t_mat2[i];
    }
    return t_mat1;
}

//* const * vector -> vector
template <typename T, typename TT>
std::vector<T> operator*(const TT& c, const std::vector<T>& t_vec) {
    std::vector<T> result;
    result.reserve(t_vec.size());
    for (const T& e : t_vec) {
        result.emplace_back(c * e);
    }
    return result;
}
template <typename T, typename TT>
std::vector<T> operator*(const std::vector<T>& t_vec, const TT& c) {
    return c * t_vec;
}
template <typename T, typename TT>
std::vector<T>& operator*=(std::vector<T>& t_vec, const TT& c) {
    for (T& e : t_vec) {
        e *= c;
    }
    return t_vec;
}

//* const * matrix -> matrix
template <typename T, typename TT>
std::vector<std::vector<T>> operator*(const TT& c, const std::vector<std::vector<T>>& t_mat) {
    const unsigned n = t_mat.size();
    const unsigned m = t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            result[i][j] = c * t_mat[i][j];
        }
    }
    return result;
}
template <typename T, typename TT>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& t_mat, const TT& c) {
    return c * t_mat;
}
template <typename T, typename TT>
std::vector<std::vector<T>>& operator*=(std::vector<std::vector<T>>& t_mat, const TT& c) {
    for (std::vector<T>& e : t_mat) {
        e *= c;
    }
    return t_mat;
}

//* vector / const -> vector
template <typename T, typename TT>
std::vector<T> operator/(const std::vector<T>& t_vec, const TT& c) {
    std::vector<T> result(t_vec.size());
    for (unsigned i = 0; i < t_vec.size(); ++i) {
        result[i] = t_vec[i] / c;
    }
    return result;
}
template <typename T, typename TT>
std::vector<T>& operator/=(std::vector<T>& t_vec, const TT& c) {
    for (T& e : t_vec) {
        e /= c;
    }
    return t_vec;
}

//* matrix / c -> matrix
template <typename T, typename TT>
std::vector<std::vector<T>> operator/(const std::vector<std::vector<T>>& t_mat, const TT& c) {
    const unsigned n = t_mat.size();
    const unsigned m = t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            result[i][j] = t_mat[i][j] / c;
        }
    }
    return result;
}
template <typename T, typename TT>
std::vector<std::vector<T>>& operator/=(std::vector<std::vector<T>>& t_mat, const TT& c) {
    for (std::vector<T>& e : t_mat) {
        e /= c;
    }
    return t_mat;
}

//* vector + c -> vector
template <typename T, typename TT>
std::vector<T> operator+(const std::vector<T>& t_vec, const TT& c) {
    const unsigned n = t_vec.size();
    std::vector<T> result(n);
    for (unsigned i = 0; i < n; ++i) {
        result[i] = t_vec[i] + c;
    }
    return result;
}
template <typename T, typename TT>
std::vector<T>& operator+=(std::vector<T>& t_vec, const TT& c) {
    for (T& e : t_vec) {
        e += c;
    }
    return t_vec;
}

//* vector - const -> vector
template <typename T, typename TT>
std::vector<T> operator-(const std::vector<T>& t_vec, const TT& c) {
    std::vector<T> result(t_vec.size());
    for (unsigned i = 0; i < t_vec.size(); ++i) {
        result[i] = t_vec[i] - c;
    }
    return result;
}
template <typename T, typename TT>
std::vector<T> operator-(const TT& c, const std::vector<T>& t_vec) {
    std::vector<T> result(t_vec.size());
    for (unsigned i = 0; i < t_vec.size(); ++i) {
        result[i] = c - t_vec[i];
    }
    return result;
}
template <typename T, typename TT>
std::vector<T>& operator-=(std::vector<T>& t_vec, const TT& c) {
    for (T& e : t_vec) {
        e -= c;
    }
    return t_vec;
}

//* const + matrix -> matrix
template <typename T, typename TT>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& t_mat, const TT& c) {
    const unsigned n = t_mat.size();
    const unsigned m = t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            result[i][j] = t_mat[i][j] + c;
        }
    }
    return result;
}
template <typename T, typename TT>
std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& t_mat, const TT& c) {
    for (T& e : t_mat) {
        e += c;
    }
    return t_mat;
}

//* matrix - const -> matrix
template <typename T, typename TT>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& t_mat, const TT& c) {
    const unsigned n = t_mat.size();
    const unsigned m = t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            result[i][j] = t_mat[i][j] - c;
        }
    }
    return result;
}
template <typename T, typename TT>
std::vector<std::vector<T>>& operator-=(std::vector<std::vector<T>>& t_mat, const TT& c) {
    for (std::vector<T>& e : t_mat) {
        e -= c;
    }
    return t_mat;
}

//* vector dot vector -> const
template <typename T>
T innerProduct(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    if (t_vec1.size() != t_vec2.size()) {
        std::cout << "In inner product, two vectors have different size" << std::endl;
        exit(1);
    }
    T result = 0;
    for (unsigned i = 0; i < t_vec1.size(); ++i) {
        result += t_vec1[i] * t_vec2[i];
    }
    return result;
}

//* vector * vector -> matrix
template <typename T>
std::vector<std::vector<T>> product(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    const unsigned n = t_vec1.size();
    const unsigned m = t_vec2.size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            result[i][j] = t_vec1[i] * t_vec2[j];
        }
    }
    return result;
}
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    return product(t_vec1, t_vec2);
}

//* matrix * vector -> vector
template <typename T>
std::vector<T> product(const std::vector<std::vector<T>>& t_mat, const std::vector<T>& t_vec) {
    if (t_mat[0].size() != t_vec.size()) {
        std::cout << "In product between matrix and vector, size is different" << std::endl;
        exit(1);
    }
    const unsigned n = t_mat.size();
    const unsigned m = t_mat[0].size();
    std::vector<T> result(n);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            result[i] += t_mat[i][j] * t_vec[j];
        }
    }
    return result;
}
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& t_mat, const std::vector<T>& t_vec) {
    return product(t_mat, t_vec);
}

//* matrix * matrix -> matrix
template <typename T>
std::vector<std::vector<T>> product(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    if (t_mat1[0].size() != t_mat2.size()) {
        std::cout << "In product between matrix and matrix, size is different" << std::endl;
        exit(1);
    }
    const unsigned n = t_mat1.size();
    const unsigned m = t_mat2[0].size();
    const unsigned k_max = t_mat2.size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));

    std::vector<std::vector<T>> tempMatrix = transpose(t_mat2);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            T temp = 0;
            for (unsigned k = 0; k < k_max; ++k) {
                temp += t_mat1[i][k] * tempMatrix[j][k];
            }
            result[i][j] = temp;
        }
    }
    return result;
}
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    return product(t_mat1, t_mat2);
}

//* Elementwise Calculation
template <typename T>
std::vector<T> elementProduct(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2) {
    if (t_vec1.size() != t_vec2.size()) {
        std::cout << "In inner product, two vectors have different size" << std::endl;
        exit(1);
    }
    const unsigned n = t_vec1.size();
    std::vector<T> result(n);
    for (unsigned i = 0; i < n; ++i) {
        result[i] = t_vec1[i] * t_vec2[i];
    }
    return result;
}

template <typename T>
std::vector<std::vector<T>> elementProduct(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2) {
    if (t_mat1.size() != t_mat2.size() || t_mat1[0].size() != t_mat2[0].size()) {
        std::cout << "In plus operation, two matrices have different size" << std::endl;
        exit(1);
    }
    const unsigned n = t_mat1.size();
    const unsigned m = t_mat1[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n, std::vector<T>(m));
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < m; ++j) {
            result[i][j] = t_mat1[i][j] * t_mat2[i][j];
        }
    }
}

template <typename T>
std::vector<T> elementTanh(const std::vector<T>& t_vec) {
    std::vector<T> result;
    result.reserve(t_vec.size());
    for (const T& e : t_vec) {
        result.emplace_back(tanh(e));
    }
    return result;
}

template <typename T, typename TT>
std::vector<T> elementPow(const std::vector<T>& t_vec, const TT& c) {
    std::vector<T> result;
    result.reserve(t_vec.size());
    for (const T& e : t_vec) {
        result.emplace_back(pow(e, c));
    }
    return result;
}
template <typename T, typename TT>
std::vector<T> elementPow(const TT& c, const std::vector<T>& t_vec) {
    std::vector<T> result;
    result.reserve(t_vec.size());
    for (const T& e : t_vec) {
        result.emplace_back(pow(c, e));
    }
    return result;
}

//* Include t_begin, t_end
std::vector<double> linspace(const double& t_begin, const double& t_end, const unsigned& t_size) {
    std::vector<double> result(t_size);
    const double delta = (t_end - t_begin) / (t_size - 1);
    for (unsigned i = 0; i < t_size; ++i) {
        result[i] = delta * i + t_begin;
    }
    return result;
}

//* Include t_begin, t_end
template <typename T>
std::vector<T> arange(const T& t_begin, const T& t_end, const T& t_delta) {
    if (t_begin > t_end) {
        std::cout << "In arange, input values are not ordered\n";
        exit(1);
    }
    const unsigned size = std::ceil((t_end - t_begin) / t_delta);
    std::vector<T> result(size);
    for (unsigned i = 0; i < size; ++i) {
        result[i] = t_begin + t_delta * i;
    }
    result.emplace_back(t_end);
    return result;
}

template <typename T>
T accumulate(const std::vector<T>& t_vec) {
    return std::accumulate(t_vec.begin(), t_vec.end(), (T)0);
}

//* Pop back of vector
template <typename T>
void popBack(std::vector<T>& t_vec) {
    t_vec.pop_back();
    return;
}

//! Map Calculation
//* map += map
template <typename T, typename TT>
std::map<T, TT> operator+=(std::map<T, TT>& t_map1, const std::map<T, TT>& t_map2) {
    for (const auto& e : t_map2) {
        t_map1[e.first] += e.second;
    }
    return t_map1;
}

//* map -= map
template <typename T, typename TT>
std::map<T, TT> operator-=(std::map<T, TT>& t_map1, const std::map<T, TT>& t_map2) {
    for (const auto& e : t_map2) {
        t_map1[e.first] -= e.second;
    }
    return t_map1;
}

//* const - map = map
template <typename T, typename TT>
std::map<T, TT> minus_first(const T& t_c, const std::map<T, TT>& t_map) {
    std::map<T, TT> result;
    for (const auto& e : t_map) {
        result[t_c - e.first] = e.second;
    }
    return result;
}

template <typename T, typename TT>
std::map<T, TT> minus_second(const TT& t_c, const std::map<T, TT>& t_map) {
    std::map<T, TT> result;
    for (const auto& e : t_map) {
        result[e.first] = t_c - e.second;
    }
    return result;
}

//* map - const = map
template <typename T, typename TT>
std::map<T, TT> minus_first(const std::map<T, TT>& t_map, const T& t_c) {
    std::map<T, TT> result;
    for (const auto& e : t_map) {
        result[e.first - t_c] = e.second;
    }
    return result;
}

template <typename T, typename TT>
std::map<T, TT> minus_second(const std::map<T, TT>& t_map, const T& t_c) {
    std::map<T, TT> result;
    for (const auto& e : t_map) {
        result[e.first] = e.second - t_c;
    }
    return result;
}

//* map * const = map
template <typename T, typename TT>
std::map<T, TT> operator*(const std::map<T, TT>& t_map, const TT& t_c) {
    std::map<T, TT> result;
    for (const auto& e : t_map) {
        result[e.first] = e.second * t_c;
    }
    return result;
}

//* const * map = map
template <typename T, typename TT>
std::map<T, TT> operator*(const TT& t_c, const std::map<T, TT>& t_map) {
    return t_map * t_c;
}

//* map *= const
template <typename T, typename TT>
std::map<T, TT> operator*=(std::map<T, TT>& t_map, const TT& t_c) {
    for (auto e : t_map) {
        e.second *= t_c;
    }
    return t_map;
}

//* map *= map
template <typename T, typename TT>
std::map<T, TT> operator*=(std::map<T, TT>& t_map1, const std::map<T, TT>& t_map2) {
    for (const auto& e : t_map2) {
        t_map1[e.first] *= e.second;
    }
    return t_map1;
}

//* map /= map
template <typename T, typename TT>
std::map<T, TT> operator/=(std::map<T, TT>& t_map1, const std::map<T, TT>& t_map2) {
    for (const auto& e : t_map2) {
        t_map1[e.first] /= e.second;
    }
    return t_map1;
}

//* map /= const
template <typename T, typename TT>
std::map<T, TT> operator/=(std::map<T, TT>& t_map, const TT& t_c) {
    for (auto& e : t_map) {
        e.second /= t_c;
    }
    return t_map;
}

//* plus 1 to the 'value' of t_map1['key'] for every 'key' inside t_map2
template <typename T, typename TT>
void sampleNum(std::map<T, int>& t_map1, const std::map<T, TT>& t_map2) {
    for (const auto& e : t_map2) {
        ++t_map1[e.first];
    }
}

//* sum all values of Map
template <typename T, typename TT>
TT accumulate(const std::map<T, TT>& t_map) {
    TT result = 0;
    for (const auto& e : t_map) {
        result += e.second;
    }
    return result;
}

//* pop back of map
template <typename T, typename TT>
void popBack(std::map<T, TT>& t_map) {
    t_map.erase(t_map.end()->first);
    return;
}

//* print
template <typename T>
void print(const std::vector<T>& t_vec, const int& t_precision = 6, const char& t_seperate = ',') {
    t_precision < 0 ? std::cout.precision(std::numeric_limits<double>::digits10 + 1) : std::cout.precision(t_precision);
    std::cout << "[";
    for (unsigned i = 0; i < t_vec.size() - 1; ++i) {
        std::cout << t_vec[i] << t_seperate;
    }
    std::cout << t_vec.back() << "]\n";
    std::cout.precision(6);
}

template <typename T>
void print(const std::vector<std::vector<T>>& t_mat, const int& t_precision = 6, const char& t_seperate = ',') {
    t_precision < 0 ? std::cout.precision(std::numeric_limits<double>::digits10 + 1) : std::cout.precision(t_precision);
    std::cout << "[";
    for (unsigned i = 0; i < t_mat.size() - 1; ++i) {
        std::cout << "[";
        for (unsigned j = 0; j < t_mat[i].size() - 1; ++j) {
            std::cout << t_mat[i][j] << t_seperate;
        }
        std::cout << t_mat[i].back() << "]\n";
    }
    for (unsigned j = 0; j < t_mat.back().size() - 1; ++j) {
        std::cout << t_mat.back()[j] << t_seperate;
    }
    std::cout << t_mat.back().back() << "]]\n";
    std::cout.precision(6);
}

template <typename T, typename TT>
void print(const std::map<T, TT>& t_map, const int& t_precision = 6, const char& t_seperate = ',') {
    t_precision < 0 ? std::cout.precision(std::numeric_limits<double>::digits10 + 1) : std::cout.precision(t_precision);
    for (const auto& e : t_map) {
        std::cout << e.first << "," << e.second << "\n";
    }
    std::cout.precision(6);
}

}  // namespace linearAlgebra