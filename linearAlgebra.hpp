#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>


namespace linearAlgebra{
template<typename T>
void print(const std::vector<T>& t_vec){
    std::cout<<"(";
    for (int i=0; i<t_vec.size(); ++i){
        std::cout<<t_vec[i]<<", ";
    }
    std::cout<<")\n";
}

template<typename T>
void print(const std::vector<std::vector<T>>& t_mat){
    for (int i=0; i<t_mat.size(); ++i){
        for (int j=0; j<t_mat[0].size(); ++j){
            std::cout<<t_mat[i][j]<<"\t";
        }
        std::cout<<"\n";
    }
}

//* transpose matrix
template<typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& t_mat){
    std::vector<std::vector<T>> result;
    result.resize(t_mat[0].size(),std::vector<T>(t_mat.size()));
    for (int i=0; i<t_mat.size(); ++i){
        for (int j=0; j<t_mat[0].size();++j){
            result[j][i]=t_mat[i][j];
        }
    }
    return result;
}

//* vector + vector -> vector
template<typename T>
std::vector<T> operator+ (const std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    if (t_vec1.size()!=t_vec2.size()){
        std::cout<<"In plus operation, two vectors have different length"<<std::endl;
        exit(1);
    }

    std::vector<T> result(t_vec1.size());
    for (int i=0; i<t_vec1.size(); ++i){
        result[i]=t_vec1[i]+t_vec2[i];
    }
    return result;
}
template<typename T>
std::vector<T>& operator+= (std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    if (t_vec1.size()!=t_vec2.size()){
        std::cout<<"In plus operation, two vectors have different length"<<std::endl;
        exit(1);
    }
    for (int i=0; i<t_vec1.size(); ++i){
        t_vec1[i]+=t_vec2[i];
    }
    return t_vec1;
}

//* vector - vector -> vector
template<typename T>
std::vector<T> operator- (const std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    if (t_vec1.size()!=t_vec2.size()){
        std::cout<<"In plus operation, two vectors have different length"<<std::endl;
        exit(1);
    }
    const int n=t_vec1.size();
    std::vector<T> result(n);
    for (int i=0; i<n; ++i){
        result[i]=t_vec1[i]-t_vec2[i];
    }
    return result;
}
template<typename T>
std::vector<T>& operator-= (std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    if (t_vec1.size()!=t_vec2.size()){
        std::cout<<"In plus operation, two vectors have different length"<<std::endl;
        exit(1);
    }
    for (int i=0; i<t_vec1.size(); ++i){
        t_vec1[i]-=t_vec2[i];
    }
    return t_vec1;
}

//* matrix + matrix -> matrix
template<typename T>
std::vector<std::vector<T>> operator+ (const std::vector<std::vector<T>>& t_mat1, 
                                       const std::vector<std::vector<T>>& t_mat2){
    if (t_mat1.size()!=t_mat2.size() || t_mat1[0].size()!=t_mat2[0].size()){
        std::cout<<"In plus operation, two matrices have different size"<<std::endl;
        exit(1);
    }
    std::vector<std::vector<T>> result;
    result.resize(t_mat1.size(),std::vector<T>(t_mat1[0].size()));
    for (int i=0; i<t_mat1.size(); ++i){
        for (int j=0; j<t_mat1[0].size(); ++j){
            result[i][j]=t_mat1[i][j]+t_mat2[i][j];
        }
    }
    return result;
}
template<typename T>
std::vector<std::vector<T>>& operator+= (std::vector<std::vector<T>>& t_mat1, 
                                       const std::vector<std::vector<T>>& t_mat2){
    if (t_mat1.size()!=t_mat2.size() || t_mat1[0].size()!=t_mat2[0].size()){
        std::cout<<"In plus operation, two matrices have different size"<<std::endl;
        exit(1);
    }
    for (int i=0; i<t_mat1.size(); ++i){
        for (int j=0; j<t_mat1[0].size(); ++j){
            t_mat1[i][j]+=t_mat2[i][j];
        }
    }
    return t_mat1;
}

//* matrix - matrix -> matrix
template<typename T>
std::vector<std::vector<T>> operator- (const std::vector<std::vector<T>>& t_mat1, 
                                       const std::vector<std::vector<T>>& t_mat2){
    if (t_mat1.size()!=t_mat2.size() || t_mat1[0].size()!=t_mat2[0].size()){
        std::cout<<"In plus operation, two matrices have different size"<<std::endl;
        exit(1);
    }

    std::vector<std::vector<T>> result;
    result.resize(t_mat1.size(),std::vector<T>(t_mat1[0].size()));
    for (int i=0; i<t_mat1.size(); ++i){
        for (int j=0; j<t_mat1[0].size(); ++j){
            result[i][j]=t_mat1[i][j]-t_mat2[i][j];
        }
    }
    return result;
}
template<typename T>
std::vector<std::vector<T>>& operator-= (std::vector<std::vector<T>>& t_mat1, 
                                       const std::vector<std::vector<T>>& t_mat2){
    if (t_mat1.size()!=t_mat2.size() || t_mat1[0].size()!=t_mat2[0].size()){
        std::cout<<"In plus operation, two matrices have different size"<<std::endl;
        exit(1);
    }
    for (int i=0; i<t_mat1.size(); ++i){
        t_mat1[i]-=t_mat2[i];
    }
    return t_mat1;
}


//* const * vector -> vector
template<typename T, typename TT>
std::vector<T> operator*(const TT& c, const std::vector<T>& t_vec){
    std::vector<T> result(t_vec.size());
    for (int i=0; i<t_vec.size(); ++i){
        result[i]=c*t_vec[i];
    }
    return result;
}
template<typename T, typename TT>
std::vector<T> operator*(const std::vector<T>& t_vec, const TT& c){
    return c*t_vec;
}
template<typename T, typename TT>
std::vector<T>& operator*=(std::vector<T>& t_vec, const TT& c){
    for (T& e : t_vec){
        e*=c;
    }
    return t_vec;
}

//* const * matrix -> matrix
template<typename T, typename TT>
std::vector<std::vector<T>> operator*(const TT& c, const std::vector<std::vector<T>>& t_mat){
    const int n=t_mat.size();
    const int m=t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            result[i][j]=c*t_mat[i][j];
        }
    }
    return result;
}
template<typename T, typename TT>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& t_mat, const TT& c){
    return c*t_mat;
}
template<typename T, typename TT>
std::vector<std::vector<T>>& operator*= (std::vector<std::vector<T>>& t_mat, const TT& c){
    for (std::vector<T>& e : t_mat){
        e*=c;
    }
    return t_mat;
}


//* vector / const -> vector
template<typename T, typename TT>
std::vector<T> operator/(const std::vector<T>& t_vec, const TT& c){
    std::vector<T> result(t_vec.size());
    for (int i=0; i<t_vec.size(); ++i){
        result[i] = t_vec[i]/c;
    }
    return result;
}
template<typename T, typename TT>
std::vector<T>& operator/=(std::vector<T>& t_vec, const TT& c){
    for (T& e : t_vec){
        e/=c;
    }
    return t_vec;
}

//* matrix / c -> matrix
template<typename T, typename TT>
std::vector<std::vector<T>> operator/(const std::vector<std::vector<T>>& t_mat, const TT& c){
    const int n=t_mat.size();
    const int m=t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            result[i][j] = t_mat[i][j]/c;
        }
    }
    return result;
}
template<typename T, typename TT>
std::vector<std::vector<T>>& operator/= (std::vector<std::vector<T>>& t_mat, const TT& c){
    for (std::vector<T>& e : t_mat){
        e /= c;
    }
    return t_mat;
}

//* vector + c -> vector
template<typename T, typename TT>
std::vector<T> operator+(const std::vector<T>& t_vec, const TT& c){
    const int n=t_vec.size();
    std::vector<T> result(n);
    for (int i=0; i<n; ++i){
        result[i]=t_vec[i]+c;
    }
    return result;
}
template<typename T, typename TT>
std::vector<T>& operator+= (std::vector<T>& t_vec, const TT& c){
    for (T& e : t_vec){
        e+=c;
    }
    return t_vec;
}

//* vector - const -> vector
template<typename T, typename TT>
std::vector<T> operator-(const std::vector<T>& t_vec, const TT& c){
    std::vector<T> result(t_vec.size());
    for (int i=0; i<t_vec.size(); ++i){
        result[i]=t_vec[i]-c;
    }
    return result;
}
template <typename T, typename TT>
std::vector<T> operator- (const TT&c, const std::vector<T>& t_vec){
    std::vector<T> result(t_vec.size());
    for (int i=0; i<t_vec.size(); ++i){
        result[i]=c-t_vec[i];
    }
    return result;
}
template<typename T, typename TT>
std::vector<T>& operator-= (std::vector<T>& t_vec, const TT& c){
    for (T& e : t_vec){
        e-=c;
    }
    return t_vec;
}

//* const + matrix -> matrix
template<typename T, typename TT>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& t_mat, const TT& c){
    const int n=t_mat.size();
    const int m=t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            result[i][j]=t_mat[i][j]+c;
        }
    }
    return result;
}
template<typename T, typename TT>
std::vector<std::vector<T>>& operator+= (std::vector<std::vector<T>>& t_mat, const TT& c){
    for (T& e : t_mat){
        e+=c;
    }
    return t_mat;
}

//* matrix - const -> matrix
template<typename T, typename TT>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& t_mat, const TT& c){
    const int n=t_mat.size();
    const int m=t_mat[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            result[i][j]=t_mat[i][j]-c;
        }
    }
    return result;
}
template<typename T, typename TT>
std::vector<std::vector<T>>& operator-= (std::vector<std::vector<T>>& t_mat, const TT& c){
    for (std::vector<T>& e : t_mat){
        e-=c;
    }
    return t_mat;
}


//* vector dot vector -> const
template<typename T>
T innerProduct(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    // Check two vectors have same length
    if (t_vec1.size()!=t_vec2.size()){
        std::cout<<"In inner product, two vectors have different size"<<std::endl;
        exit(1);
    }
    T result=0;
    for (int i=0; i<t_vec1.size(); ++i){
        result+=t_vec1[i]*t_vec2[i];
    }
    return result;
}

//* vector * vector -> matrix
template<typename T>
std::vector<std::vector<T>> product(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    const int n=t_vec1.size();
    const int m=t_vec2.size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            result[i][j]=t_vec1[i]*t_vec2[j];
        }
    }
    return result;
}
template<typename T>
std::vector<std::vector<T>> operator* (const std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    return product(t_vec1, t_vec2);
}

//* matrix * vector -> vector
template<typename T>
std::vector<T> product(const std::vector<std::vector<T>>& t_mat, const std::vector<T>& t_vec){
    if (t_mat[0].size()!=t_vec.size()){
        std::cout<<"In product between matrix and vector, size is different"<<std::endl;
        exit(1);
    }
    const int n=t_mat.size();
    const int m=t_mat[0].size();
    std::vector<T> result(n);
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            result[i]+=t_mat[i][j]*t_vec[j];
        }
    }
    return result;
}
template<typename T>
std::vector<T> operator* (const std::vector<std::vector<T>>& t_mat, const std::vector<T>& t_vec){
    return product(t_mat,t_vec);
}

//* matrix * matrix -> matrix
template<typename T>
std::vector<std::vector<T>> product(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2){
    if (t_mat1[0].size()!=t_mat2.size()){
        std::cout<<"In product between matrix and matrix, size is different"<<std::endl;
        exit(1);
    }
    const int n=t_mat1.size();
    const int m=t_mat2[0].size();
    const int k_max=t_mat2.size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));

    std::vector<std::vector<T>> tempMatrix=transpose(t_mat2);
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            T temp=0;
            for (int k=0; k<k_max; ++k){
                temp+=t_mat1[i][k]*tempMatrix[j][k];
            }
            result[i][j]=temp;
        }
    }
    return result;
}
template<typename T>
std::vector<std::vector<T>> operator* (const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2){
    return product(t_mat1, t_mat2);
}


//* Elementwise Calculation
template <typename T>
std::vector<T> elementProduct(const std::vector<T>& t_vec1, const std::vector<T>& t_vec2){
    if (t_vec1.size()!=t_vec2.size()){
        std::cout<<"In inner product, two vectors have different size"<<std::endl;
        exit(1);
    }
    const int n=t_vec1.size();
    std::vector<T> result(n);
    for (int i=0; i<n; ++i){
        result[i]=t_vec1[i]*t_vec2[i];
    }
    return result;
}

template<typename T>
std::vector<std::vector<T>> elementProduct(const std::vector<std::vector<T>>& t_mat1, const std::vector<std::vector<T>>& t_mat2){
    if (t_mat1.size()!=t_mat2.size() || t_mat1[0].size()!=t_mat2[0].size()){
        std::cout<<"In plus operation, two matrices have different size"<<std::endl;
        exit(1);
    }
    const int n=t_mat1.size();
    const int m=t_mat1[0].size();
    std::vector<std::vector<T>> result;
    result.resize(n,std::vector<T>(m));
    for (int i=0; i<n; i++){
        for (int j=0; j<m; ++j){
            result[i][j]=t_mat1[i][j]*t_mat2[i][j];
        }
    }
}

template<typename T>
std::vector<T> elementTanh(const std::vector<T>& t_vec){
    const int n=t_vec.size();
    std::vector<T> result(n);
    for (int i=0; i<n; ++i){
        result[i]=tanh(t_vec[i]);
    }
    return result;
}
template<typename T, typename TT>
std::vector<T> elementPow(const std::vector<T>& t_vec, const TT& c){
    const int n=t_vec.size();
    std::vector<T> result(n);
    for (int i=0; i<n; ++i){
        result[i]=pow(t_vec[i],c);
    }
    return result;
}

std::vector<double> linspace(const double& t_begin, const double& t_end, const int& t_size){
    std::vector<double> result(t_size);
    const double delta = (t_end-t_begin)/(t_size-1);
    for (int i=0; i<t_size; ++i){
        result[i] = delta*i + t_begin;
    }
    return result;
}

std::vector<double> arange(const double& t_begin, const double& t_end, const double& t_delta){
    int size = std::ceil((t_end-t_begin)/t_delta);
    std::vector<double> result(size);
    for (int i=0; i<size; ++i){
        result[i] = t_begin + t_delta*i;
    }
    result.emplace_back(t_end);
    return result;
}

}//* End of namespace linearAlgebra