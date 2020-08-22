# pragma once
# include <iomanip>
# include <string>

//* Change to_string with precision for main folder name
template <typename T>
std::string to_stringWithPrecision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::fixed <<std::setprecision(n) << a_value;
    return out.str();
}

template <typename T>
std::string to_stringWithExponent(const T &a_value, const int &n=1)
{
    std::ostringstream out;
    out<<std::scientific<<std::setprecision(n) << a_value;
    return out.str();
}
