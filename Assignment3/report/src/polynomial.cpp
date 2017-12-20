#include "polynomial.h"

#include <iostream>

namespace Numeric {

double Polynomial::operator()(double x) const
{
    double result = 0.0;
    double power = 1.0;
    for (auto coefficient : coefficients)
    {
        result += coefficient * power;
        power *= x;
    }
    return result;
}

bool Polynomial::operator==(const Polynomial& other) const
{
    return coefficients == other.coefficients;
}

bool Polynomial::operator!=(const Polynomial& other) const
{
    return !(*this == other);
}

Polynomial operator+(Polynomial lhs, const Polynomial& rhs)
{
    lhs += rhs;
    return lhs;
}

Polynomial& operator+=(Polynomial& lhs, const Polynomial& rhs)
{
    if (rhs.coefficients.size() > lhs.coefficients.size())
        lhs.coefficients.resize(rhs.coefficients.size());
    for (size_t i = 0; i < rhs.coefficients.size(); ++i)
        lhs.coefficients[i] += rhs.coefficients[i];
    return lhs;
}

Polynomial operator-(Polynomial lhs, const Polynomial& rhs)
{
    lhs -= rhs;
    return lhs;
}

Polynomial& operator-=(Polynomial& lhs, const Polynomial& rhs)
{
    if (rhs.coefficients.size() > lhs.coefficients.size())
        lhs.coefficients.resize(rhs.coefficients.size());
    for (size_t i = 0; i < rhs.coefficients.size(); ++i)
        lhs.coefficients[i] -= rhs.coefficients[i];
    return lhs;
}

Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs)
{
    std::vector<double> coefficients(
        lhs.coefficients.size() + rhs.coefficients.size() - 1);
    for (size_t i = 0; i < lhs.coefficients.size(); ++i)
        for (size_t j = 0; j < rhs.coefficients.size(); ++j)
            coefficients[i + j] += lhs.coefficients[i] * rhs.coefficients[j];
    return Polynomial{coefficients};
}

Polynomial& operator*=(Polynomial& lhs, const Polynomial& rhs)
{
    lhs = lhs * rhs;
    return lhs;
}

Polynomial operator/(Polynomial lhs, double rhs)
{
    lhs /= rhs;
    return lhs;
}

Polynomial& operator/=(Polynomial& lhs, double rhs)
{
    for (auto& coefficient : lhs.coefficients)
        coefficient /= rhs;
    return lhs;
}

std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial)
{
    int n = static_cast<int>(polynomial.coefficients.size());
    for (int i = n - 1; i > 0; --i)
        out << polynomial.coefficients[i] << " x^" << i << " + ";
    out << polynomial.coefficients[0];
    return out;
}

} // namespace Numeric
