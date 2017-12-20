#pragma once

#include <vector>
#include <iosfwd>

namespace Numeric {

struct Polynomial
{
    std::vector<double> coefficients;

    Polynomial() = default;

    // Contructs a polynomial from the given coefficients in ascending order.
    Polynomial(std::vector<double> coefficients) : coefficients{coefficients}
    {
    }

    // Implicit conversion from double.
    Polynomial(double a0) : coefficients{{a0}}
    {
    }

    // Evaluates the polynomial at the given value of x.
    double operator()(double x) const;

    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;
};

// Polynomial addition.
Polynomial operator+(Polynomial lhs, const Polynomial& rhs);
Polynomial& operator+=(Polynomial& lhs, const Polynomial& rhs);

// Polynomial subtraction.
Polynomial operator-(Polynomial lhs, const Polynomial& rhs);
Polynomial& operator-=(Polynomial& lhs, const Polynomial& rhs);

// Polynomial multiplication.
Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs);
Polynomial& operator*=(Polynomial& lhs, const Polynomial& rhs);

// Scalar division.
Polynomial operator/(Polynomial lhs, double rhs);
Polynomial& operator/=(Polynomial& lhs, double rhs);

// Stream I/O.
std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial);

} // namespace Numeric
