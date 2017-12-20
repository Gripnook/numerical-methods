#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>

namespace Numeric {

template <typename T>
class Matrix;

// Matrix addition.
template <typename T>
Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs);
template <typename T>
Matrix<T>& operator+=(Matrix<T>& lhs, const Matrix<T>& rhs);

// Matrix subtraction.
template <typename T>
Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs);
template <typename T>
Matrix<T>& operator-=(Matrix<T>& lhs, const Matrix<T>& rhs);

// Matrix multiplication.
template <typename T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <typename T>
Matrix<T>& operator*=(Matrix<T>& lhs, const Matrix<T>& rhs);

// Scalar multiplication.
template <typename T>
Matrix<T> operator*(Matrix<T> lhs, const T& rhs);
template <typename T>
Matrix<T> operator*(const T& lhs, Matrix<T> rhs);
template <typename T>
Matrix<T>& operator*=(Matrix<T>& lhs, const T& rhs);

// Scalar division.
template <typename T>
Matrix<T> operator/(Matrix<T> lhs, const T& rhs);
template <typename T>
Matrix<T>& operator/=(Matrix<T>& lhs, const T& rhs);

// Matrix transpose.
template <typename T>
Matrix<T> transpose(const Matrix<T>& m);

// Computes the lower and upper bandwidths of matrix m.
template <typename T>
std::pair<int, int> bandwidth(const Matrix<T>& m);

// Matrix I/O.
template <typename T>
std::istream& operator>>(std::istream& in, Matrix<T>& m);
template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& m);

template <typename T>
class Matrix
{
public:
    // Constructs a matrix of size 0.
    Matrix() = default;

    // Constructs a matrix of size [rows x cols] filled with zeros.
    Matrix(int rows, int cols)
        : m_rows{rows}, m_cols{cols}, m_elements(rows * cols)
    {
    }

    // Implicit conversion from a string.
    Matrix(const std::string& str)
    {
        std::stringstream ss{str};
        ss >> *this;
    }

    // Implicit conversion from a c string.
    Matrix(const char* str)
    {
        std::stringstream ss{str};
        ss >> *this;
    }

    // Implicit conversion to the element type for a single element matrix.
    operator T() const
    {
        if (size() != 1)
            throw std::runtime_error{"cannot convert matrix to single element"};
        return (*this)(0);
    }

    // Indexes the matrix in a row-major order.
    T& operator()(int index)
    {
        return m_elements.at(index);
    }

    // Indexes the matrix in a row-major order.
    const T& operator()(int index) const
    {
        return m_elements.at(index);
    }

    T& operator()(int row, int col)
    {
        return (*this)(row * m_cols + col);
    }

    const T& operator()(int row, int col) const
    {
        return (*this)(row * m_cols + col);
    }

    // Returns the total number of elements in the matrix.
    int size() const
    {
        return m_rows * m_cols;
    }

    int rows() const
    {
        return m_rows;
    }

    int cols() const
    {
        return m_cols;
    }

    bool operator==(const Matrix<T>& other) const
    {
        if (size() == 0 && other.size() == 0)
            return true;
        return m_rows == other.m_rows && m_cols == other.m_cols &&
            m_elements == other.m_elements;
    }

    bool operator!=(const Matrix<T>& other) const
    {
        return !(*this == other);
    }

    friend std::istream& operator>><>(std::istream& in, Matrix<T>& m);
    friend std::ostream& operator<<<>(std::ostream& out, const Matrix<T>& m);

private:
    int m_rows{};
    int m_cols{};
    std::vector<T> m_elements;

    // Constructs a matrix by interpreting the elements
    // as a matrix of size [rows x cols].
    Matrix(int rows, int cols, std::vector<T> elements)
        : m_rows{rows}, m_cols{cols}, m_elements{elements}
    {
        if (rows * cols != elements.size())
            throw std::runtime_error{"inconsistent matrix dimensions"};
    }
};

template <typename T>
Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs)
{
    return lhs += rhs;
}

template <typename T>
Matrix<T>& operator+=(Matrix<T>& lhs, const Matrix<T>& rhs)
{
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
        throw std::runtime_error{"add: inconsistent matrix dimensions: [" +
                                 std::to_string(lhs.rows()) + "x" +
                                 std::to_string(lhs.cols()) + "], [" +
                                 std::to_string(rhs.rows()) + "x" +
                                 std::to_string(rhs.cols()) + "]"};

    int rows = lhs.rows(), cols = lhs.cols();
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            lhs(i, j) += rhs(i, j);
        }
    }
    return lhs;
}

template <typename T>
Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs)
{
    return lhs -= rhs;
}

template <typename T>
Matrix<T>& operator-=(Matrix<T>& lhs, const Matrix<T>& rhs)
{
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
        throw std::runtime_error{"subtract: inconsistent matrix dimensions: [" +
                                 std::to_string(lhs.rows()) + "x" +
                                 std::to_string(lhs.cols()) + "], [" +
                                 std::to_string(rhs.rows()) + "x" +
                                 std::to_string(rhs.cols()) + "]"};

    int rows = lhs.rows(), cols = lhs.cols();
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            lhs(i, j) -= rhs(i, j);
        }
    }
    return lhs;
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
    if (lhs.cols() != rhs.rows())
        throw std::runtime_error{"multiply: inconsistent matrix dimensions: [" +
                                 std::to_string(lhs.rows()) + "x" +
                                 std::to_string(lhs.cols()) + "], [" +
                                 std::to_string(rhs.rows()) + "x" +
                                 std::to_string(rhs.cols()) + "]"};

    int rows = lhs.rows(), cols = rhs.cols(), innerSize = lhs.cols();
    Matrix<T> result{rows, cols};
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            for (int k = 0; k < innerSize; ++k)
            {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}

template <typename T>
Matrix<T>& operator*=(Matrix<T>& lhs, const Matrix<T>& rhs)
{
    lhs = lhs * rhs;
    return lhs;
}

template <typename T>
Matrix<T> operator*(Matrix<T> lhs, const T& rhs)
{
    return lhs *= rhs;
}

template <typename T>
Matrix<T> operator*(const T& lhs, Matrix<T> rhs)
{
    return rhs *= lhs;
}

template <typename T>
Matrix<T>& operator*=(Matrix<T>& lhs, const T& rhs)
{
    for (int i = 0; i < lhs.rows(); ++i)
    {
        for (int j = 0; j < lhs.cols(); ++j)
        {
            lhs(i, j) *= rhs;
        }
    }
    return lhs;
}

template <typename T>
Matrix<T> operator/(Matrix<T> lhs, const T& rhs)
{
    return lhs /= rhs;
}

template <typename T>
Matrix<T>& operator/=(Matrix<T>& lhs, const T& rhs)
{
    if (rhs == 0)
        throw std::runtime_error{"divide: division by zero"};

    for (int i = 0; i < lhs.rows(); ++i)
    {
        for (int j = 0; j < lhs.cols(); ++j)
        {
            lhs(i, j) /= rhs;
        }
    }
    return lhs;
}

template <typename T>
Matrix<T> transpose(const Matrix<T>& m)
{
    Matrix<T> result{m.cols(), m.rows()};
    for (int i = 0; i < m.rows(); ++i)
    {
        for (int j = 0; j < m.cols(); ++j)
        {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

template <typename T>
std::pair<int, int> bandwidth(const Matrix<T>& m)
{
    int lower = 0;
    int upper = 0;
    for (int i = 0; i < m.rows(); ++i)
    {
        for (int j = 0; j < m.cols(); ++j)
        {
            if (m(i, j) != 0)
            {
                auto diag = j - i;
                if (diag == 0)
                {
                    lower = std::max(lower, diag + 1);
                    upper = std::max(upper, diag + 1);
                }
                else if (diag > 0)
                {
                    upper = std::max(upper, diag + 1);
                }
                else
                {
                    lower = std::max(lower, -diag + 1);
                }
            }
        }
    }
    return {lower, upper};
}

template <typename T>
std::istream& operator>>(std::istream& in, Matrix<T>& m)
{
    int rows = 1;
    int size = 0;
    std::vector<T> elements;

    char ch;
    in >> ch;
    if (!in || ch != '[')
        throw std::runtime_error{"input: missing '['"};

    while (true)
    {
        T elem;
        in >> elem;
        if (!in)
            throw std::runtime_error{"input: invalid element"};
        elements.push_back(elem);
        ++size;

        in >> ch;
        if (!in)
            throw std::runtime_error{"input: missing ']'"};

        if (ch == ',')
        {
            // Go to next element.
        }
        else if (ch == ';')
        {
            ++rows;
        }
        else if (ch == ']')
        {
            break;
        }
        else
        {
            throw std::runtime_error{"input: invalid separator"};
        }
    }

    int cols = size / rows;
    if (rows * cols != size)
        throw std::runtime_error{"input: inconsistent matrix dimensions"};

    m = Matrix<T>{rows, cols, std::move(elements)};
    return in;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& m)
{
    out << "[" << std::endl;
    for (int i = 0; i < m.rows(); ++i)
    {
        out << "    ";
        for (int j = 0; j < m.cols(); ++j)
        {
            out << m(i, j);
            if (j != m.cols() - 1)
                out << ", ";
        }
        if (i != m.rows() - 1)
            out << ";";
        out << std::endl;
    }
    out << "]" << std::endl;
    return out;
}
}
