#include "solver.h"

#include <algorithm>
#include <cmath>

#include "matrix-solver.h"

namespace Numeric {

std::pair<double, int> newtonRaphson(
    std::function<double(double)> f,
    std::function<double(double)> fprime,
    double x0)
{
    auto error = [&](auto x) { return std::fabs(f(x) / f(x0)); };

    int iterations = 0;
    double x = x0;
    while (error(x) >= 1e-6)
    {
        ++iterations;
        x -= f(x) / fprime(x);
    }
    return {x, iterations};
}

std::pair<double, int>
    successiveSubstitution(std::function<double(double)> f, double x0)
{
    auto error = [&](auto x) { return std::fabs((f(x) - x) / (f(x0) - x0)); };

    int iterations = 0;
    double x = x0;
    while (error(x) >= 1e-6)
    {
        ++iterations;
        x = f(x);
    }
    return {x, iterations};
}

std::pair<Matrix<double>, int> newtonRaphson(
    std::function<Matrix<double>(const Matrix<double>&)> f,
    std::function<Matrix<double>(const Matrix<double>&)> jacobian,
    const Matrix<double>& x0,
    std::function<void(const Matrix<double>&, double, int)> callback)
{
    auto error = [&](const auto& x) {
        auto func = f(x);
        auto numerator = std::accumulate(
            func.begin(), func.end(), 0.0, [](auto lhs, auto rhs) {
                return lhs + std::abs(rhs);
            });
        func = f(x0);
        auto denominator = std::accumulate(
            func.begin(), func.end(), 0.0, [](auto lhs, auto rhs) {
                return lhs + std::abs(rhs);
            });
        return numerator / denominator;
    };

    int iterations = 0;
    Matrix<double> x = x0;
    while (error(x) >= 1e-6)
    {
        callback(x, error(x), iterations);
        ++iterations;
        x = solve(jacobian(x), jacobian(x) * x - f(x));
    }
    callback(x, error(x), iterations);
    return {x, iterations};
}

} // namespace Numeric
