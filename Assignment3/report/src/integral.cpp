#include "integral.h"

namespace Numeric {

double integral(std::function<double(double)> f, double a, double b, int n)
{
    double result = 0.0;
    double dx = (b - a) / n;
    for (int i = 0; i < n; ++i)
    {
        double x0 = a + i * dx;
        double x1 = a + (i + 1) * dx;
        result += (x1 - x0) * f((x0 + x1) / 2);
    }
    return result;
}

double integral(
    std::function<double(double)> f,
    double a,
    const std::vector<double>& segments)
{
    double result = 0.0;
    double x0 = a;
    for (auto segment : segments)
    {
        double x1 = x0 + segment;
        result += (x1 - x0) * f((x0 + x1) / 2);
        x0 = x1;
    }
    return result;
}

} // namespace Numeric
