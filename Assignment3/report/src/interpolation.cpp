#include "interpolation.h"

#include "polynomial.h"

namespace Numeric {

std::function<double(double)>
    lagrange(const std::vector<std::pair<double, double>>& data)
{
    int n = static_cast<int>(data.size());
    Polynomial y{std::vector<double>(n)};

    for (int i = 0; i < n; ++i)
    {
        Polynomial lagrangePolynomial{{1.0}};
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
                continue;
            lagrangePolynomial *= Polynomial{{-data[j].first, 1.0}};
            lagrangePolynomial /= data[i].first - data[j].first;
        }
        y += data[i].second * lagrangePolynomial;
    }

    return y;
}

std::function<double(double)> hermite(
    const std::vector<std::pair<double, double>>& data,
    const std::vector<double>& derivatives)
{
    int n = static_cast<int>(data.size());
    std::vector<std::pair<double, Polynomial>> y;

    for (int i = 0; i < n - 1; ++i)
    {
        Polynomial interval{std::vector<double>(2 * n - 1)};

        double diff = data[i + 1].first - data[i].first;
        Polynomial F1{{-data[i].first, 1.0}};
        Polynomial F2{{-data[i + 1].first, 1.0}};

        auto L1 = F2 / -diff;
        auto L2 = F1 / diff;
        auto L1prime = 1 / -diff;
        auto L2prime = 1 / diff;

        auto U1 = (1.0 - 2.0 * L1prime * F1) * L1 * L1;
        auto U2 = (1.0 - 2.0 * L2prime * F2) * L2 * L2;
        auto V1 = F1 * L1 * L1;
        auto V2 = F2 * L2 * L2;

        interval = data[i].second * U1 + data[i + 1].second * U2 +
            derivatives[i] * V1 + derivatives[i + 1] * V2;
        y.emplace_back(data[i + 1].first, interval);
    }

    return [=](double x) {
        for (const auto& interval : y)
            if (x < interval.first)
                return interval.second(x);
        return y.back().second(x);
    };
}

std::function<double(double)>
    pwl(const std::vector<std::pair<double, double>>& data)
{
    return [=](double x) {
        int index = static_cast<int>(data.size()) - 1;
        for (size_t i = 1; i < data.size(); ++i)
        {
            if (x < data[i].first)
            {
                index = i;
                break;
            }
        }
        return slope(data[index - 1], data[index]) * (x - data[index].first) +
            data[index].second;
    };
}

std::function<double(double)>
    pwlprime(const std::vector<std::pair<double, double>>& data)
{
    return [=](double x) {
        int index = static_cast<int>(data.size()) - 1;
        for (size_t i = 1; i < data.size(); ++i)
        {
            if (x < data[i].first)
            {
                index = i;
                break;
            }
        }
        return slope(data[index - 1], data[index]);
    };
}

std::function<double(double)>
    pwlinverse(const std::vector<std::pair<double, double>>& data)
{
    return [=](double y) {
        int index = static_cast<int>(data.size()) - 1;
        for (size_t i = 1; i < data.size(); ++i)
        {
            if (y < data[i].second)
            {
                index = i;
                break;
            }
        }
        return 1 / slope(data[index - 1], data[index]) *
            (y - data[index].second) +
            data[index].first;
    };
}

double slope(std::pair<double, double> p1, std::pair<double, double> p2)
{
    return (p2.second - p1.second) / (p2.first - p1.first);
}

std::vector<std::pair<double, double>>
    sample(std::function<double(double)> f, double x0, double x1, double dx)
{
    std::vector<std::pair<double, double>> data;
    for (double x = x0; x <= x1; x += dx)
        data.emplace_back(x, f(x));
    return data;
}

} // namespace Numeric
