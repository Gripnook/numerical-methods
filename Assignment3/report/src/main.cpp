#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <cmath>

#include "interpolation.h"
#include "solver.h"
#include "matrix.h"
#include "integral.h"

using namespace Numeric;

void question1();
void question2();
void question3();
void question4();

void printToFile(
    const std::vector<std::pair<double, double>>& data,
    const std::string& filename,
    const std::string& header);

int main()
{
    question1();
    question2();
    question3();
    question4();
    return 0;
}

void question1()
{
    std::cout << "========== Question 1 ==========" << std::endl;
    std::cout << std::endl;

    static const std::vector<std::pair<double, double>> data{{0.0, 0.0},
                                                             {0.2, 14.7},
                                                             {0.4, 36.5},
                                                             {0.6, 71.7},
                                                             {0.8, 121.4},
                                                             {1.0, 197.4},
                                                             {1.1, 256.2},
                                                             {1.2, 348.7},
                                                             {1.3, 540.6},
                                                             {1.4, 1062.8},
                                                             {1.5, 2318.0},
                                                             {1.6, 4781.9},
                                                             {1.7, 8687.4},
                                                             {1.8, 13924.3},
                                                             {1.9, 22650.2}};
    printToFile(data, "q1data.csv", "B (T),H (A/m)");

    auto f = lagrange({data.begin(), data.begin() + 6});
    auto interpol = sample(f, 0.0, 1.01, 0.01);
    printToFile(interpol, "q1a.csv", "B (T),H (A/m)");
    std::cout << "Question 1 (a) data successfully printed to q1a.csv"
              << std::endl;
    std::cout << std::endl;

    f = lagrange({data[0], data[8], data[9], data[12], data[13], data[14]});
    interpol = sample(f, 0.0, 1.91, 0.01);
    printToFile(interpol, "q1b.csv", "B (T),H (A/m)");
    std::cout << "Question 1 (b) data successfully printed to q1b.csv"
              << std::endl;
    std::cout << std::endl;

    f = hermite(
        {data[0], data[8], data[9], data[12], data[13], data[14]},
        {slope(data[0], data[1]),
         slope(data[7], data[9]),
         slope(data[8], data[10]),
         slope(data[11], data[13]),
         slope(data[12], data[14]),
         slope(data[13], data[14])});
    interpol = sample(f, 0.0, 1.91, 0.01);
    printToFile(interpol, "q1c.csv", "B (T),H (A/m)");
    std::cout << "Question 1 (c) data successfully printed to q1c.csv"
              << std::endl;
    std::cout << std::endl;
}

void question2()
{
    std::cout << "========== Question 2 ==========" << std::endl;
    std::cout << std::endl;

    static const std::vector<std::pair<double, double>> data{{0.0, 0.0},
                                                             {0.2, 14.7},
                                                             {0.4, 36.5},
                                                             {0.6, 71.7},
                                                             {0.8, 121.4},
                                                             {1.0, 197.4},
                                                             {1.1, 256.2},
                                                             {1.2, 348.7},
                                                             {1.3, 540.6},
                                                             {1.4, 1062.8},
                                                             {1.5, 2318.0},
                                                             {1.6, 4781.9},
                                                             {1.7, 8687.4},
                                                             {1.8, 13924.3},
                                                             {1.9, 22650.2}};
    static const double N = 1000;
    static const double I = 8;
    static const double S = 1e-4;
    static const double Lc = 30e-2;
    static const double Lg = 0.5e-2;
    static const double mu0 = 4 * M_PI * 1e-7;

    auto H = pwl(data);
    auto Hprime = pwlprime(data);
    auto B = pwlinverse(data);

    std::function<double(double)> f = [&](double flux) {
        return flux * Lg / (mu0 * S) + H(flux / S) * Lc - N * I;
    };
    std::function<double(double)> fprime = [&](double flux) {
        return Lg / (mu0 * S) + Lc * Hprime(flux / S) / S;
    };

    double flux;
    int iterations;
    std::tie(flux, iterations) = newtonRaphson(f, fprime, 0.0);
    std::cout << "Newton-Raphson" << std::endl;
    std::cout << "Flux: " << flux << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << std::endl;

    f = [&](double flux) {
        return S * B((N * I - Lg * flux / (mu0 * S)) / Lc);
    };

    std::tie(flux, iterations) = successiveSubstitution(f, 0.0);
    std::cout << "Successive Substitution" << std::endl;
    std::cout << "Flux: " << flux << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << std::endl;
}

void question3()
{
    std::cout << "========== Question 3 ==========" << std::endl;
    std::cout << std::endl;

    static const double E = 0.22;
    static const double R = 500;
    static const double vt = 25e-3;
    static const double IsA = 0.6e-6;
    static const double IsB = 1.2e-6;

    auto f = [&](const Matrix<double>& v) {
        Matrix<double> result{2, 1};
        result(0) = (v(0) - E) / R + IsA * (std::exp((v(0) - v(1)) / vt) - 1);
        result(1) = -IsA * (std::exp((v(0) - v(1)) / vt) - 1) +
            IsB * (std::exp(v(1) / vt) - 1);
        return result;
    };
    auto jacobian = [&](const Matrix<double>& v) {
        Matrix<double> result{2, 2};
        result(0, 0) = 1 / R + IsA / vt * std::exp((v(0) - v(1)) / vt);
        result(0, 1) = -IsA / vt * std::exp((v(0) - v(1)) / vt);
        result(1, 0) = -IsA / vt * std::exp((v(0) - v(1)) / vt);
        result(1, 1) = IsA / vt * std::exp((v(0) - v(1)) / vt) +
            IsB / vt * std::exp(v(1) / vt);
        return result;
    };

    Matrix<double> v;
    int iterations;
    std::cout << "Newton-Raphson" << std::endl;
    std::cout << "iteration,v1,v2,f1,f2,error" << std::endl;
    std::tie(v, iterations) = newtonRaphson(
        f, jacobian, "[0.0; 0.0]", [&](auto v, auto error, auto iteration) {
            auto residual = f(v);
            std::cout << iteration << "," << v(0) << "," << v(1) << ","
                      << residual(0) << "," << residual(1) << "," << error
                      << std::endl;
        });
    std::cout << std::endl;
}

void question4()
{
    std::cout << "========== Question 4 ==========" << std::endl;
    std::cout << std::endl;

    std::cout << "Integrating cos(x)" << std::endl;
    std::cout << "segments,error" << std::endl;
    for (int n = 1; n <= 20; ++n)
    {
        double result =
            integral([](double x) { return std::cos(x); }, 0.0, 1.0, n);
        double error = std::abs(std::sin(1.0) - result);
        std::cout << n << "," << error << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Integrating log(x)" << std::endl;
    std::cout << "segments,error" << std::endl;
    for (int n = 10; n <= 200; n += 10)
    {
        double result =
            integral([](double x) { return std::log(x); }, 0.0, 1.0, n);
        double error = std::abs(-1.0 - result);
        std::cout << n << "," << error << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Integrating log(x) with uneven segments" << std::endl;
    std::cout << "r,alpha,error" << std::endl;
    for (double r = 1.0; r < 2.0; r += 0.01)
    {
        double alpha =
            r == 1.0 ? 1 / 10.0 : (r - 1.0) / (std::pow(r, 10) - 1.0);
        std::vector<double> segments;
        double totalLength = 0.0;
        for (int i = 0; i < 9; ++i)
        {
            double segment = alpha * std::pow(r, i);
            segments.push_back(segment);
            totalLength += segment;
        }
        double segment = 1.0 - totalLength;
        segments.push_back(segment);
        double result =
            integral([](double x) { return std::log(x); }, 0.0, segments);
        double error = std::abs(-1.0 - result);
        std::cout << r << "," << alpha << "," << error << std::endl;
    }
    std::cout << std::endl;
}

void printToFile(
    const std::vector<std::pair<double, double>>& data,
    const std::string& filename,
    const std::string& header)
{
    std::ofstream out{filename};
    out << header << std::endl;
    for (auto p : data)
        out << p.first << "," << p.second << std::endl;
}
