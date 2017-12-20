#include "circuit-solver.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

#include "solver.h"

using namespace Numeric;

namespace Circuits {

Numeric::Matrix<double> csolve(std::istream& in)
{
    Matrix<double> A, J, Y, E;
    std::tie(A, J, Y, E) = parse(in);
    auto m = A * Y * transpose(A);
    auto b = A * (J + Y * E);
    return solve(m, b);
}

Numeric::Matrix<double> cbsolve(std::istream& in)
{
    Matrix<double> A, J, Y, E;
    std::tie(A, J, Y, E) = parse(in);
    auto m = A * Y * transpose(A);
    auto b = A * (J + Y * E);
    return bsolve(m, b);
}

std::tuple<int, int, double, double, double> parse(const std::string& line);

std::tuple<Matrix<double>, Matrix<double>, Matrix<double>, Matrix<double>>
    parse(std::istream& in)
{
    std::vector<std::tuple<int, int, double, double, double>> lines;
    std::set<int> nodes;
    for (std::string line; std::getline(in, line);)
    {
        if (line.empty() || line[0] == '#')
            continue;
        auto result = parse(line);
        lines.push_back(result);
        nodes.emplace(std::get<0>(result));
        nodes.emplace(std::get<1>(result));
    }

    if (nodes.count(0) == 0)
        throw std::runtime_error{"parse: missing ground node"};

    // Perform node mapping.
    std::unordered_map<int, int> nodeMapping;
    int index = 0;
    for (const auto& node : nodes)
    {
        if (node == 0)
            continue; // Ignore ground node.
        nodeMapping.emplace(node, index++);
    }

    int nodeCount = nodeMapping.size();
    int branchCount = lines.size();

    Matrix<double> A{nodeCount, branchCount};
    Matrix<double> J{branchCount, 1};
    Matrix<double> Y{branchCount, branchCount};
    Matrix<double> E{branchCount, 1};

    for (int i = 0; i < branchCount; ++i)
    {
        auto line = lines[i];
        int N1 = std::get<0>(line);
        if (N1 != 0)
            A(nodeMapping.find(N1)->second, i) = 1;
        int N2 = std::get<1>(line);
        if (N2 != 0)
            A(nodeMapping.find(N2)->second, i) = -1;
        J(i) = std::get<2>(line);
        Y(i, i) = std::get<3>(line);
        E(i) = std::get<4>(line);
    }

    return std::make_tuple(A, J, Y, E);
}

std::tuple<int, int, double, double, double> parse(const std::string& line)
{
    std::stringstream in{line};
    int N1, N2;
    double Jk, Rk, Ek;
    in >> N1;
    if (!in)
        throw std::runtime_error{"parse: invalid node N+"};
    in >> N2;
    if (!in)
        throw std::runtime_error{"parse: invalid node N-"};
    in >> Jk;
    if (!in)
        throw std::runtime_error{"parse: invalid branch current Jk"};
    in >> Rk;
    if (!in)
        throw std::runtime_error{"parse: invalid branch resistance Rk"};
    if (Rk == 0)
        throw std::runtime_error{"parse: branch resistance Rk cannot be zero"};
    in >> Ek;
    if (!in)
        throw std::runtime_error{"parse: invalid branch voltage Ek"};
    return std::make_tuple(N1, N2, Jk, 1.0 / Rk, Ek);
}
}
