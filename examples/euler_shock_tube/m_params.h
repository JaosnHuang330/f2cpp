#pragma once

#include <string>

namespace m_params {

// Grid
constexpr int nx = 200;
constexpr double xmin = 0.0;
constexpr double xmax = 1.0;
constexpr double dx = (xmax - xmin) / static_cast<double>(nx);

// Time stepping
constexpr double cfl = 0.8;
constexpr double t_final = 0.2;
constexpr int max_steps = 10000;

// Gas constant
constexpr double gamma_gas = 1.4;

// IO
constexpr int output_freq = 50;
const std::string output_file = "euler_output.dat";
const std::string restart_file = "restart.bin";

} // namespace m_params
