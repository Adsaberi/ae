#include <vector>
#include <cmath>
#include <iostream>
const int width = 100;
const int height = 100;
const double viscosity = 0.1;
const double dt = 0.01;
const double dx = 1.0;
const double dy = 1.0;
std::vector<std::vector<double>> u(width, std::vector<double>(height, 0.0));
std::vector<std::vector<double>> v(width, std::vector<double>(height, 0.0));
std::vector<std::vector<double>> p(width, std::vector<double>(height, 0.0));
void applyBoundaryConditions() 
void solvePressurePoisson(std::vector<std::vector<double>>& p)
void updateVelocities(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, const std::vector<std::vector<double>>& p) {
void simulateFluid() {
    for (int t = 0; t < 1000; ++t) {
        applyBoundaryConditions();
        solvePressurePoisson(p);
        updateVelocities(u, v, p);
        std::cout << "Step " << t << " completed." << std::endl;
    }
}

int main() {
    simulateFluid();
    return 0;
}
