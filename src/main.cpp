#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

#include <armadillo>
#include "solver.hpp"
#include "params.hpp"
#include "utils.hpp"

int main(int argc, char **argv)
{
    // read nanoparticle from file
    std::string filename = argv[1];
    std::ifstream file;
    file.open(filename, std::ios::in);

    size_t num_atoms;
    file >> num_atoms;
    arma::mat coordinate(3, num_atoms, arma::fill::zeros);
    double x, y, z;
    for (size_t i = 0; i < num_atoms; ++i)
    {
        file >> x >> y >> z;
        coordinate.col(i) = {x * wFQ::ang2au, y * wFQ::ang2au, z * wFQ::ang2au};
    }

    double strength;
    arma::vec3 component;
    file >> strength >> component[0] >> component[1] >> component[2];

    double scale;
    file >> scale;
    wFQ::sigma0 *= scale;
    wFQ::tau *= scale;

    double min_omega, max_omega;
    size_t num_points;
    file >> min_omega >> max_omega >> num_points;
    arma::vec freq = arma::linspace(min_omega * wFQ::ev2au, max_omega * wFQ::ev2au, num_points);

    double dt;
    size_t maxstep;
    file >> dt >> maxstep;

    file >> wFQ::t0 >> wFQ::sigma >> wFQ::w0;
    wFQ::t0 *= wFQ::fs2au;
    wFQ::sigma *= wFQ::fs2au;
    wFQ::w0 *= wFQ::ev2au;

    file.close();

    //
    size_t lastindex = filename.find_last_of(".");
    std::string dirname = filename.substr(0, lastindex);
    std::filesystem::create_directory(dirname);
    std::ofstream cross_section, polarizability;
    // std::ofstream charges;

    cross_section.open("./" + dirname + "/wFQ_cross_section.out", std::ios::out);
    cross_section << "# scale: " << std::to_string(scale) << std::endl;
    polarizability.open("./" + dirname + "/wFQ_polarizability.out", std::ios::out);
    cross_section << "# scale: " << std::to_string(scale) << std::endl;
    // charges.open("./" + dirname + "/wFQ_charge_dipole.out", std::ios::out);

    wFQ::freq_domain_solver fsolver(coordinate);

    for (auto omega : freq)
    {
        arma::cx_vec Q = fsolver.solve_charge(omega, strength, component);
        arma::cx_mat33 alpha = wFQ::polarizability(Q, coordinate, strength, component);
        arma::vec3 sigma = wFQ::cross_section(alpha, omega);

        // print cross_section
        cross_section << std::setw(8) << std::setprecision(4) << std::fixed << omega / wFQ::ev2au
                      << std::setw(15) << std::setprecision(4) << std::fixed << sigma[0] / wFQ::ang2au / wFQ::ang2au
                      << std::setw(15) << std::setprecision(4) << std::fixed << sigma[1] / wFQ::ang2au / wFQ::ang2au
                      << std::setw(15) << std::setprecision(4) << std::fixed << sigma[2] / wFQ::ang2au / wFQ::ang2au
                      << std::endl;

        // print polarizability
        polarizability << std::setw(8) << std::setprecision(4) << std::fixed << omega / wFQ::ev2au
                       << std::setw(15) << std::setprecision(4) << std::fixed << alpha(0, 0).real()
                       << std::setw(15) << std::setprecision(4) << std::fixed << alpha(0, 0).imag()
                       << std::setw(15) << std::setprecision(4) << std::fixed << alpha(1, 1).real()
                       << std::setw(15) << std::setprecision(4) << std::fixed << alpha(1, 1).imag()
                       << std::setw(15) << std::setprecision(4) << std::fixed << alpha(2, 2).real()
                       << std::setw(15) << std::setprecision(4) << std::fixed << alpha(2, 2).imag()
                       << std::endl;

        // wFQ::dump_charge_dipole(charges, Q, coordinate, std::to_string(omega / wFQ::ev2au));
    }

    std::ofstream induced_dipole;
    induced_dipole.open("./" + dirname + "/wFQ_time_propagation.out", std::ios::out);

    wFQ::time_domain_solver tsolver(coordinate);
    tsolver.propagate(induced_dipole, strength, component, dt, maxstep);
}
