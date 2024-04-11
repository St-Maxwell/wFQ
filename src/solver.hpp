#pragma once
#include <iostream>
#include <iomanip>
#include <armadillo>
#include "params.hpp"
#include "utils.hpp"

namespace wFQ
{
    inline double damping(double lmn)
    {
        return 1 / (1 + exp(-d * (lmn / (s * lij0) - 1)));
    }

    inline std::complex<double> drude_z(double omega)
    {
        constexpr std::complex<double> I{0, 1};
        return -omega * (I + omega * tau) / (2 * sigma0);
    }

    /**
     * Kbar_mn = (1 - f(l_mn)) * A_mn / l_mn
     * \param coordinate unit: au
     */
    inline arma::mat build_Kbar(const arma::mat &coordinate)
    {
        const size_t num_atoms = coordinate.n_cols;
        arma::mat Kbar(num_atoms, num_atoms, arma::fill::zeros);

        // loop lower triangular matrix
        for (size_t j = 0; j < num_atoms - 1; ++j) // loop columns
        {
            for (size_t i = j + 1; i < num_atoms; ++i) // loop rows
            {
                double lij = arma::norm(coordinate.col(i) - coordinate.col(j));
                Kbar(i, j) = (1 - damping(lij)) * Aij / lij;
            }
        }
        // symetrize
        Kbar += Kbar.t();

        return Kbar;
    }

    /**
     * build charge-charge interaction kernal
     * \param coordinate unit: au
     */
    inline arma::mat build_T(const arma::mat &coordinate)
    {
        const size_t num_atoms = coordinate.n_cols;
        arma::mat T(num_atoms, num_atoms, arma::fill::zeros);
        const double Rij = 2 / sqrt(arma::datum::pi) / eta;

        // loop lower triangular matrix
        for (size_t j = 0; j < num_atoms - 1; ++j) // loop columns
        {
            for (size_t i = j + 1; i < num_atoms; ++i) // loop rows
            {
                double lij = arma::norm(coordinate.col(i) - coordinate.col(j));
                T(i, j) = erf(lij / Rij) / lij;
            }
        }
        // symetrize
        T.diag() += eta / 2;
        T += T.t();

        return T;
    }

    inline arma::cx_mat33 polarizability(const arma::cx_vec &charge, const arma::mat &coordinate, double strength, arma::vec3 component)
    {
        const size_t num_atoms = coordinate.n_cols;

        std::complex<double> mu_x{0, 0};
        std::complex<double> mu_y{0, 0};
        std::complex<double> mu_z{0, 0};

        for (size_t i = 0; i < num_atoms; ++i)
        {
            mu_x += charge(i) * coordinate(0, i);
            mu_y += charge(i) * coordinate(1, i);
            mu_z += charge(i) * coordinate(2, i);
        }

        component /= arma::norm(component);
        double Ex = strength * component[0];
        double Ey = strength * component[1];
        double Ez = strength * component[2];

        arma::cx_mat33 alpha = {
            {mu_x / Ex, mu_x / Ey, mu_x / Ez},
            {mu_y / Ex, mu_y / Ey, mu_y / Ez},
            {mu_z / Ex, mu_z / Ey, mu_z / Ez},
        };

        return alpha;
    }

    inline arma::vec3 cross_section(arma::cx_mat33 alpha, double omega)
    {
        double sigma_x = 4 * arma::datum::pi * omega * alpha(0, 0).imag() / 3 / wFQ::c2au;
        double sigma_y = 4 * arma::datum::pi * omega * alpha(1, 1).imag() / 3 / wFQ::c2au;
        double sigma_z = 4 * arma::datum::pi * omega * alpha(2, 2).imag() / 3 / wFQ::c2au;
        return {sigma_x, sigma_y, sigma_z};
    }

    inline double gaussian_envelope(double t)
    {
        return exp(-(t - t0) * (t - t0) / (2 * sigma * sigma)) * cos(w0 * (t - t0));
    }

    class freq_domain_solver
    {
    private:
        arma::mat coordinate;
        arma::mat Kbar;
        arma::mat A;

    public:
        explicit freq_domain_solver(const arma::mat &cluster_coordinate)
            : coordinate{cluster_coordinate}
        {
            Kbar = build_Kbar(coordinate);
            arma::mat T = build_T(coordinate);

            arma::mat P = arma::diagmat(arma::sum(Kbar, 1));

            A = (Kbar - P) * T;
        }

        freq_domain_solver(const freq_domain_solver &rhs) = delete;
        ~freq_domain_solver() = default;

        arma::cx_mat solve_charge(double omega, double strength, arma::vec3 component)
        {
            const size_t num_atoms = coordinate.n_cols;

            component /= arma::norm(component);

            arma::vec Vext(num_atoms, arma::fill::zeros);
            const double factor = -eta * eta * strength / 2;

            for (size_t i = 0; i < num_atoms; ++i)
            {
                Vext(i) = factor * arma::dot(component, coordinate.col(i));
            }

            arma::cx_vec R(num_atoms, arma::fill::zeros);
            for (size_t i = 0; i < num_atoms; ++i)
            {
                for (size_t j = 0; j < num_atoms; ++j)
                {
                    R(i) += Kbar(j, i) * (Vext(i) - Vext(j));
                }
            }

            arma::cx_mat LHS = A - arma::eye(arma::size(A)) * drude_z(omega);

            return arma::solve(LHS, R);
        }
    };

    class time_domain_solver
    {
    private:
        arma::mat coordinate;
        arma::mat Kbar;
        arma::mat A;

    public:
        explicit time_domain_solver(const arma::mat &cluster_coordinate)
            : coordinate{cluster_coordinate}
        {
            Kbar = std::move(build_Kbar(coordinate));
            arma::mat T = std::move(build_T(coordinate));

            arma::mat P = arma::diagmat(arma::sum(Kbar, 1));

            A = (Kbar - P) * T;
        }

        time_domain_solver(const time_domain_solver &rhs) = delete;
        ~time_domain_solver() = default;

        void propagate(std::ostream &os, double strength, arma::vec3 component, const double dt, const size_t maxstep)
        {
            const size_t num_atoms = coordinate.n_cols;
            arma::vec R(num_atoms, arma::fill::zeros);
            arma::vec Q(num_atoms, arma::fill::zeros);
            arma::vec J(num_atoms, arma::fill::zeros);

            const double fac1 = (2 * tau - dt) / (2 * tau + dt);
            const double fac2 = (4 * sigma0 * dt) / (2 * tau + dt);
            const double factor = -eta * eta / 2;

            component /= arma::norm(component);

            os << "# step      time (au)         mux             muy             muz             Ex              Ey              Ez" << std::endl;

            for (size_t n = 0; n < maxstep; ++n)
            {
                arma::vec Vext(num_atoms, arma::fill::zeros);
                const double field = strength * gaussian_envelope(n * dt);
                for (size_t i = 0; i < num_atoms; ++i)
                {
                    Vext(i) = factor * field * arma::dot(component, coordinate.col(i));
                }

                R.fill(0);
                for (size_t i = 0; i < num_atoms; ++i)
                {
                    for (size_t j = 0; j < num_atoms; ++j)
                    {
                        R(i) += Kbar(j, i) * (Vext(i) - Vext(j));
                    }
                }

                J = fac1 * J + fac2 * (A * Q - R);
                Q = Q + J * dt;

                // induced dipole
                double mu_x{0};
                double mu_y{0};
                double mu_z{0};
                for (size_t i = 0; i < num_atoms; ++i)
                {
                    mu_x += Q(i) * coordinate(0, i);
                    mu_y += Q(i) * coordinate(1, i);
                    mu_z += Q(i) * coordinate(2, i);
                }

                // clang-format off
                os << std::setw(6) << n << std::setw(15) << std::setprecision(8) << std::fixed << n * dt
                   << std::setw(16) << std::setprecision(12) << std::fixed << mu_x
                   << std::setw(16) << std::setprecision(12) << std::fixed << mu_y
                   << std::setw(16) << std::setprecision(12) << std::fixed << mu_z
                   << std::setw(16) << std::setprecision(12) << std::fixed << component[0]*field
                   << std::setw(16) << std::setprecision(12) << std::fixed << component[1]*field
                   << std::setw(16) << std::setprecision(12) << std::fixed << component[2]*field
                   << "\n";
                // clang-format on
            }
            os << std::endl;
        }
    };
} // namespace wFQ
