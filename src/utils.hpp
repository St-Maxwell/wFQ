#pragma once

#include <iostream>
#include <string>
#include <armadillo>

namespace wFQ
{
   constexpr double ang2au = 1.8897260;
   constexpr double ev2au = 0.0367493;
   constexpr double fs2au = 41.3412764;
   constexpr double c2au = 137;

   inline void dump_charge_dipole(std::ostream &os, const arma::cx_vec &charge, const arma::mat &coordinate, const std::string &title)
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

      os << "# " + title << std::endl;
      os << "# mu_x: "
         << std::setw(15) << std::setprecision(9) << std::fixed << mu_x.real()
         << std::setw(15) << std::setprecision(9) << std::fixed << mu_x.imag()
         << std::endl;

      os << "# mu_y: "
         << std::setw(15) << std::setprecision(9) << std::fixed << mu_y.real()
         << std::setw(15) << std::setprecision(9) << std::fixed << mu_y.imag()
         << std::endl;

      os << "# mu_z: "
         << std::setw(15) << std::setprecision(9) << std::fixed << mu_z.real()
         << std::setw(15) << std::setprecision(9) << std::fixed << mu_z.imag()
         << std::endl;

      os << "# |mu|: " << std::setw(15) << std::setprecision(9) << std::fixed
         << sqrt(abs(mu_x) * abs(mu_x) + abs(mu_y) * abs(mu_y) + abs(mu_z) * abs(mu_z)) << std::endl;

      os << "#      Re(Q)          Im(Q)          Abs(Q)" << std::endl;
      for (auto Q : charge)
      {
         os << std::setw(15) << std::setprecision(9) << std::fixed << Q.real()
            << std::setw(15) << std::setprecision(9) << std::fixed << Q.imag()
            << std::setw(15) << std::setprecision(9) << std::fixed << abs(Q)
            << std::endl;
      }
      os << std::endl;
   }

} // namespace wFQ
