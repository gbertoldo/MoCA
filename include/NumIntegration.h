/**
    MoCA - Method of Characteristics for Axisymmetric flows

    Copyright (C) 2020 by

    \author: Guilherme Bertoldo
             Federal University of Technology - Parana (UTFPR)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NUMINTEGRATION_H
#define NUMINTEGRATION_H

#include <vector>

namespace NumIntegration
{
    /**

        \brief Integrates a function f by the trapezoidal integration rule.

        Code verified!

    */
    double int_trapezoidal_rule(const std::vector<double>& f, const std::vector<double>& x);


};

#endif // NUMINTEGRATION_H
