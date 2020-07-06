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

#include "NumIntegration.h"
#include <cassert>

namespace NumIntegration
{

    double int_trapezoidal_rule(const std::vector<double>& f, const std::vector<double>& x)
    {
        double S {0.0};

        size_t N = x.size();

        assert( N == f.size() );

        for (size_t i = 1; i < N; ++i)
        {
            S+= ( f[i] + f[i-1] ) * ( x[i] - x[i-1] ) * 0.5;
        }

        return S;
    }
};
