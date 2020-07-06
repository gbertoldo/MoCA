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

#include "NumInterpolation1D.h"


NumInterpolation1D::NumInterpolation1D()
{
}

void NumInterpolation1D::init(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const NumInterpolation1DOption& option)
{
    if ( isInit ) free();

    size_t N = x.size();

    acc = gsl_interp_accel_alloc();

    switch (option)
    {
    case NumInterpolation1DOption::linear :

        spline_type = gsl_spline_alloc(gsl_interp_linear, N);

        break;

    case NumInterpolation1DOption::polynomial :

        spline_type = gsl_spline_alloc(gsl_interp_polynomial, N);

        break;

    case NumInterpolation1DOption::cspline :

        spline_type = gsl_spline_alloc(gsl_interp_cspline, N);

        break;

    case NumInterpolation1DOption::cspline_periodic :

        spline_type = gsl_spline_alloc(gsl_interp_cspline_periodic, N);

        break;

    case NumInterpolation1DOption::akima :

        spline_type = gsl_spline_alloc(gsl_interp_akima, N);

        break;

    case NumInterpolation1DOption::akima_periodic :

        spline_type = gsl_spline_alloc(gsl_interp_akima_periodic, N);

        break;

    case NumInterpolation1DOption::steffen :

        spline_type = gsl_spline_alloc(gsl_interp_steffen, N);

        break;

    default:
        break;
    }

    gsl_spline_init(spline_type, x.data(), y.data(), N);

    isInit = true;
}

NumInterpolation1D::~NumInterpolation1D()
{
    if ( isInit )
    {
        free();
    }
}

double NumInterpolation1D::eval(const double& x)
{
    return gsl_spline_eval(spline_type, x, acc);
}

void NumInterpolation1D::free()
{
    isInit = false;

    gsl_interp_accel_free(acc);

    gsl_spline_free(spline_type);

    acc = nullptr;

    spline_type = nullptr;
}

