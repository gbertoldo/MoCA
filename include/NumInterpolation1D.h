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

#ifndef NUMINTERPOLATION1D_H
#define NUMINTERPOLATION1D_H

#include<vector>
#include<string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>


/**
    \brief Interpolation options
*/
enum NumInterpolation1DOption{linear,
                              polynomial,
                              cspline,
                              cspline_periodic,
                              akima,
                              akima_periodic,
                              steffen};

/**

    \brief Wrapper to GSL interpolation

*/
class NumInterpolation1D
{
    public:

        /**
            \brief Constructor just requires the vector of x and y. x must be ordered from lowest
            to greatest value. Uniform partitioning is not required.
        */
        NumInterpolation1D();


        /**
            \brief Destructor
        */
        virtual ~NumInterpolation1D();


        /**
            \brief Initializer
        */
        void init(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const NumInterpolation1DOption& option = NumInterpolation1DOption::linear);


        /**
            \brief eval performs the interpolation
        */
        double eval(const double& x);

    private:
        /**
            \brief Cleans memory
        */
        void free();

    private:
        bool                    isInit    {false}; // Initialization flag
        gsl_interp_accel          *acc  {nullptr}; // 1D Index Look-up and Acceleration
        gsl_spline        *spline_type  {nullptr}; // Type of the spline
};

#endif // NUMINTERPOLATION1D_H
