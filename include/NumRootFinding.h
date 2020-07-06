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

#ifndef NUMROOTFINDING_H
#define NUMROOTFINDING_H

#include <cmath>
#include <cstdio>


namespace NumRootFinding
{

/**

    \brief Applies the bisection method to find a root of f in the interval [xb, xe]

*/
template<typename T>
double bisection(T f, double xb, double xe, size_t Nit = 1000, double tol = 1E-13){

    double xm {0.5 * ( xb + xe)};
    double fm;
    double fb = f(xb);
    double fe = f(xe);

    if ( fb*fe <= 0.0 )
    {
        for (size_t i = 0; i < Nit; ++i)
        {
            xm = 0.5 * ( xb + xe);

            fm = f(xm);

            if ( fb * fm <= 0.0 )
            {
                xe = xm;
                fe = fm;
            }
            else
            {
                xb = xm;
                fb = fm;
            }

            if ( std::abs(xe-xb) < tol )
            {
                xm = 0.5 * ( xb + xe);

                return xm;
            }
        }

        printf("bisection: impossible to find the root with prescribed tolerance within Nit iterations.\n");
    }
    else
    {
        printf("bisection: function has the same sign at the boundaries of the interval.\n");
    }

    return xm;
}

};

#endif // NUMROOTFINDING_H
