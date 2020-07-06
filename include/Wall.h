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

#ifndef WALL_H
#define WALL_H

/*

   Interface for a generic nozzle wall contour

*/

class Wall
{
public:

    /**
      \brief Returns a linear approximation to the wall profile as r=a*x+b near x
    */
    virtual void linear_approximation(const double& x, double& a, double& b)
    {
        a=drdx(x);
        b=r(x)-drdx(x)*x;
    };


    /**
      \brief Returns the radial coord. r of the wall at x
    */
    virtual double r(const double& x) = 0;


    /**
      \brief Returns the derivative dr/dx at x
    */
    virtual double drdx(const double& x) = 0;

    /**
      \brief Returns the value of x at the nozzle end
    */
    virtual double x_end() const = 0;

};

#endif // WALL_H
