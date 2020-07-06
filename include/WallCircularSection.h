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

#ifndef WALLCIRCULARSECTION_H
#define WALLCIRCULARSECTION_H

#include "Wall.h"

class WallCircularSection: public Wall
{
public:
    /**
        \brief This constructor assumes that the circular section
        of curvature radius Rc is at the throat. So, x0=0, y0=1+Rc.
        Users must also define the inclination of the wall at the
        end of the wall thtEnd [rad].
    */
    WallCircularSection(const double& Rc, const double& thtEnd):
        Rc(Rc),
        thtEnd(thtEnd)
    {
        double tant = tan(thtEnd);
        xEnd = Rc * tant / sqrt(1.0+tant*tant);
    };


    /**
        \brief Destructor
    */
    virtual ~WallCircularSection(){};


    /**
        \brief Radius at the cross section position x
    */
    virtual double r(const double& x)
    {
        return Rc + 1.0 - sqrt(Rc*Rc-x*x);
    }

    /**
        \brief Inclination of the wall at the cross section position x
    */
    virtual double drdx(const double& x)
    {
        return x/sqrt(Rc*Rc-x*x);
    }


    /**
         \brief x position at the end of the wall
    */
    double x_end() const
    {
        return xEnd;
    }


private:

    double     Rc; // Curvature radius
    double thtEnd; // Inclination of the wall at the end of the wall
    double   xEnd; // x at the end of the wall

};

#endif // WALLCIRCULARSECTION_H
