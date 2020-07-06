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

#ifndef WALLINTERPOLATEDDIVERGENT_H
#define WALLINTERPOLATEDDIVERGENT_H

#include "Wall.h"
#include "NumInterpolation1D.h"

class WallInterpolatedDivergent: public Wall
{
public:

    WallInterpolatedDivergent(
        const NumInterpolation1DOption& opt,
        const double&                    Rc,
        const std::vector<double>&        x,
        const std::vector<double>&        r,
        const std::vector<double>&      tht);

    virtual ~WallInterpolatedDivergent();


    /**
      \brief Returns the radial coord. r of the wall at x
    */
    virtual double r(const double& x)
    {
        double aux;

        if ( x <= xc ) // circular section
        {
            aux = Rc + 1.0 - sqrt(Rc*Rc-x*x);
        }
        else if ( x <= xEnd ) // interpolated section
        {
            aux = rInterp.eval(x);
        }
        else // linear extrapolation section
        {
            aux = rEnd + drdxEnd * (x-xEnd);
        }
        return aux;
    }


    /**
      \brief Returns the derivative dr/dx at x
    */
    virtual double drdx(const double& x)
    {
        double aux;

        if ( x <= xc ) // circular section
        {
            aux = x/sqrt(Rc*Rc-x*x);
        }
        else if ( x <= xEnd ) // interpolated section
        {
            aux = tan(thtInterp.eval(x));
        }
        else // linear extrapolation section
        {
            aux = drdxEnd;
        }
        return aux;
    };

    /**
      \brief Returns the value of x at the nozzle end
    */
    virtual double x_end() const
    {
        return xEnd;
    };

private:
    double Rc;
    double xc;
    double rEnd;
    double xEnd;
    double drdxEnd;
    NumInterpolation1D   rInterp;
    NumInterpolation1D thtInterp;

};

#endif // WALLINTERPOLATEDDIVERGENT_H
