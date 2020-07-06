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

#ifndef WALLCONICALDIVERGENT_H
#define WALLCONICALDIVERGENT_H

#include <Wall.h>
#include <cmath>
#include <stdexcept>

/**
    \brief Initialization options.
*/
enum class WallConicalDivergentOption {AngleAndLength, AngleAndExitRadius, LengthAndExitRadius};

/**

    \brief Class WallConicalDivergent provides the contour of
    a nozzle conical divergent section.

*/
class WallConicalDivergent: public Wall
{
public:

    /**
        \brief Constructor:
                   Rc: curvature radius of the nozzle wall at the throat section;
        par1 and par2: parameters to define the nozzle. They depend on WallConicalDivergentOption.

        Variables are parameterized in terms of the throat radius. tht is given in radians.
    */
    WallConicalDivergent(const WallConicalDivergentOption& opt, const double& Rc, const double& par1, const double& par2)
    {
        // Setting up the radius of curvature of the nozzle at the throat section
        this->Rc = Rc;

        // Setting up the divergent angle
        switch (opt)
        {
        case WallConicalDivergentOption::AngleAndLength :
            {
                tht  = par1;
                xEnd = par2;
            }
            break;

        case WallConicalDivergentOption::AngleAndExitRadius :
            {
                tht  = par1;
                Re   = par2;
            }
            break;

        case WallConicalDivergentOption::LengthAndExitRadius :
            {
                xEnd  = par1;
                Re    = par2;

                // The following code calculates tht iteratively

                double tht0 = atan2((Re-1.0),xEnd);

                tht = tht0;

                while(true)
                {
                    xc = xIntersec(Rc, tht);

                    rc = rIntersec(Rc, xc);

                    tht = atan2((Re-rc),(xEnd-xc));

                    if ( std::abs(tht-tht0) < 1E-15 ) break;

                    tht0 = tht;
                };
            }
            break;

        default:
            throw std::invalid_argument("Unknown option.");
            break;
        }

        // x at the intersection of the circular section and the conical section
        xc = xIntersec(Rc, tht);

        // r at the intersection of the circular section and the conical section
        rc = rIntersec(Rc, xc);

        // Setting up the remaining parameters...
        switch (opt)
        {
        case WallConicalDivergentOption::AngleAndLength :
            Re   = r(xEnd);
            break;

        case WallConicalDivergentOption::AngleAndExitRadius :
            xEnd = xc + (Re-rc)/tan(tht);
            break;

        case WallConicalDivergentOption::LengthAndExitRadius :
            // Everything already defined...
            break;

        default:
            throw std::invalid_argument("Unknown option.");
            break;
        }

        // xEnd <= 0.0 => keeps only the circular section
        if (xEnd <= 0.0) xEnd = xc;

    };


    /**
        \brief Destructor.
    */
    virtual ~WallConicalDivergent(){};


    /**
      \brief Returns the radial coord. r of the wall at x
    */
    virtual double r(const double& x)
    {
        return ( x<=xc ? Rc + 1.0 - sqrt(Rc*Rc-x*x) : tan(tht) * (x-xc) + rc );
    };


    /**
      \brief Returns the derivative dr/dx at x
    */
    virtual double drdx(const double& x)
    {
        return ( x <= xc ? x/sqrt(Rc*Rc-x*x) : tan(tht) );
    };


    /**
      \brief Returns the value of x at the nozzle end
    */
    virtual double x_end() const
    {
        return xEnd;
    };

private:

    /**
        \brief Calculates the x position of intersection of the conical and the circular section.
    */
    double xIntersec(const double& Rc, const double& tht)
    {
        // Calculating the intersection of the conical section and the circular section
        double tant = tan(tht);

        // x at the intersection of the circular section and the conical section
        return Rc * tant / sqrt(1.0+tant*tant);
    };


    /**
        \brief Calculates the r position of intersection of the conical and the circular section.
    */
    double rIntersec(const double& Rc, const double& xc)
    {
        // r at the intersection of the circular section and the conical section
        return Rc + 1.0 - sqrt(Rc*Rc-xc*xc);
    };


private:

    double   Rc; // Curvature radius of the wall at the throat
    double   Re; // Radius of the wall at the exit section
    double  tht; // Cone semi-angle (rad)
    double xEnd; // x position at the end of the nozzle
    double   xc; // x position at the intersection between the cone and the circular section
    double   rc; // r position at the intersection between the cone and the circular section
};

#endif // WALLCONICALDIVERGENT_H
