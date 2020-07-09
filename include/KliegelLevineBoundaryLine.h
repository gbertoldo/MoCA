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

#ifndef KLIEGELLEVINEBOUNDARYLINE_H
#define KLIEGELLEVINEBOUNDARYLINE_H

#include <vector>

/**
    \brief \class KliegelLevineBoundaryLine provides the initial
    boundary line to apply the Method of Characteristics*.



    * Kliegel and Levine. Transonic Flow in Small Throat Radius
      of Curvature Nozzles. AIAA Journal, 1969.
*/

struct BLPoint
{
    double   z; // Axial coordinate
    double   r; // Radial coordinate
    double   M; // Mach number
    double tht; // Velocity vector angle relatively to the axial coordinate [rad]
};

class KliegelLevineBoundaryLine
{
public:

    /**
        \brief The constructor needs only the specific heat
        ratio \var gamma and the radius of curvature \var Rc at the
        throat wall (more precisely, at the throat from the subsonic
        section of the wall).
    */
    KliegelLevineBoundaryLine(
        const double& gamma,
        const double& Rc,
        const size_t& Nit = 1000,
        const double& tol = 1E-14);

    /**
        \brief Destructor
    */
    virtual ~KliegelLevineBoundaryLine();

    /**
        \brief u(z,r) = (axial component of the velocity vector)/(critical speed of sound)
                    z = (axial coordinate)/(radius of the throat)
                    r = (radial coordinate)/(radius of the throat)
    */
    double u(const double& z, const double& r);


    /**
        \brief v(z,r) = (radial component of the velocity vector)/(critical speed of sound)
                    z = (axial coordinate)/(radius of the throat)
                    r = (radial coordinate)/(radius of the throat)
    */
    double v(const double& z, const double& r);


    /**
        \brief M(z,r) = Mach number as a function of z and r
                    z = (axial coordinate)/(radius of the throat)
                    r = (radial coordinate)/(radius of the throat)
    */
    double Mzr(const double& z, const double& r);


    /**
        \brief M(u,v) = Mach number as a function of u and v
                    u = (axial component of the velocity vector)/(critical speed of sound)
                    v = (radial component of the velocity vector)/(critical speed of sound)
    */
    double Muv(const double& u, const double& v);


    /**
        \brief Returns the Mach line (line of constant Mach) that pass
        through the throat wall (z=0, r=1).

          *__   r  __* nozzle
              **|**
                |\
                | \ Mach line
                |  |
                --------------------> z

    */
    std::vector<BLPoint> MachLineThroughTheThroat(size_t N);

private:

    double   k; // Specific heat ratio
    double  Rc; // Radius of curvature at the throat of the nozzle (on the subsonic section side)
    size_t Nit; // Number of iterations for bisection method
    double tol; // Tolerance for bisection method
};

#endif // KLIEGELLEVINEBOUNDARYLINE_H
