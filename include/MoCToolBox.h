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

#ifndef MOCTOOLBOX_H
#define MOCTOOLBOX_H

#include <utility>
#include <cmath>
#include "Wall.h"

/**
    \brief MoCPoint contains elementary information for MoC
*/
struct MoCPoint
{
    double   x; // x coordinate
    double   r; // radial coordinate
    double tht; // angle of velocity relatively to x
    double  nu; // Prandtl-Meyer function
    double   M; // Mach number
    double  mu; // Mach angle
};

/**
    \brief MoCToolBox provides basic operations to perform the Method of Characteristics
*/
class MoCToolBox
{
public:

    /**
      \brief Constructor
    */
    MoCToolBox(double gamma = 1.4, int Nit = 1000, double eps = 1E-13, double mzero = 1E-13);

    /**
      \brief Destructor
    */
    virtual ~MoCToolBox();

    /**
      \brief Returns the Prandtl-Meyer function
      Code verified!
    */
    double nu(const double& M);

    /**
      \brief Returns the Mach number that is the inverse of Prandtl-Meyer function
       Code verified!
    */
    double nu_inv(const double& nu);

    /**
      \brief Returns the Mach angle
      Code verified!
    */
    double mu(const double& M);


    /**
      \brief Given a point p1 lying on the C- characteristics
             and   a point p2 lying on the C+ characteristics,
             returns the fields of a point p3 at the intersection of C- and C+

                 p1
               C- \
                   \ p3
                   /
               C+ /
                 p2

             Code verified!
    */
    MoCPoint Cminus_to_Cplus(const MoCPoint& p1, const MoCPoint& p2);


    /**
      \brief Given a point p1 lying on the C+ characteristics
             and the wall, returns the point p2 at
             the intersection of C+ and the wall

             wall    p1
               \    / C+
                \  /
                 \/ p2
                  \
                   \

             Code verified!
    */
    MoCPoint Cplus_to_wall(Wall& wall, const MoCPoint& p1);


    /**
      \brief Given a point p1 lying on the C- characteristics
             and the wall, returns the point p2 at the intersection
             of C- and the wall

              p1
               \    wall
             C- \  /
                 \/
                 / p2

             Code verified!
    */
    MoCPoint Cminus_to_wall(Wall& wall, const MoCPoint& p1);


    /**
      \brief Given a point p1 lying on the C- characteristics,
             returns the point p2 at the intersection
             of C- and the symmetry line

                 p1
                 /
             C- /
               / p2
             -------

             Code verified!
    */
    MoCPoint Cminus_to_symmetry_line(const MoCPoint& p1);


    /**
      \brief Given a point p1 lying on the C+ characteristics,
             returns the point p2 at the intersection
             of C+ and the symmetry line

              p1
               \
             C+ \
                 \ p2
               -------

             Code verified!
    */
    MoCPoint Cplus_to_symmetry_line(const MoCPoint& p1);


    /**
      \brief Given the points p1 and p2 lying on the C- characteristics
             and a point pwall on the wall (only x, r and theta are known)
             returns the point p3 at the intersection of C- and the C+
             from pwall. It also calculates the remainder variables (M,
             mu and nu) at pwall.

                   pwall
               p1   / C+
              C-\  /
                 \/
                p3\
                   \
                   p2

             Code verified!
    */
    MoCPoint wall_to_Cminus(const MoCPoint& p1, const MoCPoint& p2, MoCPoint& pwall);


    /**
      \brief Given the points p1 and p2 lying on the C+ characteristics
             and a point pwall on the wall (only x, r and theta are known)
             returns the point p3 at the intersection of C+ and the C-
             from pwall. It also calculates the remainder variables (M,
             mu and nu) at pwall.

                   p2
                   / C+
                p3/
                 /\
                /  \
               /    \ C-
             p1    pwall

             Code verified!
    */
    MoCPoint wall_to_Cplus(const MoCPoint& p1, const MoCPoint& p2, MoCPoint& pwall);


    /**
      \brief Performs linear interpolation along a straight line between p1 and p2
      using the dimensionless variable zeta [0,1]. Returns p1 for zeta=0, p2 for zeta=1.
      Code verified!
    */
    MoCPoint linear_interp_zeta(const MoCPoint& p1, const MoCPoint& p2, const double& zeta)
    {
        MoCPoint p3;

        p3.x   = (p2.x   - p1.x  ) * zeta + p1.x;
        p3.r   = (p2.r   - p1.r  ) * zeta + p1.r;
        p3.tht = (p2.tht - p1.tht) * zeta + p1.tht;
        p3.M   = (p2.M   - p1.M  ) * zeta + p1.M;

        p3.mu  = this->mu(p3.M);
        p3.nu  = this->nu(p3.M);

        return std::move(p3);
    };


    /**
      \brief Performs linear interpolation along a straight line between p1 and p2
      using the variable x. Returns p1 for x=p1.x, p2 for x=p2.x.
      Code verified!
    */
    MoCPoint linear_interp_x(const MoCPoint& p1, const MoCPoint& p2, const double& x)
    {
        double zeta = (x-p1.x)/(p2.x-p1.x);

        return linear_interp_zeta(p1, p2, zeta);
    };


    /**
      \brief Performs linear interpolation along a straight line between p1 and p2
      using the variable r. Returns p1 for r=p1.r, p2 for r=p2.r.
      Code verified!
    */
    MoCPoint linear_interp_r(const MoCPoint& p1, const MoCPoint& p2, const double& r)
    {
        double zeta = (r-p1.r)/(p2.r-p1.r);

        return linear_interp_zeta(p1, p2, zeta);
    };


    /**
      \brief Performs linear interpolation along a straight line between p1 and p2
      using the variable tht. Returns p1 for tht=p1.tht, p2 for tht=p2.tht.
      Code verified!
    */
    MoCPoint linear_interp_tht(const MoCPoint& p1, const MoCPoint& p2, const double& tht)
    {
        double zeta = (tht-p1.tht)/(p2.tht-p1.tht);

        return linear_interp_zeta(p1, p2, zeta);
    };


    /**
      \brief Performs linear interpolation along a straight line between p1 and p2
      using the variable M. Returns p1 for M=p1.M, p2 for M=p2.M.
      Code verified!
    */
    MoCPoint linear_interp_M(const MoCPoint& p1, const MoCPoint& p2, const double& M)
    {
        double zeta = (M-p1.M)/(p2.M-p1.M);

        return linear_interp_zeta(p1, p2, zeta);
    };


    /**
        \brief Returns the distance between two MoCPoints
        Code verified!
    */
    double dist(const MoCPoint& p1, const MoCPoint& p2) const
    {
        return sqrt( (p2.x-p1.x) * (p2.x-p1.x) + (p2.r-p1.r) * (p2.r-p1.r) );
    }

private:

    // Derivative of nu relatively to M
    double dnudM(const double& M);

    // f+
    double fp(const double& tht, const double& mu);

    // f-
    double fm(const double& tht, const double& mu);

    // g+
    double gp(const double& tht, const double& M, const double& r);

    // g-
    double gm(const double& tht, const double& M, const double& r);

private:

    double gamma    {1.4}; // Specific heat ratio
    int    Nit     {1000}; // Maximum number of iterations
    double eps    {1E-13}; // Tolerance for numerical error
    double mzero  {1E-13}; // Machine zero

};


/**
    \brief Prints a MoC point p
*/
void printp(const MoCPoint& p, bool header=false);

#endif // MOCTOOLBOX_H
