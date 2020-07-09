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

#include "KliegelLevineBoundaryLine.h"
#include <cmath>
#include <cstdio>
#include "NumRootFinding.h"

KliegelLevineBoundaryLine::KliegelLevineBoundaryLine(
    const double& gamma,
    const double& Rc,
    const size_t& Nit,
    const double& tol):
    k(gamma),
    Rc(Rc),
    Nit(Nit),
    tol(tol)
{
    //ctor
}

KliegelLevineBoundaryLine::~KliegelLevineBoundaryLine()
{
    //dtor
}


double KliegelLevineBoundaryLine::u(const double& z, const double& r)
{
    double zm = sqrt(2*Rc/(k+1))*z;
    return
        pow(Rc+1,-3)*(((4*pow(k,2)-57*k+27)*pow(zm,3))/1.44E+2+2*(((3-2*k)*pow(zm,2))
                      /6.0E+0+(pow(r,2)+(-5.0E+0)/8.0E+0)*zm+((2*k+9)*pow(r,4))/2.4E+1-((4*k+15)*pow(r,2))
                      /2.4E+1+(10*k+57)/2.88E+2)+(((3-7*k)*pow(r,2))/8.0E+0+(13*k-27)/4.8E+1)*pow(zm,2)
                      +(((52*pow(k,2)+51*k+327)*pow(r,4))/3.84E+2-((52*pow(k,2)+75*k+279)*pow(r,2))
                        /1.92E+2+(92*pow(k,2)+180*k+639)/1.152E+3)*zm+zm+((556*pow(k,2)+1737*k+3069)
                                *pow(r,6))/1.0368E+4-((388*pow(k,2)+1161*k+1881)*pow(r,4))/2.304E+3
                      +((304*pow(k,2)+831*k+1242)*pow(r,2))/1.728E+3+pow(r,2)/2.0E+0
                      +((-2708*pow(k,2))-7839*k-14211)/8.2944E+4+(-1.0E+0)/4.0E+0)+pow(Rc+1,-2)
        *(((3-2*k)*pow(zm,2))/6.0E+0+(pow(r,2)+(-5.0E+0)/8.0E+0)*zm+zm+((2*k+9)*pow(r,4))
          /2.4E+1-((4*k+15)*pow(r,2))/2.4E+1+pow(r,2)/2.0E+0+(10*k+57)/2.88E+2
          +(-1.0E+0)/4.0E+0)+pow(Rc+1,-1)*(zm+pow(r,2)/2.0E+0+(-1.0E+0)/4.0E+0)+1;
}


double KliegelLevineBoundaryLine::v(const double& z, const double& r)
{
    double zm = sqrt(2*Rc/(k+1))*z;

    return
        pow(2,(-1.0E+0)/2.0E+0)*pow(pow(Rc+1,-1)*(k+1),1.0E+0/2.0E+0)
        *(pow(Rc+1,-3)*((-((7*k-3)*r*pow(zm,3))/1.2E+1)+(5.0E+0*(r*pow(zm,2)
                        +(((2*k+9)*pow(r,3))/6.0E+0-((4*k+15)*r)/1.2E+1)*zm+((k+3)*pow(r,5))
                        /9.0E+0-((20*k+63)*pow(r,3))/9.6E+1+((28*k+93)*r)/2.88E+2))/2.0E+0
                        +(((52*pow(k,2)+51*k+327)*pow(r,3))/1.92E+2-((52*pow(k,2)+75*k+279)*r)
                          /1.92E+2)*pow(zm,2)+(1.5E+1*(r*zm+pow(r,3)/4.0E+0-r/4.0E+0))/8.0E+0
                        +(((556*pow(k,2)+1737*k+3069)*pow(r,5))/1.728E+3-((388*pow(k,2)+1161*k+1181)
                                *pow(r,3))/5.76E+2+((304*pow(k,2)+831*k+1242)*r)/8.64E+2)*zm
                        +((6836*pow(k,2)+23031*k+30627)*pow(r,7))/8.2944E+4-((3380*pow(k,2)
                                +11391*k+15291)*pow(r,5))/1.3824E+4+((3424*pow(k,2)+11271*k+15228)
                                        *pow(r,3))/1.3824E+4+(((-7100*pow(k,2))-22311*k-30249)*r)/8.2944E+4)
          +pow(Rc+1,-2)*(r*pow(zm,2)+(3.0E+0*(r*zm+pow(r,3)/4.0E+0-r/4.0E+0))/2.0E+0
                         +(((2*k+9)*pow(r,3))/6.0E+0-((4*k+15)*r)/1.2E+1)*zm+((k+3)*pow(r,5))
                         /9.0E+0-((20*k+63)*pow(r,3))/9.6E+1+((28*k+93)*r)/2.88E+2)+pow(Rc+1,-1)
          *(r*zm+pow(r,3)/4.0E+0-r/4.0E+0));
}


double KliegelLevineBoundaryLine::Mzr(const double& z, const double& r)
{
    return Muv(u(z,r),v(z,r));
}


double KliegelLevineBoundaryLine::Muv(const double& u, const double& v)
{
    return pow(2.0,0.5)*pow(pow(k-1.0,-1.0)*(pow(1.0-(pow(v,2.0)
                            +pow(u,2.0))*(k-1.0)*pow(k+1.0,-1.0),-1.0)-1.0),0.5);
}


std::vector<BLPoint> KliegelLevineBoundaryLine::MachLineThroughTheThroat(size_t N)
{
    // Boundary line
    std::vector<BLPoint> bline;

    // Number of partitions
    size_t Np = N-1;

    // dr
    double dr = 1.0 / Np;

    // Mach number at the throat wall
    double Mt = Mzr(0.0, 1.0);

    // Creating the first point (at the throat wall)
    {
        BLPoint p;

        p.z   = 0.0;
        p.r   = 1.0;
        p.M   =  Mt;
        p.tht = 0.0;

        bline.push_back(std::move(p));
    }

    // Creating the Mach line
    for (size_t i = 1; i <= Np; ++i)
    {
        BLPoint p;

        // Defining the next r along the Mach line
        double r = dr * (Np-i);

        p.r = r;

        // Defining the function to determine the Mach line
        auto func = [&](const double& z)
        {
            return Mzr(z,r)-Mt;
        };

        /*
            Given r, determines z such that M(z,r)=Mt.
            The solution is within the z of the last point (z0) and
            probably z0+1.
        */

        // Taking of the last point
        auto z0 = bline.back().z;

        // Calculating z of the current point by the bisection method
        p.z = NumRootFinding::bisection(func, z0, z0+0.5, Nit, tol);

        // Calculating M(z,r)
        p.M = Mzr(p.z,p.r);

        // Calculating the angle of the velocity vector relatively to the axis of symmetry
        p.tht = atan(v(p.z,p.r)/u(p.z,p.r));

        //printf("%23.16e %23.16e %23.16e %23.16e\n",p.z, p.r, p.M, p.tht);

        // Adding the point to the boundary line
        bline.push_back(std::move(p));
    }
    return bline;
}
