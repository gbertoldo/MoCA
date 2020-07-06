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

#include "RaoNozzleOptContour.h"

#include <cmath>
#include <algorithm>
#include <cstdio>


RaoNozzleOptContour::RaoNozzleOptContour(MoCToolBox& MoC, size_t Nit, double tol):
    MoC(MoC),
    Nit(Nit),
    tol(tol)
{
    //ctor
}

RaoNozzleOptContour::~RaoNozzleOptContour()
{
    //dtor
}

MoCPoint RaoNozzleOptContour::contourPoint(const MoCPoint& p1, const MoCPoint& p2, const MoCPoint& p3, bool& isInside)
{
    MoCPoint p4 {p1};

    isInside = false;

    double h;
    double den;
    double zeta;
    double zeta0;

    zeta0 = 0.5;

    for (size_t i = 0; i < Nit; ++i)
    {
        h = tan(0.5*(p1.tht+p4.tht));

        den = p3.r - p2.r - h * ( p3.x - p2.x );

        if ( std::abs(den) > 1.E-15 )
        {
            zeta = -( p2.r - p1.r - h*(p2.x-p1.x) ) / den;

            p4.tht = (p3.tht-p2.tht) * zeta + p2.tht;
        }
        else
        {
            zeta = 1000.0;
        }

        if ( std::abs(zeta) > 2.0 ) return std::move(p4);

        if ( std::abs(zeta-zeta0) < tol ) break;

        zeta0 = zeta;
    }

    if ( (0.0 <= zeta) && (zeta <= 1.0) ){

        isInside = true;

        p4 = MoC.linear_interp_zeta(p2, p3, zeta);
    }

    return std::move(p4);
}


std::vector<MoCPoint> RaoNozzleOptContour::optContour(const std::vector<MoCPoint>& CminusBD,
                                                      const std::vector<MoCPoint>& CplusDE)
{
    // Optimized wall
    std::vector<MoCPoint> optWall;

    // Creating a new characteristic net from C- BD and C+ DE
    std::vector<std::vector<MoCPoint>> cnet;

    // Creating the first C- of the net.
    {
        // Sorting C- from the greatest to the lowest value of x, in order to travel
        // from the Rao control surface to the nozzle wall
        std::vector<MoCPoint> Cm = CminusBD;

        std::sort(Cm.begin(), Cm.end(), [](const MoCPoint& a, const MoCPoint& b){return a.x>b.x;});

        cnet.push_back(std::move(Cm));
    }

    // Each point of the C+ DE is the first point of a C-. Lets add these points to cnet.
    // Skipping the first point (point D), because it was already added through the C- BD.
    for (size_t i = 1; i < CplusDE.size(); ++i)
    {
        // Creating a C-
        std::vector<MoCPoint> cminus;

        // Copying p to C-
        cminus.push_back(CplusDE[i]);

        // Copying the C- to cnet
        cnet.push_back(cminus);
    }

    auto Np = cnet.size();

    // Adding the first point of the optimized wall: point B
    optWall.push_back(cnet[0].back());

    // This external loop uses the previous C- and the first point of the current
    // C- to generate the other points along the current C-. The last C- is the
    // point E (exit of the nozzle). It has only one point. That is why j index
    // stops before Np-1.
    for (size_t j = 1; j < Np-1; ++j)
    {
        // Getting the number of points of the previous C-
        auto Nm = cnet[j-1].size();

        // For each point of the previous C- (except first),
        // generate another point of current C-
        for (size_t i = 1; i < Nm; ++i)
        {
            // Point on C-
            auto& pm = cnet[j][i-1];

            // Point on C+
            auto& pp = cnet[j-1][i];

            // Calculating the point of intersection of C- and C+
            auto p = MoC.Cminus_to_Cplus(pm, pp);

            // Reference to the last point on the wall
            auto& pWall = optWall.back();

            //bool is_pw1_OnTheWall = false;
            bool is_pw2_OnTheWall = false;

            // Checking the intersection of the wall contour with C+
            //auto pw1 = contourPoint(pWall, p, pp, is_pw1_OnTheWall);

            // Checking the intersection of the wall contour with C-
            auto pw2 = contourPoint(pWall, p, pm, is_pw2_OnTheWall);

            // Adding p to cnet
            cnet[j].push_back(std::move(p));

            // Adding pw1 to optWall
            //if ( is_pw1_OnTheWall ) optWall.push_back(std::move(pw1));

            // Adding pw2 to optWall
            if ( is_pw2_OnTheWall ) optWall.push_back(std::move(pw2));

            // If pw1 or pw2 is on the wall, break the internal loop
            //if ( is_pw1_OnTheWall || is_pw2_OnTheWall ) break;
            if ( is_pw2_OnTheWall ) break;
        }
    }

    // Adding the last point of the optimized wall: point E
    optWall.push_back(cnet.back().front());

    // Returning optimized wall
    return optWall;
}
