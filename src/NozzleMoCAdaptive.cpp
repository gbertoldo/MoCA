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

#include "NozzleMoCAdaptive.h"

#include <algorithm>
#include <cmath>

NozzleMoCAdaptive::NozzleMoCAdaptive(double                  gamma,
                                     double                    Par,
                                     MoCToolBox                moc,
                                     std::shared_ptr<Wall>    wall,
                                     std::vector<MoCPoint>   bline,
                                     double                dlminus,
                                     double                dlplus):
    NozzleMoC(gamma, Par, moc, wall, bline),
    dlminus(dlminus),
    dlplus(dlplus)
{
    //ctor
}

NozzleMoCAdaptive::~NozzleMoCAdaptive()
{
    //dtor
}


void NozzleMoCAdaptive::buildSubCNetFromBLine(const double & xBegMax, const double & xEndMax)
{

    // Sorting r from lowest to greatest value
    std::sort(bline.begin(), bline.end(), [](const MoCPoint& a, const MoCPoint& b){return a.r<=b.r;});

    // Each element of the boundary condition (bline) will become
    // the first element of a C- characteristic of the kernel of
    // characteristics
    for (auto p: bline){

        std::vector<MoCPoint> cminus;

        cminus.push_back(p);

        // Adding the C- characteristics to the net of characteristics
        cnet.push_back(std::move(cminus));
    }

    /*

        Building a sub net of characteristics from bline

    */

    // Checking the possibility to get a downstream point for each point of the initial line

    // At first, suppose that all points of bline will generate a C-
    size_t imax = cnet.size()-1;

    for (size_t i = 1; i < cnet.size(); ++i)
    {
        auto& p1 = cnet[i  ].back();
        auto& p2 = cnet[i-1].back();

        // Calculating the angle, relatively to x axis, of the line linking points p1 and p2
        double angle_between_points = atan((p2.r-p1.r)/(p2.x-p1.x));

        // Calculating the angle of C- from point 1
        double cminus_angle = p1.tht-p1.mu;

        // If the following condition is satisfied, it is not possible to obtain
        // a new point downstream p1 and p2. Save the maximum value of i and exit the cycle.
        if ( cminus_angle < angle_between_points ) {

            imax = i-1;

            break;
        }
    }

    // For each C- characteristic starting at the initial line,
    // calculate new points along C- until it reaches the axis
    // of symmetry or until it reaches the end line (line of x=xEndMax)

    for (size_t i = 1; i <= imax; ++i){

        // For each point of the previous C- (i-1) and each point of
        // current C- (i), create a new intersection point. Note: each
        // point of previous C- emits a C+ characteristic that is used
        // in the calculation of the new point

        size_t jmax = (cnet[i-1]).size()-1;

        for (size_t j = 0; j <= jmax; ++j) {

            MoCPoint& p1 = cnet[i][j];   // Point of C-
            MoCPoint& p2 = cnet[i-1][j]; // Point of C+ (of previous C-)

            // Calculates the intersection point of C- and C+
            auto p3 = moc.Cminus_to_Cplus(p1, p2);

            // Checking if the calculated point has x>xEndMax
            if ( p3.x > xEndMax )
            {
                auto p4 = moc.linear_interp_x(p1, p3, xEndMax);

                cnet[i].push_back(std::move(p4));

                break;
            }

            // Adding p3 to current C-
            cnet[i].push_back(std::move(p3));

            // Get the point at the axis of symmetry
            if ( j == jmax ) {

                // Calculating the extrapolation of C- to the symmetry line
                auto& p1 = cnet[i].back(); // p1 is the last point of the current C- characteristic
                auto  p2 = moc.Cminus_to_symmetry_line(p1);

                // Checking if the calculated point is outside the allowed region
                if ( p2.x > xEndMax )
                {
                    auto p3 = moc.linear_interp_x(p1, p2, xEndMax);

                    cnet[i].push_back(std::move(p3));

                    break;
                }

                // Adding p2 to current C-
                cnet[i].push_back(std::move(p2));
            }
        }
    }

    {
        std::vector<MoCPoint> nline;

        for (size_t i = cnet.size()-2; i > imax; --i)
        {
            // The first point results from C+ intersection to the wall
            if ( i == cnet.size()-2 ) {

                auto& p0 = cnet[i+1].back();

                auto& p1 = cnet[i].back();

                auto  p2 = moc.Cplus_to_wall(*wall, p1);

                if ( p2.x > xBegMax )
                {
                    MoCPoint pwall;

                    pwall.x = xBegMax;

                    pwall.r = wall->r(pwall.x);

                    pwall.tht = atan(wall->drdx(pwall.x));

                    moc.wall_to_Cminus(p0, p1, pwall);

                    nline.push_back(pwall);

                    break;
                }
                else
                {
                    nline.push_back(std::move(p2));
                }
            }
            else
            {
                auto &p1 = nline.back();

                auto &p2 = cnet[i].back();

                auto p3 = moc.Cminus_to_Cplus(p1, p2);

                if ( p3.x > xEndMax )
                {
                    auto p4 = moc.linear_interp_x(p1, p3, xEndMax);

                    nline.push_back(std::move(p4));

                    break;
                }

                nline.push_back(std::move(p3));
            }
        }

        if ( ! nline.empty() ) {

            auto& pline = cnet[imax];

            for (size_t j = 0; j < pline.size(); j++)
            {
                auto& p1 = nline.back();
                auto& p2 = pline[j];

                auto p3 = moc.Cminus_to_Cplus(p1,p2);

                if ( p3.x > xEndMax )
                {
                    auto p4 = moc.linear_interp_x(p1, p3, xEndMax);

                    nline.push_back(std::move(p4));

                    break;
                }

                nline.push_back(std::move(p3));

                if ( j == pline.size()-1 ) {

                    // Calculating the extrapolation of C- to the symmetry line
                    auto& p1 = nline.back(); // p1 is the last point of the current C- characteristic
                    auto  p2 = moc.Cminus_to_symmetry_line(p1);

                    // Checking if the calculated point is outside the allowed region
                    if ( p2.x > xEndMax )
                    {
                        auto p3 = moc.linear_interp_x(p1, p2, xEndMax);

                        nline.push_back(std::move(p3));

                        break;
                    }

                    // Adding p2 to current C-
                    nline.push_back(std::move(p2));
                }
            }
        }

        fillCminus(nline);

        if ( ! nline.empty() ) cnet.push_back(nline);
    }
}


void NozzleMoCAdaptive::solve(const double & xBegMax, const double & xEndMax)
{
    // Cleaning memory
    cnet.clear();

    // Building sub-cnet from bline
    buildSubCNetFromBLine(xBegMax, xEndMax);

    // Checking if the first point of the last C- characteristic
    // is less than xBegMax. If true, calculates the next C-.
    while( cnet.back().front().x < xBegMax )
    {
        // Calculating the next C-
        //auto nline = nextCminus(xBegMax, xEndMax);
        auto nline = adaptiveNextCminus(cnet.back(), xBegMax, xEndMax);

        // If the next C- is empty, break the cycle
        if ( nline.empty() )
        {
            break;
        }
        // Otherwise, add the C- to cnet
        else
        {
            // Fills Cminus such that the distance between points is less or equal to dlminus
            fillCminus(nline);

            cnet.push_back(nline);
        }
    }

    // Copying the boundary points
    extractBoundaryPoints(xEndMax);
}


std::vector<MoCPoint> NozzleMoCAdaptive::adaptiveNextCminus(const std::vector<MoCPoint>& pline,
                                                              const double& xBegMax,
                                                              const double& xEndMax)
{
   double xBeg = pline.front().x;
   double   dx = xBegMax - pline.front().x;

   while (true)
   {
      auto nline = nextCminus(pline, xBeg+dx, xEndMax);

      if (nline.empty()) return nline;

      if ( distCminus(pline, nline) > dlplus )
      {
         dx = dx/2.0;
      }
      else
      {
         return nline;
      }
   }
}


double NozzleMoCAdaptive::distCminus(const std::vector<MoCPoint>& pline, const std::vector<MoCPoint>& nline)
{

   double dmax = 0.0;

   auto N = std::min( pline.size(), nline.size() );

   for (size_t i = 0; i < N; ++i)
   {
      double d = moc.dist( pline[i], nline[i] );

      if ( dmax < d ) dmax = d;
   }
   return dmax;
}


void NozzleMoCAdaptive::fillCminus(std::vector<MoCPoint>& line)
{
    std::vector<MoCPoint> refined_line;

    for (size_t i = 0; i < line.size()-1; ++i)
    {
        auto p1 = line[i];
        auto p2 = line[i+1];

        refined_line.push_back(p1);

        double d = sqrt( pow(p1.x-p2.x,2.0) + pow(p1.r-p2.r,2.0) );

        if ( d > dlminus ) {

            std::vector<MoCPoint> subline;

            size_t Np = (size_t) ceil(d/dlminus);

            for (size_t j = 1; j < Np; ++j) {

                double zeta = ( (double) j ) / ( (double) Np );

                MoCPoint p3;

                p3.x   = (p2.x  -p1.x  ) * zeta + p1.x;
                p3.r   = (p2.r  -p1.r  ) * zeta + p1.r;
                p3.M   = (p2.M  -p1.M  ) * zeta + p1.M;
                p3.tht = (p2.tht-p1.tht) * zeta + p1.tht;
                p3.mu  = moc.mu(p3.M);
                p3.nu  = moc.nu(p3.M);
                //printf("------ %g %g %g %g %g %g\n", p3.x, p3.r, p3.M, p3.tht, p3.mu, p3.nu);

                subline.push_back(std::move(p3));
            }

            if ( ! subline.empty() ) refined_line.insert(refined_line.end(),subline.begin(),subline.end());
        }
    }

    refined_line.push_back(line.back());

    line.clear();

    line = refined_line;
}

