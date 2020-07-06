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

#include "NozzleMoC.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include "NumIntegration.h"
#include "Numerics.h"


NozzleMoC::NozzleMoC(double                  gamma,
                     double                    Par,
                     MoCToolBox                moc,
                     std::shared_ptr<Wall>    wall,
                     std::vector<MoCPoint>   bline):
    NozzleMoCInterface(),
    gamma(gamma),
    Par(Par),
    wall(wall),
    bline(bline),
    isenFlow(gamma),
    moc(moc)
{
    // Sorting r from the greatest to the lowest value
    sort(bline.begin(), bline.end(), [](const MoCPoint& a, const MoCPoint& b){return a.r>b.r;});
}


NozzleMoC::~NozzleMoC()
{
}


std::vector<MoCPoint> NozzleMoC::nextCminus(const std::vector<MoCPoint>& pline,
                                            const double               & xBegMax,
                                            const double               & xEndMax)
{
    assert( xEndMax-xBegMax >= -Numerics::eps );

    std::vector<MoCPoint> nline;

    // If the first point of pline has xBegMax<=x, returns an empty line
    if ( xBegMax <= pline.front().x ) return nline;

    // pline has N elements
    auto N = pline.size();

    // For each element of pline, starting from the second one,
    // creates another element for nline.
    for (size_t i = 1; i < N; ++i)
    {
        // Calculating the first point of nline (on the nozzle's wall)
        if ( i == 1 )
        {
            auto& p1 = pline[i];

            auto pw = moc.Cplus_to_wall(*wall, p1);

            // If pw is within the allowed interval,
            // just add it to nline.
            if ( pw.x <= xBegMax )
            {
                nline.push_back(pw);
            }
            // Otherwise, creates a point on the
            // wall with x=xBegMax
            else
            {
                MoCPoint pwall;

                pwall.x = xBegMax;

                pwall.r = wall->r(pwall.x);

                pwall.tht = atan(wall->drdx(pwall.x));

                auto& p0 = pline[0];

                moc.wall_to_Cminus(p0, p1, pwall);

                nline.push_back(pwall);

                // If xBegMax=xEndMax, the end of C- was reached,
                // just return.
                if ( std::abs(xBegMax-xEndMax) < mtol ) return nline;
            }
        }
        // Calculating a new point for nline for each
        // element of pline 2<= i <= N-1.
        else
        {
            // Taking a point along C-
            auto& p1 = nline.back();

            // A point along C+
            auto& p2 = pline[i];

            // Generating another point
            auto p3 = moc.Cminus_to_Cplus(p1,p2);

            // Is the new point within the acceptable interval?
            if ( xEndMax < p3.x) // If not within the acceptable interval
            {
                auto p4 = moc.linear_interp_x(p1, p3, xEndMax);

                nline.push_back(p4);

                break;
            }
            else  // If within the acceptable interval
            {
                // Just add p3 to nline and advance to the next point
                nline.push_back(p3);
            }
        }
    }
    // The last point to be added is an extension of the current C-
    // to the symmetry line (if possible) or to the end line
    {
        auto Nn = nline.size();

        // If nline has at least two points, then proceed...
        if ( Nn > 1 )
        {
            auto& p0 = nline[Nn-2]; // Penultimate point
            auto& p1 = nline[Nn-1]; // Last point

            // If the last point of nline is within the acceptable range
            if ( p1.x < xEndMax )
            {
                // Lets check if the last point of C- will touch the end line
                // or the symmetry line. If r at xEndMax > 0 along C-, the
                // last point will be on the end line (x=xEndMax).
                if ( (tan(p1.tht-p1.mu)*(xEndMax-p1.x)+p1.r) < 0.0 )
                {
                    auto p2 = moc.Cminus_to_symmetry_line(p1);

                    if ( p2.x <= xEndMax )
                    {
                        nline.push_back(p2);
                    }
                    else
                    {
                        auto p3 = moc.linear_interp_x(p1, p2, xEndMax);

                        nline.push_back(p3);
                    }

                }
                // Otherwise, C- will touch the end line
                else
                {
                    auto p2 = moc.linear_interp_x(p0, p1, xEndMax);

                    nline.push_back(p2);
                }
            }
        }

    }

    return nline;
}


std::vector<MoCPoint> NozzleMoC::nextCminus(const double & xBegMax,
                                            const double & xEndMax)
{
    return nextCminus(cnet.back(), xBegMax, xEndMax);
}


void NozzleMoC::solve(const double & xBegMax,
                      const double & xEndMax)
{
    // Cleaning memory
    cnet.clear();

    // Adding the first line to cnet
    cnet.push_back(bline);

    // Checking if the first point of the last C- characteristic
    // is less than xBegMax. If true, calculates the next C-.
    while( cnet.back().front().x < xBegMax )
    {

        // Calculating the next C-
        auto nline = nextCminus(xBegMax, xEndMax);

        // If the next C- is empty, break the cycle
        if ( nline.empty() )
        {
            break;
        }
        // Otherwise, add the C- to cnet
        else
        {
            cnet.push_back(nline);
        }
    }

    // Copying the boundary points
    extractBoundaryPoints(xEndMax);
}


void NozzleMoC::solve(const double & xEndMax)
{
    solve(wall->x_end(), xEndMax);
}


void NozzleMoC::solve()
{
    solve(wall->x_end(), wall->x_end());
}


double NozzleMoC::throatThrustCoeffient()
{
    std::vector<double>   xvec;
    std::vector<double>   rvec;
    std::vector<double>     f1;
    std::vector<double>     f2;

    // Sorting r from lowest to greatest value
    sort(bline.begin(), bline.end(), [](const MoCPoint& a, const MoCPoint& b){return a.r<=b.r;});

    // Extracting main data from boundary line
    for (auto& p: bline)
    {

        xvec.push_back(p.x);
        rvec.push_back(p.r);

        double r    = p.r;
        double M    = p.M;
        double sint = sin(p.tht);
        double cost = cos(p.tht);
        double pr   = isenFlow.Pr(M);

        double f1v =  pr * r * ( gamma * M * M * cost * cost + 1.0 );
        double f2v = -pr * r *   gamma * M * M * sint * cost;

        f1.push_back(f1v);
        f2.push_back(f2v);
    }

    double CT1 = NumIntegration::int_trapezoidal_rule(f1, rvec);

    double CT2 = NumIntegration::int_trapezoidal_rule(f2, xvec);

    //printf("CT1: %10g\n", CT1);
    //printf("CT2: %10g\n", CT2);

    return (CT1+CT2)*2.0;
}


double NozzleMoC::divSectionThrustCoeffient()
{
    std::vector<double>   rvec;
    std::vector<double>     f3;

    std::vector<MoCPoint> nline; // Line on the nozzle's wall

    // Obtaining data of the flow on the divergent of the nozzle's wall
    for (auto& line: cnet)
    {
        for (auto p: line)
        {
            // If the point belongs to the wall, copy it to nline
            if ( std::abs(wall->r(p.x)-p.r) < mtol )
            {
                nline.push_back(p);
            }
        }
    }

    // Sorting nline from lowest to greatest x
    std::sort(nline.begin(), nline.end(), [](const MoCPoint& a, const MoCPoint& b){return a.x <= b.x;});

    // Calculating f3
    for (auto& p: nline)
    {
        rvec.push_back(p.r);

        f3.push_back( p.r * isenFlow.Pr(p.M) );

        //printp(p);
    }

    double CT3 = NumIntegration::int_trapezoidal_rule(f3, rvec);

    //printf("CT3: %10g\n", CT3);

    return CT3*2.0;
}


double NozzleMoC::envPressureThrustCoeffient()
{
    double re  = wall->r(wall->x_end());

    double CT4 = -Par * re * re * 0.5;

    //printf("CT4: %10g\n", CT4);

    return CT4*2.0;
}


const std::vector<std::vector<MoCPoint>>& NozzleMoC::getCharacteristicNet()
{
    return cnet;
}


const std::vector<MoCPoint>& NozzleMoC::getWallMoCPoints()
{
    return wallMoCPoints;
}


const std::vector<MoCPoint>& NozzleMoC::getSymmMoCPoints()
{
    return symmMoCPoints;
}


const std::vector<MoCPoint>& NozzleMoC::getEndLMoCPoints()
{
    return endLMoCPoints;
}

void NozzleMoC::extractBoundaryPoints(const double& xEndMax)
{
    // Cleaning memory
    wallMoCPoints.clear();
    symmMoCPoints.clear();
    endLMoCPoints.clear();

    // For each line in cnet
    for (auto& line: cnet)
    {
        // For each point p in the line
        for (auto& p: line)
        {
            // Is p on the nozzle's wall?
            if ( std::abs(wall->r(p.x)-p.r) < mtol ) wallMoCPoints.push_back(p);

            // Is p on the symmetry line?
            if ( std::abs(p.r) < mtol ) symmMoCPoints.push_back(p);

            // Is p on the end line?
            if ( std::abs(p.x-xEndMax) < mtol ) endLMoCPoints.push_back(p);
        }
    }
}
