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

#include "NozzleMoCRao.h"
#include "WallInterpolatedDivergent.h"
#include "IsentropicFlowRelations.h"

NozzleMoCRao::NozzleMoCRao(
    double                                  Rc,
    std::shared_ptr<NozzleMoCInterface> solver,
    double                                  ME,
    double                               gamma,
    double                                 Par,
    std::shared_ptr<Wall>              wallExp,
    size_t        NumberOfPoints4Interpolation,
    size_t                             NpCplus,
    double                     fStartBisection,
    size_t                                 Nit,
    double                                 tol
    ):
    Rc(Rc),
    solver(solver),
    moc(gamma),
    wallExp(wallExp),
    rcs(gamma, Par, ME, this->moc, NumberOfPoints4Interpolation),
    optCon(this->moc, Nit, tol),
    NpCplus(NpCplus),
    fStartBisection(fStartBisection)
{

};

NozzleMoCRao::~NozzleMoCRao()
{
    //dtor
}



void NozzleMoCRao::solve(const double & xEndMax)
{
    // Calculating the net of characteristics from the expansion section
    solver->solve(wallExp->x_end(), xEndMax);

    // Getting a reference to the characteristics net
    auto& cnet = solver->getCharacteristicNet();

    // Now, it is necessary to evaluate the residual of the mass equation
    // for each C- in order to build the Rao's control surface. Starting
    // from the rightmost characteristic...
    double res0 = rcs.residual(cnet.back());

    // Indices of the C- that bracket the solution
    size_t ibeg;
    size_t iend;

    for (size_t i = cnet.size()-2; i >= 0; --i)
    {
        double res = rcs.residual(cnet[i]);
        //printf("%d %f\n", i, res);

        if (res*res0 < 0.0)
        {
            ibeg = i;
            iend = i+1;
            break;
        }

        res0 = res;
    }

    // x of the first point of the C- that bracket the solution
    double xbeg = cnet[ibeg].front().x;
    double xend = cnet[iend].front().x;

    // This C- will be used to generate the characteristic BD
    // that satisfies eq. 19 of Rao's paper (mass conservation)
    auto pline = cnet[ibeg];

    // Just a local function to solve the transcendental equation
    auto func = [&](const double& x)
    {
        // nextCminus generates from pline a C- from x to xEndMax
        return rcs.residual(nextCminus(pline, x, xEndMax));
    };

    // Solving the transcendental equation func
    double xB = NumRootFinding::bisection(func, xbeg*fStartBisection, xend);

    // This is the C- from xB to xEndMax applied to solve eq. (19) of Rao's paper
    auto nline = nextCminus(pline, xB, xEndMax);

    printf("Residual of Eq. (19) of Rao's paper: %E\n", rcs.residual(nline));

    // Getting the C- from B to D
    auto cminusBD = rcs.cminusBD();

    // Getting the C+ from D to E
    auto cplusDE = rcs.cplusDE(NpCplus);

    // Getting the optimized contour
    optContour = optCon.optContour(cminusBD, cplusDE);

    std::vector<double> xvec;
    std::vector<double> rvec;
    std::vector<double> thtvec;

    // Creating a new wall to the nozzle (optimized one)
    for (auto& p: optContour)
    {
        xvec.push_back(p.x);
        rvec.push_back(p.r);
        thtvec.push_back(p.tht);
    }

    wallOpt = std::make_shared<WallInterpolatedDivergent>(NumInterpolation1DOption::akima, Rc, xvec, rvec, thtvec);

    // Setting up the new wall
    solver->setWall(wallOpt);

    // Solving the entire optimized nozzle
    solver->solve(wallOpt->x_end(), wallOpt->x_end());
}
