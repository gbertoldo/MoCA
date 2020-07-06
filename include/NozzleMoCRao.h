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

#ifndef NOZZLEMOCRAO_H
#define NOZZLEMOCRAO_H

#include <memory>
#include <algorithm>
#include "IsentropicFlowRelations.h"
#include "NumIntegration.h"
#include "NozzleMoCInterface.h"
#include "RaoControlSurface.h"
#include "RaoNozzleOptContour.h"

/**

    \brief NozzleMoCRao applies the Method of Characteristics
    and the Rao's method to generate an axisymmetric nozzle
    profile that maximizes the thrust.
*/
class NozzleMoCRao: NozzleMoCInterface
{
public:
    /**
        \brief NozzleMoCRao uses an external solver to solve the
        Method of Characteristics.
    */
    NozzleMoCRao(
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
        double                                 tol);


    /**
        \brief Destructor.
    */
    virtual ~NozzleMoCRao();


    /**
        \brief Sets the initial expansion section of the nozzle's wall.
    */
    virtual void setWall(std::shared_ptr<Wall> wall)
    {
        this->wallExp = wall;
    };


    /**
        \brief Gets a shared pointer to nozzle's wall.
    */
    virtual std::shared_ptr<Wall> getWall() const
    {
        return wallOpt;
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual std::vector<MoCPoint> nextCminus(const std::vector<MoCPoint>& pline,
            const double               & xBegMax,
            const double               & xEndMax)
    {
        return solver->nextCminus(pline, xBegMax, xEndMax);
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual std::vector<MoCPoint> nextCminus(const double & xBegMax,
            const double & xEndMax)
    {
        return solver->nextCminus(xBegMax, xEndMax);
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual void solve(const double & xBegMax,
                       const double & xEndMax)
    {
        return solver->solve(xBegMax, xEndMax);
    };


    /**
        \brief The calculation of this method is overwritten by current class.
    */
    virtual void solve(const double & xEndMax);


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual void solve()
    {
        return solver->solve();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual double throatThrustCoeffient()
    {
        return solver->throatThrustCoeffient();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual double divSectionThrustCoeffient()
    {
        return solver->divSectionThrustCoeffient();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual double envPressureThrustCoeffient()
    {
        return solver->envPressureThrustCoeffient();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual const std::vector<std::vector<MoCPoint>>& getCharacteristicNet()
    {
        return solver->getCharacteristicNet();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual const std::vector<MoCPoint>& getWallMoCPoints()
    {
        return solver->getWallMoCPoints();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual const std::vector<MoCPoint>& getSymmMoCPoints()
    {
        return solver->getSymmMoCPoints();
    };


    /**
        \brief The calculation of this method is delegated to the external solver.
    */
    virtual const std::vector<MoCPoint>& getEndLMoCPoints()
    {
        return solver->getEndLMoCPoints();
    };


    /**
        \brief Returns the vector of MoC points on the optimized wall.
    */
    const std::vector<MoCPoint>& mocPointsOnOptimizedWall()
    {
        return optContour;
    };

private:

    double                                  Rc; // Curvature radius at the throat
    std::shared_ptr<NozzleMoCInterface> solver; // This is an external solver to apply the MoC to the current nozzle
    MoCToolBox                             moc; // MoC tool box
    std::shared_ptr<Wall>              wallExp; // Wall expansion section
    RaoControlSurface                      rcs; // Rao's control surface
    RaoNozzleOptContour                 optCon; // Rao's  optimized contour
    size_t                             NpCplus; // Number of points along the Cplus control surface
    double                     fStartBisection; // Fraction of the search interval to start the bisection method
    std::shared_ptr<Wall>              wallOpt; // Optimized wall
    std::vector<MoCPoint>           optContour; // MoC points along the optimized wall

};

#endif // NOZZLEMOCRAO_H
