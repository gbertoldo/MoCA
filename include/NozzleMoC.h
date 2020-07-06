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

#ifndef NOZZLEMOC_H
#define NOZZLEMOC_H

#include <memory>
#include "Wall.h"
#include "NozzleMoCInterface.h"
#include "IsentropicFlowRelations.h"
#include "MoCToolBox.h"

/**
    \brief NozzleMoC is an implementation class that provides
    methods to solve the isentropic flow through a nozzle
    by the method of characteristics (MoC).
*/
class NozzleMoC: public NozzleMoCInterface
{
public:

    /**
        \brief Constructor requires the ratio of specific heats gamma,
        the ratio of environment pressure to the stagnation pressure Par,
        the nozzle wall and
        the line of boundary conditions bline.
    */
    NozzleMoC(double                  gamma,
              double                    Par,
              MoCToolBox                moc,
              std::shared_ptr<Wall>    wall,
              std::vector<MoCPoint>   bline);

    /**
        \brief Destructor.
    */
    virtual ~NozzleMoC();


    /**
        \brief Sets the nozzle's wall.
    */
    virtual void setWall(std::shared_ptr<Wall> wall){ this->wall = wall;};


    /**
        \brief Gets a shared pointer to nozzle's wall.
    */
    virtual std::shared_ptr<Wall> getWall() const { return wall;};


    /**
        \brief Implements NozzleMoCInterface::nextCminus.
    */
    virtual std::vector<MoCPoint> nextCminus(const std::vector<MoCPoint>& pline,
                                     const double                       & xBegMax,
                                     const double                       & xEndMax);


    /**
        \brief Implements NozzleMoCInterface::nextCminus.
    */
    virtual std::vector<MoCPoint> nextCminus(const double & xBegMax,
                                             const double & xEndMax);


    /**
        \brief Implements NozzleMoCInterface::solve.
    */
    virtual void solve(const double & xBegMax,
               const double & xEndMax);


    /**
        \brief Implements NozzleMoCInterface::solve.
    */
    virtual void solve(const double & xEndMax);


    /**
        \brief Implements NozzleMoCInterface::solve.
    */
    virtual void solve();


    /**
        \brief Returns the thrust coefficient due to the momentum
        transport across the nozzle throat.
    */
    virtual double throatThrustCoeffient();


    /**
        \brief Returns the thrust coefficient due to the pressure
        over the divergent section of the nozzle.
    */
    virtual double divSectionThrustCoeffient();


    /**
        \brief Returns the thrust coefficient due to the pressure
        difference between outlet pressure and environment pressure.
    */
    virtual double envPressureThrustCoeffient();


    /**
        \brief Returns a const reference to the characteristic net.
    */
    virtual const std::vector<std::vector<MoCPoint>>& getCharacteristicNet();


    /**
        \brief Returns a const reference to the MoCPoints on the nozzle wall.
    */
    virtual const std::vector<MoCPoint>& getWallMoCPoints();


    /**
        \brief Returns a const reference to the MoCPoints on the symmetry line.
    */
    virtual const std::vector<MoCPoint>& getSymmMoCPoints();


    /**
        \brief Returns a const reference to the MoCPoints on the end line.
    */
    virtual const std::vector<MoCPoint>& getEndLMoCPoints();

protected:

    /**
        \brief Extracts (copy) boundary MoCPoints from characteristic net.
        (along the nozzle's wall, symmetry line and end line)
    */
    void extractBoundaryPoints(const double& xEndMax);


protected:

    double                      gamma  {1.4}; // Specific heat ratio
    double                        Par  {0.0}; // Ratio of environment pressure to stagnation pressure
    std::shared_ptr<Wall>     wall {nullptr}; // Nozzle's wall
    // Line of boundary conditions from the nozzle's throat to the axis of symmetry
    std::vector<MoCPoint>              bline; // Boundary line
    IsentropicFlowRelations         isenFlow; // Isentropic flow relations
    MoCToolBox                           moc; // Procedures to apply the method of characteristics
    std::vector<std::vector<MoCPoint>>  cnet; // Net of characteristics
    double                      mtol {1e-13}; // Machine tolerance
    std::vector<MoCPoint>      wallMoCPoints; // MoCPoints along the nozzle's wall
    std::vector<MoCPoint>      symmMoCPoints; // MoCPoints along the symmetry line
    std::vector<MoCPoint>      endLMoCPoints; // MoCPoints along the end line
};

#endif // NOZZLEMOC_H
