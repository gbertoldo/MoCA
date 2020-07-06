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

#ifndef NOZZLEMOCINTERFACE_H
#define NOZZLEMOCINTERFACE_H

#include <vector>
#include "MoCToolBox.h"

/**
    \brief NozzleMoCInterface is an abstract class that provides
    methods to solve the isentropic flow through a nozzle
    by the method of characteristics (MoC).
*/
class NozzleMoCInterface
{
public:

    /**
        \brief Destructor.
    */
    virtual ~NozzleMoCInterface(){};


    /**
        \brief Sets the nozzle's wall.
    */
    virtual void setWall(std::shared_ptr<Wall>) = 0;


    /**
        \brief Gets a shared pointer to nozzle's wall.
    */
    virtual std::shared_ptr<Wall> getWall() const = 0;


    /**
        \brief Given the previous C- (pline), nextCminus calculates and
        returns the next C-. The next C- starts on the nozzle wall. The
        first point of the next C- on the wall must have x <= xBegMax.
        The next C- ends on the symmetry line or on the end line. The end
        line is defined by x=xEndMax. By default, the next C- is added
        to the characteristic net, that is a vector of C-. This behavior
        is disable is addToNet is set as false.
    */
    virtual std::vector<MoCPoint> nextCminus(const std::vector<MoCPoint>& pline,
                                             const double               & xBegMax,
                                             const double               & xEndMax) = 0;


    /**
        \brief This version of nextCminus takes the previous C- in the
        characteristic net to calculate the next C-.
    */
    virtual std::vector<MoCPoint> nextCminus(const double & xBegMax,
                                             const double & xEndMax) = 0;


    /**
        \brief This method applies the nextCminus method until the
        x position of the C- characteristics on the nozzle wall is
        less or equal to xBegMax. The C- characteristics end on the
        symmetry line or on the end line (line of x=xEndMax).
    */
    virtual void solve(const double & xBegMax,
               const double & xEndMax) = 0;


    /**
        \brief This method overloads solve by setting xBegMax equals
        the x of the nozzle outlet.
    */
    virtual void solve(const double & xEndMax) = 0;


    /**
        \brief This method overloads solve by setting xBegMax and xEndMax
        equal the x of the nozzle outlet.
    */
    virtual void solve() = 0;


    /**
        \brief Returns the thrust coefficient due to the momentum
        transport across the nozzle throat.
    */
    virtual double throatThrustCoeffient() = 0;


    /**
        \brief Returns the thrust coefficient due to the pressure
        over the divergent section of the nozzle.
    */
    virtual double divSectionThrustCoeffient() = 0;


    /**
        \brief Returns the thrust coefficient due to the pressure
        difference between outlet pressure and environment pressure.
    */
    virtual double envPressureThrustCoeffient() = 0;


    /**
        \brief Returns a const reference to the characteristic net.
    */
    virtual const std::vector<std::vector<MoCPoint>>& getCharacteristicNet() = 0;


    /**
        \brief Returns a const reference to the MoCPoints on the nozzle wall.
    */
    virtual const std::vector<MoCPoint>& getWallMoCPoints() = 0;


    /**
        \brief Returns a const reference to the MoCPoints on the symmetry line.
    */
    virtual const std::vector<MoCPoint>& getSymmMoCPoints() = 0;


    /**
        \brief Returns a const reference to the MoCPoints on the end line.
    */
    virtual const std::vector<MoCPoint>& getEndLMoCPoints() = 0;

};

#endif // NOZZLEMOCINTERFACE_H
