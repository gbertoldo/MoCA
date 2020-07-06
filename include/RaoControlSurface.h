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

#ifndef RAOCONTROLSURFACE_H
#define RAOCONTROLSURFACE_H

#include <vector>
#include "MoCToolBox.h"
#include "IsentropicFlowRelations.h"
#include "NumRootFinding.h"
#include "NumInterpolation1D.h"


/**

    \brief RaoControlSurface provides methods to calculate the
    control surface of the Rao's optimized nozzle

*/
class RaoControlSurface
{
public:

    /**
        \brief This constructor does not build MoCToolBox.
        Instead, it uses reference to the object of
        this classes.
    */
    RaoControlSurface(const double& gamma,
                      const double& Par,
                      const double& ME,
                      MoCToolBox& MoC,
                      const size_t N = 1000);


    /**
        \brief Destructor
    */
    virtual ~RaoControlSurface();


    /**
        \brief Returns the residual of the mass conservation (Eq. (19) of Rao's paper).
        This method also calculates the coordinates of the point D. If the residual
        is within the desired tolerance, user may call cplusDE method to get the control
        surface and cminusBD to get the C- from B to D.
        Code verified!
    */
    double residual(const std::vector<MoCPoint>& cminus);



    /**
        \brief Returns the control surface DE with N points
        The control surface is a C+ characteristics. Hence,
        this method returns a vector of MoCPoints.

        Caution: only call this method if the residual() is sufficiently small!
        Code verified!
    */
    std::vector<MoCPoint> cplusDE(size_t N);



    /**
        \brief Returns the C- characteristics from B to D.
        This method returns a vector of MoCPoints.

        Caution: only call this method if the residual() is sufficiently small!
        Code verified!
    */
    std::vector<MoCPoint> cminusBD(){ return CmBD; };

private:

    /**
        \brief Initializes the class
        Code verified!
    */
    void init();


    /**
        \brief Calculates M(theta) in accordance to Eq. (17) of Rao's paper.
        Code verified!
    */
    double Mt(const double & tht, size_t Nit = 1000, double tol = 1E-14);


    /**
        \brief Calculates eta(theta,M) in accordance to Eq. (18) of Rao's paper.
        Code verified!
    */
    double eta(const double & tht, const double & M);


    /**
        \brief Calculates g(theta, M). This is an auxiliary function to obtain eta(tht,M).
        Code verified!
    */
    double g(const double & tht, const double & M);


    /**
        \brief isThereACommonIntervalOfTheta returns true if there is a common
        interval o theta between C- from B to D and the control surface in order
        to search for the point D. This method also returns the common interval:
        [tlower, tupper].
        Code verified!
    */
    bool isThereACommonIntervalOfTheta(const double& tht1,
                                       const double& tht2,
                                       double& tlower,
                                       double& tupper)
    {
        // Getting the upper and lower value of tht of the current interval
        // of C-
        double tl = std::min(tht1,tht2);
        double tu = std::max(tht1,tht2);

        tlower = std::max(tl, thtmin);
        tupper = std::min(tu, thtmax);

        return tlower < tupper;
    };


private:

    double                            gamma   {1.4}; // Specific heat ratio
    double                              Par   {0.0}; // Ratio of environment pressure to stagnation pressure
    size_t                                N;         // Number of points to perform interpolation
    double                               ME;         // Exit Mach number at nozzle wall
    double                             thtE;         // Angle of the wall at nozzle exit
    double                              muE;         // Mach angle at the exit of the wall of the nozzle
    double                               xE;         // Axial position at the exit of the wall of the nozzle
    double                               rE;         // Radius at the exit of the wall of the nozzle
    MoCToolBox                          MoC;         // MoC tool box
    IsentropicFlowRelations        isenFlow;         // Isentropic flow relations
    NumInterpolation1D            thtInterp;         // tht(eta) interpolator
    NumInterpolation1D             I1Interp;         // I1(eta)  interpolator
    NumInterpolation1D             I2Interp;         // I2(eta)  interpolator

    std::vector<MoCPoint>              CmBD;         // C- from point B to point D

    double                           etamin;         // Minimum value of eta
    double                           thtmin;         // Minimum value of theta
    double                           thtmax;         // Maximum value of theta
    double                             Mmin;         // Minimum value of M
    double                             Mmax;         // Maximum value of M
};

#endif // RAOCONTROLSURFACE_H
