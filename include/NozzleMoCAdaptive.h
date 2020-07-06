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

#ifndef NOZZLEMOCADAPTIVE_H
#define NOZZLEMOCADAPTIVE_H

#include "NozzleMoC.h"

/**
    \brief NozzleMoCAdaptive provides adaptive features for characteristics
    construction.
*/
class NozzleMoCAdaptive: public NozzleMoC
{
    public:

        /**
            \brief Constructor
        */
        NozzleMoCAdaptive(double                  gamma,
                          double                    Par,
                          MoCToolBox                moc,
                          std::shared_ptr<Wall>    wall,
                          std::vector<MoCPoint>   bline,
                          double                dlminus,
                          double                dlplus);

        /**
            \brief Destructor
        */
        virtual ~NozzleMoCAdaptive();


        /**
            \brief Overrides NozzleMoC::solve.
        */
        virtual void solve(const double & xBegMax, const double & xEndMax);


    protected:

        /**
            \brief Calculates a subnet of characteristics from boundary line.
        */
        void buildSubCNetFromBLine(const double & xBegMax, const double & xEndMax);


        /**
            \brief Calculates a subnet of characteristics from the boundary line.
        */
        std::vector<MoCPoint> adaptiveNextCminus(const std::vector<MoCPoint>& pline,
                                                   const double& xBegMax,
                                                   const double& xEndMax);

        /**
            \brief Calculates the maximum distance between corresponding points of two
            C- characteristics pline and nline.
        */
        double distCminus(const std::vector<MoCPoint>& pline, const std::vector<MoCPoint>& nline);


        /**
            \brief fillCminus inserts points in the Cminus 'line'
            by interpolation such that the distance between points
            is less or equal to 'dlminus'.
        */
        void fillCminus(std::vector<MoCPoint>& line);

    protected:

        double  dlminus  {100.0}; // Maximum distance between points along a Cminus
        double  dlplus   {100.0}; // Maximum distance between two Cminus

};

#endif // NOZZLEMOCADAPTIVE_H
