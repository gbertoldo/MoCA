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

#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

#include "MoCToolBox.h"
#include <vector>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "json.hpp"

/**
    \brief Manages input and output
*/
namespace io
{


    /**
        \brief Stores the configuration parameters read from the configuration file
    */
    extern nlohmann::json configParameters;


    /**
        \brief Loads the configuration file
    */
    void loadConfigFile(const std::string& filename);


    /**
        \brief Prints characteristic net to output file filename
    */
    void print_cnet(const std::vector<std::vector<MoCPoint>>& cnet, std::string filename="cnet.txt");


    /**
        \brief Prints characteristic line to output file filename
    */
    void print_cline(const std::vector<MoCPoint>& cnet, std::string filename="cline.txt");


    /**
        \brief Reads the initial line, near the nozzle throat, as the
        boundary condition to MoC
    */
    std::vector<MoCPoint> loadInitialLine(MoCToolBox& MoC);
}


#endif // IO_H_INCLUDED
