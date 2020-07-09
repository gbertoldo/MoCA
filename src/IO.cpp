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

#include <iostream>
#include <exception>
#include "IO.h"
#include "CSVHandler.h"
#include "KliegelLevineBoundaryLine.h"

namespace io
{

nlohmann::json configParameters;


void loadConfigFile(const std::string& filename){

    std::ifstream configFile(filename.c_str());

    if ( configFile.good() )
    {
        configFile >> configParameters;
    } else
    {
        throw std::invalid_argument("File "+filename+" is not good.");
    }

    configFile.close();

 }


void print_cnet(const std::vector<std::vector<MoCPoint>>& cnet, std::string filename)
{
    FILE* ofile;

    ofile = fopen(filename.c_str(),"w");

    fprintf(ofile,"%23s, %23s, %23s, %23s, %23s, %23s \n","# x", "y", "M", "tht", "nu", "mu");
    for(auto& vec: cnet)
    {
        for (auto& p: vec)
        {
            fprintf(ofile,"%23.15g, %23.15g, %23.15g, %23.15g, %23.15g, %23.15g \n",p.x, p.r, p.M, p.tht, p.nu, p.mu);
        }

        //auto& p = vec.back();
        //printf("%23.15g, %23.15g, %23.15g, %23.15g, %23.15g, %23.15g \n",p.x, p.r, p.M, p.tht, p.nu, p.mu);

        fprintf(ofile,"\n");

    }

    fclose(ofile);
};


void print_cline(const std::vector<MoCPoint>& cline, std::string filename)
{
    FILE* ofile;

    ofile = fopen(filename.c_str(),"w");

    fprintf(ofile,"%23s, %23s, %23s, %23s, %23s, %23s \n","# x", "y", "M", "tht", "nu", "mu");
    for (auto& p: cline)
    {
        fprintf(ofile,"%23.15g, %23.15g, %23.15g, %23.15g, %23.15g, %23.15g \n",p.x, p.r, p.M, p.tht, p.nu, p.mu);
    }
    fclose(ofile);
};


std::vector<MoCPoint> loadInitialLine(MoCToolBox& MoC)
{
    // Boundary line
    std::vector<MoCPoint> bline;

    bool          loadFromFile  = configParameters["NozzleInitialLine"]["LoadFromFile"];

    // Loading data from file
    if ( loadFromFile )
    {
        std::string       filename  = configParameters["NozzleInitialLine"]["InputFile"];
        std::string fieldDelimiter  = configParameters["NozzleInitialLine"]["FieldDelimiter"];

        csv::CSVFile csvfile(filename.c_str(), fieldDelimiter[0]);

        csvfile.load({});

        auto x   = csv::cast_to_vector<double>(csvfile.submatrix(1,csvfile.nrows()-1,0,0));
        auto y   = csv::cast_to_vector<double>(csvfile.submatrix(1,csvfile.nrows()-1,1,1));
        auto M   = csv::cast_to_vector<double>(csvfile.submatrix(1,csvfile.nrows()-1,2,2));
        auto tht = csv::cast_to_vector<double>(csvfile.submatrix(1,csvfile.nrows()-1,3,3));

        for (size_t i=0; i<x.size(); ++i)
        {
            MoCPoint p;

            p.x   = x[i];
            p.r   = y[i];
            p.M   = M[i];
            p.tht = tht[i];

            // Calculating mu and nu
            p.mu = MoC.mu(p.M);
            p.nu = MoC.nu(p.M);

            bline.push_back(std::move(p));
        }
    }
    // Generating initial line by the Kliegel-Levine approach
    else
    {
        // Reading parameters
        double gamma  = configParameters["Gas"]["SpecificHeatRatio"];
        double    Rc  = configParameters["NozzleInitialLine"]["CurvRadiusAtThroatLeft"];
        size_t    Np  = configParameters["NozzleInitialLine"]["NumOfPoints"];
        size_t   Nit  = configParameters["NozzleInitialLine"]["MaxIter"];
        double   tol  = configParameters["NozzleInitialLine"]["Tolerance"];

        // Creating the Kliegel-Levine calculator
        KliegelLevineBoundaryLine kl(gamma, Rc, Nit, tol);

        // Generating a Mach line that passes through the nozzle throat and contains Np points
        auto line = kl.MachLineThroughTheThroat(Np);

        // Generating MoCPoints from the Mach line
        for (size_t i=0; i<line.size(); ++i)
        {
            MoCPoint p;

            p.x   = (line[i]).z;
            p.r   = (line[i]).r;
            p.M   = (line[i]).M;
            p.tht = (line[i]).tht;

            // Calculating mu and nu
            p.mu = MoC.mu(p.M);
            p.nu = MoC.nu(p.M);

            bline.push_back(std::move(p));
        }
    }

    return bline;
}

}
