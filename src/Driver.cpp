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

#include "Driver.h"
#include <exception>

namespace Driver{

void run()
{
    std::string MoCSolver {io::configParameters["MoCSolver"]};

    if ( MoCSolver.compare("NozzleMoCRao") == 0 )
    {
        runNozzleMoCRao();
    }
    else if ( MoCSolver.compare("NozzleMoCAdaptive") == 0 )
    {
        runNozzleMoCAdaptive();
    }
    else if ( MoCSolver.compare("NozzleMoC") == 0 )
    {
        runNozzleMoC();
    }
    else
    {
        throw std::invalid_argument("Unknown option.");
    }

}

void runNozzleMoC()
{

}


void runNozzleMoCAdaptive()
{

        // Parameters
    double gamma  = io::configParameters["Gas"]["SpecificHeatRatio"];
    double Par    = io::configParameters["Environment"]["PressureRatio"];

    size_t Nit    = io::configParameters["MoCToolBox"]["MaxIter"];
    double tol    = io::configParameters["MoCToolBox"]["Tolerance"];
    double mzero  = io::configParameters["MoCToolBox"]["MachineZero"];


    /*
        Building the MoCToolBox
    */
    printf("Building the MoCToolBox...");

    MoCToolBox moc(gamma, Nit, tol, mzero);

    printf(" Ok\n");


    /*
        Building boundary line
    */
    printf("Loading the initial line...");

    auto bline = io::loadInitialLine(moc);

    printf(" Ok\n");


    /*
        Building the nozzle wall
    */
    printf("Building the nozzle wall...");

    std::shared_ptr<Wall> wall;

    std::string wallOption = io::configParameters["NozzleWall"];

    if ( wallOption.compare("NozzleWallConicalDivergent") == 0 )
    {

        double Rc = io::configParameters["NozzleWallConicalDivergent"]["CurvRadiusAtThroatRight"];

        std::string buildOption = io::configParameters["NozzleWallConicalDivergent"]["BuildOption"];

        double par1;
        double par2;

        if ( buildOption.compare("LengthAndExitRadius") == 0 )
        {
            par1 =  io::configParameters["NozzleWallConicalDivergent"]["DivergentLength"];
            par2 =  io::configParameters["NozzleWallConicalDivergent"]["DivergentExitRadius"];

            wall = std::make_shared<WallConicalDivergent>(WallConicalDivergentOption::LengthAndExitRadius, Rc, par1, par2);
        }
        else if ( buildOption.compare("AngleAndExitRadius") == 0 )
        {
            par1 =  io::configParameters["NozzleWallConicalDivergent"]["DivergentAngleDeg"];
            par2 =  io::configParameters["NozzleWallConicalDivergent"]["DivergentExitRadius"];

            // Converting angle from deg to rad
            par1 = par1 * acos(-1.0) / 180.0;

            wall = std::make_shared<WallConicalDivergent>(WallConicalDivergentOption::AngleAndExitRadius, Rc, par1, par2);
        }
        else if ( buildOption.compare("AngleAndLength") == 0 )
        {
            par1 =  io::configParameters["NozzleWallConicalDivergent"]["DivergentAngleDeg"];
            par2 =  io::configParameters["NozzleWallConicalDivergent"]["DivergentLength"];

            // Converting angle from deg to rad
            par1 = par1 * acos(-1.0) / 180.0;

            wall = std::make_shared<WallConicalDivergent>(WallConicalDivergentOption::AngleAndLength, Rc, par1, par2);
        }
        else throw std::runtime_error("Unknown option for building NozzleWallConicalDivergent.");
    }
    else if ( wallOption.compare("NozzleWallInterpolatedDivergent") == 0 )
    {
        double                       Rc = io::configParameters["NozzleWallInterpolatedDivergent"]["CurvRadiusAtThroatRight"];
        std::string interpolationOption = io::configParameters["NozzleWallInterpolatedDivergent"]["InterpolationOption"];
        std::string           inputFile = io::configParameters["NozzleWallInterpolatedDivergent"]["InputFile"];
        std::string      fieldDelimiter = io::configParameters["NozzleWallInterpolatedDivergent"]["FieldDelimiter"];
        std::string           angleUnit = io::configParameters["NozzleWallInterpolatedDivergent"]["AngleUnit"];

        NumInterpolation1DOption interpOption;

        // Checking interpolation option
        if ( interpolationOption.compare("linear") == 0 )
        {
            interpOption = NumInterpolation1DOption::linear;
        }
        else if ( interpolationOption.compare("cspline") == 0 )
        {
            interpOption = NumInterpolation1DOption::cspline;
        }
        else if ( interpolationOption.compare("akima") == 0 )
        {
            interpOption = NumInterpolation1DOption::akima;
        }
        else if ( interpolationOption.compare("steffen") == 0 )
        {
            interpOption = NumInterpolation1DOption::steffen;
        }
        else throw std::runtime_error("Invalid interpolation option. Options are: linear, cspline, akima and steffen.");

        // Reading the input file
        csv::CSVFile csvfile(inputFile, fieldDelimiter[0]);

        // Load file ignoring lines starting with '#'
        csvfile.load({'#'});

        auto x   = csv::cast_to_vector<double>(csvfile.col(0));
        auto r   = csv::cast_to_vector<double>(csvfile.col(1));
        auto tht = csv::cast_to_vector<double>(csvfile.col(2));

        // Check angle units
        if ( angleUnit.compare("deg") == 0 )
        {
            // Converting tht from deg to rad
            for (auto& angle: tht)
            {
                angle *= acos(-1.0)/180.0;
            }
        }
        else if ( angleUnit.compare("rad") == 0 )
        {
            // Do nothing
        }
        else throw std::runtime_error("Invalid angle unit option. Options are: deg, rad");

        wall = std::make_shared<WallInterpolatedDivergent>(interpOption, Rc, x, r, tht);

    }
    else throw std::runtime_error("Unknown option for building NozzleWall.");

    printf(" Ok\n");


     /*
        Building the MoC solver
    */
    printf("Building the MoC solver...");
    double dlminus = io::configParameters["NozzleMoCAdaptive"]["dlminus"];
    double dlplus  = io::configParameters["NozzleMoCAdaptive"]["dlplus"];

    std::shared_ptr<NozzleMoCInterface> nozzle = std::make_shared<NozzleMoCAdaptive>(gamma, Par, moc, wall, bline, dlminus, dlplus);

    printf(" Ok\n");

    /*
        Solving...
    */
    printf("Solving...\n");

    nozzle->solve();

    printf("Saving the solution... ");

    io::print_cnet(nozzle->getCharacteristicNet(),        "NozzleMoCAdaptive_cnet.txt");
    io::print_cline(nozzle->getWallMoCPoints(),     "NozzleMoCAdaptive_wallPoints.txt");
    io::print_cline(nozzle->getSymmMoCPoints(),     "NozzleMoCAdaptive_symmPoints.txt");
    io::print_cline(nozzle->getEndLMoCPoints(),     "NozzleMoCAdaptive_endLPoints.txt");

    printf(" Ok\n");


    double Ct1 = nozzle->throatThrustCoeffient();
    double Ct2 = nozzle->divSectionThrustCoeffient();
    double Ct3 = nozzle->envPressureThrustCoeffient();
    double Ct  = Ct1+Ct2+Ct3;

    FILE* fout;

    fout = fopen("NozzleMoCAdaptive_Ct.txt","w");

    printf("\n\n Thrust coefficient Ct: \n");
    printf("Ct due to momentum transfer through the throat: %23.16f\n", Ct1);
    printf("Ct due to pressure on the divergent section:    %23.16f\n", Ct2);
    printf("Ct due to environment pressure:                 %23.16f\n", Ct3);
    printf("Ct (total):                                     %23.16f\n", Ct);

    fprintf(fout, "\n\n Thrust coefficient Ct: \n");
    fprintf(fout, "Ct due to momentum transfer through the throat: %23.16f\n", Ct1);
    fprintf(fout, "Ct due to pressure on the divergent section:    %23.16f\n", Ct2);
    fprintf(fout, "Ct due to environment pressure:                 %23.16f\n", Ct3);
    fprintf(fout, "Ct (total):                                     %23.16f\n", Ct);

    fclose(fout);
}


void runNozzleMoCRao()
{

    // Parameters
    double gamma  = io::configParameters["Gas"]["SpecificHeatRatio"];
    double Par    = io::configParameters["Environment"]["PressureRatio"];
    double ME     = io::configParameters["NozzleWallRao"]["ExitMachNumberAtLip"];
    double Rc     = io::configParameters["NozzleWallRao"]["CurvRadiusAtThroatRight"];
    double thtMax = io::configParameters["NozzleWallRao"]["ThtMaxOnCircularExpansionDeg"];
    // Converting deg to rad
    thtMax = thtMax * acos(-1.0) / 180.0;

    size_t Nit    = io::configParameters["MoCToolBox"]["MaxIter"];
    double tol    = io::configParameters["MoCToolBox"]["Tolerance"];
    double mzero  = io::configParameters["MoCToolBox"]["MachineZero"];


    /*
        Building the MoCToolBox
    */
    printf("Building the MoCToolBox...");

    MoCToolBox moc(gamma, Nit, tol, mzero);

    printf(" Ok\n");


    /*
        Building boundary line
    */
    printf("Loading the initial line...");

    auto bline = io::loadInitialLine(moc);

    printf(" Ok\n");

    /*
        Initial expansion section of the wall of the nozzle
    */
    printf("Building the initial expansion section...");

    auto wallExpSec = std::make_shared<WallCircularSection>(Rc, thtMax);

    printf(" Ok\n");


    /*
        Building the MoC solver
    */
    printf("Building the MoC solver for nozzles...");

    std::shared_ptr<NozzleMoCInterface> solver;

    double          upToX = io::configParameters["NozzleMoCRao"]["solveUpToX"];
    std::string MoCSolver = io::configParameters["NozzleMoCRao"]["MoCSolver"];

    if ( MoCSolver.compare("NozzleMoC") == 0 )
    {
        solver = std::make_shared<NozzleMoC>(gamma, Par, moc, wallExpSec, bline);
    }
    if ( MoCSolver.compare("NozzleMoCAdaptive") == 0 )
    {
        double dlminus = io::configParameters["NozzleMoCAdaptive"]["dlminus"];
        double dlplus  = io::configParameters["NozzleMoCAdaptive"]["dlplus"];

        solver = std::make_shared<NozzleMoCAdaptive>(gamma, Par, moc, wallExpSec, bline, dlminus, dlplus);
    }
    else
    {
        throw std::invalid_argument("Unknown option.");
    }

    printf(" Ok\n");


    /*
        Building the Rao's nozzle
    */
    printf("Building the Rao's MoC solver for nozzles...");

    size_t NumberOfPoints4Interpolation = io::configParameters["NozzleMoCRao"]["RaoControlSurface"]["NumberOfPoints4Interpolation"];
    size_t                      NpCplus = io::configParameters["NozzleMoCRao"]["NumberOfPointsAlongDE"];
    double              fStartBisection = io::configParameters["NozzleMoCRao"]["fStartBisection"];
    size_t                    NitOptCtr = io::configParameters["NozzleMoCRao"]["RaoNozzleOptContour"]["MaxIter"];
    size_t                    tolOptCtr = io::configParameters["NozzleMoCRao"]["RaoNozzleOptContour"]["Tolerance"];

    NozzleMoCRao nozzle(
        Rc,
        solver,
        ME,
        gamma,
        Par,
        wallExpSec,
        NumberOfPoints4Interpolation,
        NpCplus,
        fStartBisection,
        NitOptCtr,
        tolOptCtr
        );

    printf(" Ok\n");


    /*
        Calculating the optimized nozzle
    */
    printf("Calculating the optimized nozzle...\n");

    nozzle.solve(upToX);


    /*
        Saving the solution
    */
    printf("Saving the solution...");

    io::print_cnet(nozzle.getCharacteristicNet(),             "NozzleMoCRao_cnet.txt");
    io::print_cline(nozzle.mocPointsOnOptimizedWall(),"NozzleMoCRao_optWallPoints.txt");
    io::print_cline(nozzle.getWallMoCPoints(),           "NozzleMoCRao_wallPoints.txt");
    io::print_cline(nozzle.getSymmMoCPoints(),           "NozzleMoCRao_symmPoints.txt");
    io::print_cline(nozzle.getEndLMoCPoints(),           "NozzleMoCRao_endLPoints.txt");

    printf(" Ok\n");


    double Ct1 = nozzle.throatThrustCoeffient();
    double Ct2 = nozzle.divSectionThrustCoeffient();
    double Ct3 = nozzle.envPressureThrustCoeffient();
    double Ct  = Ct1+Ct2+Ct3;

    FILE* fout;

    fout = fopen("NozzleMoCRao_Ct.txt","w");

    printf("\n\n Thrust coefficient Ct: \n");
    printf("Ct due to momentum transfer through the throat: %23.16f\n", Ct1);
    printf("Ct due to pressure on the divergent section:    %23.16f\n", Ct2);
    printf("Ct due to environment pressure:                 %23.16f\n", Ct3);
    printf("Ct (total):                                     %23.16f\n", Ct);

    fprintf(fout, "\n\n Thrust coefficient Ct: \n");
    fprintf(fout, "Ct due to momentum transfer through the throat: %23.16f\n", Ct1);
    fprintf(fout, "Ct due to pressure on the divergent section:    %23.16f\n", Ct2);
    fprintf(fout, "Ct due to environment pressure:                 %23.16f\n", Ct3);
    fprintf(fout, "Ct (total):                                     %23.16f\n", Ct);

    fclose(fout);
}

}
