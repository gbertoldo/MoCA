/**

    \brief CSVHandler is a C++ utility to handle csv files.
    Handling includes reading, manipulating and writing csv files.

    \version 1.1.1
    \date 25/06/2020

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

#include "CSVHandler.h"

namespace csv{

const std::string about {"\n CSVHandler.\n Version 1.1.1.\n Last modified on 25/06/2020. \n\n Guilherme Bertoldo. \n\n"};

/**
    \brief Parse a string 'str' into a standard vector of strings using delim as the field delimiter.
    By default, merge delimiter is activated. To disable it, set 'merge_delim' as false.
*/
std::vector<std::string> ParseCSVRow(std::string str, char delim, bool merge_delim)
{
    // Stores a cell
    std::string cell;

    // Converting str to a stringstream in order to parse with getline command
    std::stringstream strstream(str);

    std::vector<std::string> row;

    while (std::getline(strstream,cell,delim))
    {
        if ( ! ( merge_delim && cell.empty() ) ) row.push_back(cell);
    }

    // This checks for a trailing delimiter with no data after it.
    if (!strstream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        row.push_back(cell);
    }

    return row;
}


/**
    \brief Overloads ostream << operator for a CSVMatrix
*/
std::ostream& operator<<(std::ostream& os, csv::CSVMatrix const& csvmatrix)
{
    for ( std::size_t i=0; i<csvmatrix.nrows(); ++i)
    {
        for ( std::size_t j=0; j<csvmatrix.ncols(); ++j)
        {
            os << csvmatrix(i,j) << "\t";
        }
        os << std::endl;
    }
    return os;
};
}
