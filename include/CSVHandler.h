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

#ifndef CSVHANDLER_H_INCLUDED
#define CSVHANDLER_H_INCLUDED

#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include "StringManip.h"

namespace csv
{

// About this code
extern const std::string about;


/**

   \brief CSVMatrix defines a matrix storage for CSV data.
   All elements of the matrix are strings (std::string).

*/

class CSVMatrix
{

public:

    /**
        \brief Creates an empty CSVMatrix
    */
    CSVMatrix() {};


    /**
        \brief Creates a CSVMatrix with rows and cols filled with string str
    */
    CSVMatrix(std::size_t rows, std::size_t cols, std::string str)
    {
        for (std::size_t i=0; i<rows; ++i)
        {
            std::vector<std::string> row;

            for (std::size_t j=0; j<cols; ++j)
            {
                row.push_back(str);
            }

            csvmatrix.push_back(row);
        }
    };


    /**
        \brief Copy constructor
    */
    CSVMatrix(const CSVMatrix& other): csvmatrix(other.csvmatrix){};


    /**
        \brief Move constructor
    */
    CSVMatrix(CSVMatrix&& other): csvmatrix(std::move(other.csvmatrix)){};


    /**
        \brief Copy assignment
    */
    CSVMatrix& operator=(const CSVMatrix& other)
    {
        csvmatrix = other.csvmatrix;

        return *this;
    }


    /**
        \brief Move assignment deleted
    */
    CSVMatrix& operator=(CSVMatrix&& other) = delete;


    /**
        \brief Clear the data matrix
    */
    void clear()
    {
        csvmatrix.clear();
    }

    /**
        \brief Returns the number of rows of the CSVMatrix
    */
    std::size_t nrows() const
    {
        return csvmatrix.size();
    }

    /**
        \brief Returns the number of columns of the CSVMatrix
    */
    std::size_t ncols() const
    {
        return ( 0 < nrows() ? csvmatrix[0].size() : 0);
    }

    /**
        \brief Index operator which allow modification of the data matrix
    */
    std::string& operator()(std::size_t row, std::size_t col)
    {
        assert(row<nrows() && col<ncols());

        return (csvmatrix[row])[col];
    }

    /**
        \brief Index operator which DO NOT allow modification of the data matrix
    */
    std::string const& operator()(std::size_t row, std::size_t col) const
    {
        assert(row<nrows() && col<ncols());

        return (csvmatrix[row])[col];
    }

    /**
        \brief Adds a row of the same size of CSVMatrix number of columns
        containing strings str at row position i. If the CSVMatrix
        is empty, just add a row of one element.
    */
    void add_row(std::size_t i, std::string str = "")
    {

        assert( i <= nrows() );

        // If the csvmatrix is not empty
        if ( 0 < nrows() )
        {
            std::vector<std::string> row(ncols(), str);

            csvmatrix.insert(csvmatrix.begin()+i,row);
        }
        else
        {
            std::vector<std::string> row(1, str);

            csvmatrix.push_back(row);
        }
    }

    /**
        \brief Adds a vector of strings vrow at row position i.
        If the vector size is lower then the matrix number of columns,
        the size of the vector is adjusted to match the size of the matrix.
        The new elements are empty strings. If the reverse is true, the matrix
        number of columns is adjusted to match the size of the vector.
    */
    void add_row(std::size_t i, std::vector<std::string> vrow)
    {

        assert( i <= nrows() );

        // If the number of columns of the added matrix is less than the
        // number of columns of the current matrix, just add empty cells
        // to match sizes
        while ( vrow.size() < ncols() ) vrow.push_back("");

        // Otherwise, add columns to the current matrix if it is not empty
        while ( ncols() < vrow.size() && 0 < ncols() ) add_col(ncols(),"");

        csvmatrix.insert(csvmatrix.begin()+i,move(vrow));
    }


    /**
        \brief Adds a column of the same size of CSVMatrix number of rows
        containing strings str at column position j. If the CSVMatrix
        is empty, just add a column of one element.
    */
    void add_col(std::size_t j, std::string str = "")
    {

        assert( j <= ncols() );

        // If the csvmatrix is not empty
        if ( 0 < nrows() )
        {
            for (auto& row: csvmatrix)
            {
                row.insert(row.begin()+j,str);
            }
        }
        else // Otherwise
        {
            std::vector<std::string> row;

            row.push_back(str);

            csvmatrix.push_back(row);
        }
    }


    /**
        \brief Adds a vector of strings vcol at column position j.
        If the vector size is lower then the matrix number of rows,
        the size of the vector is adjusted to match the size of the matrix.
        The new elements are empty strings. If the reverse is true, the matrix
        number of rows is adjusted to match the size of the vector.
    */
    void add_col(std::size_t j, std::vector<std::string> vcol)
    {

        assert( j <= ncols() );

        // If the number of rows of the added vector is less than the
        // number of rows of the current matrix, just add empty cells
        // to the vector to match sizes
        while ( vcol.size() < nrows() ) vcol.push_back("");

        // Otherwise, add empty rows to the current matrix if it is not already empty
        while ( nrows() < vcol.size() && 0 < nrows() ) add_row(nrows(),"");

        // If csvmatrix is not empty
        if ( 0 < nrows() )
        {
            size_t i {0};

            for (auto& row: csvmatrix)
            {
                row.insert(row.begin()+j,vcol[i]);
                i++;
            }
        }
        else // Otherwise, if the csvmatrix is empty
        {
            for (size_t i = 0;  i < vcol.size(); ++i)
            {
                std::vector<std::string> vec;

                vec.push_back(vcol[i]);

                csvmatrix.push_back(vec);
            }
        }
    }


    /**
        \brief Removes the row at position i
    */
    void remove_row(std::size_t i)
    {
        assert( i < nrows() );

        csvmatrix.erase(csvmatrix.begin()+i);
    }


    /**
        \brief Removes the column at position j
    */
    void remove_col(std::size_t j)
    {
        assert( j < ncols() );

        for (auto& row: csvmatrix)
        {
            row.erase(row.begin()+j);
        }
    }


    /**
        \brief Adds a CSVMatrix m as the last row of the current CSVMatrix.
        If the number of columns of m is lower then the current matrix number of columns,
        the number of columns of m is adjusted to match the number of columns of
        the current matrix. The new elements are empty strings. If the reverse is true, the
        current matrix number of columns is adjusted to match the number of columns of m.
    */
    void add_rows(CSVMatrix m)
    {
        // If the number of columns of the added matrix is less than the
        // number of columns of the current matrix, just add empty cells
        // to match sizes
        while ( m.ncols() < ncols() ) m.add_col(m.ncols(),"");

        // Otherwise, add columns to the current matrix if it is not empty
        while ( ncols() < m.ncols() && 0 < ncols() ) add_col(ncols(),"");


        for (std::size_t i=0; i<m.nrows(); ++i)
        {
            std::vector<std::string> row;

            for (std::size_t j=0; j<m.ncols(); ++j)
            {
                row.push_back(move(m(i,j)));
            }
            csvmatrix.push_back(move(row));
        }
    }


     /**
        \brief Adds a CSVMatrix m as the last column of the current CSVMatrix.
        If the number of rows of m is lower then the current matrix number of rows,
        the number of rows of m is adjusted to match the number of rows of
        the current matrix. The new elements are empty strings. If the reverse is true, the
        current matrix number of rows is adjusted to match the number of rows of m.
    */
    void add_cols(CSVMatrix m)
    {
        for (std::size_t i=0; i<m.nrows(); ++i)
        {
            for (std::size_t j=0; j<m.ncols(); ++j)
            {
                csvmatrix[i].push_back(m(i,j));
            }
        }
    }


    /**
        \brief Returns a copy of the row i of the current CSVMatrix as a CSVMatrix.
    */
    CSVMatrix row(size_t i) const
    {
        return submatrix(i,i,0,ncols()-1);
    }


    /**
        \brief Returns a copy of the column j of the current CSVMatrix as a CSVMatrix.
    */
    CSVMatrix col(size_t j) const
    {
        return submatrix(0,nrows()-1,j,j);
    }

    /**
        \brief Returns a copy of column 'label' as a CSVMatrix. The 'label' is mapped by
        the labels in the header. By default, the header is assumed to be in the first
        row. If this is not true, user may indicate the row position of the header by index i.
        The returned column, by default, does not include rows before label. This behavior may
        be changed by marking 'remove_label' as false.
    */
    CSVMatrix col(std::string label, bool remove_label=true, size_t i = 0) const
    {
        CSVMatrix m;

        for ( size_t j = 0; j < ncols(); ++j)
        {
            std::string str {csvmatrix[i][j]};

            strmanip::trim(str);

            if ( label.compare(str) == 0 ) return (remove_label ? submatrix(i+1,nrows()-1,j,j) : submatrix(0,nrows()-1,j,j) );
        }

        class ColumnNotFound {};

        throw ColumnNotFound();

        return m;
    }

    /**
        \brief Returns a sub-matrix (slice) of the current CSVMatrix as a CSVMatrix.
        User must provide the row interval [ibeg, iend] and column interval [jbeg,jend]
        of the current CSVMatrix that will produce the sub-matrix.
    */
    CSVMatrix submatrix(size_t ibeg, size_t iend, size_t jbeg, size_t jend) const
    {
        assert( 0 <= ibeg && ibeg <= iend    );
        assert( 0 <= iend && iend <  nrows() );
        assert( 0 <= jbeg && jbeg <= jend    );
        assert( 0 <= jend && jend <  ncols() );

        CSVMatrix submatrix;

        for (size_t i = ibeg; i <= iend; ++i)
        {
            std::vector<std::string> row;

            for (size_t j = jbeg; j <= jend; ++j)
            {
                row.push_back(csvmatrix[i][j]);
            }
            submatrix.add_row(submatrix.nrows(),row);
        }

        return submatrix;
    }


protected:

    std::vector<std::vector<std::string>> csvmatrix; //! CSVMatrix basic data is a vector of a vector of strings

};






/**
    \brief Casts a CSVMatrix m to a standard vector of type T
*/
template <typename T>
std::vector<T> cast_to_vector(CSVMatrix m)
{
    T aux;

    std::vector<T> vec;

    for (size_t i = 0; i < m.nrows(); ++i)
    {
        for (size_t j = 0; j < m.ncols(); ++j)
        {
            std::istringstream ( m(i,j) ) >> aux;

            vec.push_back(aux);
        }
    }

    return vec;
}


/**
    \brief Casts a CSVMatrix m to a matrix (standard vector of vectors) of type T
*/
template <typename T>
std::vector<std::vector<T>> cast_to_matrix(CSVMatrix m)
{
    T aux;

    std::vector<std::vector<T>> matrix;

    for (size_t i = 0; i < m.nrows(); ++i)
    {
        std::vector<T> row;

        for (size_t j = 0; j < m.ncols(); ++j)
        {
            std::istringstream ( m(i,j) ) >> aux;

            row.push_back(aux);
        }

        matrix.push_back(row);
    }

    return matrix;
}


/**
    \brief Casts a standard vector of type T, vecT, to a standard vector of strings
*/
template <typename T>
std::vector<std::string> cast_vector_to_string(std::vector<T> vecT)
{
    std::vector<std::string> vec;

    for (auto& elem: vecT)
    {
        std::ostringstream oss;

        oss << elem;

        vec.push_back(oss.str());
    }

    return vec;
}


/**
    \brief Parse a string 'str' into a standard vector of strings using delim as the field delimiter.
    By default, merge delimiter is activated. To disable it, set 'merge_delim' as false.
*/
std::vector<std::string> ParseCSVRow(std::string str, char delim, bool merge_delim = true);


/**
    \brief Overloads ostream << operator for a CSVMatrix
*/
std::ostream& operator<<(std::ostream& os, csv::CSVMatrix const& csvmatrix);



/**

    \brief Class CSVFile extends CSVMatrix to read, manipulate and write csv files.

*/

class CSVFile: public CSVMatrix
{
public:
    /**
        \brief CSVFile constructor
        \param filename: name of the file to be read or to write in.
        \param field_delimiter: field delimiter
    */
    CSVFile(std::string filename, char field_delimiter):
        filename{filename}, field_delimiter{field_delimiter} {};


    /**
        \brief CSVFile copy constructor
    */
    CSVFile(const CSVFile& other):
        CSVMatrix(other),
        filename(other.filename),
        field_delimiter(other.field_delimiter)
        {};


    /**
        \brief CSVFile move constructor
    */
    CSVFile(CSVFile&& other):
        CSVMatrix(std::move(other)),
        filename(std::move(other.filename)),
        field_delimiter(std::move(other.field_delimiter))
        {};


    /**
        \brief CSVFile copy assignment operator
    */
    CSVFile& operator=(const CSVFile& other)
    {
        csvmatrix = other.csvmatrix;

        filename = other.filename;

        field_delimiter = other.field_delimiter;

        return *this;
    }


    /**
        \brief CSVFile move assignment operator (deleted)
    */
    CSVFile& operator=(CSVFile&& other) = delete;


    /**
        \brief Loads file from disk
        \param comments: a standard vector of characters. Lines that starts
                         with one of these characters are treated as comments.
                         They are not loaded. True by default.
        \param skip_empty_lines: skips empty lines. True by default.
        \param merge_delim: ignores empty cells. True by default.
    */
    void load(std::vector<char> comments = {}, bool skip_empty_lines = true, bool merge_delim = true)
    {
        std::fstream file(filename);

        std::string line;

        while( std::getline(file,line) )
        {

            std::vector<std::string> row;

            std::string str = line;

            strmanip::trim(str);

            // Check for blank lines
            if ( ! (skip_empty_lines && str.empty()) ){

                // Variables to check for skipping commented lines
                bool skip_commented_line {false};

                // Check for comments. If a line starts with a character of 'comments', it is ignored.
                if ( !str.empty() && !comments.empty()) {

                    char first_char = str[0];

                    for (auto c: comments){
                        if (first_char == c) {

                            skip_commented_line = true;

                            break;
                        }
                    }
                }

                if ( ! skip_commented_line )
                {
                    // Converting a line string to a CSVMatrix (row matrix)
                    row = ParseCSVRow(line,field_delimiter, merge_delim);

                    add_row(nrows(),row);
                }
            }
        }
        file.close();
    };



    /**
        \brief Saves the file to disk.
    */
    void save()
    {
        std::ofstream file(filename);

        for (std::size_t i=0; i<nrows(); ++i)
        {
            for (std::size_t j=0; j<ncols()-1; ++j)
            {
                file << (*this)(i,j) << field_delimiter;
            }

            file << (*this)(i,ncols()-1) << std::endl;
        }

        file.close();
    };

    /**
        \brief Saves the file to disk, and closes file.
    */
    void close()
    {
        save();
    };

public:

    std::string filename;
    char field_delimiter;

};

}


#endif // CSVHANDLER_H_INCLUDED
