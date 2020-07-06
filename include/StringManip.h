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

#ifndef STRINGMANIP_H_INCLUDED
#define STRINGMANIP_H_INCLUDED

#include <algorithm>

namespace strmanip
{

   // trim from start (in place)
   static inline void ltrim(std::string &s) {
       s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
           return !std::isspace(ch);
       }));
   }

   // trim from end (in place)
   static inline void rtrim(std::string &s) {
       s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
           return !std::isspace(ch);
       }).base(), s.end());
   }

   // trim from both ends (in place)
   static inline void trim(std::string &s) {
       ltrim(s);
       rtrim(s);
   }

   // trim from start (copying)
   static inline std::string ltrim_copy(std::string s) {
       ltrim(s);
       return s;
   }

   // trim from end (copying)
   static inline std::string rtrim_copy(std::string s) {
       rtrim(s);
       return s;
   }

   // trim from both ends (copying)
   static inline std::string trim_copy(std::string s) {
       trim(s);
       return s;
   }

}


#endif // STRINGMANIP_H_INCLUDED
