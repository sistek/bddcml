/* BDDCML - Multilevel BDDC
 *
 * This program is a free software.
 * You can redistribute it and/or modify it under the terms of 
 * the GNU Lesser General Public License 
 * as published by the Free Software Foundation, 
 * either version 3 of the license, 
 * or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details
 * <http://www.gnu.org/copyleft/lesser.html>.
 *_______________________________________________________________*/

/*****************************************
* Macro for producing Fortran symbols based on 
* useful for calling Fortran/C interoperability
* Behaviour defined by preprocesor variable:
* UPPER - selects the upper case name
* Add_  - appends "_" at the end of the  lower case name (most common gfortran, ifort, etc.)
* Add__ - appends "__" at the end of the  lower case name
* if none is defined, lower case name is used
* inspired by MUMPS package
******************************************/

#ifndef F_SYMBOL

# if defined(UPPER) 
#  define F_SYMBOL(lower_case,upper_case) upper_case
# elif defined(Add_)
#  define F_SYMBOL(lower_case,upper_case) lower_case##_
# elif defined(Add__)
#  define F_SYMBOL(lower_case,upper_case) lower_case##__
# else
#  define F_SYMBOL(lower_case,upper_case) lower_case
# endif

#endif
