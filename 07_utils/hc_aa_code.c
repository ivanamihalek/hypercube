/*
This source code is part of specs, an application for
protein residue specialization and conservation scoring.
Written by Ivana Mihalek.
Copyright (C) 2010-2013 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "hypercube.h"

char single_letter ( char code[]){

    switch ( code[0] ) {
    case 'A':
	switch ( code [1]) {
	case 'L':
	    return 'A';
	    break;
	case 'R':
	    return 'R';
	    break;
	case 'S':
	    switch ( code[2] ) {
	    case 'N':
		return 'N';
		break;
	    case 'P':
		return  'D';
		break;
	    }
	    break;
	}
	break;
    case 'C':
	return 'C'; 
	break;
    case 'G':
	/* the second letter is always L */ 
	switch ( code[2] ) {
	case 'U':
	    return 'E';
	    break;
	case 'N':
	    return  'Q';
	    break;
	case 'Y':
	    return 'G';
	    break;
	}
	break;
    case 'H':
	return  'H';
	break;
    case 'I':
	return  'I';
	break;
    case 'L':
	switch ( code [1]) {
	case 'E':
	    return 'L';
	    break;
	case 'Y':
	    return 'K';
	    break;
	}
	break;
    case 'M':
	return 'M';
	break;
    case 'P':
	switch ( code [1]) {
	case 'H':
	    return 'F';
	    break;
	case 'R':
	    return 'P';
	    break;
	}
	break;
    case 'S':
	return 'S';
	break;
    case 'T':
	switch ( code [1]) {
	case 'H':
	    return 'T';
	    break;
	case 'R':
	    return 'W';
	    break;
	case 'Y':
	    return 'Y';
	    break;
	}
	break;
    case 'V':
	return 'V';
	break;
	
    }


    fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
    return 0;
}

