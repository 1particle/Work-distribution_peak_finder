//============================================================================
//
// 	Work-distribution_peak_finder.cpp
// 	Copyright 2014  © Mostafa Nategholeslam
//
//	This file is part of Work-distribution_peak_finder.
//
//    Work-distribution_peak_finder is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Work-distribution_peak_finder is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
 //============================================================================

#include <iostream>
using namespace std;

#include <stdio.h>
#include <gsl/gsl_multifit.h>
using namespace std;
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
# include "functions.h"


int  main (int argc, char **argv) {
	double xInitial, xFinal;
	int numberOfBins, order;
	vector <xBin> forward, reverse, PMF;
	cout<<endl<<"Initial value of the range of reaction coordinate (in Angstroms) = ";
	cin>>xInitial;
	cout<<"Final value of the range of reaction coordinate (in Angstroms) = ";
	cin>>xFinal;
	cout<<"Number of bins along this reaction path = ";
	cin>>numberOfBins;
	cout<<"Order of polynomial to which the data should be fitted = ";
		cin>>order;
	xBinVectorBuilder	(xInitial, xFinal,numberOfBins, forward);
	xBinVectorBuilder	(xInitial, xFinal,numberOfBins, reverse);
	xBinVectorBuilder	(xInitial, xFinal,numberOfBins, PMF);
	solver (forward, reverse, PMF, order);
	outputter(forward, reverse, PMF, order);

	return 0;
}
