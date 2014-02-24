//============================================================================
//
// 	functions.h
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


#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

class workCount {
public:
	double w, count;
	workCount ()	{
		w = count = 0;
	}
};

class xBin {
public:
	double xInitial, xFinal , binSize;
	double w;
	double SDw;
	double R2;
	vector <workCount>  W;
	xBin (double x1, double x2)	{
		xInitial = x1;
		xFinal = x2;
		binSize = x2 - x1;
		w = SDw = R2 = 0;
		vector <workCount>  W;
	}
};


int readFile (int a, vector<xBin> & forward, vector<xBin> & reverse );
void wMatrix_builder (vector <workCount> & forward, vector <workCount> & reverse,  gsl_matrix* wForward,
		gsl_vector * countForward, gsl_vector * cForward, gsl_matrix * covForward,	 gsl_matrix * wReverse,
		gsl_vector * countReverse, gsl_vector * cReverse, gsl_matrix * covReverse, int order);
void solver (vector<xBin> & forward , vector<xBin> & reverse, vector<xBin> & PMF , int order);
void xBinVectorBuilder	(double xInitial, double xFinal, int numberOfBins, vector<xBin> & X);
void outputter (vector<xBin> & forward, vector<xBin> & reverse, vector<xBin> & PMF , int order);
double workPeakFinder 	(vector <workCount> & W, gsl_vector * c, int order);

#endif /* FUNCTIONS_H_ */
