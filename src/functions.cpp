//============================================================================
//
// 	functions.cpp
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
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <cstring>
# include "functions.h"

#include <iomanip>
#include <locale>
#include <sstream>
#include <sys/stat.h>




int readFile (int a, vector<xBin> & forward, vector<xBin> & reverse )	{
	string line;
	char* st;
	ifstream data;
	data.precision(18);
	int i = 0;

	ostringstream x1, x2;
	x1 <<  forward[a].xInitial;
	x2 <<  forward[a].xFinal;
	string fileName = "Bin_"+ x1.str() + "_to_"+ x2.str() + "_Angstrom.txt";
	data.open(fileName.c_str());
	while ( data.peek() != EOF ) 	{
		getline (data,line);
		//		cout<<i<<"  "<<line<<endl;
		if (i > 0  ) {
		workCount *wf 	= new workCount();
		workCount *wr 	= new workCount();
		wf->w 			= strtod(line.c_str(),&st);
		wf->count 		= strtod(st,&st);
		wr->w 			= strtod(st,&st);
		wr->count 		= strtod(st,&st);
		forward[a].W.push_back(*wf);
		reverse[a].W.push_back(*wr);
					}
			i++;
		}

	i--;  // That annoying newline character!
	data.close();
	return i;
}

void wMatrix_builder (vector <workCount> & forwardW, vector <workCount> & reverseW,  gsl_matrix* wForward,
		gsl_vector * countForward, gsl_vector * cForward, gsl_matrix * covForward,	 gsl_matrix * wReverse,
		gsl_vector * countReverse, gsl_vector * cReverse, gsl_matrix * covReverse, int order)	{

		int n = (int) forwardW.size();  // Number of work bins.

		for (int i = 0; i<n ; i++)	{  // Loop over work bins
			for (int j=0 ; j<= order ; j++)	{
		gsl_matrix_set (wForward, i, j, pow(forwardW[i].w,j));
		gsl_matrix_set (wReverse, i, j, pow(reverseW[i].w,j));
			}

		gsl_vector_set (countForward, i, forwardW[i].count);
		gsl_vector_set (countReverse, i, reverseW[i].count);
		}
}

void solver (vector<xBin> & forward , vector<xBin> & reverse, vector<xBin> & PMF , int order)	{

	for (int i = 0 ; i <(int)forward.size() ; i++)	{
		int n= readFile (i, forward, reverse);

			double chisqForward, chisqReverse;
			gsl_matrix *wForward, *wReverse, *covForward, *covReverse;
			gsl_vector *countForward, *countReverse, *cForward, *cReverse;

			wForward = gsl_matrix_alloc (n, order+1);
			wReverse = gsl_matrix_alloc (n, order+1);
			countForward = gsl_vector_alloc (n);
			countReverse = gsl_vector_alloc (n);

			cForward = gsl_vector_alloc (order+1);
			cReverse = gsl_vector_alloc (order+1);
			covForward = gsl_matrix_alloc (order+1, order+1);
			covReverse = gsl_matrix_alloc (order+1, order+1);


			wMatrix_builder ( forward[i].W, reverse[i].W, wForward, countForward, cForward, covForward,
								wReverse, countReverse, cReverse, covReverse, order);

			{
			gsl_multifit_linear_workspace * workForward = gsl_multifit_linear_alloc (n, order+1);
			gsl_multifit_linear_workspace * workReverse = gsl_multifit_linear_alloc (n, order+1);
			gsl_multifit_linear (wForward, countForward, cForward, covForward,  &chisqForward, workForward);
			gsl_multifit_linear (wReverse, countReverse, cReverse, covReverse,  &chisqReverse, workReverse);
			gsl_multifit_linear_free (workForward);
			gsl_multifit_linear_free (workReverse);
			}

			double * wCount = new double [n];

			for (int a = 0; a<n ; a++)	wCount[a]= forward[i].W[a].count;
			forward[i].R2 = 1 - chisqForward / gsl_stats_tss (wCount, 1, n );

			for (int a = 0; a<n ; a++)	wCount[a]= reverse[i].W[a].count;
			reverse[i].R2 = 1 - chisqReverse / gsl_stats_tss (wCount, 1, n );


			if (order ==2)	{
			double c1, c2, dc1, dc2;
			c1 = gsl_vector_get(cForward,1);   						c2 = gsl_vector_get(cForward,2);
			dc1 = sqrt( gsl_matrix_get (covForward, 1, 1) );		dc2 = sqrt( gsl_matrix_get (covForward, 2, 2) );

			forward[i].w = -c1/(2*c2);
			forward[i].SDw = fabs(-c1/(2*c2)) * (sqrt(pow(dc1/c1, 2) + pow(dc2/c2, 2)));

			c1 = gsl_vector_get(cReverse,1);   						c2 = gsl_vector_get(cReverse,2);
			dc1 = sqrt( gsl_matrix_get (covReverse, 1, 1) );		dc2 = sqrt( gsl_matrix_get (covReverse, 2, 2) );
			reverse[i].w = -c1/(2*c2);
			reverse[i].SDw = fabs(-c1/(2*c2)) * (sqrt(pow(dc1/c1, 2) + pow(dc2/c2, 2)));

			PMF[i].SDw = sqrt (forward[i].SDw * forward[i].SDw + reverse[i].SDw * reverse[i].SDw);
			}
			else	{
				forward[i].w = workPeakFinder 	(forward[i].W, cForward , order);
				reverse[i].w = workPeakFinder 	(reverse[i].W, cReverse , order);
			}
			PMF[i].w = (forward[i].w - reverse[i].w) / 2;

			cout<<i<<"  "<<PMF[i].w <<endl;

			gsl_matrix_free (wForward);
			gsl_vector_free (countForward);
			gsl_vector_free (cForward);
			gsl_matrix_free (covForward);

			gsl_matrix_free (wReverse);
			gsl_vector_free (countReverse);
			gsl_vector_free (cReverse);
			gsl_matrix_free (covReverse);
	}

}


void xBinVectorBuilder	(double xInitial, double xFinal, int numberOfBins, vector<xBin> & X)	{
	double binSize = (xFinal - xInitial)/(double)numberOfBins;
	for (int i = 0 ; i<numberOfBins ; i ++)	{
		xBin *x = new xBin (xInitial + (double)i*binSize , xInitial + (double)(i+1)*binSize );
		X.push_back(*x);
	}

}

void outputter (vector<xBin> & forward, vector<xBin> & reverse, vector<xBin> & PMF, int order)	{

	ostringstream ord;
	ord <<  order;
	string fileName;
	if (order == 2)
		fileName = "works_"+ ord.str() + "nd_order.txt";
	else if (order == 3)
		fileName = "works_"+ ord.str() + "rd_order.txt";
	else
		fileName = "works_"+ ord.str() + "th_order.txt";

	ofstream works;
	works.open(fileName.c_str());
	works.precision(18);
	(order == 2)  ?	works<<"x  PMF SD_PMF PMF_difference wForward SDwForward wReverse SDwReverse R2Forward R2Reverse" : works<<"x  PMF PMF_difference wForward wReverse R2Forward R2Reverse" ;

	double * error_temp = new double [(int)PMF.size()];
	double * PMF_difference = new double [(int)PMF.size()];

	for (int i = 0 ; i <(int)PMF.size() ; i++)
		PMF_difference[i] = PMF[i].w;

	for (int i = 1 ; i <(int)PMF.size() ; i++)	{
	              PMF[i].w += PMF[i-1].w;
	              PMF[i].SDw = sqrt(pow(PMF[i].SDw,2)+pow(PMF[i-1].SDw,2));
	            }
	            // Getting the errors from the opposite size and combining the errors with appropriate weights
	            for (int i=(int)PMF.size()-1; i >= 0 ; i --)         {
	                error_temp[i] = sqrt (pow(forward[i].SDw,2)+pow(reverse[i].SDw,2));
	                if (i < (int)PMF.size()-1)
	                      error_temp[i] = sqrt( pow(error_temp[i], 2) + pow(error_temp[i+1], 2) );
	                PMF[i].SDw = (PMF[i].SDw * error_temp[i]) / sqrt( pow(PMF[i].SDw, 2) + pow(error_temp[i],2) );

	            }


	            // Finally, writing to the output file:
	for (int i = 0 ; i <(int)forward.size() ; i++)	{
		if (i == 0)
			works<<endl<<forward[0].xInitial<< "   0   0  0";
		if (order == 2)
			works<<endl<<forward[i].xFinal<< "   "<<PMF[i].w<<"  "<<PMF[i].SDw<<"   "<<PMF_difference[i]<<"  "<<forward[i].w<<"  "<<forward[i].SDw<<"  "
			<<reverse[i].w<<"  "<<reverse[i].SDw<<"  "<<forward[i].R2<<"  "<<reverse[i].R2;
		else
			works<<endl<<forward[i].xFinal<< "   "<<PMF[i].w<<"   "<<PMF_difference[i]<<"  "<<forward[i].w<<"  "<<reverse[i].w<<"  "<<forward[i].R2<<"  "<<reverse[i].R2;
	}
	works.close();
}



double workPeakFinder 	(vector <workCount> & W, gsl_vector * c, int order)	{
	double dw = 0.0001;
	double peakWork = 0, countMax = 0;
	int nWorkBins = (int)W.size();

		int nWorks = (W[nWorkBins-1].w - W[0].w) / dw;

		for (int b=0 ; b < nWorks ; b++ )	{
			double work = W[0].w + (double) b * dw;
			double count = 0;
			for (int d = 0; d <= order ; d++)
				count += gsl_vector_get(c,(d)) * pow(work, d);
			if (countMax < count)	{	countMax = count; 	peakWork = work; }

		}
		return peakWork;

}
