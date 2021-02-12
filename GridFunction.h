#pragma once

#include <vector>
#include <cmath>

using Pair = std::pair<double, double>;
using Vector_Pair = std::vector<std::pair<double, double>>;

// This object allows us to get the grids for functions used in the TG schema.
// It needs an interval which would be the domain of the functions, and a number of points to discretize this interval.
class gridFunction
{
public:
	gridFunction(Pair interval, int number_points);

	Vector_Pair getGridFunctionR();
	Vector_Pair getGridFunctionMu();
	Vector_Pair getGridFunctionSigma();

	double functionR(double psi, Vector_Pair rGrid);
	double functionMu(double psi, Vector_Pair muGrid);
	double functionSigma(double psi, Vector_Pair sigmaGrid);


protected:

	int _number_points;
	Pair _interval;
	
};