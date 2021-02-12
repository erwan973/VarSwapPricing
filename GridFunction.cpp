#include <iostream>

#include "GridFunction.h"
#include "brent.h"
#define _USE_MATH_DEFINES
#include <math.h>


double densityGaussian(double x) {
	return 1. / sqrt(2. * M_PI) * exp(-x * x * 0.5);
}

double normalCDF(double x) {
	return 0.5 * erfc(-x * M_SQRT1_2);
}



gridFunction::gridFunction(Pair interval, int number_points):
	_interval(interval), _number_points(number_points)
{
}

Vector_Pair gridFunction::getGridFunctionR()
{
	Vector_Pair rGrid;
	double psi_bounded_min = _interval.first;
	double psi_bounded_max = _interval.second;

	for (int i = 0; i < _number_points; i++) {
		double psi = psi_bounded_min + (double)i * (psi_bounded_max - psi_bounded_min) / ((double)_number_points - 1);
		auto funct = [psi](double r) {
			return r * densityGaussian(r) + normalCDF(r) * (1. + r * r) - (1. + psi) *
				(densityGaussian(r) + r * normalCDF(r)) * (densityGaussian(r) + r * normalCDF(r));
		};
		rGrid.push_back({ psi , (ocl::brent_zero(psi_bounded_min - 10., psi_bounded_max + 10., 1e-6, 1e-4, funct)) });
	}

	return rGrid;
}

Vector_Pair gridFunction::getGridFunctionMu()
{
	Vector_Pair fMuGrid;
	Vector_Pair rGrid = getGridFunctionR();
	double psi_bounded_min = _interval.first;
	double psi_bounded_max = _interval.second;


	for (int i = 0; i < _number_points; i++) {
		double psi = rGrid[i].first;
		double functionR = rGrid[i].second;
		fMuGrid.push_back({ psi , functionR / (densityGaussian(functionR) + functionR * normalCDF(functionR)) });
	}

	return fMuGrid;
}

Vector_Pair gridFunction::getGridFunctionSigma()
{
	Vector_Pair fSigmaGrid;
	Vector_Pair rGrid = getGridFunctionR();
	double psi_bounded_min = _interval.first;
	double psi_bounded_max = _interval.second;


	for (int i = 0; i < _number_points; i++) {
		double psi = rGrid[i].first;
		double functionR = rGrid[i].second;
		fSigmaGrid.push_back({ psi , 1. / (sqrt(psi) * (densityGaussian(functionR) + functionR * normalCDF(functionR)) ) });
	}

	return fSigmaGrid;
}

double gridFunction::functionR(double psi, Vector_Pair rGrid)
{
	if (rGrid[0].first >= psi) {
		return rGrid[0].second;
	}
	for (int i = 1; i < rGrid.size(); i++) {
		if (rGrid[i].first >= psi) {
			double x_b = rGrid[i].first;
			double f_x_b = rGrid[i].second;
			double x_a = rGrid[i - 1].first;
			double f_x_a = rGrid[i - 1].second;
			return f_x_a + (psi - x_a) * (f_x_b - f_x_a) / (x_b - x_a);
		}
	}
	return rGrid[rGrid.size() - 1].second;
}

double gridFunction::functionMu(double psi, Vector_Pair muGrid)
{
	if (muGrid[0].first >= psi) {
		return muGrid[0].second;
	}
	for (int i = 1; i < muGrid.size(); i++) {
		if (muGrid[i].first >= psi) {
			double x_b = muGrid[i].first;
			double f_x_b = muGrid[i].second;
			double x_a = muGrid[i - 1].first;
			double f_x_a = muGrid[i - 1].second;
			return f_x_a + (psi - x_a) * (f_x_b - f_x_a) / (x_b - x_a);
		}
	}
	return muGrid[muGrid.size() - 1].second;
}

double gridFunction::functionSigma(double psi, Vector_Pair sigmaGrid)
{
	if (sigmaGrid[0].first >= psi) {
		return sigmaGrid[0].second;
	}
	for (int i = 1; i < sigmaGrid.size(); i++) {
		if (sigmaGrid[i].first >= psi) {
			double x_b = sigmaGrid[i].first;
			double f_x_b = sigmaGrid[i].second;
			double x_a = sigmaGrid[i - 1].first;
			double f_x_a = sigmaGrid[i - 1].second;
			return f_x_a + (psi - x_a) * (f_x_b - f_x_a) / (x_b - x_a);
		}
	}
	return sigmaGrid[sigmaGrid.size() - 1].second;
}