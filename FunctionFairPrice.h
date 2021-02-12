#pragma once

#include <complex>
#include <stdlib.h>
#include <cmath>
#include "Schema.h"

using namespace std::complex_literals;

class FairPriceFunction
{
public:
	FairPriceFunction(double h, double rate, schema& schema);
	// Pointer in member variable, so copy, assignment and destructor needed
	FairPriceFunction(const FairPriceFunction& fair_price_function);
	FairPriceFunction& operator=(const FairPriceFunction& fair_price_function);
	~FairPriceFunction();

	double getFairPrice();
private:

	double _h;
	double _rate;
	const schema* _schema;

	std::complex<double> CFunction(double tau, double omega);
	std::complex<double> CFunctionPrime(double tau, double omega);
	std::complex<double> CFunctionPrimePrime(double tau, double omega);
	std::complex<double> DFunction(double tau, double omega);
	std::complex<double> DFunctionPrime(double tau, double omega);
	std::complex<double> DFunctionPrimePrime(double tau, double omega);
	double GFunction(double tau, int index);

	double getFairPriceIndex(int index);

};