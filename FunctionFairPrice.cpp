#include "FunctionFairPrice.h"

FairPriceFunction::FairPriceFunction(double h, double rate, schema& schema)
{
	_schema = schema.clone();
	_h = h;
	_rate = rate;
}

FairPriceFunction::FairPriceFunction(const FairPriceFunction& fair_price_function)
{
	_schema = fair_price_function._schema->clone();
	_h = fair_price_function._h;
	_rate = fair_price_function._rate;
}

FairPriceFunction& FairPriceFunction::operator=(const FairPriceFunction& fair_price_function)
{
	// check for "self assignment" and do nothing in that case
	if (!(this == &fair_price_function)) {
		delete _schema;
		_schema = fair_price_function._schema->clone();
		_h = fair_price_function._h;
		_rate = fair_price_function._rate;
	}
	return *this;
}

FairPriceFunction::~FairPriceFunction()
{
	delete _schema;
}

std::complex<double> FairPriceFunction::CFunction(double tau, double omega)
{
	const Model2D* model = _schema->getModel();
	std::complex<double> a = model->get_mean_reversion_speed() - model->get_correlation() * model->get_vol_of_vol() * 1.0i * omega;
	std::complex<double> b = sqrt(a * a + model->get_vol_of_vol() * model->get_vol_of_vol() * (1.0i * omega + omega * omega));
	std::complex<double> g = (a - b) / (a + b);
	return (model->get_drift() * 1.0i * omega - _rate) * tau + model->get_mean_reversion_speed() * model->get_mean_reversion_level()
		/ (model->get_vol_of_vol() * model->get_vol_of_vol()) * ((a - b) * tau - 2. * log((1. - g * exp(-b * tau)) / (1. - g)));
}

std::complex<double> FairPriceFunction::CFunctionPrime(double tau, double omega)
{
	return (CFunction(tau, omega + _h) - CFunction(tau, omega)) / _h;
}

std::complex<double> FairPriceFunction::CFunctionPrimePrime(double tau, double omega)
{
	return (CFunction(tau, omega + _h) - 2. * CFunction(tau, omega) + CFunction(tau, omega - _h)) / (_h * _h);
}

std::complex<double> FairPriceFunction::DFunction(double tau, double omega)
{
	const Model2D* model = _schema->getModel();
	std::complex<double> a = model->get_mean_reversion_speed() - model->get_correlation() * model->get_vol_of_vol() * 1.0i * omega;
	std::complex<double> b = sqrt(a * a + model->get_vol_of_vol() * model->get_vol_of_vol() * (1.0i * omega + omega * omega));
	std::complex<double> g = (a - b) / (a + b);
	return (a - b) / (model->get_vol_of_vol() * model->get_vol_of_vol()) * ((1. - exp(-b * tau)) / (1. - g * exp(-b * tau)));
}

std::complex<double> FairPriceFunction::DFunctionPrime(double tau, double omega)
{
	return (DFunction(tau, omega + _h) - DFunction(tau, omega)) / _h;
}

std::complex<double> FairPriceFunction::DFunctionPrimePrime(double tau, double omega)
{
	return (DFunction(tau, omega + _h) - 2. * DFunction(tau, omega) + DFunction(tau, omega - _h)) / (_h * _h);
}

// We added the index parameter, because this function needs t_i and t_(i-1), so it's implicit and we had to add it for completeness
double FairPriceFunction::GFunction(double tau, int index)
{
	if (index == 0) return 0;
	double delta = _schema->getTimePoints()[index] - _schema->getTimePoints()[index - 1];
	std::complex<double> firstPart = DFunctionPrime(delta, 0.) * DFunctionPrime(delta, 0.) * tau * tau;
	std::complex<double> secondPart = (2. * CFunctionPrime(delta, 0.) * DFunctionPrime(delta, 0.) + DFunctionPrimePrime(delta, 0.)) * tau;
	std::complex<double> thirdPart = CFunctionPrime(delta, 0.) * CFunctionPrime(delta, 0.) + CFunctionPrimePrime(delta, 0.);
	return (firstPart + secondPart + thirdPart).real();
}

// Function to compute the expectation of the log return at each step, to make sure the computation is right at each step
double FairPriceFunction::getFairPriceIndex(int index) {
	const Model2D* model = _schema->getModel();
	double s0 = _schema->getInitialFactors().first;
	double v0 = _schema->getInitialFactors().second;
	Vector timePoints = _schema->getTimePoints();
	if (index < 1) return 0.;
	else if (index == 1) return GFunction(v0, 1);
	else {
		double T_I = timePoints[index];
		double T_IMinusOne = timePoints[index - 1];
		double delta = T_I - T_IMinusOne;
		std::complex<double> qTilde = 2. * model->get_mean_reversion_speed() * model->get_mean_reversion_level()
			/ (model->get_vol_of_vol() * model->get_vol_of_vol());
		std::complex<double> ci = 2. * model->get_mean_reversion_speed() / (model->get_vol_of_vol() * model->get_vol_of_vol()
			* (1. - exp(-model->get_mean_reversion_speed() * T_IMinusOne)));
		std::complex<double> Wi = ci * v0 * exp(-model->get_mean_reversion_speed() * T_IMinusOne);
		std::complex<double> price = -DFunctionPrime(delta, 0.) * DFunctionPrime(delta, 0.)
			* (qTilde + 2. * Wi + (qTilde + Wi) * (qTilde + Wi)) / (ci * ci)
			- (2. * CFunctionPrime(delta, 0.) * DFunctionPrime(delta, 0.)
				+ DFunctionPrimePrime(delta, 0.)) * (qTilde + Wi) / ci
			- (CFunctionPrime(delta, 0.) * CFunctionPrime(delta, 0.)
				+ CFunctionPrimePrime(delta, 0.));
		return price.real();
	}
}

// Final function to calculate the fair price
double FairPriceFunction::getFairPrice()
{
	Vector timePoints = _schema->getTimePoints();
	double strikeCalc = 0;
	for (int i = 0; i < timePoints.size(); i++) {
		strikeCalc += getFairPriceIndex(i);
	}
	return strikeCalc / timePoints[timePoints.size() - 1];
}
