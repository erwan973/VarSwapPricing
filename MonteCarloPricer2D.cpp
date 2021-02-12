#include "MonteCarloPricer2D.h"

MonteCarloPricer2D::MonteCarloPricer2D(const PathSimulator2D & path_simulator, size_t number_of_simulations, double discount_rate)
	: _path_simulator(new PathSimulator2D(path_simulator)), _number_of_simulations(number_of_simulations), _discount_rate(discount_rate)
{
}

MonteCarloPricer2D::MonteCarloPricer2D(const MonteCarloPricer2D & pricer)
	: _path_simulator(new PathSimulator2D(*(pricer._path_simulator))), _number_of_simulations(pricer._number_of_simulations), _discount_rate(pricer._discount_rate)
{
}

MonteCarloPricer2D & MonteCarloPricer2D::operator=(const MonteCarloPricer2D & pricer)
{
	if (this == &pricer)
		return *this;
	else
	{
		delete _path_simulator;												// free the storage pointed to by _model
		_path_simulator = new PathSimulator2D(*(pricer._path_simulator));		// allocate new memory for the pointer

		// assignment for other fields
		_number_of_simulations = pricer._number_of_simulations;
		_discount_rate = pricer._discount_rate;
	}
	return *this;
}

MonteCarloPricer2D::~MonteCarloPricer2D()
{
	delete _path_simulator; 
}

double MonteCarloPricer2D::price() const
{
	double price = 0.;

	for (size_t simulation_index = 0; simulation_index < _number_of_simulations; ++simulation_index)
	{
		Vector_Pair path = _path_simulator->path();
		price += path_price(path);
	}
	price /= _number_of_simulations;
	return price;
}


MonteCarloVarianceSwapPricer2D::MonteCarloVarianceSwapPricer2D(const PathSimulator2D& path_simulator, size_t number_of_simulations, double discount_rate, double strike, bool is_call)
	: MonteCarloPricer2D(path_simulator, number_of_simulations, discount_rate), _strike(strike), _is_call(is_call)
{}

double MonteCarloVarianceSwapPricer2D::path_price(const Vector_Pair& path) const
{
	// payoff for this specific path scenario
	double volatility_at_maturity = path.at(path.size() - 1).second;
	double path_payoff = (_is_call ? volatility_at_maturity - _strike : _strike - volatility_at_maturity);

	// Discounted payoff = PV
	double maturity = _path_simulator->getTimePoints().at(_path_simulator->getTimePoints().size() - 1);
	double path_price = std::exp(-_discount_rate * maturity) * path_payoff;

	return path_price;
}
