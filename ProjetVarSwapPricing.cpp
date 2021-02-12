#include <iostream>
#include <vector>

#include "MonteCarloPricer2D.h"
#include "Schema.h"
#include "FunctionFairPrice.h"
#include <time.h>

using Vector = std::vector<double>;
using Pair = std::pair<double, double>;


std::vector<double> create_discretization_time_points()
{
	Vector time_points;
	size_t number_time_points = 365.;
	double maturity = 1.;

	for (size_t time_index = 0; time_index < number_time_points; ++time_index)
	{
		double time_index_double = (double)time_index;
		double number_time_points_double = (double)number_time_points;
		time_points.push_back(time_index_double * maturity / (number_time_points_double - 1.));
	}

	return time_points;
}


HestonModel create_heston_model() {
	double drift = 0.0;
	double mean_reversion_speed = 0.5;
	double mean_reversion_level = 0.04;
	double vol_of_vol = 1.;
	double correlation = 0.5;

	HestonModel model_Heston(correlation, drift, mean_reversion_speed, mean_reversion_level, vol_of_vol);
	return model_Heston;
}


PathSimulator2D create_pathsimulator_heston_schemaQE()
{
	Pair initial_factors(10., 0.04);
	double psiC = 1.5;
	HestonModel model_Heston = create_heston_model();

	// Defines the discretization of time space
	Vector time_points = create_discretization_time_points();

	// Defines the path simulator
	schemaQE schemaqe(initial_factors,
		time_points,
		psiC,
		model_Heston);

	PathSimulator2D path_simulator_Heston(initial_factors, time_points, model_Heston, schemaqe);
	return path_simulator_Heston;
}


PathSimulator2D create_pathsimulator_heston_schemaTG()
{
	Pair initial_factors(10., 0.04);
	HestonModel model_Heston = create_heston_model();

	double alpha = 5.;
	double psi_bounded_min = 1. / (alpha * alpha);
	double psi_bounded_max = model_Heston.get_vol_of_vol() * model_Heston.get_vol_of_vol()
		/ (2. * model_Heston.get_mean_reversion_speed() * model_Heston.get_mean_reversion_level());

	Pair interval(psi_bounded_min, psi_bounded_max);
	int number_points = 500;

	// Defines the discretization of time space
	Vector time_points = create_discretization_time_points();

	// Defines the path simulator
	schemaTG schematg(initial_factors,
		time_points,
		model_Heston,
		interval,
		number_points);

	PathSimulator2D path_simulator_Heston(initial_factors, time_points, model_Heston, schematg);
	return path_simulator_Heston;
}


double get_fair_strike(schema* schema, double h, double rate) {
	const Model2D* model = schema->getModel();

	Pair initial_factors = schema->getInitialFactors();
	double drift = model->get_drift();
	double mean_reversion_speed = model->get_mean_reversion_speed();
	double mean_reversion_level = model->get_mean_reversion_level();
	double vol_of_vol = model->get_vol_of_vol();
	double correlation = model->get_correlation();

	FairPriceFunction fair_price(h, rate, *schema);
	return fair_price.getFairPrice();
}


void testing_pricer_2D()
{
	size_t number_of_simulations = 2E2;
	double h = 1E-3;
	double rate = 0.;

	PathSimulator2D path_simulator_Heston_QE = create_pathsimulator_heston_schemaQE();
	PathSimulator2D path_simulator_Heston_TG = create_pathsimulator_heston_schemaTG();

	schema* schemaQE_get = path_simulator_Heston_QE.getSchema();
	schema* schemaTG_get = path_simulator_Heston_TG.getSchema();

	double strike = get_fair_strike(schemaQE_get, h, rate);
	bool isCall = true;

	MonteCarloVarianceSwapPricer2D* pricer_Heston_SchemaQE = new MonteCarloVarianceSwapPricer2D(path_simulator_Heston_QE, number_of_simulations, rate, strike, isCall);
	MonteCarloVarianceSwapPricer2D* pricer_Heston_SchemaTG = new MonteCarloVarianceSwapPricer2D(path_simulator_Heston_TG, number_of_simulations, rate, strike, isCall);

	// Testing convergence ...
	size_t number_of_tests = 5;
	HestonModel model = create_heston_model();

	std::cout << "--------- Heston Model parameters: ---------\n";
	std::cout << "Drift: " << model.get_drift() << "\n";
	std::cout << "Mean Reversion Speed: " << model.get_mean_reversion_speed() << "\n";
	std::cout << "Mean Reversion Level: " << model.get_mean_reversion_level() << "\n";
	std::cout << "Vol of Vol: " << model.get_vol_of_vol() << "\n";
	std::cout << "Correlation: " << model.get_correlation() << "\n\n";
	std::cout << "Initial Price: " << schemaQE_get->getInitialFactors().first << "\n";
	std::cout << "Initial Volatility: " << schemaQE_get->getInitialFactors().second << "\n";
	std::cout << "--------------------------------------------\n\n";

	std::cout << "The Variance Strike calculated with the Analytical Formula is: " << strike << "\n\n";


	for (int test_index = 0; test_index < number_of_tests; ++test_index)
	{
		double pv_varianceSwapQE = pricer_Heston_SchemaQE->price();
		std::cout << "Variance Swap Return with Heston model and schema QE for test number " << test_index << " is " << pv_varianceSwapQE << "\n";
	}

	std::cout << "\n";
	for (int test_index = 0; test_index < number_of_tests; ++test_index)
	{
		double pv_varianceSwapTG = pricer_Heston_SchemaTG->price();
		std::cout << "Variance Swap with Heston model and schema TG for test number " << test_index << " is " << pv_varianceSwapTG << "\n";
	}

	std::cout << "\n";

}



int main() {
	srand(time(NULL));
	testing_pricer_2D();

	return 0;
}