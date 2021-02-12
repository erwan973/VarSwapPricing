#ifndef MONTECARLOPRICER2D_H
#define MONTECARLOPRICER2D_H

#ifndef PATHSIMULATOR2D_H
#include "PathSimulator2D.h"
#endif 

class MonteCarloPricer2D
{
public:
	MonteCarloPricer2D(const PathSimulator2D& path_simulator, size_t number_of_simulations, double discount_rate);

	// Copy constructor, Assignement operator and Destructor are NEEDED because one of the member variable is a POINTER
	MonteCarloPricer2D(const MonteCarloPricer2D& pricer);
	MonteCarloPricer2D& operator=(const MonteCarloPricer2D& pricer);
	virtual ~MonteCarloPricer2D();

	virtual double path_price(const Vector_Pair& path) const = 0;
	double price() const;


protected:
	const PathSimulator2D* _path_simulator;
	size_t _number_of_simulations;
	double _discount_rate;
};

// abstract as well
class MonteCarloVarianceSwapPricer2D : public MonteCarloPricer2D
{
public:
	MonteCarloVarianceSwapPricer2D(const PathSimulator2D& path_simulator, size_t number_of_simulations, double discount_rate, double strike, bool is_call);

	double path_price(const Vector_Pair& path) const override;
protected:

	double _strike;
	bool _is_call;
};
#endif