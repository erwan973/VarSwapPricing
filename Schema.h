#ifndef MODEL2D_H
#include "Model2D.h"
#endif
#pragma once

#include <vector>
#include <cmath>
#include "GridFunction.h"

using Vector = std::vector<double>;
using Pair = std::pair<double, double>;
using Vector_Pair = std::vector<std::pair<double, double>>;

class schema
{
public:
	schema(Pair initial_factors,
		const Vector& time_points,
		const Model2D& model);

	// Copy constructor, Assignement operator and Destructor are NEEDED because one of the member variable is a POINTER
	schema(const schema& schema);
	virtual schema* clone() const = 0;
	schema& operator=(const schema& schema);

	virtual ~schema();

	Pair getInitialFactors() const;
	Vector getTimePoints() const;
	const Model2D* getModel() const;
	virtual double nextStepVolatility(int current_index, Pair current_factors) const = 0;
	double nextStepSpot(double v_delta, int current_index, Pair current_factors) const;
protected:

	Pair _initial_factors;
	Vector _time_points;
	const Model2D* _model;

};

class schemaQE final : public schema
{
public:
	schemaQE(Pair initial_factors,
		const Vector& time_points,
		const double psiC,
		const Model2D& model);

	schemaQE* clone() const override;
	double getPsiC() const;
	double nextStepVolatility(int current_index, Pair current_factors) const override;

private:
	const double _psiC;
};

class schemaTG final : public schema
{
public:
	schemaTG(Pair initial_factors,
		const Vector& time_points,
		const Model2D& model,
		Pair interval,
		int number_points);

	schemaTG* clone() const override;
	double nextStepVolatility(int current_index, Pair current_factors) const override;

private:

	const Pair _interval;
	const int _number_points;
	Vector_Pair _gridR;
	Vector_Pair _gridMu;
	Vector_Pair _gridSigma;

};
