#ifndef PATHSIMULATOR2D_H
#define PATHSIMULATOR2D_H

#ifndef MODEL2D_H
#include "Model2D.h"
#endif


#ifndef SCHEMA_H
#include "Schema.h"
#endif

#include <vector>
#include <cmath>

using Vector = std::vector<double>;
using Vector_Pair = std::vector<std::pair<double, double> >;

class PathSimulator2D final
{
public:
	PathSimulator2D(Pair initial_factors, 
                    const Vector& time_points, 
                    const Model2D& model,
					const schema& schema);
	// Copy constructor, Assignement operator and Destructor are NEEDED because one of the member variable is a POINTER 
	PathSimulator2D(const PathSimulator2D& path_simulator);
	PathSimulator2D& operator=(const PathSimulator2D& path_simulator);
	~PathSimulator2D();


	Vector_Pair path() const;
	Vector getTimePoints() const;
	schema* getSchema() const;
	const Model2D* getModel() const;

private:
	// This method is internal to the class, not needed outside it, so we set it as being private
	Pair nextStep(int current_index, Pair current_factors) const; 

	Pair _initial_factors;
	Vector _time_points;
	const Model2D* _model;
	schema* _schema;

};

#endif