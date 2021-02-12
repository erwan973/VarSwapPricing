#include "PathSimulator2D.h"
#include "RandomNormalGenerator.h"

PathSimulator2D::PathSimulator2D(Pair initial_factors, 
                const Vector& time_points, 
                const Model2D& model,
                const schema& schema):
    _initial_factors(initial_factors), _time_points(time_points), _model(model.clone()), _schema(schema.clone())
{
}

PathSimulator2D::PathSimulator2D(const PathSimulator2D& path_simulator):
    _initial_factors(path_simulator._initial_factors), _time_points(path_simulator._time_points),
    _model(path_simulator._model->clone()), _schema(path_simulator._schema->clone())
{}

// P2 = P1 equivalent to P2.operator=(P1)
PathSimulator2D& PathSimulator2D::operator=(const PathSimulator2D& path_simulator){
    // check for "self assignment" and do nothing in that case
	if (!(this == &path_simulator)){
		delete _model;								// free the storage pointed to by _model
		_model = path_simulator._model->clone();	// allocate new memory for the pointer

		// assignment for other fields
		_initial_factors = path_simulator._initial_factors;
		_time_points = path_simulator._time_points;
    }
    return *this;								 // return this PathSimulator2D
}

PathSimulator2D::~PathSimulator2D(){
    delete _model;
    delete _schema;
}

Vector_Pair PathSimulator2D::path() const
{
	Vector_Pair path2D{ _initial_factors };

	for (int index = 0; index < _time_points.size() - 1; ++index)
	{
		path2D.push_back(nextStep(index, path2D[index]));
	}

	return path2D;
}

Vector PathSimulator2D::getTimePoints() const
{
	return _time_points;
}

schema* PathSimulator2D::getSchema() const
{
	return _schema;
}

const Model2D* PathSimulator2D::getModel() const
{
    return _model;
}

Pair PathSimulator2D::nextStep(int current_index,
    Pair current_factors) const {

    double cur_time = _time_points[current_index];
    double time_gap = _time_points[current_index + 1] - cur_time;

    Pair nextStep;

    nextStep.second = _schema->nextStepVolatility(current_index, current_factors);
    nextStep.first = _schema->nextStepSpot(nextStep.second, current_index, current_factors);

    return nextStep;


}