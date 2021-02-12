#include "Schema.h"
#include "RandomNormalGenerator.h"
#include "GridFunction.h"

schema::schema(Pair initial_factors,
    const Vector& time_points,
    const Model2D& model) :
    _initial_factors(initial_factors), _time_points(time_points), _model(model.clone())
{
}

schema::schema(const schema& schema_ex) :
    _initial_factors(schema_ex._initial_factors), _time_points(schema_ex._time_points),
    _model(schema_ex._model->clone())
{
}

// P2 = P1 equivalent to P2.operator=(P1)
schema& schema::operator=(const schema& schema_ex) {
    // check for "self assignment" and do nothing in that case
    if (!(this == &schema_ex)) {
        delete _model;
        _model = schema_ex._model->clone();	// allocate new memory for the pointer

        // assignment for other fields
        _initial_factors = schema_ex._initial_factors;
        _time_points = schema_ex._time_points;
    }
    return *this;
}

schema::~schema() {
    delete _model;
}

Pair schema::getInitialFactors() const
{
    return _initial_factors;
}

Vector schema::getTimePoints() const
{
    return _time_points;
}

const Model2D* schema::getModel() const
{
    return _model;
}

// TODO: Enhance the method (trapeze method ?)
double schema::nextStepSpot(double v_delta, int current_index,
    Pair current_factors) const {
    double randomNormal = RandomNormalGenerator::normalRandom();
    double cur_time = _time_points[current_index];
    double time_gap = _time_points[current_index + 1] - cur_time;
    double v = current_factors.second;
    double log_spot = log(current_factors.first);


    double intApproximationTime = time_gap * (v + v_delta) * 0.5;
    double intApproximationBrown = sqrt(time_gap * (v + v_delta) * 0.5) * randomNormal;
    double log_spot_delta = log_spot + _model->get_correlation() / _model->get_vol_of_vol()
        * (v_delta - v - _model->get_mean_reversion_speed() * _model->get_mean_reversion_level() * time_gap) +
        (_model->get_mean_reversion_speed() * _model->get_correlation() / _model->get_vol_of_vol() - 0.5) * intApproximationTime
        + sqrt(1. - _model->get_correlation() * _model->get_correlation()) * intApproximationBrown;
    return exp(log_spot_delta);
}

schemaQE::schemaQE(Pair initial_factors,
    const Vector& time_points,
    const double psiC,
    const Model2D& model) :
    schema(initial_factors, time_points, model), _psiC(psiC)
{
}

schemaQE* schemaQE::clone() const
{
    return new schemaQE(*this);
}

double schemaQE::getPsiC() const
{
	return _psiC;
}

double schemaQE::nextStepVolatility(int current_index, Pair current_factors) const {
    double cur_time = _time_points[current_index];
    double time_gap = _time_points[current_index + 1] - cur_time;
    double v_hat = current_factors.second;

    //We have two independent normal random variables N(0,1)
    double randomNormal = RandomNormalGenerator::normalRandom();
   
    double m = _model->get_mean_reversion_level() + (v_hat - _model->get_mean_reversion_level()) * exp(-_model->get_mean_reversion_speed() * time_gap);
    double s_square = ((v_hat * _model->get_vol_of_vol() * _model->get_vol_of_vol() * exp(-_model->get_mean_reversion_speed() * time_gap)) / _model->get_mean_reversion_speed())
        * (1. - exp(-_model->get_mean_reversion_speed() * time_gap)) +
        (_model->get_mean_reversion_level() * _model->get_vol_of_vol() * _model->get_vol_of_vol() *
            (1. - exp(-_model->get_mean_reversion_speed() * time_gap)) * (1 - exp(-_model->get_mean_reversion_speed() * time_gap)))
        / (2. * _model->get_mean_reversion_speed());

    double psi = s_square / (m * m);
    double psiInv = 1. / psi;
    double uV = RandomNormalGenerator::uniformRandom();
    double nextStep;

    if (psi <= _psiC) {
        double b_square = 2. * psiInv - 1. + sqrt(2. * psiInv) * sqrt(2. * psiInv - 1.);
        double b = sqrt(b_square);
        double a = m / (1. + b_square);
        nextStep = a * (b + randomNormal) * (b + randomNormal);
        return nextStep;
    }
    else {
        double p = (psi - 1.) / (psi + 1.);
        double beta = (1. - p) / m;

        if (p >= uV && uV >= 0.) {
            nextStep = 0.;
            return nextStep;
        }
        nextStep = (1. / beta) * log((1. - p) / (1. - uV));
        return nextStep;
    }
}

schemaTG::schemaTG(Pair initial_factors, const Vector& time_points, const Model2D& model, Pair interval, int number_points) :
    schema(initial_factors, time_points, model), _interval(interval), _number_points(number_points)
{
    gridFunction gridObj(interval, number_points);
    _gridR = gridObj.getGridFunctionR();
    _gridMu = gridObj.getGridFunctionMu();
    _gridSigma = gridObj.getGridFunctionSigma();
}

schemaTG* schemaTG::clone() const
{
    return new schemaTG(*this);
}

double schemaTG::nextStepVolatility(int current_index, Pair current_factors) const
{
    double cur_time = _time_points[current_index];
    double time_gap = _time_points[current_index + 1] - cur_time;
    double randomNormal = RandomNormalGenerator::normalRandom();
    double v_hat = current_factors.second;
    gridFunction gridFunc(_interval, _number_points);

    double m = _model->get_mean_reversion_level() + (v_hat - _model->get_mean_reversion_level()) * exp(-_model->get_mean_reversion_speed() * time_gap);
    double s_square = ((v_hat * _model->get_vol_of_vol() * _model->get_vol_of_vol() * exp(-_model->get_mean_reversion_speed() * time_gap)) / _model->get_mean_reversion_speed())
        * (1. - exp(-_model->get_mean_reversion_speed() * time_gap)) +
        (_model->get_mean_reversion_level() * _model->get_vol_of_vol() * _model->get_vol_of_vol() *
        (1. - exp(-_model->get_mean_reversion_speed() * time_gap)) * (1 - exp(-_model->get_mean_reversion_speed() * time_gap)))
        / (2. * _model->get_mean_reversion_speed());

    double psi = s_square / (m * m);

    double fMu = gridFunc.functionMu(psi, _gridMu);
    double fSigma = gridFunc.functionSigma(psi, _gridSigma);

    double mu = fMu * m;
    double sigma = fSigma * sqrt(s_square);
    double v_hat_delta = mu + sigma * randomNormal;
    if (v_hat_delta < 0.) v_hat_delta = 0.;
    return v_hat_delta;
}

// TODO: Cache exponentials to gain speed, and check the speed of the TG Schema