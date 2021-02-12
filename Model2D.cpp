#include "Model2D.h"

Model2D::Model2D(double correlation) : _correlation(correlation) {}

double Model2D::get_correlation() const {
    return _correlation;
}

HestonModel::HestonModel(double correlation,
    double drift,
    double mean_reversion_speed,
    double mean_reversion_level,
    double vol_of_vol) :
    Model2D(correlation), _drift(drift), _mean_reversion_speed(mean_reversion_speed),
    _mean_reversion_level(mean_reversion_level), _vol_of_vol(vol_of_vol)
{}

HestonModel* HestonModel::clone() const {
    return new HestonModel(*this);
}

double HestonModel::get_drift() const {
    return _drift;
}

double HestonModel::get_mean_reversion_speed()const {
    return _mean_reversion_speed;
}

double HestonModel::get_mean_reversion_level()const {
    return _mean_reversion_level;
}

double HestonModel::get_vol_of_vol()const {
    return _vol_of_vol;
}

Pair HestonModel::drift_term(double time, Pair factors) const {
    Pair drift_pair;
    drift_pair.first = _drift * factors.first;
    drift_pair.second = _mean_reversion_speed * (_mean_reversion_level - factors.second);
    return drift_pair;
}

Pair HestonModel::volatility_term(double time, Pair factors) const {
    Pair volatility_pair;
    volatility_pair.first = sqrt(factors.second) * factors.first;
    volatility_pair.second = _vol_of_vol * sqrt(factors.second);
    return volatility_pair;
}