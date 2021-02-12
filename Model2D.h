#ifndef MODEL2D_H
#define MODEL2D_H

#include <utility>
using Pair = std::pair<double, double>;

// Abstract class = it contains at least one pure virtual method [virtual return_type method_name(arguments) = 0;]
class Model2D
{
public:
	// Public methods
	// Since our base class does have data, we need to declare the destructor virtual
	virtual ~Model2D() = default;

	Model2D(double correlation);
	virtual Model2D* clone() const = 0;

	virtual Pair drift_term(double time, Pair factors) const = 0;
	virtual Pair volatility_term(double time, Pair factors) const = 0;

	//Parameters Heston Model
	virtual double get_drift() const = 0;
	virtual double get_mean_reversion_speed() const = 0;
	virtual double get_mean_reversion_level() const = 0;
	virtual double get_vol_of_vol() const = 0;


	double get_correlation() const;
protected:
	double _correlation;
};

class HestonModel final : public Model2D
{
public:
	HestonModel(double correlation,
		double drift,
		double mean_reversion_speed,
		double mean_reversion_level,
		double vol_of_vol);
	HestonModel* clone() const override;
	Pair drift_term(double time, Pair factors) const override;
	Pair volatility_term(double time, Pair factors) const override;
	double get_drift() const override;
	double get_mean_reversion_speed() const override;
	double get_mean_reversion_level() const override;
	double get_vol_of_vol() const override;
private:
	double _drift;
	double _mean_reversion_speed;
	double _mean_reversion_level;
	double _vol_of_vol;
};

#endif