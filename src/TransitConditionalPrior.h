#ifndef DNest4_TransitConditionalPrior
#define DNest4_TransitConditionalPrior

#include "RNG.h"
#include "RJObject/ConditionalPriors/ConditionalPrior.h"

// whether the model includes hyper-priors 
// for the orbital period and semi-amplitude
extern const bool hyperpriors;

class TransitConditionalPrior:public DNest4::ConditionalPrior
{
	private:
		// Parameters of bi-exponential hyper-distribution for log-periods
		double center, width;

		// Mean of exponential hyper-distribution for semi-amplitudes
		double muK;

		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		TransitConditionalPrior();

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
		static const int weight_parameter = 1;

};

#endif

