#include "TransitConditionalPrior.h"
#include "DNest4.h"
#include "Utils.h"
#include <cmath>
#include <typeinfo>

using namespace std;
using namespace DNest4;


extern ContinuousDistribution *Pprior; // Orbital period
extern ContinuousDistribution *Rratprior; // radius ratio
extern ContinuousDistribution *aRprior; // a over Rstar
extern ContinuousDistribution *PhiTprior; // Transit time in Phase


TransitConditionalPrior::TransitConditionalPrior()
{
}


void TransitConditionalPrior::from_prior(RNG& rng)
{}

double TransitConditionalPrior::perturb_hyperparameters(RNG& rng)
{
    double logH = 0.;
    return logH;
}

// vec[0] = period
// vec[1] = Rp/R*
// vec[2] = a/R*
// vec[3] = phi (phase of transit)

double TransitConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    return Pprior->log_pdf(vec[0]) +
           Rratprior->log_pdf(vec[1]) +
           aRprior->log_pdf(vec[2]) +
           PhiTprior->log_pdf(vec[3]);
}

void TransitConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    vec[0] = Pprior->cdf_inverse(vec[0]);
    vec[1] = Rratprior->cdf_inverse(vec[1]);
    vec[2] = aRprior->cdf_inverse(vec[2]);
    vec[3] = PhiTprior->cdf_inverse(vec[3]);
}

void TransitConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    vec[0] = Pprior->cdf(vec[0]);
    vec[1] = Rratprior->cdf(vec[1]);
    vec[2] = aRprior->cdf(vec[2]);
    vec[3] = PhiTprior->cdf(vec[3]);
}

void TransitConditionalPrior::print(std::ostream& out) const
{}
