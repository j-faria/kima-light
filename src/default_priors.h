#include "DNest4.h"

Uniform *Cprior = new Uniform(0.9, 1.1); // normalized out-of-transit level (without unit)
Uniform *Jprior = new Uniform(0.0, 0.001); // additional white noise, m/s

Uniform *Rratprior = new Uniform(0, 1); // radius ratio (without unit)
Uniform *aRprior = new Uniform(1, 50); // a over Rstar (without unit)
Uniform *PhiTprior = new Uniform(0, 1); // Transit time in Phase (without unit)

LogUniform *Pprior = new LogUniform(1, 10); // orbital period days
