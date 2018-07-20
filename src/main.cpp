#include <iostream>
#include <typeinfo>
#include "DNest4.h"
#include "Data.h"
#include "TransitModel.h"
#include "TransitConditionalPrior.h"

using namespace std;
using namespace DNest4;

/* edit from here on */

#include "default_priors.h"

const bool obs_after_HARPS_fibers = false;
const bool GP = true;
const bool hyperpriors = false;
const bool trend = false;

TransitModel::TransitModel()
:planets(5, 0, true, TransitConditionalPrior())
,mu(Data::get_instance().get_t().size())
,C(Data::get_instance().get_t().size(), Data::get_instance().get_t().size())
{
    auto data = Data::get_instance();
    // Cprior = new Uniform(ymin, ymax);

    // save the current model for further analysis
    // save_setup();
}

int main(int argc, char** argv)
{
    /* set the RV data file */
    char* datafile = "data/ktwo228801451c102_lpd_LC.txt";

    /* load the file (RVs are in km/s) */
    /* don't skip any lines in the header */
	Data::get_instance().load(datafile, "ms", 7);

    // set the sampler and run it!
	Sampler<TransitModel> sampler = setup<TransitModel>(argc, argv);
	// sampler.run();

	return 0;
}
