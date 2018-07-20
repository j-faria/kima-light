#ifndef DNest4_TransitModel
#define DNest4_TransitModel

#include <vector>
#include "TransitConditionalPrior.h"
#include "RJObject/RJObject.h"
#include "RNG.h"
#include "Data.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include "celerite/celerite.h"

// whether the model includes a GP component
extern const bool GP;

// whether there are observations after the change in HARPS fibers
extern const bool obs_after_HARPS_fibers;

// whether the model includes a linear trend
extern const bool trend;


/// From pysyzygy
// Models
#define QUADRATIC               0
#define KIPPING                 1
#define NONLINEAR               2
#define RIEMANN                 5
#define TRAPEZOID               6
#define SMARTINT                7
#define SLOWINT                 8
#define MDFAST                  9
#define NEWTON                  10

// Errors
#define ERR_NONE                0                                                     // We're good!
#define ERR_NOT_IMPLEMENTED     1                                                     // Function/option not yet implemented
#define ERR_MAX_PTS             2                                                     // Maximum number of points exceeded in transit. Increase settings.maxpts.
#define ERR_KEPLER              3                                                     // Error in the Kepler solver; probably didn't converge
#define ERR_NO_TRANSIT          4                                                     // The planet doesn't transit the star
#define ERR_BAD_ECC             5                                                     // Bad value for eccentricity
#define ERR_RC                  6                                                     // Error in rc() function
#define ERR_RJ                  7                                                     // Error in rj() function
#define ERR_RF                  8                                                     // Error in rf() function
#define ERR_RADIUS              9                                                     // Bad input radius
#define ERR_EXP_PTS             10                                                    // The number of exposure points cannot be odd
#define ERR_NOT_COMPUTED        11                                                    // User attempted to bin before computing
#define ERR_STAR_CROSS          12                                                    // Star-crossing orbit
#define ERR_PER                 13                                                    // Bad period
#define ERR_RHOS_ARS            14                                                    // Must specify either rhos or aRs!
#define ERR_RHOS                15                                                    // Bad rhos
#define ERR_ECC_W               16                                                    // Bad eccentricity/omega
#define ERR_LD                  17                                                    // Bad limb darkening coeffs
#define ERR_T0                  18                                                    // Bad t0

// Arrays
#define ARR_FLUX                0
#define ARR_BFLX                1
#define ARR_M                   2
#define ARR_E                   3
#define ARR_F                   4
#define ARR_R                   5
#define ARR_X                   6
#define ARR_Y                   7
#define ARR_Z                   8
#define ARR_B                   9

// Numerical
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
                  (dmaxarg1) : (dmaxarg2))
static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
                  (dminarg1) : (dminarg2))
#define RC_ERRTOL 0.04   
#define RC_TINY 1.69e-38   
#define RC_SQRTNY 1.3e-19   
#define RC_BIG 3.e37   
#define RC_TNBG (RC_TINY*RC_BIG)   
#define RC_COMP1 (2.236/RC_SQRTNY)   
#define RC_COMP2 (RC_TNBG*RC_TNBG/25.0)   
#define RC_THIRD (1.0/3.0)   
#define RC_C1 0.3   
#define RC_C2 (1.0/7.0)   
#define RC_C3 0.375   
#define RC_C4 (9.0/22.0)
#define RJ_ERRTOL 0.05   
#define RJ_TINY 2.5e-13   
#define RJ_BIG 9.0e11   
#define RJ_C1 (3.0/14.0)   
#define RJ_C2 (1.0/3.0)   
#define RJ_C3 (3.0/22.0)   
#define RJ_C4 (3.0/26.0)   
#define RJ_C5 (0.75*RJ_C3)   
#define RJ_C6 (1.5*RJ_C4)   
#define RJ_C7 (0.5*RJ_C2)   
#define RJ_C8 (RJ_C3+RJ_C3)   
#define RF_ERRTOL 0.08   
#define RF_TINY 1.5e-38   
#define RF_BIG 3.0e37   
#define RF_THIRD (1.0/3.0)   
#define RF_C1 (1.0/24.0)   
#define RF_C2 0.1   
#define RF_C3 (3.0/44.0)   
#define RF_C4 (1.0/14.0) 

// Constants
#define PI                      acos(-1.)
#define G                       6.672e-8
#define DAYSEC                  86400.
#define KEPLONGEXP              (1765.5/86400.)
#define KEPLONGCAD              (1800./86400.)
#define KEPSHRTEXP              (58.89/86400.)
#define KEPSHRTCAD              (60./86400.)
#define MAXTRANSITS             500

// Structs
typedef struct {
  double bcirc;
  double rhos;
  double MpMs;
  double esw;
  double ecw;
  double per;
  double RpRs;
  double t0;
  double ecc;
  double w;
  double aRs;
  int ntrans;
  double tN[MAXTRANSITS];
} TRANSIT;

typedef struct {
  int ldmodel;
  double u1;
  double u2;
  double q1;
  double q2;
  double c1;
  double c2;
  double c3;
  double c4;
} LIMBDARK;

typedef struct {
  int nstart;
  int nend;
  int ipts;
  int calloc;
  int balloc;
  int ialloc;
  double *time;
  double *flux;
  double *bflx;
  double *M;
  double *E;
  double *f;
  double *r;
  double *x;
  double *y;
  double *z;
  double *b;
  double *iarr;  
} ARRAYS;

typedef struct {
  double exptime;
  double keptol;
  int fullorbit;
  int maxpts;
  int exppts;
  int binmethod;
  int intmethod;
  int maxkepiter;
  int computed;
  int binned;
  int kepsolver;
} SETTINGS;


class TransitModel
{
    private:
        // Fix the number of planets? (by default, yes)
        bool fix {true};
        // Maximum number of planets
        int npmax {1};

        DNest4::RJObject<TransitConditionalPrior> planets =
            DNest4::RJObject<TransitConditionalPrior>(5, npmax, fix, TransitConditionalPrior());

        // set the LD model
        LIMBDARK ld = {
            0,    // int ldmodel,
            0.0,    // double u1,
            0.0,    // double u2,
            0.0,    // double q1,
            0.0,    // double q2,
            0.0,    // double c1,
            0.0,    // double c2,
            0.0,    // double c3,
            0.0,    // double c4,
        };

        // set the settings
        SETTINGS settings = {
            KEPLONGEXP,   // double exptime;
            1.e-15,       // double keptol;
            0,            // int fullorbit;
            10000,        // int maxpts;
            50,           // int exppts;
            RIEMANN,      // int binmethod;
            SMARTINT,     // int intmethod;
            100,          // int maxkepiter;
            0,            // int computed;
            0,            // int binned;
            NEWTON        // int kepsolver;
        };


        double background;
        //std::vector<double> offsets;
        double slope, quad;
        double fiber_offset;

        double extra_sigma;

        // Parameters for the quasi-periodic extra noise
        double eta1, eta2, eta3, eta4, eta5;
        double log_eta1, log_eta2, log_eta3, log_eta4, log_eta5;
        double a,b,c,P;

        celerite::solver::CholeskySolver<double> solver;
        Eigen::VectorXd alpha_real,
                 beta_real,
                 alpha_complex_real,
                 alpha_complex_imag,
                 beta_complex_real,
                 beta_complex_imag;

        // The signal
        std::vector<long double> mu = 
                            std::vector<long double>(Data::get_instance().N());
        void calculate_mu();

        // From pysyzygy
        // Functions
        double ellec(double k);
        double ellk(double k);
        double rc(double x, double y, int *err);
        double rj(double x, double y, double z, double p, int *err);
        double rf(double x, double y, double z, int *err);
        double sgn(double x);
        double TrueAnomaly(double E, double ecc);
        double EccentricAnomalyFast(double dMeanA, double dEcc, double tol, int maxiter);
        double EccentricAnomaly(double dMeanA, double dEcc, double tol, int maxiter);
        int Compute(TRANSIT *transit, LIMBDARK *limbdark, SETTINGS *settings, ARRAYS *arr);
        int Bin(TRANSIT *transit, LIMBDARK *limbdark, SETTINGS *settings, ARRAYS *arr);
        int Interpolate(double *t, int ipts, int array, TRANSIT *transit, LIMBDARK *limbdark, SETTINGS *settings, ARRAYS *arr);
        void dbl_free(double *ptr);


        // eccentric and true anomalies
        double ecc_anomaly(double time, double prd, double ecc, double peri_pass);
        double eps3(double e, double M, double x);
        double keplerstart3(double e, double M);
        double true_anomaly(double time, double prd, double ecc, double peri_pass);

        // The covariance matrix for the data
        Eigen::MatrixXd C {Data::get_instance().N(), Data::get_instance().N()};
        void calculate_C();

        //QPkernel *kernel;
        //HODLR_Tree<QPkernel> *A;

        unsigned int staleness;

    public:
        TransitModel();

        void save_setup();

        // Generate the point from the prior
        void from_prior(DNest4::RNG& rng);

        // Metropolis-Hastings proposals
        double perturb(DNest4::RNG& rng);

        // Likelihood function
        double log_likelihood() const;

        // Print to stream
        void print(std::ostream& out) const;

        // Return string with column information
        std::string description() const;

};

#endif

