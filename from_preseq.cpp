/*
  This file contains functions from preseq preseq_v2.0.3.tar.bz2
  From this paper: https://www.nature.com/articles/nmeth.2375
  Daley, T., Smith, A. Predicting the molecular complexity of sequencing libraries. Nat Methods 10, 325–327 (2013). https://doi.org/10.1038/nmeth.2375
  
  It is being called from main function in superduper that removes cluster duplicates.
  Credit and copyright for from_preseq.* and continued_fraction should go the the above authors.

  Thorfinn 20may 2020 copenhagen, denmark.
 */

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>
#include <cassert>
#include "continued_fraction.hpp"

using namespace std;


static void
write_predicted_complexity_curve(const string outfile,
                                 const double c_level, const double step_size,
                                 const vector<double> &yield_estimates,
                                 const vector<double> &yield_lower_ci_lognormal,
                                 const vector<double> &yield_upper_ci_lognormal) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << "TOTAL_READS\tEXPECTED_DISTINCT\t"
      << "LOWER_" << c_level << "CI\t"
      << "UPPER_" << c_level << "CI" << endl;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(1);

  out << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << endl;
  for (size_t i = 0; i < yield_estimates.size(); ++i)
    out << (i + 1)*step_size << '\t'
        << yield_estimates[i] << '\t'
        << yield_lower_ci_lognormal[i] << '\t'
        << yield_upper_ci_lognormal[i] << endl;
}


static void
median_and_ci(const vector<double> &estimates,
              const double ci_level,
              double &median_estimate,
              double &lower_ci_estimate,
              double &upper_ci_estimate){
  assert(!estimates.empty());
  const double alpha = 1.0 - ci_level;
  const size_t n_est = estimates.size();
  vector<double> sorted_estimates(estimates);
  sort(sorted_estimates.begin(), sorted_estimates.end());
  median_estimate =
    gsl_stats_median_from_sorted_data(&sorted_estimates[0], 
                                      1, n_est);

  lower_ci_estimate = 
    gsl_stats_quantile_from_sorted_data(&sorted_estimates[0],
					1, n_est, alpha/2);
  upper_ci_estimate = 
    gsl_stats_quantile_from_sorted_data(&sorted_estimates[0],
					1, n_est, 1.0 - alpha/2);

}

static void
vector_median_and_ci(const vector<vector<double> > &bootstrap_estimates,
                     const double ci_level, 
                     vector<double> &yield_estimates,
                     vector<double> &lower_ci_lognormal,
                     vector<double> &upper_ci_lognormal) {

  yield_estimates.clear();
  lower_ci_lognormal.clear();
  upper_ci_lognormal.clear();
  assert(!bootstrap_estimates.empty());

  const size_t n_est = bootstrap_estimates.size();
  vector<double> estimates_row(bootstrap_estimates.size(), 0.0);
  for (size_t i = 0; i < bootstrap_estimates[0].size(); i++) {

    // estimates is in wrong order, work locally on const val
    for (size_t k = 0; k < n_est; ++k)
      estimates_row[k] = bootstrap_estimates[k][i];

    double median_estimate, lower_ci_estimate, upper_ci_estimate;
    median_and_ci(estimates_row, ci_level, median_estimate,
                  lower_ci_estimate, upper_ci_estimate);
    sort(estimates_row.begin(), estimates_row.end());

    yield_estimates.push_back(median_estimate);
    lower_ci_lognormal.push_back(lower_ci_estimate);
    upper_ci_lognormal.push_back(upper_ci_estimate);
  }
}

static bool
check_yield_estimates(const vector<double> &estimates) {

  if (estimates.empty())
    return false;

  // make sure that the estimate is increasing in the time_step and is
  // below the initial distinct per step_size
  if (!isfinite(accumulate(estimates.begin(), estimates.end(), 0.0)))
    return false;

  for (size_t i = 1; i < estimates.size(); ++i)
    if ((estimates[i] < estimates[i - 1]) ||
        (i >= 2 && (estimates[i] - estimates[i - 1] >
                    estimates[i - 1] - estimates[i - 2])) ||
        (estimates[i] < 0.0))
      return false;

  return true;
}

static double
interpolate_distinct(vector<double> &hist, size_t N,
                      size_t S, const size_t n) {
  double denom = gsl_sf_lngamma(N + 1) - gsl_sf_lngamma(n + 1) - gsl_sf_lngamma(N - n + 1);
  vector<double> numer(hist.size(), 0); 
  for (size_t i = 1; i < hist.size(); i++) {
	// N - i -n + 1 should be greater than 0
	if (N < i + n) {
	  numer[i] = 0;
	} else {
	  numer[i] = gsl_sf_lngamma(N - i + 1) - gsl_sf_lngamma(n + 1) - gsl_sf_lngamma(N - i - n + 1);
	  numer[i] = exp(numer[i] - denom) * hist[i];
	}
  }
  return S - accumulate(numer.begin(), numer.end(), 0);
}

void
resample_hist(const gsl_rng *rng, const vector<size_t> &vals_hist_distinct_counts,
              const vector<double> &distinct_counts_hist,
              vector<double> &out_hist) {

  vector<unsigned int> sample_distinct_counts_hist(distinct_counts_hist.size(), 0);

  const unsigned int distinct =
    static_cast<unsigned int>(accumulate(distinct_counts_hist.begin(),
                                         distinct_counts_hist.end(), 0.0));

  gsl_ran_multinomial(rng, distinct_counts_hist.size(), distinct,
                      &distinct_counts_hist.front(),
                      &sample_distinct_counts_hist.front());

  out_hist.clear();
  out_hist.resize(vals_hist_distinct_counts.back() + 1, 0.0);
  for(size_t i = 0; i < sample_distinct_counts_hist.size(); i++)
    out_hist[vals_hist_distinct_counts[i]] =
      static_cast<double>(sample_distinct_counts_hist[i]);
}


void
extrap_bootstrap(const bool VERBOSE, const bool DEFECTS,
		 const unsigned long int seed,
		 const vector<double> &orig_hist,
                 const size_t bootstraps, const size_t orig_max_terms,
                 const int diagonal, const double bin_step_size,
                 const double max_extrapolation, const size_t max_iter,
                 vector< vector<double> > &bootstrap_estimates) {
  // clear returning vectors
  bootstrap_estimates.clear();

  //setup rng
  srand(time(0) + getpid());
  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, seed);

  double vals_sum = 0.0;
  for(size_t i = 0; i < orig_hist.size(); i++)
    vals_sum += orig_hist[i]*i;

  const double initial_distinct 
    = accumulate(orig_hist.begin(), orig_hist.end(), 0.0);


  vector<size_t> orig_hist_distinct_counts;
  vector<double> distinct_orig_hist;
  for (size_t i = 0; i < orig_hist.size(); i++){
    if (orig_hist[i] > 0) {
      orig_hist_distinct_counts.push_back(i);
      distinct_orig_hist.push_back(orig_hist[i]);
    }
  }
  
  for (size_t iter = 0;
       (iter < max_iter && bootstrap_estimates.size() < bootstraps);
       ++iter) {

    vector<double> yield_vector;
    vector<double> hist;
    resample_hist(rng, orig_hist_distinct_counts, distinct_orig_hist, hist);

    double sample_vals_sum = 0.0;
    for(size_t i = 0; i < hist.size(); i++)
      sample_vals_sum += i*hist[i];

    //resize boot_hist to remove excess zeros
    while (hist.back() == 0)
      hist.pop_back();

    // compute complexity curve by random sampling w/out replacement
    const size_t upper_limit = static_cast<size_t>(sample_vals_sum);
	const size_t distinct = static_cast<size_t>(accumulate(hist.begin(), hist.end(), 0.0));
    const size_t step = static_cast<size_t>(bin_step_size);
    size_t sample = step;
    while(sample < upper_limit){
      yield_vector.push_back(interpolate_distinct(hist, upper_limit, distinct, sample));
      sample += step;
    }

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < hist.size() &&
           hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;
    
    size_t max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
    // refit curve for lower bound (degree of approx is 1 less than
    // max_terms)
    max_terms = max_terms - (max_terms % 2 == 1);
    
    // defect mode, simple extrapolation
    if(DEFECTS){
      vector<double> ps_coeffs;
      for (size_t j = 1; j <= max_terms; j++)
	ps_coeffs.push_back(hist[j]*std::pow((double)(-1), (int)(j + 1)) );
    
      const ContinuedFraction
	defect_cf(ps_coeffs, diagonal, max_terms);

      double sample_size = static_cast<double>(sample);
      while(sample_size < max_extrapolation){
	double t = (sample_size - sample_vals_sum)/sample_vals_sum;
	assert(t >= 0.0);
	yield_vector.push_back(initial_distinct + t*defect_cf(t));
	sample_size += bin_step_size;
      }
      // no checking of curve in defect mode
      bootstrap_estimates.push_back(yield_vector);
      if (VERBOSE) cerr << '.';
    }
    else{
      //refit curve for lower bound
      const ContinuedFractionApproximation
	lower_cfa(diagonal, max_terms);

      const ContinuedFraction
	lower_cf(lower_cfa.optimal_cont_frac_distinct(hist));

      //extrapolate the curve start
      if (lower_cf.is_valid()){
	double sample_size = static_cast<double>(sample);
	while(sample_size < max_extrapolation){
	  double t = (sample_size - sample_vals_sum)/sample_vals_sum;
	  assert(t >= 0.0);
	  yield_vector.push_back(initial_distinct + t*lower_cf(t));
	  sample_size += bin_step_size;
	}

	// SANITY CHECK
	if (check_yield_estimates(yield_vector)) {
	  bootstrap_estimates.push_back(yield_vector);
	  if (VERBOSE) cerr << '.';
	}
	else if (VERBOSE){
	  cerr << "_";
	}
      }
      else if (VERBOSE){
	cerr << "_";
      }

    }
  }
  if (VERBOSE)
    cerr << endl;
  if (bootstrap_estimates.size() < bootstraps)
    fprintf(stderr,"too many defects in the approximation, consider running in defect mode");
}


static double
GoodToulmin2xExtrap(const vector<double> &counts_hist){
  double two_fold_extrap = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    two_fold_extrap += pow(-1.0, i + 1)*counts_hist[i];

  return two_fold_extrap;
}


int lc_extrap(vector<double> &counts_hist,char *nam,double max_extrapolation, double step_size, size_t bootstraps, double c_level,size_t orig_max_terms, int DEFECTS,int VERBOSE, unsigned long int seed) {
  #if 0
    for(uint jj=0;jj<counts_hist.size();jj++)
      fprintf(stderr,"jj:%d :%f\n",jj,counts_hist[jj]);
    fprintf(stderr,"counts_hist.size(): %lu\n",counts_hist.size());
    exit(0);
#endif
  std::string outfile=std::string(nam);
  const size_t MIN_REQUIRED_COUNTS = 4;
  int diagonal = 0;


  
  // if seed is not set, make it random
  if(seed == 0){
    seed = rand();
  }

 
    size_t n_reads = 0;

 
    const size_t max_observed_count = counts_hist.size() - 1;
    const double distinct_reads = accumulate(counts_hist.begin(),
                                             counts_hist.end(), 0.0);

    // ENSURE THAT THE MAX TERMS ARE ACCEPTABLE
    size_t counts_before_first_zero = 1;
    while (counts_before_first_zero < counts_hist.size() &&
           counts_hist[counts_before_first_zero] > 0)
      ++counts_before_first_zero;

    orig_max_terms = std::min(orig_max_terms, counts_before_first_zero - 1);
    orig_max_terms = orig_max_terms - (orig_max_terms % 2 == 1);


    const size_t distinct_counts =
      static_cast<size_t>(std::count_if(counts_hist.begin(), counts_hist.end(),
                                        bind2nd(std::greater<double>(), 0.0)));
    if (VERBOSE)
      std::cerr << "TOTAL READS     = " << n_reads << endl
           << "DISTINCT READS  = " << distinct_reads << endl
           << "DISTINCT COUNTS = " << distinct_counts << endl
           << "MAX COUNT       = " << max_observed_count << endl
           << "COUNTS OF 1     = " << counts_hist[1] << endl
           << "MAX TERMS       = " << orig_max_terms << endl;

    if (VERBOSE) {
      // OUTPUT THE ORIGINAL HISTOGRAM
      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
        if (counts_hist[i] > 0)
          cerr << i << '\t' << static_cast<size_t>(counts_hist[i]) << endl;
      cerr << endl;
    }

    // check to make sure library is not overly saturated
    const double two_fold_extrap = GoodToulmin2xExtrap(counts_hist);
    if(two_fold_extrap < 0.0)
      fprintf(stderr,"Library expected to saturate in doubling of "
                              "size, unable to extrapolate");


    size_t total_reads = 0;
    for(size_t i = 0; i < counts_hist.size(); i++){
      total_reads += i*counts_hist[i];
    }
    //assert(total_reads == n_reads);

    // catch if all reads are distinct
    if (orig_max_terms < MIN_REQUIRED_COUNTS)
      fprintf(stderr,"max count before zero is les than min required "
                              "count (4), sample not sufficiently deep or "
                              "duplicates removed");

    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // ESTIMATE COMPLEXITY CURVE

    if(VERBOSE)
      cerr << "[ESTIMATING YIELD CURVE]" << endl;
    vector<double> yield_estimates;


    if (VERBOSE)
      cerr << "[BOOTSTRAPPING HISTOGRAM]" << endl;
    
    const size_t max_iter = 10*bootstraps;
    
    vector<vector <double> > bootstrap_estimates;
    extrap_bootstrap(VERBOSE, DEFECTS, seed, counts_hist, bootstraps, 
		     orig_max_terms, diagonal, step_size, max_extrapolation, 
		     max_iter, bootstrap_estimates);
    
    
    /////////////////////////////////////////////////////////////////////
    if (VERBOSE)
      cerr << "[COMPUTING CONFIDENCE INTERVALS]" << endl;
    
    // yield ci
    vector<double> yield_upper_ci_lognormal, yield_lower_ci_lognormal;
    
    vector_median_and_ci(bootstrap_estimates, c_level, yield_estimates,
			 yield_lower_ci_lognormal, yield_upper_ci_lognormal);
    
    /////////////////////////////////////////////////////////////////////
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
    
    write_predicted_complexity_curve(outfile, c_level, step_size,
				     yield_estimates, yield_lower_ci_lognormal,
				     yield_upper_ci_lognormal);
  
    fprintf(stderr,"\t-> Output in: \'%s\'\n",nam);
  return EXIT_SUCCESS;
}

