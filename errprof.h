#ifndef RB3_ERRPROF_H
#define RB3_ERRPROF_H

// Empirical sequencing-error profile (refmap --err-max-events). Pure statistics on a histogram
// of per-read (or per-pair) mismatch "events" -- no index, no I/O beyond the sidecar writer -- so
// it can be unit tested in isolation. The events themselves are counted in search.c
// (refmap_err_events: greedy edit-tolerant FM-index extension of the read's exact core).
//
// Model: events ~ pi0 * Poisson(mu0) + (1 - pi0) * Poisson(mu1), mu0 < mu1. Component 0 is the
// sequencing-error floor (a read matching its carrier, plus its own errors); component 1 is the
// divergence tail (reads from sequence no indexed assembly carries exactly). The last histogram
// bin is censored ("cap or more"). The cutoff is the smallest K with P(Poisson(mu0) > K) <= alpha:
// a read/pair with more than K events is unlikely to be explained by sequencing error alone.

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RB3_ERR_CAP 20 // events counted per read (and per pair); the histogram's last bin is ">= cap"

typedef struct {
	int64_t n;       // observations fitted
	double pi0;      // weight of the sequencing-error component
	double mu0, mu1; // Poisson means: sequencing error, divergence tail
	int32_t iter;    // EM iterations run
	int ok;          // 1 if the fit converged on a usable histogram (n > 0)
} rb3_errfit_t;

// Poisson pmf and upper tail P(X > k).
double rb3_err_pois_pmf(double mu, int32_t k);
double rb3_err_pois_sf(double mu, int32_t k);

// Fit the two-component mixture to hist[0..cap] (hist[cap] = count of ">= cap") by EM.
void rb3_err_fit(const int64_t *hist, int32_t cap, rb3_errfit_t *fit);

// Smallest K in [0, cap) with P(Poisson(mu) > K) <= alpha (cap-1 if none).
int32_t rb3_err_cutoff(double mu, double alpha, int32_t cap);

// Sidecar TSV block for one histogram: a '#' header with the fit and cutoff, then one row per
// bin: unit k count frac cum_frac err_expected post_err (err_expected = n*pi0*Pois(k;mu0);
// post_err = the posterior probability that a count of k came from the error component).
void rb3_err_write_hist(FILE *fp, const char *unit, const int64_t *hist, int32_t cap,
						const rb3_errfit_t *fit, double alpha, int32_t cutoff);

#ifdef __cplusplus
}
#endif

#endif
