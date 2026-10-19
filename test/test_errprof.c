// Unit tests for errprof.c (error-profile statistics), no index needed.
#include <stdio.h>
#include <math.h>
#include "errprof.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(cond, msg) do { if (cond) ++n_pass; else { ++n_fail; fprintf(stderr, "FAIL: %s\n", msg); } } while (0)
#define NEAR(a, b, tol) (fabs((a) - (b)) <= (tol))

// expected counts of n draws from pi0*Pois(mu0) + (1-pi0)*Pois(mu1), last bin censored
static void make_hist(int64_t n, double pi0, double mu0, double mu1, int64_t *h)
{
	int k;
	double tail0 = rb3_err_pois_sf(mu0, RB3_ERR_CAP - 1), tail1 = rb3_err_pois_sf(mu1, RB3_ERR_CAP - 1);
	for (k = 0; k < RB3_ERR_CAP; ++k)
		h[k] = (int64_t)(n * (pi0 * rb3_err_pois_pmf(mu0, k) + (1.0 - pi0) * rb3_err_pois_pmf(mu1, k)) + .5);
	h[RB3_ERR_CAP] = (int64_t)(n * (pi0 * tail0 + (1.0 - pi0) * tail1) + .5);
}

int main(void)
{
	int64_t h[RB3_ERR_CAP + 1] = {0};
	rb3_errfit_t fit;
	int k;
	// Poisson basics (hand-computed: e^-0.3 = 0.740818, 0.3 e^-0.3 = 0.222245)
	CHECK(NEAR(rb3_err_pois_pmf(0.3, 0), 0.740818, 1e-6), "pmf(0.3, 0)");
	CHECK(NEAR(rb3_err_pois_pmf(0.3, 1), 0.222245, 1e-6), "pmf(0.3, 1)");
	CHECK(NEAR(rb3_err_pois_sf(0.3, 1), 1.0 - 0.740818 - 0.222245, 1e-6), "sf(0.3, 1) = P(X > 1)");
	CHECK(rb3_err_pois_sf(0.3, -1) == 1.0, "sf(mu, -1) = 1");
	// cutoffs: P(Pois(0.3) > 2) = 0.0036, > 3 = 0.00027 -> alpha 1e-3 gives K=3, alpha 1e-2 gives K=2
	CHECK(rb3_err_cutoff(0.3, 1e-3, RB3_ERR_CAP) == 3, "cutoff(mu=0.3, alpha=1e-3) = 3");
	CHECK(rb3_err_cutoff(0.3, 1e-2, RB3_ERR_CAP) == 2, "cutoff(mu=0.3, alpha=1e-2) = 2");
	CHECK(rb3_err_cutoff(1.5, 1e-3, RB3_ERR_CAP) == 6, "cutoff(mu=1.5, alpha=1e-3) = 6 (noisier library, higher K)");
	CHECK(rb3_err_cutoff(0.0, 1e-3, RB3_ERR_CAP) == 0, "cutoff(mu=0) = 0");
	// empty histogram
	rb3_err_fit(h, RB3_ERR_CAP, &fit);
	CHECK(fit.ok == 0 && fit.n == 0, "empty histogram: no fit");
	// pure error component (in-panel reads): mu0 recovered, pi0 ~ 1
	make_hist(1000000, 1.0, 0.15, 5.0, h);
	rb3_err_fit(h, RB3_ERR_CAP, &fit);
	// (a single Poisson is not identifiable as a mixture: EM may split it into two near-equal
	// components; what matters is that the error component and the cutoff are right)
	CHECK(fit.ok && NEAR(fit.mu0, 0.15, 0.01), "pure Poisson(0.15): mu0 recovered");
	CHECK(rb3_err_cutoff(fit.mu0, 1e-3, RB3_ERR_CAP) == rb3_err_cutoff(0.15, 1e-3, RB3_ERR_CAP), "pure Poisson(0.15): same cutoff as the true mu");
	// mixture with a divergence tail (held-out reads): 85% Pois(0.3) + 15% Pois(4)
	make_hist(1000000, 0.85, 0.3, 4.0, h);
	rb3_err_fit(h, RB3_ERR_CAP, &fit);
	CHECK(NEAR(fit.mu0, 0.3, 0.01) && NEAR(fit.mu1, 4.0, 0.1) && NEAR(fit.pi0, 0.85, 0.01), "85/15 mixture recovered");
	CHECK(rb3_err_cutoff(fit.mu0, 1e-3, RB3_ERR_CAP) == 3, "85/15 mixture: auto cutoff 3");
	// a heavy censored tail (10% of reads at >= cap) must not drag the error component
	make_hist(1000000, 0.9, 0.3, 30.0, h);
	rb3_err_fit(h, RB3_ERR_CAP, &fit);
	CHECK(NEAR(fit.mu0, 0.3, 0.01) && NEAR(fit.pi0, 0.9, 0.01), "censored heavy tail: error component unaffected");
	CHECK(fit.mu1 > (double)RB3_ERR_CAP, "censored heavy tail: tail mean beyond the cap");
	// sidecar writer: header + cap+1 rows, cumulative fraction ends at 1
	{
		FILE *fp = tmpfile();
		char line[512];
		int n_rows = 0, n_hdr = 0;
		double last = 0;
		rb3_err_write_hist(fp, "pair", h, RB3_ERR_CAP, &fit, 1e-3, 3);
		rewind(fp);
		while (fgets(line, sizeof(line), fp)) {
			if (line[0] == '#') ++n_hdr;
			else { char u[16], kk[8]; long c; double fr, cf; ++n_rows; if (sscanf(line, "%15s %7s %ld %lf %lf", u, kk, &c, &fr, &cf) == 5) last = cf; }
		}
		fclose(fp);
		CHECK(n_hdr == 1 && n_rows == RB3_ERR_CAP + 1, "write_hist: 1 header, cap+1 rows");
		CHECK(NEAR(last, 1.0, 1e-6), "write_hist: cum_frac ends at 1");
	}
	for (k = 0; k <= RB3_ERR_CAP; ++k) h[k] = 0;
	printf("test_errprof: %d passed, %d failed\n", n_pass, n_fail);
	return n_fail? 1 : 0;
}
