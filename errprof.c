#include <math.h>
#include "errprof.h"

double rb3_err_pois_pmf(double mu, int32_t k)
{
	if (k < 0) return 0.0;
	if (mu <= 0.0) return k == 0? 1.0 : 0.0;
	return exp(k * log(mu) - mu - lgamma(k + 1.0));
}

double rb3_err_pois_sf(double mu, int32_t k) // P(X > k)
{
	double cdf = 0.0;
	int32_t i;
	if (k < 0) return 1.0;
	for (i = 0; i <= k; ++i) cdf += rb3_err_pois_pmf(mu, i);
	return cdf >= 1.0? 0.0 : 1.0 - cdf;
}

// P(X >= c) and E[X; X >= c] = mu * P(X >= c-1) for X ~ Poisson(mu)
static inline double pois_ge(double mu, int32_t c) { return c <= 0? 1.0 : rb3_err_pois_sf(mu, c - 1); }

void rb3_err_fit(const int64_t *hist, int32_t cap, rb3_errfit_t *fit)
{
	int64_t n = 0;
	int32_t k, it;
	double s = 0.0, pi0, mu0, mu1, ll_prev = -1e300;
	fit->n = 0, fit->pi0 = 1.0, fit->mu0 = fit->mu1 = 0.0, fit->iter = 0, fit->ok = 0;
	for (k = 0; k <= cap; ++k) n += hist[k], s += (double)hist[k] * k;
	if (n == 0) return;
	fit->n = n;
	// start: error component from the 0/1 ratio, tail at the mean of the >=1 counts
	mu0 = hist[0] > 0? (double)hist[1] / hist[0] : 1.0;
	if (mu0 < 1e-4) mu0 = 1e-4;
	if (mu0 > 2.0) mu0 = 2.0;
	mu1 = n > hist[0]? s / (double)(n - hist[0]) : mu0 * 4.0;
	if (mu1 < mu0 * 2.0 + 0.5) mu1 = mu0 * 2.0 + 0.5;
	pi0 = 0.9;
	for (it = 0; it < 2000; ++it) {
		double w0 = 0.0, x0 = 0.0, x1 = 0.0, ll = 0.0;
		for (k = 0; k <= cap; ++k) {
			double f0, f1, e0, e1, r0, m;
			if (hist[k] == 0) continue;
			if (k < cap) f0 = rb3_err_pois_pmf(mu0, k), f1 = rb3_err_pois_pmf(mu1, k), e0 = k, e1 = k;
			else { // censored bin: tail probabilities and conditional means
				f0 = pois_ge(mu0, cap), f1 = pois_ge(mu1, cap);
				e0 = f0 > 0.0? mu0 * pois_ge(mu0, cap - 1) / f0 : cap;
				e1 = f1 > 0.0? mu1 * pois_ge(mu1, cap - 1) / f1 : cap;
			}
			m = pi0 * f0 + (1.0 - pi0) * f1;
			if (m <= 0.0) { r0 = mu0 < mu1? 0.0 : 1.0; m = 1e-300; } else r0 = pi0 * f0 / m;
			ll += hist[k] * log(m);
			w0 += hist[k] * r0, x0 += hist[k] * r0 * e0, x1 += hist[k] * (1.0 - r0) * e1;
		}
		pi0 = w0 / n;
		if (w0 > 0.0) mu0 = x0 / w0;
		if (n - w0 > 0.0) mu1 = x1 / (n - w0);
		if (pi0 > 1.0 - 1e-12) pi0 = 1.0 - 1e-12;
		if (mu0 < 1e-9) mu0 = 1e-9;
		fit->iter = it + 1;
		if (fabs(ll - ll_prev) < 1e-9 * (fabs(ll) + 1.0)) break;
		ll_prev = ll;
	}
	if (mu0 > mu1) { double t = mu0; mu0 = mu1, mu1 = t, pi0 = 1.0 - pi0; } // keep component 0 the lower one
	fit->pi0 = pi0, fit->mu0 = mu0, fit->mu1 = mu1, fit->ok = 1;
}

int32_t rb3_err_cutoff(double mu, double alpha, int32_t cap)
{
	int32_t k;
	for (k = 0; k < cap; ++k)
		if (rb3_err_pois_sf(mu, k) <= alpha) return k;
	return cap - 1;
}

void rb3_err_write_hist(FILE *fp, const char *unit, const int64_t *hist, int32_t cap,
						const rb3_errfit_t *fit, double alpha, int32_t cutoff)
{
	int64_t n = 0, cum = 0;
	int32_t k;
	for (k = 0; k <= cap; ++k) n += hist[k];
	fprintf(fp, "#unit=%s\tn=%ld\tpi0=%.6f\tmu0=%.6f\tmu1=%.6f\tem_iter=%d\talpha=%g\tcutoff=%d\n",
			unit, (long)n, fit->pi0, fit->mu0, fit->mu1, fit->iter, alpha, cutoff);
	for (k = 0; k <= cap; ++k) {
		double f0 = k < cap? rb3_err_pois_pmf(fit->mu0, k) : pois_ge(fit->mu0, cap);
		double f1 = k < cap? rb3_err_pois_pmf(fit->mu1, k) : pois_ge(fit->mu1, cap);
		double m = fit->pi0 * f0 + (1.0 - fit->pi0) * f1;
		cum += hist[k];
		fprintf(fp, "%s\t%d%s\t%ld\t%.6f\t%.6f\t%.1f\t%.6f\n", unit, k, k == cap? "+" : "", (long)hist[k],
				n? (double)hist[k] / n : 0.0, n? (double)cum / n : 0.0, n * fit->pi0 * f0, m > 0.0? fit->pi0 * f0 / m : 0.0);
	}
}
