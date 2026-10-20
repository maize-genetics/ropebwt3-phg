#include <stdlib.h>
#include <string.h>
#include "edstats.h"
#include "errprof.h"

int32_t rb3_ed_summary(const uint8_t *v, int32_t n, int32_t gbest, int32_t *n_best, int32_t *wrong)
{
	int32_t a, best = gbest >= 0 && gbest < RB3_ED_NA? gbest : RB3_ED_NA, nb = 0, w = RB3_ED_NA;
	for (a = 0; a < n; ++a)
		if (v[a] != RB3_ED_NA && v[a] < best) best = v[a];
	if (best != RB3_ED_NA)
		for (a = 0; a < n; ++a) {
			if (v[a] == best) ++nb;
			else if (v[a] != RB3_ED_NA && v[a] > best + 1 && v[a] < w) w = v[a];
		}
	if (n_best) *n_best = nb;
	if (wrong) *wrong = w;
	return best;
}

void rb3_ed_pair(const uint8_t *a, const uint8_t *b, int32_t n, uint8_t *out)
{
	int32_t k;
	for (k = 0; k < n; ++k)
		out[k] = a[k] == RB3_ED_NA || b[k] == RB3_ED_NA? RB3_ED_NA : a[k] + b[k];
}

rb3_edhist_t *rb3_edhist_init(int32_t n_asm, int32_t cap)
{
	rb3_edhist_t *h = (rb3_edhist_t*)calloc(1, sizeof(rb3_edhist_t));
	int u;
	h->n_asm = n_asm, h->cap[0] = cap, h->cap[1] = 2 * cap;
	for (u = 0; u < 2; ++u) {
		h->best[u] = (int64_t*)calloc(h->cap[u] + 2, sizeof(int64_t));
		h->asm_[u] = (int64_t*)calloc(h->cap[u] + 2, sizeof(int64_t));
		h->wrong[u] = (int64_t*)calloc(h->cap[u] + 2, sizeof(int64_t));
		h->nbest[u] = (int64_t*)calloc(n_asm + 1, sizeof(int64_t));
	}
	return h;
}

void rb3_edhist_destroy(rb3_edhist_t *h)
{
	int u;
	if (h == 0) return;
	for (u = 0; u < 2; ++u)
		free(h->best[u]), free(h->asm_[u]), free(h->wrong[u]), free(h->nbest[u]);
	free(h);
}

static inline int32_t ed_bin(int32_t x, int32_t cap) { return x == RB3_ED_NA || x > cap? cap + 1 : x; }

void rb3_edhist_add(rb3_edhist_t *h, int unit, const uint8_t *v, int32_t gbest, int flags)
{
	int32_t a, nb, w, best, cap = h->cap[unit];
	best = rb3_ed_summary(v, h->n_asm, gbest, &nb, &w);
	++h->n[unit];
	++h->best[unit][ed_bin(best, cap)];
	++h->wrong[unit][ed_bin(w, cap)];
	++h->nbest[unit][nb];
	for (a = 0; a < h->n_asm; ++a) ++h->asm_[unit][ed_bin(v[a], cap)];
	if (flags & RB3_ED_F_TRUNC) ++h->n_trunc[unit];
	if (flags & RB3_ED_F_REP) ++h->n_rep[unit];
}

static void ed_write_block(FILE *fp, const char *unit, const char *measure, const int64_t *x, int32_t nbin,
						   int na_last, const rb3_errfit_t *fit)
{
	int64_t n = 0, cum = 0;
	int32_t k, cap = nbin - 1; // fit: the last bin is censored (">= cap")
	for (k = 0; k < nbin; ++k) n += x[k];
	for (k = 0; k < nbin; ++k) {
		cum += x[k];
		fprintf(fp, "%s\t%s\t", unit, measure);
		if (na_last && k == nbin - 1) fputs("NA", fp);
		else fprintf(fp, "%d", k);
		fprintf(fp, "\t%ld\t%.6f\t%.6f", (long)x[k], n? (double)x[k] / n : 0.0, n? (double)cum / n : 0.0);
		if (fit && fit->ok) {
			double f0 = k < cap? rb3_err_pois_pmf(fit->mu0, k) : rb3_err_pois_sf(fit->mu0, cap - 1);
			double f1 = k < cap? rb3_err_pois_pmf(fit->mu1, k) : rb3_err_pois_sf(fit->mu1, cap - 1);
			double m = fit->pi0 * f0 + (1.0 - fit->pi0) * f1;
			fprintf(fp, "\t%.1f\t%.6f\n", n * fit->pi0 * f0, m > 0.0? fit->pi0 * f0 / m : 0.0);
		} else fputs("\t.\t.\n", fp);
	}
}

void rb3_edhist_write(FILE *fp, const rb3_edhist_t *h, const char *const *asm_name, double alpha)
{
	static const char *uname[2] = { "read", "pair" };
	int u;
	int32_t a;
	fprintf(fp, "#ropebwt3 refmap --ed-hist: best edit distance (substitutions + indel bases + unaligned read bases)\n");
	fprintf(fp, "#best = best over the panel; asm = best within each assembly (read x assembly); wrong = closest distinct\n");
	fprintf(fp, "#wrong copy (best among assemblies > best+1); nbest = assemblies at the best value; NA = none within the cap\n");
	if (asm_name) {
		fprintf(fp, "#assemblies=%d\t", h->n_asm);
		for (a = 0; a < h->n_asm; ++a) fprintf(fp, "%s%s", a? "," : "", asm_name[a]);
		fputc('\n', fp);
	}
	for (u = 0; u < 2; ++u) {
		rb3_errfit_t fit;
		if (h->n[u] == 0) continue;
		rb3_err_fit(h->best[u], h->cap[u] + 1, &fit);
		fprintf(fp, "#unit=%s\tn=%ld\tcap=%d\ttruncated=%ld\trepetitive_tile=%ld\tpi0=%.6f\tmu0=%.6f\tmu1=%.6f\tem_iter=%d\talpha=%g\tcutoff=%d\n",
				uname[u], (long)h->n[u], h->cap[u], (long)h->n_trunc[u], (long)h->n_rep[u], fit.pi0, fit.mu0, fit.mu1, fit.iter,
				alpha, fit.ok? rb3_err_cutoff(fit.mu0, alpha, h->cap[u] + 1) : -1);
	}
	fprintf(fp, "unit\tmeasure\tk\tcount\tfrac\tcum_frac\terr_expected\tpost_err\n");
	for (u = 0; u < 2; ++u) {
		rb3_errfit_t fit;
		if (h->n[u] == 0) continue;
		rb3_err_fit(h->best[u], h->cap[u] + 1, &fit);
		ed_write_block(fp, uname[u], "best", h->best[u], h->cap[u] + 2, 1, &fit);
		ed_write_block(fp, uname[u], "asm", h->asm_[u], h->cap[u] + 2, 1, 0);
		ed_write_block(fp, uname[u], "wrong", h->wrong[u], h->cap[u] + 2, 1, 0);
		ed_write_block(fp, uname[u], "nbest", h->nbest[u], h->n_asm + 1, 0, 0);
	}
}

void rb3_ed_npy_header(FILE *fp, int64_t rows, int64_t cols)
{
	char hdr[RB3_ED_NPY_HDR];
	int l;
	memset(hdr, ' ', sizeof(hdr));
	memcpy(hdr, "\x93NUMPY\x01\x00", 8);
	hdr[8] = (RB3_ED_NPY_HDR - 10) & 0xff, hdr[9] = (RB3_ED_NPY_HDR - 10) >> 8;
	l = snprintf(hdr + 10, RB3_ED_NPY_HDR - 10, "{'descr': '|u1', 'fortran_order': False, 'shape': (%ld, %ld), }", (long)rows, (long)cols);
	hdr[10 + l] = ' '; // snprintf's NUL
	hdr[RB3_ED_NPY_HDR - 1] = '\n';
	fwrite(hdr, 1, RB3_ED_NPY_HDR, fp);
}
