/** \file  stats.h  
\brief Statistics utilities for glylib.

Begun in the spring of 2008 by BLFoley and modified by heaven only
knows who since then. */

#if !defined(GLYLIB_STATS)
#define GLYLIB_STATS

/** \addtogroup ANALYSIS
 * @{
 */
typedef struct {
  char t; // type, population (p) or sample (s)
  int n; // number in sample/population
  double *d; // data points (n of these)
} statsarray;

typedef struct {
	char t; // type population (p) or sample (s)
	int n; // number of units in sample
	double m; // mean
	double v; // variance
	double s; // standard deviation
} meanvar;

typedef struct {
	int k; // k intervals used for autocorrelation function
	double *a; // a(k) the (estimate of) the autocorrelation function
} autocorr;

meanvar get_meanvar_array(statsarray S);
meanvar zero_meanvar();
statsarray zero_statsarray(); // if only one allocation
statsarray init_statsarray(); // for dynamic allocations
autocorr zero_autocorr(); // if only one allocation
autocorr init_autocorr(); // for dynamic allocations
autocorr get_autocorr_est_array(statsarray S,meanvar M);
/** @}*/

#endif
