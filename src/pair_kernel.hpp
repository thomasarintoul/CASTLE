/// AI GENERATED CONTENT - DO NOT USE - TAKE INSPIRATION ONLY



#pragma once
#include <cmath>

// When nvcc compiles this, functions must be callable from the GPU too.
// When g++ compiles it, the macro is empty.
#if defined(__CUDACC__)
  #define SF_HD __host__ __device__
#else
  #define SF_HD
#endif

constexpr int MAX_FIELDS = 16;
constexpr int MAX_OBS    = 64;

enum Stat { STAT_MAG = 0, STAT_PARR = 1, STAT_PERP = 2 };

// One thing to compute, e.g. {field 0 (velocity), STAT_PARR, order 2}
struct Observable { int field; int stat; int order; };

// Everything the maths needs, bundled up so it can be passed to a GPU kernel in one go.
struct Problem {
    long long n_points;
    int n_bins;
    int n_obs;
    const double* edges;   // n_bins + 1 bin edges
    const double* data;    // one big array, laid out below
    long long field_offset[MAX_FIELDS];
    int       n_components[MAX_FIELDS];   // 1 = scalar, 3 = vector
    Observable obs[MAX_OBS];
};
// data layout:  x[0..N) y[0..N) z[0..N)  then field 0, then field 1, ...
// A 3-vector field is stored as its own x[N] y[N] z[N]. A scalar field is just [N].
//
// Histogram layout (the output): slot 0 = pair counts, slot o+1 = observable o.
// Each slot is n_bins long, so the total is (n_obs + 1) * n_bins numbers.

SF_HD inline int find_bin(double r, const double* edges, int n_bins) {
    if (r < edges[0] || r >= edges[n_bins]) return -1;      // outside the range
    int lo = 0, hi = n_bins;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (r >= edges[mid]) lo = mid; else hi = mid;
    }
    return lo;
}

SF_HD inline double int_pow(double x, int p) {
    double r = 1.0;
    for (int k = 0; k < p; ++k) r *= x;
    return r;
}

// The value of one observable for the pair (i, j).
SF_HD inline double observable_value(const Problem& p, const Observable& o,
                                     long long i, long long j,
                                     double rx, double ry, double rz, double r) {
    const long long N = p.n_points;
    const double* f = p.data + p.field_offset[o.field];

    if (p.n_components[o.field] == 1)                        // scalar field
        return int_pow(fabs(f[i] - f[j]), o.order);

    double dx = f[i] - f[j];                                 // vector field
    double dy = f[N + i] - f[N + j];
    double dz = f[2 * N + i] - f[2 * N + j];
    double mag2 = dx*dx + dy*dy + dz*dz;
    double along = (dx*rx + dy*ry + dz*rz) / r;              // component along the separation

    double v;
    if      (o.stat == STAT_MAG)  v = sqrt(mag2);
    else if (o.stat == STAT_PARR) v = fabs(along);
    else                          v = sqrt(fmax(mag2 - along*along, 0.0));
    return int_pow(v, o.order);
}

// Process one pair. `add` says how to add to the histogram (see below).
template <class AddFn>
SF_HD inline void process_pair(const Problem& p, long long i, long long j,
                               double* hist, AddFn add) {
    const long long N = p.n_points;
    double rx = p.data[i] - p.data[j];
    double ry = p.data[N + i] - p.data[N + j];
    double rz = p.data[2 * N + i] - p.data[2 * N + j];
    double r  = sqrt(rx*rx + ry*ry + rz*rz);

    int b = find_bin(r, p.edges, p.n_bins);
    if (b < 0) return;

    add(&hist[b], 1.0);                                      // slot 0: counts
    for (int o = 0; o < p.n_obs; ++o)
        add(&hist[(long long)(o + 1) * p.n_bins + b],
            observable_value(p, p.obs[o], i, j, rx, ry, rz, r));
}