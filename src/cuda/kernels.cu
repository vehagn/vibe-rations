/* The time step's kernels for an NVIDIA GPU: CUDA C, compiled when a run starts (NVRTC). They
 * are written in what C and CUDA C++ share, so the same source also compiles for a CPU
 * (host.c), where a test runs it a thread at a time against the solver.
 *
 * One thread updates one node, doing what the solver's rows do (src/solver.rs, src/fd.rs,
 * src/cpml.rs) in the same order with the same single-precision operations, so the results
 * are the solver's bit for bit as long as no product is fused into a sum (NVRTC's
 * --fmad=false; a compiler's -ffp-contract=off). Keep it so: when the solver's arithmetic
 * changes, this must follow.
 *
 * Arrays are the solver's padded ones, x fastest. STENCIL_L, the stencil's half-length, is
 * defined when compiling for a GPU, so that the loops over it unroll.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Everything a kernel reads, passed by value. Eight-byte members first, then an even number
 * of four-byte ones, so that it has no padding: Args in src/cuda.rs mirrors it. */
typedef struct Args {
    float *f[9];        /* vx vy vz sxx syy szz sxy sxz syz */
    float *m[8];        /* bx by bz lam l2m mxy mxz myz, each times dt */
    float *psi_v[3];    /* CPML memory of the velocity update along x, y, z */
    float *psi_s[3];    /* and of the stress update */
    float *profile;     /* CPML coefficients: per axis, at nodes then half nodes, a, b, 1/kappa */
    float *sponge;      /* the sponge's factors along x, then y, then z */
    int *term_field;    /* point operators' terms: the field, */
    int *term_index;    /* the padded index, */
    float *term_weight; /* the weight, */
    int *term_source;   /* and the source whose value it takes */
    float *values;      /* the sources' values: a row per step */
    float *gathered;    /* what the receivers' terms read: a row per sample */
    float *peak;        /* the largest velocity of each row of nodes */
    /* A gradient's (cuda/replay.rs): */
    float *v[3];        /* the velocities strain differentiates, or copy copies into */
    float *e[6];        /* the strain rates strain writes or adds to; what inject adds to with which 1 */
    float *ea[6];       /* the adjoint strain correlate reads */
    float *before[3];   /* the forward velocities before the step correlated, */
    float *after[3];    /* and after it */
    double *sums;       /* the correlations' sums, eight a padded node (gradient::Sums) */
    double *lit;        /* the illumination, one an interior node */
    double volume;      /* a cell's volume */
    int n[3];           /* interior nodes */
    int halo[3];
    int padded[3];
    int lo[3];          /* the layers' widths at the start */
    int hi[3];          /* and the end of each axis */
    int profile_at[3];  /* where each axis's coefficients start in profile */
    int l;              /* the stencil's half-length */
    int free_surface;
    int layers;         /* whether the CPML acts */
    int which;          /* mirror: 0 the velocities, 1 the stresses; inject: 0 into the fields, 1
                         * into e; strain: 0 writes, 1 adds */
    int from;           /* the first term to inject or gather, */
    int to;             /* and the one after the last */
    int row;            /* the step whose values to inject; the sample to gather into */
    int sources;        /* values in a row */
    int check;          /* layout: comes back if the structure arrived whole */
    int skip;           /* inject: leaves out terms whose value is zero, as the adjoint does */
    int pad;            /* keeps the four-byte members even */
    float sign;         /* +1 a step, -1 a step taken back; the mirrors' sign */
    float c[3][8];      /* the stencil's coefficients over the spacing along x, y, z */
} Args;

#ifndef STENCIL_L
#define STENCIL_L a.l
#endif

/* The derivative half a node forward of r along stride s (Stencil::forward). */
static __forceinline__ __device__ float forward(const float *f, int r, int s, const float *c, int l)
{
    float d = c[0] * (f[r + s] - f[r]);
    for (int k = 1; k < l; k++)
        d += c[k] * (f[r + (k + 1) * s] - f[r - k * s]);
    return d;
}

/* The derivative half a node backward of r (Stencil::backward). */
static __forceinline__ __device__ float backward(const float *f, int r, int s, const float *c, int l)
{
    float d = c[0] * (f[r] - f[r - s]);
    for (int k = 1; k < l; k++)
        d += c[k] * (f[r + k * s] - f[r - (k + 1) * s]);
    return d;
}

/* The slab index of node p of n along an axis with layers lo and hi wide, or -1 outside
 * them (CpmlAxis::slab). */
static __forceinline__ __device__ int slab(int p, int n, int lo, int hi)
{
    if (p < lo)
        return p;
    if (p >= n - hi)
        return lo + p - (n - hi);
    return -1;
}

/* A derivative d damped by a layer (cpml::update): p holds a, b and 1/kappa for w slab
 * nodes each, psi its memory. */
static __forceinline__ __device__ float damped(float d, float *psi, const float *p, int w, int s)
{
    float memory = p[w + s] * *psi + p[s] * d;
    *psi = memory;
    return d * p[2 * w + s] + memory;
}

/* The node (i, j, k) a thread of a sweep over the interior updates, and its padded index r;
 * a thread past the last node has none. */
#define NODE \
    int i = blockIdx.x * blockDim.x + threadIdx.x; \
    int j = blockIdx.y * blockDim.y + threadIdx.y; \
    int k = blockIdx.z * blockDim.z + threadIdx.z; \
    if (i >= a.n[0] || j >= a.n[1] || k >= a.n[2]) \
        return; \
    int r = ((k + a.halo[2]) * a.padded[1] + j + a.halo[1]) * a.padded[0] + i + a.halo[0];

/* The distance in memory between neighbours along y and z. */
#define STRIDES \
    const int sy = a.padded[0]; \
    const int sz = a.padded[0] * a.padded[1];

/* Adds sign times the velocity increment (VelocityUpdate::row). */
__global__ void velocity(Args a)
{
    NODE
    STRIDES
    const float *sxx = a.f[3], *syy = a.f[4], *szz = a.f[5];
    const float *sxy = a.f[6], *sxz = a.f[7], *syz = a.f[8];
    const int l = STENCIL_L;
    float d0 = forward(sxx, r, 1, a.c[0], l);   /* dx sxx at i+1/2 */
    float d2 = backward(sxz, r, sz, a.c[2], l); /* dz sxz at k */
    float d3 = backward(sxy, r, 1, a.c[0], l);  /* dx sxy at i */
    float d5 = backward(syz, r, sz, a.c[2], l); /* dz syz at k */
    float d6 = backward(sxz, r, 1, a.c[0], l);  /* dx sxz at i */
    float d8 = forward(szz, r, sz, a.c[2], l);  /* dz szz at k+1/2 */
    float d1 = 0.0f, d4 = 0.0f, d7 = 0.0f;
    if (a.n[1] > 1) {
        d1 = backward(sxy, r, sy, a.c[1], l); /* dy sxy at j */
        d4 = forward(syy, r, sy, a.c[1], l);  /* dy syy at j+1/2 */
        d7 = backward(syz, r, sy, a.c[1], l); /* dy syz at j */
    }
    /* The layers' memory along x holds, per plane and row, each term's slab; along y, per
     * plane and slab row, each term's row; along z, per slab plane and row, the same
     * (Outputs::row). At half nodes: dx sxx, dy syy and dz szz; the rest at nodes. */
    if (a.layers) {
        int w, s;
        const float *p;
        w = a.lo[0] + a.hi[0];
        s = slab(i, a.n[0], a.lo[0], a.hi[0]);
        if (s >= 0) {
            float *q = a.psi_v[0] + (k * a.n[1] + j) * 3 * w + s;
            p = a.profile + a.profile_at[0];
            d0 = damped(d0, q, p + 3 * w, w, s);
            d3 = damped(d3, q + w, p, w, s);
            d6 = damped(d6, q + 2 * w, p, w, s);
        }
        w = a.lo[1] + a.hi[1];
        s = slab(j, a.n[1], a.lo[1], a.hi[1]);
        if (s >= 0) {
            float *q = a.psi_v[1] + (k * w + s) * 3 * a.n[0] + i;
            p = a.profile + a.profile_at[1];
            d1 = damped(d1, q, p, w, s);
            d4 = damped(d4, q + a.n[0], p + 3 * w, w, s);
            d7 = damped(d7, q + 2 * a.n[0], p, w, s);
        }
        w = a.lo[2] + a.hi[2];
        s = slab(k, a.n[2], a.lo[2], a.hi[2]);
        if (s >= 0) {
            float *q = a.psi_v[2] + (s * a.n[1] + j) * 3 * a.n[0] + i;
            p = a.profile + a.profile_at[2];
            d2 = damped(d2, q, p, w, s);
            d5 = damped(d5, q + a.n[0], p, w, s);
            d8 = damped(d8, q + 2 * a.n[0], p + 3 * w, w, s);
        }
    }
    a.f[0][r] += a.sign * (a.m[0][r] * (d0 + d1 + d2));
    a.f[1][r] += a.sign * (a.m[1][r] * (d3 + d4 + d5));
    a.f[2][r] += a.sign * (a.m[2][r] * (d6 + d7 + d8));
}

/* Adds sign times the stress increment (StressUpdate::row). */
__global__ void stress(Args a)
{
    NODE
    STRIDES
    const float *vx = a.f[0], *vy = a.f[1], *vz = a.f[2];
    const int l = STENCIL_L;
    const int surface = a.free_surface && k == 0;
    float d0 = backward(vx, r, 1, a.c[0], l); /* dx vx at i */
    /* On the free surface the condition szz = 0 replaces dz vz. */
    float d2 = surface ? 0.0f : backward(vz, r, sz, a.c[2], l); /* dz vz at k */
    float d4 = forward(vy, r, 1, a.c[0], l);  /* dx vy at i+1/2 */
    float d5 = forward(vx, r, sz, a.c[2], l); /* dz vx at k+1/2 */
    float d7 = forward(vy, r, sz, a.c[2], l); /* dz vy at k+1/2 */
    float d6 = forward(vz, r, 1, a.c[0], l);  /* dx vz at i+1/2 */
    float d1 = 0.0f, d3 = 0.0f, d8 = 0.0f;
    if (a.n[1] > 1) {
        d1 = backward(vy, r, sy, a.c[1], l); /* dy vy at j */
        d3 = forward(vx, r, sy, a.c[1], l);  /* dy vx at j+1/2 */
        d8 = forward(vz, r, sy, a.c[1], l);  /* dy vz at j+1/2 */
    }
    /* At nodes: dx vx, dy vy and dz vz (the first along each axis); the rest at half nodes. */
    if (a.layers) {
        int w, s;
        const float *p;
        w = a.lo[0] + a.hi[0];
        s = slab(i, a.n[0], a.lo[0], a.hi[0]);
        if (s >= 0) {
            float *q = a.psi_s[0] + (k * a.n[1] + j) * 3 * w + s;
            p = a.profile + a.profile_at[0];
            d0 = damped(d0, q, p, w, s);
            d4 = damped(d4, q + w, p + 3 * w, w, s);
            d6 = damped(d6, q + 2 * w, p + 3 * w, w, s);
        }
        w = a.lo[1] + a.hi[1];
        s = slab(j, a.n[1], a.lo[1], a.hi[1]);
        if (s >= 0) {
            float *q = a.psi_s[1] + (k * w + s) * 3 * a.n[0] + i;
            p = a.profile + a.profile_at[1];
            d1 = damped(d1, q, p, w, s);
            d3 = damped(d3, q + a.n[0], p + 3 * w, w, s);
            d8 = damped(d8, q + 2 * a.n[0], p + 3 * w, w, s);
        }
        w = a.lo[2] + a.hi[2];
        s = slab(k, a.n[2], a.lo[2], a.hi[2]);
        if (s >= 0) {
            float *q = a.psi_s[2] + (s * a.n[1] + j) * 3 * a.n[0] + i;
            p = a.profile + a.profile_at[2];
            d2 = damped(d2, q, p, w, s);
            d5 = damped(d5, q + a.n[0], p + 3 * w, w, s);
            d7 = damped(d7, q + 2 * a.n[0], p + 3 * w, w, s);
        }
    }
    const float lam = a.m[3][r], l2m = a.m[4][r];
    if (surface) {
        /* With szz = 0, dz vz = -lam/(lam + 2 mu) (dx vx + dy vy). */
        const float q = l2m > 0.0f ? lam * lam / l2m : 0.0f;
        const float u = l2m - q, v = lam - q;
        a.f[3][r] += a.sign * (u * d0 + v * d1);
        a.f[4][r] += a.sign * (v * d0 + u * d1);
        a.f[5][r] = 0.0f;
    } else {
        a.f[3][r] += a.sign * (l2m * d0 + lam * (d1 + d2));
        a.f[4][r] += a.sign * (l2m * d1 + lam * (d0 + d2));
        a.f[5][r] += a.sign * (l2m * d2 + lam * (d0 + d1));
    }
    a.f[6][r] += a.sign * (a.m[5][r] * (d3 + d4));
    a.f[7][r] += a.sign * (a.m[6][r] * (d5 + d6));
    a.f[8][r] += a.sign * (a.m[7][r] * (d7 + d8));
}

/* Damps every field inside the side layers by the sponge's factors (Solver::damp_layers). */
__global__ void damp(Args a)
{
    NODE
    const float *fx = a.sponge, *fy = fx + a.n[0], *fz = fy + a.n[1];
    const float yz = fy[j] * fz[k];
    float factor;
    if (yz < 1.0f)
        factor = yz * fx[i];
    else if (i < a.lo[0] || i >= a.n[0] - a.hi[0])
        factor = fx[i];
    else
        return;
    for (int q = 0; q < 9; q++)
        a.f[q][r] *= factor;
}

/* Mirrors the velocities (which = 0), or szz, sxz and syz (1), across the free surface with
 * sign (Solver::mirror_velocity, Solver::mirror_stress): a thread for each value of a padded
 * plane and each plane above the surface. */
__global__ void mirror(Args a)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    const int m = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const int len = a.padded[0] * a.padded[1];
    const int s = a.halo[2]; /* the surface's plane */
    if (p >= len || m > s)
        return;
    /* A plane above takes the plane as far below; what sits half a node down, the one before. */
    const int above = (s - m) * len + p, on = (s + m) * len + p, half = on - len;
    if (a.which == 0) {
        a.f[0][above] = a.sign * a.f[0][on];
        a.f[1][above] = a.sign * a.f[1][on];
        a.f[2][above] = a.sign * a.f[2][half];
    } else {
        a.f[5][above] = a.sign * a.f[5][on];
        a.f[7][above] = a.sign * a.f[7][half];
        a.f[8][above] = a.sign * a.f[8][half];
    }
}

/* Adds the sources' terms from..to for step `row`, one after another as the solver does
 * (PointOp::inject): one thread, since terms may share a value. Into the fields, or with
 * which 1 into e (the adjoint strain, as gradient::PressureResiduals::add does); with skip,
 * a term whose value is zero adds nothing, as in Adjoint's injections. */
__global__ void inject(Args a)
{
    if (blockIdx.x != 0 || threadIdx.x != 0)
        return;
    float **to = a.which ? a.e : a.f;
    const float *value = a.values + a.row * a.sources;
    for (int t = a.from; t < a.to; t++) {
        const float x = value[a.term_source[t]];
        if (a.skip && x == 0.0f)
            continue;
        to[a.term_field[t]][a.term_index[t]] += a.term_weight[t] * x;
    }
}

/* Copies what the receivers' terms from..to read into row `row` of gathered; the weighted
 * sums are taken where the traces are kept. */
__global__ void gather(Args a)
{
    const int t = a.from + blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= a.to)
        return;
    a.gathered[a.row * (a.to - a.from) + t - a.from] = a.f[a.term_field[t]][a.term_index[t]];
}

/* The largest particle-velocity component of each row of nodes, or infinity if any is not
 * finite (Solver::peak_velocity): a thread per row. */
__global__ void peak(Args a)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int k = blockIdx.y * blockDim.y + threadIdx.y;
    if (j >= a.n[1] || k >= a.n[2])
        return;
    const int r = ((k + a.halo[2]) * a.padded[1] + j + a.halo[1]) * a.padded[0] + a.halo[0];
    const float largest = 3.4028235e38f;
    float best = 0.0f;
    for (int q = 0; q < 3; q++) {
        for (int i = 0; i < a.n[0]; i++) {
            float x = a.f[q][r + i];
            x = x < 0.0f ? -x : x;
            if (!(x <= largest))
                x = largest * 2.0f;
            if (x > best)
                best = x;
        }
    }
    a.peak[k * a.n[1] + j] = best;
}

/* The strain rates of the velocities v (Solver::strain_rate): dx vx, dy vy, dz vz,
 * dy vx + dx vy, dz vx + dx vz and dz vy + dy vz, each where the stress it changes sits,
 * written into e, or with which 1 added to it (Solver::add_strain_rate), each sum of two
 * derivatives rounded before it is added. Without the layers; on the free surface dz vz is
 * zero. */
__global__ void strain(Args a)
{
    NODE
    STRIDES
    const float *vx = a.v[0], *vy = a.v[1], *vz = a.v[2];
    const int l = STENCIL_L;
    float d0 = backward(vx, r, 1, a.c[0], l);
    float d2 = a.free_surface && k == 0 ? 0.0f : backward(vz, r, sz, a.c[2], l);
    float d4 = forward(vy, r, 1, a.c[0], l);
    float d5 = forward(vx, r, sz, a.c[2], l);
    float d6 = forward(vz, r, 1, a.c[0], l);
    float d7 = forward(vy, r, sz, a.c[2], l);
    float d1 = 0.0f, d3 = 0.0f, d8 = 0.0f;
    if (a.n[1] > 1) {
        d1 = backward(vy, r, sy, a.c[1], l);
        d3 = forward(vx, r, sy, a.c[1], l);
        d8 = forward(vz, r, sy, a.c[1], l);
    }
    if (a.which) {
        a.e[0][r] += d0;
        a.e[1][r] += d1;
        a.e[2][r] += d2;
        a.e[3][r] += d3 + d4;
        a.e[4][r] += d5 + d6;
        a.e[5][r] += d7 + d8;
    } else {
        a.e[0][r] = d0;
        a.e[1][r] = d1;
        a.e[2][r] = d2;
        a.e[3][r] = d3 + d4;
        a.e[4][r] = d5 + d6;
        a.e[5][r] = d7 + d8;
    }
}

/* Adds a forward step's correlations with the adjoint step that pairs with it
 * (gradient::correlate), in double precision and in its order: -W ea (x) e for the stiffness,
 * the adjoint velocities (the fields' first three) times the forward velocities' change for
 * the density, and |after|^2 to the illumination. */
__global__ void correlate(Args a)
{
    NODE
    const int surface = a.free_surface && k == 0;
    /* A cell's share is half on the free surface for what sits on it (Solver::cell_share). */
    const double w_normal = a.volume * (surface ? 0.5 : 1.0);
    const double w_xy = w_normal, w_xz = a.volume, w_yz = a.volume;
    const double xx = -w_normal * (double)a.ea[0][r];
    const double yy = -w_normal * (double)a.ea[1][r];
    const double zz = -w_normal * (double)a.ea[2][r];
    const double exx = a.e[0][r], eyy = a.e[1][r], ezz = a.e[2][r];
    double *sum = a.sums + 8LL * r;
    if (surface) {
        sum[3] += xx * exx + yy * eyy;
        sum[4] += xx * eyy + yy * exx;
    } else {
        sum[3] += xx * exx + yy * eyy + zz * ezz;
        sum[4] += xx * (eyy + ezz) + yy * (exx + ezz) + zz * (exx + eyy);
    }
    sum[5] += -w_xy * (double)a.ea[3][r] * (double)a.e[3][r];
    sum[6] += -w_xz * (double)a.ea[4][r] * (double)a.e[4][r];
    sum[7] += -w_yz * (double)a.ea[5][r] * (double)a.e[5][r];
    for (int c = 0; c < 3; c++)
        sum[c] += (double)a.f[c][r] * (double)(a.after[c][r] - a.before[c][r]);
    double lit = 0.0;
    for (int c = 0; c < 3; c++) {
        const double x = a.after[c][r];
        lit = c == 0 ? x * x : lit + x * x;
    }
    a.lit[((long long)k * a.n[1] + j) * a.n[0] + i] += lit;
}

/* Copies the velocities, the whole padded arrays, into v: a thread a value. */
__global__ void copy(Args a)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (long long)a.padded[0] * a.padded[1] * a.padded[2])
        return;
    for (int c = 0; c < 3; c++)
        a.v[c][t] = a.f[c][t];
}

/* Tells whether the structure arrived as it was sent: its size and its last members. */
__global__ void layout(Args a)
{
    if (blockIdx.x != 0 || threadIdx.x != 0)
        return;
    a.peak[0] = (float)sizeof(Args);
    a.peak[1] = (float)a.check;
    a.peak[2] = a.sign;
    a.peak[3] = a.c[2][7];
}

#ifdef __cplusplus
}
#endif
