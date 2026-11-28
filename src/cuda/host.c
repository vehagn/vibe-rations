/* The GPU's kernels (kernels.cu) compiled for a CPU, to run a thread at a time: what
 * `Host` (src/cuda/host.rs) loads, so that the kernels and what drives them can be tested
 * against the solver without a GPU. Compile with -ffp-contract=off, or the compiler may fuse
 * products into sums, which the solver never does. */

#define __global__
#define __device__
#define __forceinline__ inline

typedef struct Dim {
    unsigned x, y, z;
} Dim;

/* Each thread that emulates a launch has its own. */
static __thread Dim blockIdx, threadIdx, blockDim;

#include "kernels.cu"

/* Runs every thread of every block of a launch, one after another. */
#define EMULATE(kernel) \
    void emulate_##kernel(const Args *a, const unsigned *grid, const unsigned *block) \
    { \
        blockDim.x = block[0]; \
        blockDim.y = block[1]; \
        blockDim.z = block[2]; \
        for (blockIdx.z = 0; blockIdx.z < grid[2]; blockIdx.z++) \
            for (blockIdx.y = 0; blockIdx.y < grid[1]; blockIdx.y++) \
                for (blockIdx.x = 0; blockIdx.x < grid[0]; blockIdx.x++) \
                    for (threadIdx.z = 0; threadIdx.z < block[2]; threadIdx.z++) \
                        for (threadIdx.y = 0; threadIdx.y < block[1]; threadIdx.y++) \
                            for (threadIdx.x = 0; threadIdx.x < block[0]; threadIdx.x++) \
                                kernel(*a); \
    }

EMULATE(velocity)
EMULATE(stress)
EMULATE(damp)
EMULATE(mirror)
EMULATE(inject)
EMULATE(gather)
EMULATE(peak)
EMULATE(layout)
