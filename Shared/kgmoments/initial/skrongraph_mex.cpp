/**
 * @file skronggraph_mex.cpp
 * Implement stochastic Kronecker generation via coin-flipping.
 * @author David F. Gleich
 */

// compile with
// mex -O skrongraph_mex.cpp CXXFLAGS="\$CXXFLAGS  -Wall" -largeArrayDims

#include "mex.h"
#include "math.h"

#include <vector>

#include "krongraph_coin_flip.cpp"

#define assert(a) assert1(a,__LINE__)
#define assert1(a,b) assert2(a,b)
#define assert2(a,b) mxAssert((a),"assert failed in " __FILE__ " at " #b );

/*struct thread_input {
    unsigned int seed1;
    unsigned int seed2;
    size_t start_row;
    size_t end_row; // run in [start_row,end_row]    
};*/

/* 
 * This is the gateway routine into the function.
 * We take in a sparse matrix and an ordering of the vertices.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* argument check */
    
    if (nrhs != 3) {
        mexErrMsgTxt("Three inputs required.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments");
    }
    
    
   
    if (mxGetM(prhs[0])!=mxGetN(prhs[0]) ) {
        mexErrMsgTxt("Invalid initiator matrix.");
    }
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgTxt("Invalid recursion level.");
    }
    if (mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgTxt("Invalid thread count.");
    }
    
    double *T = mxGetPr(prhs[0]);
    size_t r = (size_t)mxGetScalar(prhs[1]);
    size_t nthreads = (size_t)mxGetScalar(prhs[2]);
    
    KronGraphCoinFlip kg(T, (size_t)mxGetM(prhs[0]), r);
    
    
    // seed the random number generator;
    mxArray *num_seeds[2];
    mxArray *rand_seeds;
    double *rseeds;
    num_seeds[0] = mxCreateScalarDouble(2);
    num_seeds[1] = mxCreateScalarDouble(1);
    mexCallMATLAB(1, &rand_seeds, 2, num_seeds, "rand");
    rseeds = mxGetPr(rand_seeds);
    kg.rand_seed((unsigned int)floor(rseeds[0]*4294967296.0),
        (unsigned int)floor(rseeds[1]*4294967296.0));
    
    kg.generate_edges();
    
    // copy data out
    mwIndex dims[2];
    dims[0] = kg.edges.size();
    dims[1] = 2;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *edata_src = (double*)mxGetData(plhs[0]);
    double *edata_dst = edata_src + kg.edges.size();
    
    for (size_t i=0; i<kg.edges.size(); ++i) {
        edata_src[i] = (double)kg.edges[i].first+1.;
        edata_dst[i] = (double)kg.edges[i].second+1.;
    }
    

}

