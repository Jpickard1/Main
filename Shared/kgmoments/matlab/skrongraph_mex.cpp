/**
 * @file skronggraph_mex.cpp
 * Implement stochastic Kronecker generation via coin-flipping.
 * @author David F. Gleich
 */

// compile with
// mex -O skrongraph_mex.cpp CXXFLAGS="\$CXXFLAGS  -Wall" -largeArrayDims -I../skcoin

/** History
 * :2011-04-10: Added threaded computation
 */

#include "mex.h"
#include "math.h"

#include <vector>

#include <pthread.h>


#include "krongraph_coin_flip.cpp"

#define assert(a) assert1(a,__LINE__)
#define assert1(a,b) assert2(a,b)
#define assert2(a,b) mxAssert((a),"assert failed in " __FILE__ " at " #b );

struct kron_thread_input {
    KronGraphCoinFlip *kg;
    size_t start_row;
    size_t end_row; // run in [start_row,end_row]    
};


void* kron_edges_func(void* input) {
    kron_thread_input* args = (kron_thread_input*)input;
    
    size_t start = args->start_row;
    size_t end = args->end_row;
    
    for (size_t i=start; i < end; ++i) {
        args->kg->generate_row_edges(i);
    }
    
    return NULL;
}

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
    
    std::vector<KronGraphCoinFlip> threadkrons(nthreads);
    
    for (size_t t=0; t<nthreads; ++t) {
        // initialize the generator in each thread
        threadkrons[t].init(T, (size_t)mxGetM(prhs[0]), r);

        // randomly seed it.
        mxArray *num_seeds[2];
        mxArray *rand_seeds;
        double *rseeds;
        num_seeds[0] = mxCreateScalarDouble(2);
        num_seeds[1] = mxCreateScalarDouble(1);
        mexCallMATLAB(1, &rand_seeds, 2, num_seeds, "rand");
        rseeds = mxGetPr(rand_seeds);
        threadkrons[t].rand_seed(
            (unsigned int)floor(rseeds[0]*4294967296.0),
            (unsigned int)floor(rseeds[1]*4294967296.0));
    }
    
    if (nthreads == 1) {
       threadkrons[0].generate_edges();
    } else {
        // need to ship these off to the thread functions!
        std::vector<pthread_t> tids(nthreads);
        std::vector<kron_thread_input> tinputs(nthreads);
        
        size_t N = threadkrons[0].N;
        
        size_t currow = 0;
        
        for (size_t i=0;i<nthreads;i++){
            tinputs[i].kg = &threadkrons[i];
            tinputs[i].start_row = currow;
            // N/threads = floor(N/nthreads)
            // N%nthreads = remainder we have to allocate
            // i<N%threads will only be true for exactly N%nthreads indices i
            // e.g. remainder is 4, nthreads = 5, then it will be true unless
            // i==4, the last one.
            tinputs[i].end_row = currow + (N/nthreads + (i<N%nthreads));
            
            currow = tinputs[i].end_row;
            
            //mexPrintf("Thread %i generating rows %Zi-%Zi of %Zi\n", i+1, 
                        //tinputs[i].start_row, tinputs[i].end_row, N);
            
            pthread_create(&tids[i],NULL, kron_edges_func, (void*)&tinputs[i]);
        }
        
        // wait for threads to finish!
        assert(currow == N);
        for (size_t i=0;i<nthreads;i++){
            pthread_join(tids[i],NULL);
        }
    }
    
    size_t totaledges = 0;
    for (size_t t=0; t<nthreads; ++t) {
        //mexPrintf("Thread %i generated %Zi edges\n", t+1, threadkrons[t].edges.size());
        totaledges += threadkrons[t].edges.size();
    }
    
    // copy data out
    mwIndex dims[2];
    dims[0] = totaledges;
    dims[1] = 2;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *edata_src = (double*)mxGetData(plhs[0]);
    double *edata_dst = edata_src + totaledges;
    
    size_t curedge = 0;
    for (size_t t=0; t<nthreads; ++t) {
        for (size_t i=0; i<threadkrons[t].edges.size(); ++i) {
            assert(curedge < totaledges);
            edata_src[curedge] = (double)threadkrons[t].edges[i].first+1.;
            edata_dst[curedge] = (double)threadkrons[t].edges[i].second+1.;
            curedge ++;
        }
    }
    

}

