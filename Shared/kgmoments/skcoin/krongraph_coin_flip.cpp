/** 
 * @file krongraph_coin_flip.cc
 * @author David F. Gleich
 * Implement a Stochastic Kronecker graph generator by coin-flipping.
 */
 
#include <vector>
#include <math.h> 
#include <assert.h>
         
class KronGraphCoinFlip
{
private:
    
/* return a random float >= 0 and < 1 */
#define rand_float          ((double)myrand() / 4294967296.0)
#define rand_float_log          (log(((double)myrand() / 4294967296.0)))


    unsigned int myrand ()
    {
    /* Use any pair of non-equal numbers from this list for "a" and "b"
        18000 18030 18273 18513 18879 19074 19098 19164 19215 19584       
        19599 19950 20088 20508 20544 20664 20814 20970 21153 21243       
        21423 21723 21954 22125 22188 22293 22860 22938 22965 22974       
        23109 23124 23163 23208 23508 23520 23553 23658 23865 24114       
        24219 24660 24699 24864 24948 25023 25308 25443 26004 26088       
        26154 26550 26679 26838 27183 27258 27753 27795 27810 27834       
        27960 28320 28380 28689 28710 28794 28854 28959 28980 29013       
        29379 29889 30135 30345 30459 30714 30903 30963 31059 31083
    */
       static unsigned int rand_a = 18000, rand_b = 30903;

       SEED_X = rand_a*(SEED_X&65535) + (SEED_X>>16);
       SEED_Y = rand_b*(SEED_Y&65535) + (SEED_Y>>16);

       return ((SEED_X<<16) + (SEED_Y&65535));
    }

    unsigned int SEED_X, SEED_Y;
    void default_seed() {
        SEED_X = 521288629;
        SEED_Y = 362436069;
    }

public:
    std::vector<std::pair<int,int> > edges;
    size_t n; // size of the initiator matrix
    size_t r; // number of levels
    size_t N; // size of the full matrix N = n^r
    double total_prob; 
    double self_prob;
    std::vector<double> Tdata; // the data for the initator matrix
    double* T;
    std::vector<size_t> npowers;
    
    void rand_seed( unsigned int seed1, unsigned int seed2 )
    {
       if (seed1) SEED_X = seed1;   /* use default seeds if parameter is 0 */
       if (seed2) SEED_Y = seed2;
    }
    
    KronGraphCoinFlip() 
    : n(0), r(0), N(0), Tdata(0), npowers(0) {
        default_seed();
    }
    
    /** Reset the initialization */
    void init(const double* T_, size_t n_, size_t r_) {
        default_seed();
        
        n = n_;
        r = r_;
        
        Tdata.resize(n*n);
        std::copy(T_, T_+n*n, Tdata.begin());
        T = &Tdata[0];
        
        npowers.resize(r);
        
        // compute the big size
        N = 1;
        for (size_t i=0; i<r; i++) {
            N *= n;
        }
        
        // take logs of the data in T
        total_prob = 0.; // sum of all entries in T
        self_prob = 0.; // sum of diagonal entries in T
        for (size_t j=0; j<n; ++j) {
            for (size_t i=0; i<n; ++i) {
                total_prob += T[i+j*n];
                if (i==j) {
                    self_prob += T[i+j*n];
                }
                T[i+j*n] = log(T[i+j*n]);
            }
        }
        
        // setup some data for fast index computations
        npowers[0] = 1;
        for (size_t i=1; i<r; i++) {
            npowers[i] = npowers[i-1]*n;
        }
    }
    
    KronGraphCoinFlip(const double* T_, size_t n_, size_t r_)
    : n(0), r(0), N(0), Tdata(0), npowers(0) {
        init(T_, n_, r_);
    }
        
    
    size_t expected_edges() {
        double nedges = 0.;
        double total_edges_flipped = 1.;
        double self_edges_flipped = 1.;
        
        for (size_t i=0; i<r; i++) {
            total_edges_flipped *= total_prob;
            self_edges_flipped *= self_prob;
        }
        nedges = total_prob + self_prob;
        return (size_t)(nedges + 0.5);
    }
    
    /** Compute the log-probability of a value in an i,j cell
     * See the implementation of ind2sub in matlab
     * for a description of what's going on here.
     * @param i the row index
     * @param j the column index
     */
    double cell_probability(size_t i, size_t j) {
        double logprob = 0.;
        size_t ni = i;
        size_t nj = j;
        for (size_t ind = r; ind>0; --ind) {
            size_t k = ind-1;
            size_t v,vi,vj;
            // first index
            v = ni%npowers[k];
            vi = (ni - v)/npowers[k];
            ni = v;
            // second index
            v = nj%npowers[k];
            vj = (nj - v)/npowers[k];
            nj = v;
            assert(vi < n);
            assert(vj < n);
            logprob += T[vi+n*vj];
        }
        return logprob;
    }
    
    void generate_edges() {
        // reserve space for the edges
        edges.resize(0);
        edges.reserve(expected_edges());
        
        for (size_t i=0; i<N; ++i) {
            generate_row_edges(i);
        }
    }
    
    void generate_row_edges(size_t i) {
        for (size_t j=0; j<N; ++j) {
            // determine the coordinates of the i,j position
            // in terms of the matrix T
            double logprob = cell_probability(i,j);
            double rv = rand_float_log;
            if (rv < logprob) {
                edges.push_back(std::make_pair(i,j));
            }
        }
    }
    
    
};

