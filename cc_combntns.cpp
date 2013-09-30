#include <iostream>
using std::ios;
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::flush;

#include <iomanip>
using std::setw;
using std::setfill;
using std::setprecision;
using std::setbase;
using std::hex;
using std::dec;

#include "mex.h"
#include "matrix.h"
#include "math.h"

/*stamatiad.st@gmail.com*/

typedef  unsigned long         ulong;

#define WITH_COMP
/*#define TIMING*/

#ifdef WITH_COMP
ulong *S, *P;
ulong np, kp;
#endif

ulong ct = 0;

class combination_pref
// Combinations via prefix shifts ("cool-lex" order) as delta sets.
//.
//  Algorithm as in
//  Frank Ruskey, Aaron Williams:
//    "Generating combinations by prefix shifts"
//    Lecture Notes in Computer Science, vol.3595, 2005.
//    Extended Abstract for COCOON 2005.
{
public:
    ulong *b_;  // data as delta set
    ulong s_, t_, n_;  // combination (n choose k) where n=s+t, k=t.
private:
    ulong x, y;  // aux
    
private:  // have pointer data
    combination_pref(const combination_pref&);  // forbidden
    combination_pref & operator = (const combination_pref&);  // forbidden
    
public:
    explicit combination_pref(ulong n, ulong k)
    // Must have: n>=2, k>=1  (i.e. s!=0 and t!=0)
    {
        s_ = n - k;
        t_ = k;
        n_ = s_ + t_;
        b_ = new ulong[n_];
        first();
    }
    
    virtual ~combination_pref()
    {
        delete [] b_;
    }
    
    void first()
    {
        for (ulong j=0; j<n_; ++j)  b_[j] = 0;
        for (ulong j=0; j<t_; ++j)  b_[j] = 1;
        x = 0;  y = 0;
    }
    
    bool next()
    {
        if ( x==0 )  { x=1;  b_[t_]=1;  b_[0]=0;  return true; }
        else
        {
            if ( x>=n_-1 )  return false;
            else
            {
                b_[x] = 0; ++x;  b_[y] = 1; ++y;  // X(s,t)
                if ( b_[x]==0 )
                {
                    b_[x] = 1;  b_[0] = 0;       // Y(s,t)
                    if ( y>1 )  x = 1;           // Z(s,t)
                    y = 0;
                }
                return true;
            }
        }
    }
    
    const ulong * data()  const  { return b_; }
    ulong size()  const  { return n_; }
};

void
        print_vec(const char *bla, const ulong *x, ulong n, bool dfz/*=false*/)
// Print x[0,..,n-1] as vector, n is the number of elements in the set.
// If dfz is true then Dots are printed For Zeros.
{
    if ( bla )  cout << bla;
    
    cout << "[ ";
    for (ulong k=0; k<n; ++k)
    {
        ulong t = x[k];
        if ( t!=0 )  cout << t;
        else         cout << (dfz?'.':'0');
        if ( k<n-1 )  cout << " ";
    }
    cout << " ]";
}

template <typename Type>
        static inline void  swap2(Type &x, Type &y)
// swap values
{ Type t(x);  x = y;  y = t; }

template <typename Type>
        static inline void  swap0(Type &x, Type &y)
// swap() for y known to be zero
{ y = x;  x = 0; }

template <typename Type>
        inline void reverse(Type *f, ulong n)
// Reverse order of array f.
{
    if ( n>=2 )
    {
        for (ulong k=0, i=n-1;  k<i;  ++k, --i)
            swap2(f[k], f[i]);
    }
}
// -------------------------


template <typename Type>
        inline void reverse_0(Type *f, ulong n)
// Reverse array around index zero.
{
    if ( n>2 )  reverse(f+1, n-1);
}

inline void comp2comb_nk(ulong n, ulong k, ulong &N, ulong &K)
// A composition P(n,k) of n into (at most) k parts corresponds to
// a combination B(N,K) of K=n parts from N=n+k-1 elements:
//   P(n, k)  <-->  B(N, K) == B(n+k-1, n)
{
    N = n + k - 1;
    K = n;
}
// -------------------------

inline void comb2comp_nk(ulong N, ulong K, ulong &n, ulong &k)
// A combination B(N,K) of K elements out of N
// corresponds to a composition P(n,k) of n into (at most) k parts
// where k=N-K+1 and n=K:
//   B(N, K)  <-->  P(n, k) == P(K, N-K+1)
{
    n = K;
    k = N - K + 1;
}
// -------------------------


inline void comp2comb(const ulong *p, ulong k, ulong *b)
// Convert composition P(*, k) in p[] to combination in b[]
{
    for (ulong j=0, i=0, z=0; j<k; ++j)
    {
        ulong pj = p[j];
        for (ulong w=0; w<pj; ++w)   b[i++] = z++;
        ++z;
    }
}
// -------------------------

inline void comb2comp(const ulong *b, ulong N, ulong K, ulong *p)
// Convert combination B(N, K) in b[] to composition P(*,k) in p[]
// Must have: K>=1
{
    ulong k = N-K+1;
    for (ulong z=0; z<k; ++z)  p[z] = 0;
    --k;
    ulong c1 = N;
    while ( K-- )
    {
        ulong c0 = b[K];
        ulong d = c1 - c0;
        k -= (d-1);
        ++p[k];
        c1 = c0;
    }
}
// -------------------------


inline void reverse_combination(ulong *b, ulong N, ulong K)
// Reverse order and complement values in combination b[]
// Equivalent to order reversal of the corresponding composition:
//   B <--> P  ==>  reverse_combination(B) <--> reverse(P)
{
    for (ulong z=0; z<K; ++z)  b[z] = N-1-b[z];
    reverse(b, K);
}

void visit(const combination_pref &C)
{
    const ulong n = C.size();
    const ulong *x = C.data();
    cout << setw(4) << ct << ":";
    /*print_deltaset("    ", x, n);
     * print_deltaset_as_set("    ", x, n);*/
#ifdef WITH_COMP
    const ulong k = C.t_;
    for (ulong i=0, j=0;  j<k;  ++i)  if ( x[i] )  { S[j]=i; ++j; }  // delta set to set
//    print_set("    ", S, k);
    comb2comp(S, n, k, P);
    print_vec("    ", P, kp, true);
#endif
    
    cout << endl;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     double *out;
     int i,j, outLen;
     ulong n,k;
     
     n = (ulong)(*(mxGetPr(prhs[0])));
     k = (ulong)(*(mxGetPr(prhs[1])));
      outLen = (int)(*(mxGetPr(prhs[2])));
     
    printf("Starting...\n");
    ct = 0;
    
#ifdef WITH_COMP
    S = new ulong[k];
    comb2comp_nk(n, k, np, kp);
    P = new ulong[kp];
    
    plhs[0] = mxCreateDoubleMatrix(kp,outLen,mxREAL);
    out = mxGetPr(plhs[0]);
#endif
    
    combination_pref C(n, k);
    do
    {
        ++ct;
#ifndef TIMING
        visit( C );
        for(j=0;j<kp;j++){
            out[(ct-1)*kp+j] = P[j];
        }
#endif
    }
    while ( C.next() );
    
    printf("ct = %d\n",ct);

#ifdef WITH_COMP
    delete [] P;
#endif
    
    
    return ;
}

