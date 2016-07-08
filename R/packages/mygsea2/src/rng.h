/*
rng.h - random number generators

Sun Jan  9 18:34:13 CET 2005
AW: We use modifications of the Mersenne Twister RNG implementation
by the authors (see below).

The changes are:
- make the state into separate data structure; generating the
  numbers require the struct of the state.
- use C99 stdint.h to ensure portability (hopefully)
- create this header, with inlined definition of most functions.

*/

/*=========== original copyright notice ===============*/
/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#ifndef _MT_H_
#define _MT_H_

#include <stdint.h>
#include <math.h>

/* Period parameters */  
#define rngN 624
#define rngM 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

#if 0
extern uint32_t mt[];
extern int mti;
#endif

typedef struct
  {
  uint32_t mt[rngN];
  int mti;
  }
  rngs;

/* create new (if R=NULL) or reseed */
extern rngs*
init_genrand( rngs *R, uint32_t s);

extern rngs*
init_by_array( rngs *R, uint32_t init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
extern uint32_t 
genrand_int32(rngs* R);

/* generates a random number on [0,0x7fffffff]-interval */
static inline int32_t
genrand_int31(rngs* R)
{
    return (long)(genrand_int32(R)>>1);
}

/* generates a random number on [0,1]-real-interval */
static inline double 
genrand_real1(rngs* R)
{
    return genrand_int32(R)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
static inline double 
genrand_real2(rngs* R)
{
    return genrand_int32(R)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
static inline double 
genrand_real3(rngs* R)
{
    return (((double)genrand_int32(R)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
static inline double 
genrand_res53(rngs* R) 
{ 
    uint32_t a=genrand_int32(R)>>5, b=genrand_int32(R)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/*============== end of modified mt19937 code =========*/


/*
 * The following are more convenient idioms.
 *
 */


/* this is just an alias */
static inline rngs*
rng ( uint32_t seed ) { return init_genrand(NULL,seed); }

static inline void
rngseed ( rngs *R, uint32_t seed ) { init_genrand(R, seed); }

/* random uniform in (0,1). 
We use the open interval to make it consistent
with *the* Uniform distribution (a special case of beta dist) */
static inline double
randu ( rngs* R ) { return genrand_real3(R); }

/* generates an integer between 0 and n-1. 
Fixed 20050413, following Thierry Sengstag's comment on the use of "% n".
This is now done according to Matsumoto's suggestion.
*/
static inline int 
randi ( rngs* R, uint32_t n ) { return floor( genrand_real2(R) * n ); }
//randi ( rngs* R, uint32_t n ) { return floor(n*(double)rand() /((double)RAND_MAX+1.0) ); }

/* random permutation of an int array */
static inline void
randprmi ( rngs* R, int n, int* i )
{
  while ( n > 1 )
    {
    int k = randi(R, n);
    n--;
    int t = i[k]; i[k] = i[n]; i[n] = t;
    }
}

/* random permutation of an array of number */
/* note: 
    This is rarely used in practice. randprmi() is more commonly
    used to randomly permut indices.
*/
static inline void
randprm ( rngs* R, int n, double *x )
{
  while ( n > 1 )
    {
    int k = randi(R, n);
    n--;
    double t = x[k]; x[k] = x[n]; x[n] = t;
    }
}

/* random normal */

static inline void
box_mueller ( rngs* R, double *pu, double *pv )
{
  double u, v, rr, z;
  do
    {
    u = randu(R) * 2 - 1;
    v = randu(R) * 2 - 1;
    rr = u*u + v*v;
    }
  while ( rr >= 1 || rr == 0 );
  z = sqrt( -2 * log(rr)/rr );
  *pu = u*z; *pv = v*z; 
  return;
}

static inline void
randnorm ( rngs *R, int n, double *x, double mu, double sigma )
{
  int i; 
  double u, v;
  for ( i = 0; i < n/2; i++ )
    {
    box_mueller ( R, &u, &v );
    x[2*i] = sigma*u + mu;
    x[2*i+1] = sigma*v + mu;
    }
    
  if( n % 2 )
    { box_mueller ( R, &u, & v ); x[2*i] = sigma*u+mu; }
  return;
}

static inline void
randexp( rngs *R, int n, double *x, double lambda )  
{
  for(int i = 0; i < n; i++ )
    x[i] = -lambda * log(randu(R));
  return;
}

#endif /* _MT_H_ */
