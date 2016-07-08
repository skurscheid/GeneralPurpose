/*
rng.c
*/

#include <stdlib.h>  /* for malloc() */
#include <stdio.h>   /* for perror() */
#include "rng.h"


/* initializes mt[rngN] with a seed */
rngs*
init_genrand( rngs *RR, uint32_t s)
{
rngs *R;
if( RR == NULL )
  {
  R = (rngs*) malloc( sizeof(*R) );
  if( !R ) { perror(__func__); return 0; }
  }
else R = RR;
R->mti = rngN+1;
uint32_t *mt = R->mt; int mti = R->mti;
  
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<rngN; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
    
R->mti = mti;
return R;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
rngs*
init_by_array( rngs* RR, uint32_t init_key[], int key_length)
{
rngs* R = init_genrand(RR,19650218UL);
uint32_t *mt = R->mt; int mti = R->mti;
    
    int i, j, k;
    i=1; j=0;
    k = (rngN>key_length ? rngN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=rngN) { mt[0] = mt[rngN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=rngN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=rngN) { mt[0] = mt[rngN-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 

R->mti = mti;
return R;
}

/* generates a random number on [0,0xffffffff]-interval */
uint32_t 
genrand_int32(rngs* R)
{
uint32_t *mt = R->mt; int mti = R->mti;

    uint32_t y;
    const uint32_t mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= rngN) { /* generate N words at one time */
        int kk;

/* pjw: the following is unnecessary in rng.[hc], since creating rngs struct
        always done via rng(). we keep it, just in case. */
        
        if (mti == rngN+1)   /* if init_genrand() has not been called, */
            init_genrand(R,5489UL); /* a default initial seed is used */


        for (kk=0;kk<rngN-rngM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+rngM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<rngN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(rngM-rngN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[rngN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[rngN-1] = mt[rngM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

R->mti = mti;

    return y;
}

