
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "rng.h"
// #include "mygsea2.h"

/*==============================================================================
                               sorting routines
  ===============================================================================
  Routines from Asa Wirapati
  Output of C preprocessor on ~pwirapat/src/n/srt.[ch]
*/

/*
inline int isrt_iad_ ( int n, int *x , const double *y ) { int i, j; for( i = 1;
 i < n; i++ ) { int temp = x[i]; for( j = i; j > 0 && ((y[temp]) < (y[x[j-1]]));
 j-- ) x[j] = x[j-1]; x[j] = temp; } return n; }

int qsrt_iad ( int n, int *x , const double *y ) { int s[32]; int k = 0; int i,
 j, a, b; int t; a = 0; b = n-1; for (;;) { while ( b - a >= 32 ) { int c = (a+b
)/2; if( ((y[x[b]]) < (y[x[a]])) ) { t = x[a]; x[a] = x[b]; x[b] = t; } if( ((y[
x[c]]) < (y[x[b]])) ) { t = x[c]; x[c] = x[b]; if( ((y[t]) < (y[x[a]])) ) { x[b]
 = x[a]; x[a] = t; } else x[b] = t; } i = a-1; j = b; for (;;) { while ( ((y[x[+
+i]]) < (y[x[b]])) ) ; while ( ((y[x[b]]) < (y[x[--j]])) ) ; if( i >= j ) break;
 t = x[i]; x[i] = x[j]; x[j] = t; } t = x[i]; x[i] = x[b]; x[b] = t; if( i-a > b
-i ) { s[k++] = a; s[k++] = i-1; a = i+1; } else { s[k++] = i+1; s[k++] = b; b =
 i-1; } } if( k == 0 ) break; b = s[--k]; a = s[--k]; } isrt_iad_ (n,x , y );
return n; }
*/

/*==============================================================================
                               insertion sort
  ===============================================================================
  Recoded by TS
*/

void dirty_sort(int n, int *x , const double *y, int *xx ) {

  //  int *xx=(int *)malloc(n*sizeof(int));
  if (xx==NULL) {
    printf("NULL pointer in dirty_sort\n");
    exit(1);
  }
  //  int xx[10000];

  xx[0]=x[0];
  
  for (int i=1; i<n; i++) {
    int j=0;
    while ( y[x[i]] > y[xx[j]] ) {
      //      printf("In dirty_sort : i:%d j:%d \t %d %f \t %d %f\n",i,j,x[i],y[x[i]],xx[j],y[xx[j]]);
      j++;
      if (j==i) break;
    }
    if (j==i) {
      xx[j]=x[i];
    } else {
      //      if (i==n-1) printf("XXXX  In dirty_sort : i:%d j:%d \t %d \t %d\n",i,j,x[i],y[x[i]],xx[j],y[xx[j]]);
      for (int k=i; k>j; k--) {
	xx[k]=xx[k-1];
      }
      xx[j]=x[i];
    }
  }
  for (int i=0; i<n; i++) {
    x[i]=xx[i];
  }
  //  free(xx);
  return;
}

/* random permutation of an int array
   only m first indices are randomized,
   based on randprmi in rgn.h */
static inline void
randprmi_mfirst ( rngs* R, int n, int* i, int m )
{
  int n2,k2,nn=n;
  while ( n > nn-m )
    {
    int k = randi(R, n);
    n2=nn-n;
    k2=nn-k-1;
    n--;
    int t = i[k2]; i[k2] = i[n2]; i[n2] = t;
    }
}

/*==============================================================================
                                  weigthed_KS
 ===============================================================================

  Original version by Asa Wirapati

  m : number of genes in small list
  g : vector of indices of small-list genes in big list
      (indexing starting at 0)
  r : rank of genes in big list
  w : weigths of all genes in big list
  ng_gene: number of genes in big list

  TS changes:
  - Initialization of KS scores to zero (avoid a few f.p. operations)
  - Compute KS scores on the positive and negative side
  - Return a vector of KS scores (for absolute/positive/negative)
*/

//double *weighted_KS ( int m, int *g, double *r, int n_gene, double *w ) {
//resKS weighted_KS ( int m, int *g, double *r, int n_gene, double *w ) {
void weighted_KS ( int m, int *g, double *r, int n_gene, double *w, double *result,
		   int *iwrk) {

  /*
  printf("m: %d\n",m);
  printf("n_gene: %d\n",n_gene);
  printf("g: %d %d ...\n",g[0],g[1]);
  printf("r: %f %f ...\n",r[0],r[1]);
  printf("w: %f %f ...\n",w[0],w[1]);
  */

  //  resKS res;
  //  res.pos_score=0.0;
  //  res.neg_score=0.0;

  double pos_score=0.0;
  double neg_score=0.0;

  //  double *res=(double *)malloc(2*sizeof(double));
  result[0]=0.0;
  result[1]=0.0;

  if(m == 0 ) return;

  //  printf("Before sort\n");
  //  qsrt_iad ( m, g, r );
  /*
  for (int i=0; i<m; i++) {
    printf("before sort: %d %d %f\n",i,g[i],r[g[i]]);
  }
  */
  dirty_sort ( m, g, r, iwrk );
  /*
  for (int i=0; i<m; i++) {
    printf("after sort: %d %d %f\n",i,g[i],r[g[i]]);
  }
  */
  //  printf("After sort\n");
  
  double sum_weight = 0;
  for(int i = 0; i < m; i++ )
    sum_weight += w[g[i]];

  /*
    for (int i=0; i<m;i++) {
      printf("%d ",g[i]);
    }
    printf("\n");
  */
  //  printf("Sum weight: %f\n",sum_weight);
  
  // double d = -(r[g[0]]-1)/(n_gene-m);
  // double d = 0.0;
  double P = 0.0, Q, Qprev;
  double d,dprev;
  for(int i = 1; i <= m; i++ )
    {
      if (g[i-1]>0) {
	Qprev=(r[g[i-1]-1]-(i-1))/(n_gene-m);// -(i-1), because we have to skip (i-1) hits
      } else {
	Qprev=0.0;
      }
      dprev=P-Qprev;
      //      printf("%d  %d  r:%f  %f  %f  %f\n",i,g[i-1],-(r[g[i-1]]-(i-1)),P,Qprev,P-Qprev);

      P += w[g[i-1]]/sum_weight;
      Q = (r[g[i-1]]-i)/(n_gene-m);  // -i, because we have to skip i hits
      //      printf("%d  %d  r:%f  %f  %f  %f\n",i,g[i-1],-(r[g[i-1]]-i),P,Q,P-Q);

      d=P-Q;
      //      if(fabs(P-Q) > fabs(d) ) d = P-Q;
      //      if(fabs(P-Q) > fabs(res[0]) ) res[0] = d;
      if ( d     > pos_score ) pos_score = d;
      if ( dprev < neg_score ) neg_score = dprev;
    }
  result[0]=pos_score;
  result[1]=neg_score;
  return;
  //  return res;
  //  return -d;
}

/*==============================================================================
                                  permute_wKS
 ===============================================================================*/
void permute_wKS( int m, int *g, double *r, int n_gene, double *w,
		  int nperm, uint32_t seed, double *resks, double *resperm,
		  int *iwrk) {

  int *geneind=(int *)malloc(n_gene*sizeof(int));

  //  resKS res;
  //   resPerm resperm;

  // Compute initial KS score
  weighted_KS(m, g, r, n_gene, w, resks, iwrk);
  double KS_pos=resks[0];
  double KS_neg=resks[1];

  /*
  printf("%.15lf\t%.15lf",resks[0],resks[1]);
  for (int j=0; j<m; j++) {
    printf("\t%d",g[j]);
  }
  printf("\n");
  */

  //  double true_KS=weighted_KS(m, g, r, n_gene, w);

  double curKS_pos,curKS_neg;
  int p_count_pos=0;
  int p_count_neg=0;

  rngs *R = rng(seed);
  for(int i=0; i<n_gene; i++) {
    geneind[i]=i;
  }
  for(int i=0; i<nperm; i++) {
    // Permute genes
    // -   randprmi_mfirst is a variant of ranprmi in which only the first m positions are random
    //randprmi( R, n_gene, geneind );
    randprmi_mfirst( R, n_gene, geneind, m );
    weighted_KS(m, geneind, r, n_gene, w, resks, iwrk);
    /*
    printf("%.15lf\t%.15lf",resks[0],resks[1]);
    for (int j=0; j<m; j++) {
      printf("\t%d",geneind[j]);
    }
    printf("\n");
    */

    curKS_pos=resks[0];
    curKS_neg=resks[1];
    //    printf("Test: %f %f\n",curKS_pos,curKS_neg);
    if (curKS_pos > KS_pos) {
      /*
      printf("HIT: %f > %f ",curKS_pos,KS_pos);
      for(int j=0; j<m; j++) {
	printf("%d ",geneind[j]);
      }
      printf("\n");
      */
      p_count_pos++;
    }
    if (curKS_neg < KS_neg) { p_count_neg++; }
  }

  /* 
  printf("P+ %d %d\n",p_count_pos,nperm);
  printf("P- %d %d\n",p_count_neg,nperm);
  */
  resperm[0]=p_count_pos/(double)nperm;
  resperm[1]=p_count_neg/(double)nperm;
  resks[0]=KS_pos;
  resks[1]=KS_neg;
  //  double p_value=p_count/(double)nperm;

  free(geneind);
  return;
}
void Rpermute_wKS(int *m, int *g, double *r, int *n_gene, double *w,
		  int *nperm, int *seed, double *resks, double *resperm,
		  int *iwrk) {
  /*
  printf("Into...\n");
  printf("m     : %d\n",*m);
  printf("n_gene: %d\n",*n_gene);
  printf("nperm : %d\n",*nperm);
  printf("seed  : %d\n",*seed);
  printf("g: %d %d %d ...\n",g[0],g[1],g[2]);
  printf("r: %f %f ...\n",r[0],r[1]);
  if (w!=NULL) {
    printf("w: %f %f ...\n",w[0],w[1]);
  } else {
    printf("No weights\n");
  }
  */
  permute_wKS(*m, g, r, *n_gene, w, *nperm, *seed, resks, resperm, iwrk);
  /*
  printf("KS pos: %f\n",resks[0]);
  printf("KS neg: %f\n",resks[1]);
  printf("p.val pos: %f\n",resperm[0]);
  printf("p.val neg: %f\n",resperm[1]);
  printf("Outo...\n");
  */
  return;
}


/*==============================================================================
                                     test1
 ===============================================================================*/

void test1() {
  int n_gene=40000;
  int m=50;
  //  printf("time:%d\n",time(NULL));
  rngs *RR = rng(time(NULL));//1180636537);
  int *gg=(int *)malloc(n_gene*sizeof(double));
  //printf("\n");
  for(int j=1; j<10000;j++) {
    for (int i=0; i<n_gene; i++) {
      gg[i]=i;
    }
    randprmi_mfirst(RR, n_gene, gg, m);
    //randprmi(RR, n_gene, gg);
    for (int i=0;i<m;i++) {
      //printf("%d\t",gg[i]);
    }
    //printf("\n");
  }
  exit(0);
}

/*==============================================================================
                                     main
 ===============================================================================*/

int main() {

  // test1();

  int n_gene=12;
  int m=3;


  double *r=(double *)malloc(n_gene*sizeof(double));
  double *w=(double *)malloc(n_gene*sizeof(double));
  int *g=(int *)malloc(m*sizeof(int));
  int *iwrk=(int *)malloc(n_gene*sizeof(int));

  double *res=(double *)malloc(2*sizeof(double));
  double *rp =(double *)malloc(2*sizeof(double));

  //  resKS res;

  //  printf("Here 1\n");
  for(int i=0; i<n_gene; i++) {
    r[i]=i+1;
    w[i]=abs(i-(n_gene-1)/2)/(double)(n_gene-1);
    w[i]=1.0;
    //    printf("%d\t%f\t%f\n",i,r[i],w[i]);
  }
  //  printf("Here 2\n");
  for(int i=0;i<m;i++) {
    g[i]=600*i+4;
    //    printf("%d ",g[i]);
  }

  g[0]=0;g[1]=2;g[2]=3;
  g[0]=1;g[1]=3;g[2]=5;
  //  printf("\n");
  //  printf("Here 3\n");
  weighted_KS(m, g, r, n_gene, w, res, iwrk);
  //  printf("\n%f %f\n",res[0],res[1]);

  printf("KS pos: %f\n",res[0]);
  printf("KS neg: %f\n",res[1]);

  printf("=======================\n");

  //  exit(0);

  int nperm=10000;
  int seed=time(NULL);
  seed=1;

  printf("======== standard interface ==========\n");
  permute_wKS(m, g, r, n_gene, w, nperm, seed, res, rp, iwrk);
  printf("KS pos: %f\n",res[0]);
  printf("KS neg: %f\n",res[1]);

  printf("p.value pos: %f\n",rp[0]);
  printf("p.value neg: %f\n",rp[1]);

  printf("======== R interface ==========\n");
  Rpermute_wKS(&m, g, r, &n_gene, w, &nperm, &seed, res, rp, iwrk);

  printf("KS pos: %f\n",res[0]);
  printf("KS neg: %f\n",res[1]);

  printf("p.value pos: %f\n",rp[0]);
  printf("p.value neg: %f\n",rp[1]);

  free(iwrk);
  return 0;
}
