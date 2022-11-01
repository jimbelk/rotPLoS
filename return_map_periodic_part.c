/* return_map_periodic_part.c */

/* to compile using gcc type:
   gcc return_map_periodic_part.c -lgmp
   (add other options such as -o as desired --- see gcc documentation. */

/* This code was developed by Justin Moore, copyright 10/31/2022 */
/* It is subject to the creative commons legal code included in this posting. */
/* You may contact Justin via justin@math.cornell.edu */

/* It was used to perform a computation of the periodic parts of continued
fraction expanstions of rot(f_q,r) as detailed in the paper
"A piecewise linear homeomorphim of the circle which is periodic under
renormalization" by James Belk, James Hyde, and Justin Tatch Moore */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include<gmp.h> 

const unsigned char N=4;
const long D=5000;

struct fraction{
	unsigned long m;
	unsigned long n;
	};

struct point{
	mpq_t x;
	mpq_t y;
	};

/* PL homemorphisms of S1 = [0,1]/{0,1} are specified by a list of points (which themselves are structures
consisting of a pair of mpq's. The first point should have x-coordinate 0 and the last should have 
x-coordinate 1.  Note that a point with an x-coordinate of 1 is treated as a treminating symbol when
specifying the list of points.  The x-coordinates should occurr in increasing order.  */

mpq_t zero,one;


/* routine which computes gcd of two long integers via the Euclidean algorithm */
/* if one of the arguments is 0, it returns 0 */

long gcd(long x, long y)
{
  if(x < 0)
    x=-x;	
  if(y < 0)
    y=-y;
  if(x==0)
    return y;
  if(y%x == 0)
    return x;
  else return gcd(y%x,x);
}


/* Computes m = (y1-y0)/(x1-x0) */
/* All arguments should be initialized prior to calling the routine. */

void slope(mpq_t m, mpq_t x0, mpq_t y0, mpq_t x1, mpq_t y1)
{
  mpq_t dx,dy;
  mpq_init(dx);
  mpq_init(dy);
  mpq_sub(dx,x1,x0);
  mpq_sub(dy,y1,y0);
  mpq_div(m,dy,dx);
  mpq_clear(dx);
  mpq_clear(dy);
}

/* Computes y = [(y1-y0)/(x1-x0)]*(s-x0)+y0 */
/* i.e. the value y of the affine function mapping xi->yi at s */
/* all arguments should be initialized prior to calling the routine */

void L_eval(mpq_t y, mpq_t s,  mpq_t x0, mpq_t y0, mpq_t x1, mpq_t y1)
{
  mpq_t m,dx,dy,X0,X1,Y0,Y1;
  mpq_init(m);
  mpq_init(dx);
  mpq_init(dy);
  mpq_init(X0);
  mpq_init(Y0);
  mpq_init(X1);
  mpq_init(Y1);
  if(mpq_equal(x0,one)==0)
    mpq_set(X0,x0);
  if(mpq_equal(y0,one)==0)
    mpq_set(Y0,y0);
  if(mpq_equal(x1,zero)!=0)
    mpq_set(X1,one);
  else
    mpq_set(X1,x1);
  if(mpq_equal(y1,zero)!=0)
    mpq_set(Y1,one);
  else
    mpq_set(Y1,y1);
  slope(m,X0,Y0,X1,Y1);
  mpq_sub(dx,s,X0);
  mpq_mul(dy,m,dx);
  mpq_add(y,Y0,dy);
  mpq_clear(m);
  mpq_clear(dx);
  mpq_clear(dy);
  mpq_clear(X0);
  mpq_clear(Y0);
  mpq_clear(X1);
  mpq_clear(Y1);
}



/* Tests whether the PL function specified by f has any fixed points. */

/* (See comment just after the declaration of the point datatype for how PL functions are handled.) */

bool PL_fixedpoint(struct point *f, mpq_t t)
{
  int i;
  struct point P,Q;
  mpq_t m,tmp1,tmp2,tmp3;
  char *str1,*str2;
  mpq_init(P.x);
  mpq_init(P.y);
  mpq_init(Q.x);
  mpq_init(Q.y);
  mpq_init(m);
  mpq_init(tmp1);
  mpq_init(tmp2);
  mpq_init(tmp3);
  for(i=0;mpq_cmp(f[i].x,one)<0;i++)
    {
      if(mpq_equal(f[i].x,one)!=0)
	mpq_set(P.x,zero);
      else
	mpq_set(P.x,f[i].x);
      if(mpq_equal(f[i].y,one)!=0)
	mpq_set(P.y,zero);
      else
	mpq_set(P.y,f[i].y);
      if(mpq_equal(f[i+1].x,zero)!=0)
	mpq_set(Q.x,one);
      else
	mpq_set(Q.x,f[i+1].x);
      if(mpq_equal(f[i+1].y,zero)!=0)
	mpq_set(Q.y,one);
      else	
	mpq_set(Q.y,f[i+1].y);
      if((mpq_cmp(P.x,P.y)*mpq_cmp(Q.x,Q.y))<=0)
	{
	  slope(m,P.x,P.y,Q.x,Q.y);
	  mpq_mul(tmp1,m,P.x);
	  mpq_sub(tmp2,tmp1,P.y);
	  mpq_sub(tmp3,m,one);
	  if(mpq_equal(tmp3,zero)==0)
	    mpq_div(t,tmp2,tmp3);
	  else
	    mpq_set(t,P.x);
	  mpq_clear(P.x);
	  mpq_clear(P.y);
	  mpq_clear(Q.x);
	  mpq_clear(Q.y);
	  mpq_clear(m);
	  mpq_clear(tmp1);
	  mpq_clear(tmp2);
	  mpq_clear(tmp3);
	  return 1;
	}
    }
  mpq_clear(P.x);
  mpq_clear(P.y);
  mpq_clear(Q.x);
  mpq_clear(Q.y);
  mpq_clear(m);
  mpq_clear(tmp1);
  mpq_clear(tmp2);
  mpq_clear(tmp3);
  return 0;
}


void PL_eval(mpq_t z, mpq_t t, struct point *f)
/* Makes the computation z=f(t) for a PL function f. */
/* (See comment just after the declaration of the point datatype for how PL functions are handled.) */
{
  int i;
  for(i=0;mpq_cmp(f[i+1].x,t)<0;i++);
  L_eval(z,t,f[i].x,f[i].y,f[i+1].x,f[i+1].y);
}


void PLinv_eval(mpq_t s, mpq_t t, struct point *f)
/* Makes the computation s=f^(-1)(t) for a PL function f. */
/* (See how PL functions are coded in comment just after the declaration of the point datatype.) */
{
  int i;
  for(i=0;
      ((mpq_cmp(t,f[i].y)<0)&&(mpq_equal(f[i].y,one)==0))
	||
	((mpq_cmp(f[i+1].y,t)<0)&&(mpq_equal(f[i+1].y,zero)==0)) ;
      i++);
  L_eval(s,t,f[i].y,f[i].x,f[i+1].y,f[i+1].x);
}


int mpq_sort(mpq_t *L,int n)
/* Sort the first n entries in the array L of mpq's in strictly increasing order. */
/* Return the number of repeated entries. */
{
  int i,j,m=0;
  for(i=0;i<n-m-1;)
    {
      if(mpq_equal(L[i+1],L[i])!=0)
	{
	  mpq_swap(L[n-m-1],L[i+1]);
	  m++;
	}
      if(mpq_cmp(L[i+1],L[i])<0)
	{
	  mpq_swap(L[i],L[i+1]);
	  if(i > 0)
	    i--;
	  else
	    i++;
	}
      else
	i++;
    }
  return m;
}



int PL_return(struct point *f, mpq_t z)
/* This is the main routine.  It calculates f* given a PL homemorphism f of S1. */
/* Here f* is the renormalized return map for f^-1 on the interval [0,f(0)) --- renormalized by */
/* rescaling by 1/f(0). */
/* The routine returns the least m > 0 such that f^-m (0) is in [0,f(0)). */
/* Notice that this is the first digit of the continued fraction expansion of rot(f). */
/* If f has a fixed point t (which occurrs iff no such m exists), then the routine returns -1 and sets z equal t a fixed point. */

{
  mpq_t D[N],R[N],s,t,r;
  int i=0,j,n,m,l;
  if(PL_fixedpoint(f,z))
      return -1;
  for(l=1;mpq_cmp(f[l].x,one)<0;l++);
  for(i=0;i<=l;i++)
    {
      mpq_init(D[i]);
      mpq_init(R[i]);
    }
  mpq_init(r);
  mpq_init(s);
  mpq_init(t);
  mpq_set(r,f[0].y);
  mpq_set(t,one);
  for(m=0;mpq_cmp(r,t)<=0;m++)
    {
    PLinv_eval(s,t,f);
    mpq_set(t,s);
    }
  j=l;
  for(i=0;i<l;i++)  /* This loop copies the breakpoints of f^-1 into the array D and sets j 
		       equal to the location of 0 */
    {
      mpq_set(D[i],f[i].y);
      if(j==l)
	{
	  if(mpq_equal(D[i],one)!=0)
	    {
	      mpq_set(D[i],zero);
	      j=i;
	    }
	  else if(mpq_equal(D[i],zero)!=0)
	    j=i;
	}
    }
  for(i=0;i<l ;i++)
    if(i!=j) /* Skip if D[i]=0 */
      {
	while(mpq_cmp(r,D[i])<=0) /* This loop applies f to D[i] until it is in [0,r) */
	  {
	    PL_eval(t,D[i],f);
	    if(mpq_equal(t,one)!=0)
	      mpq_set(D[i],zero);
	    else
	      mpq_set(D[i],t);
	  }
      }
  n=l-mpq_sort(D,i); /* Sort D in increasing order */
  for(i=0;i<n;i++)
    {
      PLinv_eval(s,D[i],f);
      if(mpq_equal(s,zero)!=0)
	mpq_set(s,one); 
      while(mpq_cmp(r,s)<0)
	{
	  PLinv_eval(t,s,f);
	  if(mpq_equal(t,zero)!=0)
	    mpq_set(s,one);
	  else
	    mpq_set(s,t);
	}
      mpq_set(R[i],s);
    }
  for(i=0;i<n;i++)
    {
      mpq_div(f[i].x,D[i],r);
      mpq_div(f[i].y,R[i],r);
      if((mpq_equal(f[i].x,zero)!=0)&&(i > 0))
	mpq_set(f[i].x,one);
      if(mpq_equal(f[i].y,zero)!=0)
	mpq_set(f[i].y,one);
    }
  mpq_set(f[n].x,one);
  mpq_set(f[n].y,f[0].y);
  for(i=0;i<=l;i++)
    {
      mpq_clear(D[i]);
      mpq_clear(R[i]);
    }
  mpq_clear(s);
  mpq_clear(t);
  mpq_clear(r);
  return m;
}

void print_PL(struct point *f)
/* This routine outputs the details of the PL function f to the standard output. */
{
  int i;
  char *str1,*str2,*str3,*str4;
  mpq_t m,b,s,y0;
  mpq_init(m);
  mpq_init(b);
  mpq_init(s);
  mpq_init(y0);
  for(i=0;mpq_cmp(f[i].x,one)<0;i++)
    {
      if(mpq_equal(f[i].y,one)!=0)
	mpq_set(y0,zero);
      else
	mpq_set(y0,f[i].y);
      slope(m,f[i].x,y0,f[i+1].x,f[i+1].y);
      mpq_mul(s,m,f[i].x);
      mpq_sub(b,y0,s);
      str1=mpq_get_str(NULL,10,f[i].x);
      str2=mpq_get_str(NULL,10,f[i+1].x);
      str3=mpq_get_str(NULL,10,y0);
      str4=mpq_get_str(NULL,10,f[i+1].y);
      printf("[%s,%s] -> [%s,%s]  ",str1,str2,str3,str4);
      free(str1);
      free(str2);
      free(str3);
      free(str4);
      if(mpq_equal(m,one)==0)
	{
	  str1=mpq_get_str(NULL,10,m);
	  printf("%s",str1);
	  free(str1);
	}
      putchar('t');
      if(mpq_cmp(zero,b)<0)
	putchar('+');
      if(mpq_equal(b,zero)==0)
	{
	  str1=mpq_get_str(NULL,10,b);   					
	  printf("%s\n",str1);
	  free(str1);  
	}
      else
	putchar('\n');
    }
  putchar('\n'); 
  mpq_clear(m);
  mpq_clear(b);
  mpq_clear(s);
  mpq_clear(y0);
}


void init_f_qr(mpq_t q, mpq_t r, struct point *f)
/* This function makes the assignment f = f_q,r where f_q,r is 
the composition of the interval exchange [0,1/(q+1)) <-> [1/(q+1),1) and
rotation by r. */

/* NOTE: this routine does not allocate the pointer *f and does not initialize
   the arbitrary precision coordinates which are listed by f. */
{
  mpq_t u,t,s,ratio,r1;
  char cmp;
  mpq_init(t);
  mpq_init(s);
  mpq_init(u);
  mpq_init(ratio);
  mpq_add(t,q,one); /* t=q+1 */
  mpq_div(s,one,t); /* s=1/(q+1) */
  mpq_div(u,q,t);  /* u = q/(q+1) */
  cmp=mpq_cmp(r,u);
  if(cmp < 0)
    {
      mpq_set(f[0].x,zero);
      mpq_add(f[0].y,s,r);
      mpq_div(ratio,r,q);
      mpq_sub(f[1].x,s,ratio);
      mpq_set(f[1].y,one);
      mpq_set(f[2].x,s);
      mpq_set(f[2].y,r);
      mpq_set(f[3].x,one);
      mpq_set(f[3].y,f[0].y);
    }
  else if(cmp > 0)
    {
      mpq_init(r1);
      mpq_sub(r1,r,one);
      mpq_set(f[0].x,zero);
      mpq_add(f[0].y,s,r1);     
      mpq_set(f[1].x,s);
      mpq_set(f[1].y,r);
      mpq_mul(t,q,r1); /* t= q*(r-1) */
      mpq_sub(f[2].x,s,t);
      mpq_set(f[2].y,one);
      mpq_set(f[3].x,one);
      mpq_set(f[3].y,f[0].y);
      mpq_clear(r1);
    }
  else if(cmp ==0)
    {
      mpq_set(f[0].x,zero);
      mpq_set(f[0].y,zero);
      mpq_set(f[1].x,s);
      mpq_set(f[1].y,u);
      mpq_set(f[2].x,one);
      mpq_set(f[2].y,one);     
    }
  mpq_clear(t);
  mpq_clear(s);
  mpq_clear(u);
  mpq_clear(ratio);
}


int PL_duplicate(struct point *f)
/* This routine creates a copy of the PL function pointed to by f. */
/* Here f is viewed as a sequence of points with final point the first point having
   x-coordinate 1.  The return value l is the number of points in this PL function. 
   The duplicate of f is written at the address f+l --- i.e. starting just after the terminating symbol.*/

/* CAUTION: this routine does not check whether the block of memory pointed to by f is sufficient
   for the duplication --- care should be taken only to call this routine when there is enough remaining
   allocated space in f. */
{
  int i,l;
  for(l=1;mpq_equal(f[l].x,one)==0;l++);
  l++;
  for(i=0;i<=l;i++)
    {
      mpq_set(f[i+l].x,f[i].x);
      mpq_set(f[i+l].y,f[i].y);
    }
  return l;
}



int PL_equal(struct point *f, struct point *g)
/* This routine tests whether two PL functions f and g are equal by checking whehter
   they are the same list of points. */

/* CAUTION: This routine does not check to see if points listed in either f or g are redundant - for instance
   if two consecutive segments have the same slope.  As a result, it may yield false negative results. */
/* In this program, it is used to certify periodic behavior in the map f->f*, in which case false negatives
   are not problematic in practice. */
{
  int i;
  for(i=0;mpq_equal(f[i].x,one)==0;i++)
    if((mpq_equal(f[i].x,g[i].x)*mpq_equal(f[i].y,g[i].y))==0)
      return 0;
  if((mpq_equal(f[i].x,g[i].x)*mpq_equal(f[i].y,g[i].y))==0)
    return 0;
  return 1;
}

/* The next routines change the text color and/or face written to the standard output. */

void black()
{
  printf("\033[0;30m");
}

void red()
{
  printf("\033[0;31m");
}

void green()
{
  printf("\033[0;32m");
}

void yellow()
{
  printf("\033[0;33m");
}

void blue()
{
  printf("\033[0;34m");
}

void purple()
{
  printf("\033[0;35m");
}

void cyan()
{
  printf("\033[0;36m");
}

void white()
{
  printf("\033[0;37m");
}

void black_bf()
{
  printf("\033[1;30m");
}

void red_bf()
{
  printf("\033[1;31m");
}

void green_bf()
{
  printf("\033[1;32m");
}

void yellow_bf()
{
  printf("\033[1;33m");
}

void blue_bf()
{
  printf("\033[1;34m");
}

void purple_bf()
{
  printf("\033[1;35m");
}

void cyan_bf()
{
  printf("\033[1;36m");
}

void white_bf()
{
  printf("\033[1;37m");
}

void bold()
{
  printf("\033[1m");
}


void reset()
/* Call this routine to resent the text color/face back to normal. */
{
  printf("\033[0m");
}


void mpz_excise_factors(mpz_t z, long n)
{
  mpz_t g,x,y;
  bool done=false;
  mpz_init(g);
  mpz_init(x);
  mpz_init(y);
  mpz_set_si(x,n);
  while(!done)
    {
      mpz_gcd(g,x,z);
      mpz_remove(y,z,g);
      if(mpz_cmp(y,z)!=0)
	mpz_set(z,y);
      else
	done=true;
    }
  mpz_clear(g);
  mpz_clear(x);
  mpz_clear(y);
}


long mpq_set_cf(mpq_t q, long *n)
/* This routine assigns the continued fraction [0;n0,n1,...] to q. */
{
  long i,l;
  mpq_t s,t;
  mpq_set(q,zero);
  mpq_init(s);
  mpq_init(t);
  for(l=0;n[l]!=-1;l++);
  for(i=l-1;0<=i;i--)
    {
      mpq_set_si(t,n[i],1L);
      mpq_add(s,t,q);
      mpq_div(q,one,s);
    }
  mpq_canonicalize(q);
  mpq_clear(s);
  mpq_clear(t);
  return l;
}


bool cyclic_equal(long *A, long *B,int m,int n)
/* This routine determines when the array A of length m is equal to the array B of length n up to 
   cyclic permutation. */
  
{
  int i,j;
  if(m!=n)
    return false;
  if(n==0)
    return true;
  for(i=0;i<n;i++)
    {
      for(j=0;(j<n)&&(A[j]==B[(i+j)%n]);j++);
      if(j==n)
	return true;
    }
  return false;
}

/* This routine finds the lexicographically minimal cyclic permutation of A, determines
the shortest period, and writes it over the initial values of A.
The length of the period is the return value. */

int cyclic_canonize(long *A, int n)
{
  int i=0,j,k;
  long *B;
  if(n==0)
    return 0;
  if(n==1)
    return 1;
  for(j=i+1;j<n;j++)
    {
      for(k=0;(k<n)&&(A[(i+k)%n]==A[(j+k)%n]);k++);
      if(A[(i+k)%n] > A[(j+k)%n])
	i=j;
      else if(k==n)
	{
	for(k=0;k < j-i; k++)
	  A[k]=A[i+k];
	return j-i;
	}
    }
   B=malloc(sizeof(long)*n);
   for(k=0;k<n;k++)
     B[k]=A[(i+k)%n];
   for(k=0;k<n;k++)
     A[k]=B[k];
   free(B);
   return n;       
}

/* Here is some sample code making use of the above routines. */
/* This implementation lists the periodic parts in continued faction expansions of rot(f_q,r) when it is irrational. */


int main()
{
  struct fraction Q; /* This variable contains the current value of q for f_q,r. */
  char report_rational=0,report_irrational=0,colored_output=1; /* These flags determine whether examples with rational or
								  irrational rotation numbers are reported and also if the output
								  should be color-coded. */
  long i,j,k,l,m,n;
  long K=1000; /* This parameter is the largest value of a denominaor for r which is considered. */
  long theta[D]; /* contains a continued fraction expansion */
  long loc[D]; /* Keeps track of the locations within the array f of PL functions.  loc[0] is the number of entries. */
  long R[D]; /* contains the periodic parts of eventually periodic continued fraction expansions of rotation numbers. */
  long rloc[D]; /* Contains the starting points of entries of R. */
  int rc=0; /* Number of entries in R. */
  long worst_rat,worst_irr;
  char *str1,*str2,*r_str,repeat,rational,report;
  mpq_t q,*r,s,t,tmp;
  struct point *f; /* this array will contain several PL functions.  The displacements to individual functions is specified by loc */
  bool flag;
  clock_t start_time,last_time,current_time;
  start_time=clock();
  last_time=start_time;
  if(colored_output)
    reset();
  mpq_init(zero);
  mpq_init(one);
  mpq_init(q);
  mpq_init(s);
  mpq_init(t);
  mpq_init(tmp);
  mpq_set_ui(one,1,1);
  mpq_set_ui(zero,0,1);
  f=(struct point *)malloc(D*N*sizeof(struct point));
  for(i=0;i<D*N;i++)
    {
      mpq_init(f[i].x); 	
      mpq_init(f[i].y);
    }
  /* This next block of code creates a list of values of r to consider. */
  n=0;
  for(k=2;k<=K;k++)
    for(i=1;i<k;i++)
      if(gcd(i,k)==1L)
	n++;
  r=(mpq_t *)malloc(n*sizeof(mpq_t));
  n=0;
  for(k=2;k<=K;k++)
    for(i=1;i < k ;i++)
      if(gcd(i,k)==1L)
	{
	  mpq_init(r[n]);
	  mpq_set_ui(r[n++],i,k);
	}
  /* If K <= 100, it orders
     the list (for large values of K this step is omitted to save time). */
  if(K<=100)
    mpq_sort(r,n);
  
  for(Q.n=3L;Q.n<=9L;Q.n++)
    for(Q.m=2L;Q.m<Q.n;Q.m++)
      if(gcd(Q.m,Q.n)==1)
	{
	  rc=0;
	  worst_rat=0;
	  worst_irr=0;
	  rloc[0]=0;
	  printf("q=%li/%li:\n",Q.m,Q.n);
	  /* This next block of code determines how much space to allocate for the PL functions
	     and then allocates and initializes the space. */
	  mpq_set_ui(q,Q.m,Q.n);
	  
	  /* The array f lists several PL functions.  The first function will be the PL function on which we are
	     running the algorithm -- in this case f_q,r.  Subsequent functions in the list will be obtained by applying
	     the renomalization proceedure f -> f*.
	     The array loc keeps track of the locations of each PL function in this list.
	     Keeping track of the list of all iterates of the renomalization proceedure
	     allows us to test for repetitions.*/
	  for(i=0;i<n;i++)
	    {
	      init_f_qr(q,r[i],f);
	      loc[0]=0; 
	      theta[0]=0;
	      repeat=0; /* repeat and rational are flags indicating either there has been a repeat in the
			   renormalization iteration f -> f* or else that a fixed point has been encountered,
			   indicating a rational rotation number. */
	      rational=0;
	      for(m=1;(m<D-1)&&(repeat==0)&&(rational==0);m++)
		{
		  loc[m]=loc[m-1]+PL_duplicate(f+loc[m-1]);
		  theta[m-1]=PL_return(f+loc[m],s);
		  if(theta[m-1]==-1)
		    {
		      rational=1;
		      if((m > 1)&&(mpq_equal(s,zero)==0)&&(mpq_equal(s,one)==0))
			{
			  mpq_mul(t,s,f[loc[m-2]].y);
			  PLinv_eval(s,t,f+loc[m-2]);			
			  for(theta[m-2]=1;mpq_equal(s,t)==0;theta[m-2]++)
			    {
			      PLinv_eval(tmp,s,f+loc[m-2]);
			      mpq_set(s,tmp);
			    }
			}
		    }
		  else
		    for(k=0;(k<m)&&(repeat==0); k++)
		      {
			repeat=PL_equal(f+loc[m],f+loc[k]);
		      }
		}
	      k--; /* after deciment, k is the location of the start of the first period if one occurrs. */
	      m--; /* after decriment, m denotes the number of digits of the continued fraction expansion before the end of the first period
		      if a repetition occurres, the end of the expansion if the rotation number is rational,
		      of D-2 if the cap was reached. */
	      
	      /* The next blocks of code generate output depending on how reporting flags are set. */
	      if(rational==0)
		{
		  if(repeat==1)
		    {
		      if(worst_irr < m)
			worst_irr=m;
		      if(report_irrational==1)
			{
			  r_str=mpq_get_str(NULL,10,r[i]);
			  printf("q=%li/%li r=%s : ",Q.m,Q.n,r_str); 
			  free(r_str);
			}
		      if(report_irrational==1)
			{
			  for(j=0;j<k;j++)
			    printf(" %li",theta[j]);
			  printf(" (");
			}
		      for(j=k;j<m;j++)
			{
			  if(report_irrational==1)
			    {
			      printf(" %li",theta[j]);
			    }
			  R[rloc[rc]+j-k]=theta[j];
			}
		      rloc[rc+1]=rloc[rc]+cyclic_canonize(R+rloc[rc],m-k);
		      for(j=0;(j<rc)&&(!cyclic_equal(R+rloc[j],R+rloc[rc],rloc[j+1]-rloc[j],rloc[rc+1]-rloc[rc]));j++);
		      if(j==rc)
			rc++;
		      if(report_irrational==1)
			printf(" )\n");
		    }
		  else if(repeat!=1)
		    {
		      /* If no repetition or fixed point has been found before a threshold, then we find outselves here.
			 Being in this case suggests a candidate for a function with rotation number
			 which is neither rational nor quadratic irrational. The code halts. */
		      r_str=mpq_get_str(NULL,10,r[i]);
		      if(colored_output)
			red();
		      printf("q=%li/%li r=%s : ",Q.m,Q.n,r_str); 
		      free(r_str);
		      for(j=1;j<=m;j++)
			printf(" %li",theta[j]);
		      printf("...\n");
		      reset();
		      exit(0);
		    }
		}
	      else if(rational==1)
		{
		  if(worst_rat < m)
		    worst_rat=m;
		  if(report_rational==1)
		    {
		      if(colored_output)
			cyan();
		      r_str=mpq_get_str(NULL,10,r[i]);
		      printf("q=%li/%li r=%s : ",Q.m,Q.n,r_str); 
		      free(r_str);
		      j=mpq_set_cf(tmp,theta);
		      str1=mpq_get_str(NULL,10,tmp);
		      printf("  %s ; continued fraction had %li digits.\n",str1,j);
		      free(str1);
		      if(colored_output)
			reset();
		    }
		}
	    }
	  printf("Periodic parts:\n");
	  for(i=0;i<rc;i++)
	    {
	      printf("( ");
	      for(j=rloc[i];j<rloc[i+1];j++)
		printf("%li ",R[j]);
	      printf(")\n");
	    }
	  printf("Longest continued fraction of rt(f_q,r) which was rational had %li digits.\n",worst_rat);
	  printf("Longest number of digits before peridicicty was %li.\n",worst_irr);	  
	  current_time=clock();
	  printf("Time elapsed: %1.5f seconds.\n\n",((double)current_time - (double)last_time)/CLOCKS_PER_SEC);
	  last_time=current_time;
	}
  /* clear arbitrary precision variables and free up allocated memory. */
  mpq_clear(zero);
  mpq_clear(one);
  mpq_clear(s);
  mpq_clear(q);
  mpq_clear(t);
  mpq_clear(tmp);
  for(i=0;i<100*N;i++)
    {
      mpq_clear(f[i].x); 	
      mpq_clear(f[i].y);
    }
  for(i=0;i<n;i++)
    {
      mpq_clear(r[i]);
    }
  free(r);
  free(f);
  if(colored_output)
    reset();
  current_time=clock();
  printf("Total time elapsed: %1.5f seconds.\n\n",((double)current_time - (double)start_time)/CLOCKS_PER_SEC);
  printf("Done.  No counterexamples found.\n");
}


