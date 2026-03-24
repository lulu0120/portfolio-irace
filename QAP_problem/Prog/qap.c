/*****************************************************************************/
/*                                                                           */
/*      Version:  1.00   Date: 19/04/96   File: main-trial.c                 */
/* Last Version:                          File:                              */
/* Changes:                                                                  */
/* 19/04/96 Created                                                          */
/*                                                                           */
/* Purpose:                                                                  */
/*                                                                           */
/*                                                                           */
/* Author:  Thomas Stuetzle                                                  */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/*===========================================================================*/

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "timer.h"
#include "stat.h"

#define LINE_BUF_LEN  512
#define PROG_ID_STR "Basic input routine and local search for QAP, V0.1\n"
#define CALL_SYNTAX_STR "lsqap <instance>\n" 
#define VERBOSE(x) 
#define VERYVERBOSE(x)
 
#define XOR(x,y) ((x && !y) || (!x && y))
#define MIN(x,y)        ((x)<=(y)?(x):(y))
#define MAX(x,y)        ((x)>=(y)?(x):(y))

#define DEBUG( x ) x

#define FALSE 0
#define TRUE  1

#define E     2.7182818

/* --- variables for main program ---------------------------------------- */

/* constants for a pseudo-random number generator, details see Numerical Recipes in C book */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* seed for the random number generator */

long int seed = -1; /* MAO: initialized to -1, meaning "uninitialized"... */

long int opt;

long int **d;  /* first matrix read from input file, typically distance matrix */
long int **f;  /* second matrix read from input file, typically flow matrix */

long int  *p;  /* array containing the current solution. if d is the distance matrix 
		  and f the flow matrix, p[i] is the item assigned to location i
	       */

long int  n;                    /* instance size of the QAP instance */
long int  best_known;           /* best objective function value of QAP instance */

long int  best_found;           /* best solution found in local search */

long int  best_since_restart;   /* best solution since temperature was increased again */

long int  re_heat_start_it;     /* indicates in which iteration restart or reheating took place */

long int null_diagonal_flag  = FALSE;  /* at least one matrix has zero diagonal: TRUE */
long int d_symmetric_flag    = FALSE;  /* if first (d) matrix is symmetric: TRUE */
long int f_symmetric_flag    = FALSE;  /* if second (f) matrix is symmetric: TRUE */
long int make_symmetric_flag = FALSE;  /* convert asymmetric instance into symmetric 
					  instance: TRUE
				       */

long int first_it_flag;         /* first iteration of local search: yes / no */

FILE    *output_file;           /* output file for developement of solution quality */
FILE    *output_summary_file;   /* output file for developement of solution quality */
char name_buf[LINE_BUF_LEN];    /* input file name */

long int   **move_values;       /* move values in best improvement local search */
long int   **tabu_values;       /* entries give iteration up to which an attribute is forbidden */ 

long int  *s;             /* current solution */
long int  *shadow;        /* shadow solution in ILS */
long int  *best;          /* best solution found so far */

long int  *dlb;           /* don't look bits, just one vector needed */
long int  dlb_flag;       /* set to one if don't look bits should be used */

long int  local_search_flag  = 0;  /* indicates which local search should be run
				      0: first improvement, don't look bits
				      1: first improvement, no don't look bits
				      2: best improvement 
				      3: tabu_search short = 1 * n 
				      4: tabu_search medium = 4 * n 
				      5: tabu_search medium = 10 * n
				   */

long int  s_value;       /* objective function value of current solution */
long int  shadow_value;  /* objective function value of shadow solution */
long int  best_value;    /* objective function value of best solution */

long int  rchosen, schosen;    /* locations involved in move */
long int  tabu_search_length;  /* how many iterations the tabu search is run, is a multiple of instance size */
long int  tabu_list_length;    
long int  iterations_before_aspiration;  /* robust tabu search specific, multiple of instance size */
long int  aspirating_iteration;
long int  max_trials;          /* max. number of repetitions of robust tabu search */
long int  trials;              /* Number trials */


double    T_current;                 /* Temperature for large step marcov chains */
double    T_init;                    /* Initial Temperature for LSMC */
double    T_end;                     /* End Temperature for LSMC */
long int  T_steps;                   /* No iterations at each temp for LSMC */
double    T_constant;                /* Temperature for constant temp schemes */
double    alpha;                     /* factor for geometric cooling */
double    acceptance_ratio;          /* ratio of accepted worsening moves */
long int  no_accepted;

long int  k_max;                     /* parameter for VNS */
long int  k_min;                     /* parameter for VNS */
long int  k_var;                     /* parameter for VNS */
long int  k_exact;                   /* for fixed perturbation size */
long int  perturbation_flag=INT_MAX; /* indicates what type of perturbation 
					0 - 4: fixed perturbation size
					5 - 9: VNS variation of perturbation strength
					0: n / 16, fixed
					1: n / 8, fixed
					2: n / 4, fixed
					3: n / 2, fixed
					4: 3 n / 4, fixed
					5: n / 16, variable
					6: n / 8, variable
					7: n / 4, variable
					8: n / 2, variable
					9: 3 n / 4, variable
				     */

long int  acceptance_flag=INT_MAX;   /* indicates what type of acceptance criterion 
					0 - 4: ConstTemp
					5 - 9: Variations (LSMC, Restart)
					0: ConstTemp, T = 0 (Better)
					1: ConstTemp, T = 0.5% (low-med)
					2: ConstTemp, T = 1.5% (med-high)
					3: ConstTemp, T = INFTY (RandomWalk)
					4: Restart (T = 0), quick
					5: Restart (T = 0), med
					6: Restart (T = 0), slow
					7: LSMC (as is now)
					8: n / 2, variable
					9: 3 n / 4, variable
				     */


long int  optimal;                   /* Stop when hitting a solution with that or better value */

long int  iterations;                /* Number of ILS iterations */
long int  max_iterations;            /* Max. Number of ILS iterations */
long int  iteration_best_found;      /* iteration at which best solution is found */

long int  last_improvement;          /* last time an improved solution is found since restart */

long int  restart_iterations;        /* restart fom new initial solution after .. iterations*/

double    time_limit;               /* time limit */
double    time_best_found = 0.0;    /* stores the time at which best solution is found 
					in each trial */


long int  *best_in_try;       /* store some summary data */
long int  *best_found_at;     /* store some summary data */
double    *t_best_found;   /* store some summary data */



double ran01( long *idum ) {
/* 
      FUNCTION:      returns a pseudo-random number
      INPUT:         a pointer to the seed variable 
      OUTPUT:        a pseudo-random number uniformly distributed in [0,1]
      (SIDE)EFFECTS: changes the value of seed
*/
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

void read_problem_size( FILE *input ) {
/* 
      FUNCTION:      read the dimension of the QAP instanced/
      INPUT:         a pointer to the input file 
      OUTPUT:        none
      (SIDE)EFFECTS: assigns n the instance dimension
*/
  if( fscanf(input, "%ld", &n) <= 0) {
    fprintf(stderr, "error reading qap size value in data file\n");
    exit(1);	
  }
  VERYVERBOSE ( printf("QAP instance size %ld\n",n); )
}

void read_best_known_value( FILE *input ) {
/* 
      FUNCTION:      read the best known objective function value of the QAP instance
      INPUT:         a pointer to the input file 
      OUTPUT:        none
      (SIDE)EFFECTS: assigns best_known the read value
*/
  if( fscanf(input, "%ld ", &best_known) < 0 ) {
    fprintf(stderr, "error reading best known solution in data file\n");
    exit(1);	
  }
  VERYVERBOSE ( printf(" Best known solution %ld\n",best_known); )
}

long int ** read_matrix( FILE *input, long int size )
/*    
      FUNCTION:      read in a QAP instance matrix
      INPUT:         Pointer to input file, size of QAP instance
      OUTPUT:        pointer to matrix, has to be freed before program stops
      (SIDE)EFFECTS: allocates a memory of appropriate size for the matrix
*/
{
  long int     i, j;
  long int     **matrix;

  if((matrix = malloc(sizeof(long int) * size * size +
		      sizeof(long int *) * size	 )) == NULL){
/*      printf("Out of memory, exit."); */
    exit(1);
  }
  for ( i = 0 ; i < size ; i++ ) {
    matrix[i] = (long int *)(matrix + size) + i*size;
    for ( j = 0  ; j < size ; j++ ) {
      if( fscanf(input, "%ld", &matrix[i][j]) < 0) {
	fprintf(stderr, "error reading (%ld, %ld) qap_distance in data file\n", i, j);
	exit(1);	
      }
    }
  }
  return matrix;
}

long int check_null_diagonal ( long int **matrix, long int size ) {
/* 
      FUNCTION:      check whether the Matrix matrix has a zero diagonal
      INPUT:         pointer to the matrix
      OUTPUT:        TRUE if null diagonal, otherwise FALSE
      (SIDE)EFFECTS: none
*/
  long int   i;
  
  for ( i = 0 ; i < size ; i++ ) {
    if( matrix[i][i] != 0 ) {
      return FALSE;
    }
  }
  return TRUE;
}

long int check_symmetry ( long int **matrix, long int size ) {
/* 
      FUNCTION:      check whether the Matrix matrix is symmetric
      INPUT:         pointer to the matrix
      OUTPUT:        TRUE if symmetric, otherwise FALSE
      (SIDE)EFFECTS: none
*/
  long int   i, j;
  
  for ( i = 0 ; i < size - 1 ; i++ ) {
    for ( j = i + 1 ; j < size ; j++ ) {
      if( matrix[i][j] != matrix[j][i] )
	return FALSE;
    }
  }
  return TRUE;
}
       
void make_matrix_symmetric( long int **matrix, long int size ) {
/* 
      FUNCTION:      makes an asymmetric matrix symmetric (calculates M = M + M-transpose)
      INPUT:         pointer to the matrix
      OUTPUT:        none
      (SIDE)EFFECTS: makes the Matrix matrix symmetric 
*/
  long int  i, j;  /* index variables */
  long int  help;
  
  for ( i = 0 ; i < size ; i++ ) {
    for ( j = 0 ; j < i ; j++ ) {
      help = matrix[i][j] + matrix[j][i];
      matrix[i][j] = help;
      matrix[j][i] = help;
    }
  } 
}

long int ** init_move_values( )
/*    
      FUNCTION:      allocate and initialize the tabu values
      INPUT:         instance size
      OUTPUT:        pointer to matrix, has to be freed before next trial
      (SIDE)EFFECTS: current assignment is modified
*/
{
  long int i, j;
  long int     **matrix;

  VERYVERBOSE ( printf("initialze matrix of move values\n"); ) 

  if((matrix = malloc(sizeof(long int) * n * n +
		      sizeof(long int *) * n	 )) == NULL){
/*      printf("Out of memory, exit."); */
    exit(1);
  }

  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (long int *)(matrix + n) + i*n;
    for (j = 0 ; j < n ; j++ ) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

void print_solution( long int *p ) {
/*    
      FUNCTION:      prints solution p
      INPUT:         pointer to the solution
      OUTPUT:        none
      (SIDE)EFFECTS: none
*/
  int i;

  printf("Assignment: \n");
  for ( i = 0 ; i < n ; i++ ) {
    printf(" %ld ",p[i]);
  }
  printf("\n");
}

void print_matrix ( long int **m ) {
/*    
      FUNCTION:      prints a matrix m
      INPUT:         pointer to the matrix 
      OUTPUT:        none
      (SIDE)EFFECTS: none
*/
  long int i, j;

  printf("\n");
  for ( i = 0 ; i < n ; i++ ) {
    for ( j = 0 ; j < n ; j++ ) {
      printf(" %ld ", m[i][j]);
    }
    printf("\n");
  }
}    

void print_header ( ) {
/*    
      FUNCTION:      prints the header of the output file
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: none
*/

  fprintf(output_file,"Machine Pentium III, 700 MHz, Cache unknown\n");
  fprintf(output_file,"Compiler gcc  egcs-2.91.66 19990314 Flags -O -ansi -Wall\n");
  fprintf(output_file,"OS SuSe Linux 6.1\n\n");
  fprintf(output_file,"Algorithm ILS (LSMC variant, DLB local search) QAPs, V1.0\n\n");
  fprintf(output_file,"Parameter Settings\n\n");
  fprintf(output_file,"max-tries %ld\n", max_trials);
  fprintf(output_file,"tabu-length %ld\n", tabu_search_length);
  fprintf(output_file,"time %f\n\n", time_limit);

  fprintf(output_file,"Initialization time %f\n\n", elapsed_time ( VIRTUAL ) );

  fprintf(output_file,"begin problem %s\n",name_buf);
}

long int compute_evaluation_function( long int *p ) {
/*    
      FUNCTION:      computes the objective function value for the QAP
      INPUT:         pointer to a solution 
      OUTPUT:        none
      (SIDE)EFFECTS: none
      COMMENTS:      Division by 2 has to be done if we have a n asymmetric instance which 
                     has been converted into a symmetric one (indicated by make_symmetric_flag). 
		     This is due to the particular way of doing this conversion.
*/
   long int  i, j;
   unsigned long  obj_f_value; /* unsigned, because with only "long int" we have an overflow
				  on some few instances, for example, tai100b. This is because 
				  of making this instance symmetric (see make_matrix_symmetric) 
			       */
   obj_f_value = 0;
   for ( i = 0 ; i < n ; i++ ) {
     for ( j = 0 ; j < n ; j++ ) {
       obj_f_value += d[i][j] * f[p[i]][p[j]];
     }
   }
   if ( make_symmetric_flag ) 
     obj_f_value /= 2;
   VERYVERBOSE ( printf("objective function value = %lu \n\n",obj_f_value); )
   return obj_f_value;
}

long int * generate_random_vector( long int size )
/*    
      FUNCTION:      generates a random vector, quick and dirty
      INPUT:         vector dimension
      OUTPUT:        returns pointer to vector, free memory after using the vector
      (SIDE)EFFECTS: none
*/
{
   long int  i, j, help;
   long int  *v;

   v = malloc( size * sizeof(long int) );

   for ( i = 0 ; i < size; i++ ) 
     v[i] = i;

   for ( i = 0 ; i < size-1 ; i++) {
     j = (long int) ( ran01( &seed ) * (size - i)); 
     assert( i + j < size );
     help = v[i];
     v[i] = v[i+j];
     v[i+j] = help;
   }
   VERYVERBOSE ( printf("Random vector:\n");
   for (i = 0 ; i < size ; i++ ) 
     printf(" %ld ",v[i]);
   printf("\n"); )
   return v;
}


void copy_from_to ( long int *v, long int *w )
{
  int i;

  for ( i = 0 ; i < n ; i++ ) {
    w[i] = v[i];
  }
}


long int termination_condition( void )
{

/*    printf("it %ld, maxit %ld, time %f tl %f, best %ld, opt %ld\n",iterations,max_iterations,elapsed_time( VIRTUAL ), time_limit, best_value, optimal); */
  return ( ((iterations >= max_iterations) && (elapsed_time( VIRTUAL ) > time_limit)) || 
	  (best_value <= optimal)); 

}
 

void swap( long int pos_1, long int pos_2, long int *q ) {   
/*    
      FUNCTION:      swap items at positions pos_1 and pos_2
      INPUT:         positions 1 and 2, pointer to current assignment
      OUTPUT:        none
      (SIDE)EFFECTS: current assignment is modified
*/
  long int  help;
  
  help     = q[pos_1];
  q[pos_1] = q[pos_2];
  q[pos_2] = help;
}

void best_2_opt_asymmetric ( long int * q ) {
/*    
      FUNCTION:      best improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      fast evaluation of moves with additional matrix move_values
                     the first local search iteration is slow 
                     local search looks for best neighboring solution in each iteration
*/

  long int   improvement = TRUE;
  long int   u, v, k;
  long int   tmp;
  long int   original_symmetric_factor; /* = 2: original symmetric instance
					   = 1: original asymmetric instance 
					*/
  long int   **move_values;             /* matrix of move values in previous iteration
					   allows for fast evaluation of neighbourhood
					*/
  long int   first_it_flag = TRUE;      /* first iteration of local search: TRUE */
  long int   max_decrease;              /* largest decrease found so far in neighbourhood scan */
  long int   rchosen = n, schosen = n;  /* memorize which is best move in current iteration */
  long int   r, s;                      /* memorize which is best move in previous iteration */

  VERBOSE ( printf("best imp, asymmetric case\n"); )
  if ( make_symmetric_flag )
    original_symmetric_factor = 1;
  else
    original_symmetric_factor = 2;

  best_found  = compute_evaluation_function( q );

  if((move_values = malloc(sizeof(long int) * n * n +
		      sizeof(long int *) * n )) == NULL){
/*      printf("Out of memory, exit."); */
    exit(1);
  }
  for ( k = 0 ; k < n ; k++ ) {
    move_values[k] = (long int *)(move_values + n) + k*n;
  }

/*    printf("best 2-opt\n"); */

  r = rchosen;
  s = schosen;

  while ( improvement ) {
    improvement = FALSE;
    max_decrease = LONG_MAX;
    /* in the first local search iteration the full neighborhood has to be evaluated */
    if (first_it_flag) {
      first_it_flag = FALSE;
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  tmp = 0;
	  for ( k = 0 ; k < n ; k++ ) {
	    if ( (k != u) && (k != v) ) {
	      tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
		d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
		d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
		d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	    }    
	  }
	  tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	    d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	    d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
	    d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	  move_values[u][v] = tmp;
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
	    rchosen = u;
	    schosen = v;
	  }
	}
      }
    } else {
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  if (u == r || v == s || u == s || v == r) {
	    tmp = 0;
	    for ( k = 0 ; k < n ; k++ ) {
	      if ( (k != u) && (k != v) ) {
		tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
		  d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
		  d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
		  d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	      }    
	    }
	    tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	      d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	      d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
	      d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	    move_values[u][v] = tmp;
	    if (tmp < max_decrease) {
	      max_decrease = tmp;
	      rchosen = u;
	      schosen = v;
	      /*   	    printf(" max-decr = %ld\n",tmp); */
	    }
	  } else { /* change derived from move_values */
	    tmp = ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	      ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] )
	      + ( d[u][r] - d[v][r] + d[v][s] - d[u][s] ) *
	      ( f[q[u]][q[s]] - f[q[v]][q[s]] + f[q[v]][q[r]] - f[q[u]][q[r]] );
 	    tmp += move_values[u][v];
	    move_values[u][v] = tmp;
	  }
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
	    rchosen = u;
	    schosen = v;
	  }	   
	}
      }
    }
    if ( max_decrease < 0 ) {      /* Obj. function value can be improved */
      assert (rchosen < schosen);
      improvement = TRUE;
      best_found += max_decrease;
      swap(rchosen,schosen,q);    
      r = rchosen; /* memorize previously done move */
      s = schosen;/* memorize previously done move */
      VERYVERBOSE ( printf("improvement %ld, best_found %ld, exchange %ld and %ld\n",max_decrease, best_found,rchosen, schosen); )
    }
  }
  free ( move_values );
}

void best_2_opt_symmetric ( long int * q ) {
/*    
      FUNCTION:      best improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      faster evaluation of moves with additional matrix move_values
*/

  long int   improvement = TRUE;
  long int   u, v, k;
  long int   tmp;
  long int   original_symmetric_factor;
  long int   **move_values;
  long int   first_it_flag = TRUE;
  long int   max_decrease;
  long int   rchosen = n, schosen = n;
  long int   r, s;

  VERYVERBOSE ( printf("best imp, symmetric case\n"); )
  if ( make_symmetric_flag )
    original_symmetric_factor = 1;
  else
    original_symmetric_factor = 2;

  best_found  = compute_evaluation_function( q );

  /* allocate and prepare matrix with move_values */
  if((move_values = malloc(sizeof(long int) * n * n +
		      sizeof(long int *) * n )) == NULL){
/*      printf("Out of memory, exit."); */
    exit(1);
  }
  for ( k = 0 ; k < n ; k++ ) {
    move_values[k] = (long int *)(move_values + n) + k*n;
  }
  r = rchosen;
  s = schosen;

/*    printf("best 2-opt, value before %ld\n",compute_evaluation_function(q)); */

  while ( improvement ) {
    improvement = FALSE;
    max_decrease = LONG_MAX;
    /* in the first iteration the full neighborhood has to be evaluated */
    if (first_it_flag) {
      first_it_flag = FALSE;
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  tmp = 0;
	  for ( k = 0 ; k < n ; k++ ) {
	    if ( (k != u) && (k != v) ) {
	      tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	    }    
	  }
	  tmp *= original_symmetric_factor;
	  move_values[u][v] = tmp;
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
	    rchosen = u;
	    schosen = v;
/*    	    printf(" max-decr = %ld\n",tmp); */
	  }
	}
      }
    } else {
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  if (u == r || v == s || u == s || v == r) {
	    tmp = 0.;
	    for (k = 0 ; k < n ; k++) {
	      if ( (k != u) && (k != v) ) {
		tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	      }    
	    }
	    tmp *= original_symmetric_factor;
	    move_values[u][v] = tmp;
	  } else { /* change derived from prev iteration, u and v differ from rchosen or schosen */
	    tmp = original_symmetric_factor * ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	      ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] );
	    tmp += move_values[u][v];
	    move_values[u][v] = tmp;
	  }
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
/*    	    printf(" max-decr = %ld\n",tmp); */
	    rchosen = u; /* memorize move */
	    schosen = v; /* memorize move */
	  }	   
	}
      }
    }
    if ( max_decrease < 0 ) {      /* Obj. function value can be improved */
      assert (rchosen < schosen);
      improvement = TRUE;
      best_found += max_decrease;
      swap(rchosen,schosen,q);    
/*        printf(" max-decr = %ld\n",max_decrease); */
      r = rchosen; /* memorize previously done move */
      s = schosen; /* memorize previously done move */

      VERYVERBOSE ( printf("improvement %ld, best_found %ld, exchange %ld and %ld\n",max_decrease, best_found,rchosen, schosen); )
    }
  }
/*    printf("best 2-opt, value after %ld\n",compute_evaluation_function(q)); */

  free ( move_values );
}

void best_2_opt_asymmetric_tabu ( long int * q ) {
/*    
      FUNCTION:      best improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      fast evaluation of moves with additional matrix move_values
                     the first local search iteration is slow 
                     local search looks for best neighboring solution in each iteration
*/

  long int   u, v, k;
  long int   r, s;                      /* memorize which is best move in previous iteration */
  long int   tmp;
  long int   original_symmetric_factor; /* = 2: original symmetric instance
					   = 1: original asymmetric instance 
					*/

  VERYVERBOSE ( printf("best imp, asymmetric case\n"); )
  if ( make_symmetric_flag )
    original_symmetric_factor = 1;
  else
    original_symmetric_factor = 2;

  r = rchosen;
  s = schosen;

  /* in the first local search iteration the full neighborhood has to be evaluated */
  if (first_it_flag) {
    first_it_flag = FALSE;
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	tmp = 0;
	for ( k = 0 ; k < n ; k++ ) {
	  if ( (k != u) && (k != v) ) {
	    tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
	      d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
	      d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
	      d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	  }    
	}
	tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	  d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	  d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
	  d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	move_values[u][v] = tmp;
      }
    }
  } else {
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	if (u == r || v == s || u == s || v == r) {
	  tmp = 0;
	  for ( k = 0 ; k < n ; k++ ) {
	    if ( (k != u) && (k != v) ) {
	      tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
		d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
		d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
		d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	    }    
	  }
	  tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	    d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	    d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
	    d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	  move_values[u][v] = tmp;
	} else { /* change derived from move_values */
	  tmp = ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	    ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] )
	      + ( d[u][r] - d[v][r] + d[v][s] - d[u][s] ) *
	    ( f[q[u]][q[s]] - f[q[v]][q[s]] + f[q[v]][q[r]] - f[q[u]][q[r]] );
	  tmp += move_values[u][v];
	  move_values[u][v] = tmp;
	}
      }
    }
  }
}

void best_2_opt_symmetric_tabu ( long int * q ) {
/*    
      FUNCTION:      best improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      faster evaluation of moves with additional matrix move_values
*/

  long int   u, v, k;
  long int   r, s;
  long int   tmp;
  long int   original_symmetric_factor;

  VERYVERBOSE ( printf("best imp, symmetric case\n"); )
  if ( make_symmetric_flag )
    original_symmetric_factor = 1;
  else
    original_symmetric_factor = 2;

  r = rchosen;
  s = schosen;

  /* in the first iteration the full neighborhood has to be evaluated */
  if (first_it_flag) {
    first_it_flag = FALSE;
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	tmp = 0;
	for ( k = 0 ; k < n ; k++ ) {
	  if ( (k != u) && (k != v) ) {
	    tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	  }    
	}
	move_values[u][v] = tmp * original_symmetric_factor;
      }
    }
  } else {
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	if (u == r || v == s || u == s || v == r) {
	  tmp = 0.;
	  for (k = 0 ; k < n ; k++) {
	    if ( (k != u) && (k != v) ) {
	      tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	    }    
	  }
	  move_values[u][v] = tmp * original_symmetric_factor;
	} else { /* change derived from prev iteration, u and v differ from rchosen or schosen */
	  tmp = original_symmetric_factor * ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	    ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] );
	  tmp += move_values[u][v];
	  move_values[u][v] = tmp;
	}
      }
    }
  }
}

long int choose_tabu_length( )
/*    
      FUNCTION:      choose a new tabu list length
      INPUT:         instance size
      OUTPUT:        none
      (SIDE)EFFECTS: changes the value of variable tabu__list_length
*/
{
  double  help, min, max;

  min      = 0.5  * n;
  max      = 1.1  * n;
  help     = min + ran01( &seed ) * (max - min);
  return   MAX(2,((long int) help));
}

long int ** init_tabu( )
/*    
      FUNCTION:      allocate and initialize the tabu values
      INPUT:         instance size
      OUTPUT:        pointer to matrix, has to be freed before next trial
      (SIDE)EFFECTS: current assignment is modified
*/
{
  long int i, j;
  long int     **matrix;

  VERYVERBOSE ( printf("initialze matrix of tabu values\n"); ) 

  if((matrix = malloc(sizeof(long int) * n * n +
		      sizeof(long int *) * n )) == NULL){
/*      printf("Out of memory, exit."); */
    exit(1);
  }

  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (long int *)(matrix + n) + i*n;
    for (j = 0 ; j < n ; j++ ) {
      matrix[i][j] = 0;
    }
  }
  /* print_matrix(matrix); */
  return matrix;
}

void  make_tabu( long int * q, long int iteration, long int r, long int s )
/*    
      FUNCTION:      make an invers move tabu (forbids a location for an item)
      INPUT:         pointer to some assignment, iteration number, two locations involved in move
      OUTPUT:        none
      (SIDE)EFFECTS: matrix of tabu_values is modified
*/
{
  tabu_values[r][q[r]] = iteration + tabu_list_length;
  tabu_values[s][q[s]] = iteration + tabu_list_length;
}

void select_move( long int *q, long int current, long int iteration, long int best_so_far )
/*    
      FUNCTION:      select the best move which is not tabu or is aspired
      INPUT:         pointer to some assignment, obj_f_val of current solution, iteration number
      OUTPUT:        none
      (SIDE)EFFECTS: global variables rchosen and schosen are assigned the locations involved in the move
*/
{
  long int   i, j;
  long int   max_decrease;
  long int   taboo, aspired;
  
  max_decrease = LONG_MAX;
  rchosen = n; schosen = n;
  for ( i = 0 ; i < n - 1 ; i++) {
    for (j = (i+1); j < n ; j++) {
      if ( (tabu_values[i][q[j]] > iteration ) && (tabu_values[j][q[i]] > iteration ) )
	taboo = TRUE;
      else 
	taboo = FALSE;
      if ( ((current + move_values[i][j]) < best_so_far ) && taboo )
	aspired = TRUE;
      else if ( tabu_values[i][q[j]] < aspirating_iteration && 
	      tabu_values[j][q[i]] < aspirating_iteration )
	aspired = TRUE;
      else 
	aspired = FALSE;
      if ( (move_values[i][j] < max_decrease && !taboo) || aspired ) {
	rchosen = i;
	schosen = j;
	if ( aspired && !taboo ) {
	  max_decrease = LONG_MIN;
	}
	else 
	  max_decrease = move_values[i][j];
      }
    }
  }
  /* print_matrix (tabu_values); */
  assert ( rchosen >= 0 && rchosen < n);
  assert ( schosen >= 0 && schosen < n);
}

void tabu_search( long int *s )
/*    
      FUNCTION:      perform robust tabu search
      INPUT:         pointer to initial solution, instance size
      OUTPUT:        pointer to best solution
      (SIDE)EFFECTS: 
*/
{
  long int  i;
  long int  iterations, obj_f_value, *b, *q;
  long int  best_so_far = 0;

  VERYVERBOSE ( printf("tabu search\n"); )

  move_values = init_move_values( );
  tabu_values = init_tabu( );


  /* b = generate_random_vector( n ); */

  b = malloc( n * sizeof(long int) );
  q = malloc( n * sizeof(long int) );
  for ( i = 0 ; i < n ; i++ ) {
    b[i] = s[i];
    q[i] = s[i];
  }
  /*  q = malloc( n * sizeof(long int) );
  for ( i = 0 ; i < n ; i++ ) 
    q[i] = b[i];
  */
  obj_f_value       = compute_evaluation_function( b );
  best_so_far       = obj_f_value;
  first_it_flag     = TRUE;  
  rchosen           = schosen      = 0;
  /* init_tabu(); */
  iterations        = 1;
  tabu_list_length  = choose_tabu_length();

  while ( iterations < tabu_search_length ) {

    aspirating_iteration = iterations - iterations_before_aspiration;

    if ( d_symmetric_flag && f_symmetric_flag && null_diagonal_flag )
      best_2_opt_symmetric_tabu ( q );
    else if ( make_symmetric_flag )
      best_2_opt_symmetric_tabu ( q );
    else
      best_2_opt_asymmetric_tabu ( q );

    select_move( q , obj_f_value, iterations , best_so_far );

    make_tabu( q, iterations, rchosen, schosen ); /* make_tabu has to go before swap */

    swap( rchosen, schosen, q );

    obj_f_value += move_values[rchosen][schosen];

    if ( obj_f_value < best_so_far ) {
      best_so_far = obj_f_value;
      /* printf("Best %ld\t time %f\n",best_found,elapsed_time( VIRTUAL )); */
      /* fprintf(output_file,"best %ld\tit %ld\tsteps %ld\ttime %f\n",best_found,iterations,iterations,elapsed_time( VIRTUAL )); */
      time_best_found = elapsed_time( VIRTUAL );
      /* iteration_best_found = iterations; */
      for ( i = 0 ; i < n ; i++ ) 
	b[i] = q[i];
    }

    if ( !(iterations % ( long int )(2.2 * n + 4))) {
      tabu_list_length = choose_tabu_length();
      VERYVERBOSE ( printf(" tabu_length = %ld, iteration %ld\n",tabu_list_length, iterations); )
    }
    iterations++;

  }
/*    print_solution ( b ); */
  /* fprintf(output_summary_file,"best %ld\titeration %ld\ttime to best %.2f\ttotal time %.2f\n",best_found, iteration_best_found, time_best_found,elapsed_time( VIRTUAL )); */

  if ( compute_evaluation_function( b ) != best_so_far )
    fprintf(stderr,"Some error must have occurred in local search routine,\n values do not match\n");

  
  for ( i = 0 ; i < n ; i++ ) 
    s[i] = b[i];
  
  free ( b );
  free ( q );
  free ( move_values );
  free ( tabu_values );

}

void first_2_opt_asymmetric ( long int * q ) {
/*    
      FUNCTION:      first improvement 2-opt local search for asymmetric instances
      INPUT:         pointer to initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENT:       neighborhood is scanned in random order, this often is helpful if local 
                     search is used many times from similar starting solutions
*/

  long int   improvement = TRUE;
  long int   improve_item = FALSE;
  long int   i, j, u, v, k;
  long int   tmp;
  long int   *x;            /* random vector to scan neighborhood in random order */

  VERBOSE ( printf("first imp, asymmetric case\n"); )

  x = generate_random_vector( n );
  while ( improvement ) {
    improvement = FALSE;
    for ( i = 0 ; i < n ; i++) {
      u = x[i];
      if ( dlb_flag && dlb[u] )
 	continue;
      improve_item = FALSE;
      for ( j = 0 ; j < n ; j++) {
	v = x[j];
	if (u == v)
	  continue;
	tmp = 0;
	for ( k = 0 ; k < n ; k++ ) {
	  if ( (k != u) && (k != v) ) {
	    tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
	           d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
	           d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
	           d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	  }    
	}
	tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	       d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	       d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
	       d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	if (tmp < 0) {
	  improvement = TRUE;
	  improve_item = TRUE;	  
	  dlb[u] = FALSE;
	  dlb[v] = FALSE;	  
	  swap( u, v, q);
	  VERYVERBOSE ( printf("improvement %ld, best_known %ld\n",tmp, best_value); )
	}
      }
      if ( !improve_item )
	dlb[u] = TRUE;
    }
  }
  free ( x );
}

void first_2_opt_symmetric ( long int *q ) {
/*    
      FUNCTION:      first improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENT:       neighborhood is scanned in random order
*/

  long int   improvement = TRUE;
  long int   improve_item = FALSE;
  register long int i, j, u, v, k;
  long int   tmp;
  long int   *x;  /* scan neighborhood in random order */
  long int   original_symmetric_factor; /* = 2: original symmetric instance
					   = 1: original asymmetric instance 
					*/

  
  VERBOSE (  printf("first imp, symmetric case\n"); )
  if ( make_symmetric_flag )
    original_symmetric_factor = 1; /* compensation because of not dividing matrix by 2 */
  else
    original_symmetric_factor = 2;
  improvement   = TRUE;
  x = generate_random_vector( n );
  while ( improvement ) {
    improvement = FALSE;
    for ( i = 0 ; i < n ; i++ ) {
      u = x[i];
      if ( dlb_flag && dlb[u] )
 	continue;
      improve_item = FALSE;
      for ( j = 0 ; j < n ; j++ ) {
	v = x[j];
	if (u == v)
	  continue;
	tmp = 0;
	for ( k = 0 ; k < n ; k++ ) {
	  if ( (k != u) && (k != v) ) {
	    tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	  }    
	}
   	tmp *= original_symmetric_factor;
	if (tmp < 0) {
	  improvement = TRUE;
	  improve_item = TRUE;	  
	  dlb[u] = FALSE;
	  dlb[v] = FALSE;	  
	  swap( u, v, q);
	  VERYVERBOSE ( printf("improvement %ld, best_value %ld\n",tmp, best_value); )
	}
      }
      if ( improve_item )
	dlb[u] = FALSE;
      else
	dlb[u] = TRUE;	
    }
  }
  free ( x );
}

void acceptance_criterion()
{
  long int  delta;

  delta = shadow_value - s_value;

  if ( acceptance_flag < 7 ) {
    T_current = T_constant * s_value;
/*      printf("current Temp %f\n",T_current); */
  }

/*    printf("current Temp %f\n",T_current); */


  if ( delta < 0 ) {     /* shadow is better than s */ 

    VERYVERBOSE( printf("Improvement, accept new solution%ld, old %ld\n",shadow_value,s_value); )

    copy_from_to( shadow, s);
    s_value = shadow_value;
    k_var = k_min;   

  }
  else if ( delta == 0 ) {

    copy_from_to( shadow, s); /* accept new sol if same value, plateaus possible */
    s_value = shadow_value;
    k_var++;

  }
  else if ( delta > 0 ) {      

    VERYVERBOSE( printf("delta = %ld, current = %ld, k-var = %ld\n",delta,shadow_value,k_var); )

    k_var++; /* always increase k_var if worse solution */

    if ( ran01( &seed ) <= pow(E,( - (double)delta / T_current ))) {
      /* Accept worse solution */
      k_var = k_min;
      copy_from_to( shadow, s);
      s_value = shadow_value;
      no_accepted++;
      VERBOSE ( printf(" Accept new worse solution !!!\n"); )
    }
/*      else */
/*        printf(" Rejected worse solution !!!\n"); */
  }

  if ( k_var > k_max ) {
    k_var = k_min;
  }

}

void    perturbation( long int *s )
{
  long int  i, j, h, unit, help;
  long int  *used, mem_first, mem_first_pos, mem_last, mem_last_pos, tot_assigned;

  VERBOSE ( printf("perturbation size %ld, k_exact %ld\n",k_var,k_exact); )

  used = malloc( n * sizeof(long int) );
  for ( i = 0 ; i < n ; i++ )
    used[i] = FALSE;
  i     = (long int) (n * ran01( &seed ));
/*    printf(" i %ld",i); */
  used[i] = TRUE;
  mem_first = s[i];
  mem_first_pos = i;
  mem_last  = mem_first;
  mem_last_pos = i;
  tot_assigned = 1;
  DEBUG ( assert (k_var <= k_max ); )
  DEBUG ( assert (k_var < n ); )
  for ( h = 1 ; h < k_var ; h++ ) {
    j      = ran01( &seed ) * (n - tot_assigned); 
/*      printf(" j %ld\n",j); */
    unit   = 0;
    while (used[unit]) {
      unit++;
    } /* 0-th free location */
    while (j-- > 0) {
      do {
	unit++;
      } while (used[unit]);
    } /* j-th free location */
    /* exchange the items */
/*      printf("exchange items\n"); */
    help = s[unit];
    s[unit] = mem_last;
    dlb[unit] = FALSE;
    mem_last_pos = unit;
    mem_last = help;
    used[unit] = TRUE;
    tot_assigned++;
  }
  s[mem_first_pos] = mem_last;
  free ( used );
}


void local_search( long int *s ) {

  VERYVERBOSE ( printf("local opt solution\n"); )  
    if ( d_symmetric_flag && f_symmetric_flag && null_diagonal_flag ) {
      if ( local_search_flag == 2 )
	best_2_opt_symmetric( s );
      else if ( local_search_flag == 0 || local_search_flag == 1 )
	first_2_opt_symmetric ( s );
      else if ( local_search_flag >= 3 &&  local_search_flag <= 5 )
	tabu_search( s );
    }
    else if ( make_symmetric_flag ) {
      if ( local_search_flag == 2 )
	best_2_opt_symmetric( s );
      else if ( local_search_flag == 0 || local_search_flag == 1 )
	first_2_opt_symmetric ( s );
      else if ( local_search_flag >= 3 &&  local_search_flag <= 5 )
	tabu_search( s );
    }
    else {
      if ( local_search_flag == 2 )
	best_2_opt_asymmetric( s );
      else if ( local_search_flag == 0 || local_search_flag == 1 )
	first_2_opt_asymmetric ( s );
      else if ( local_search_flag >= 3 &&  local_search_flag <= 5 )
	tabu_search( s );
    }

  VERYVERBOSE ( printf(".. done\n"); )  
  
}

void ils()
{
  long int h;

  VERYVERBOSE ( printf("ils\n"); )
/*      for ( h = 0 ; h < n ; h++ ) */
/*        dlb[h] = TRUE; */

/*    if ( (re_heat_start_it + k_max - iterations >=  3) ) { */
/*      k_min = MAX(3, re_heat_start_it + k_max - iterations);  */
/*    } */
/*    if ( k_min > 3 ) { */
/*      printf(" k_min = %ld\n",k_min); */
/*    }  */


  copy_from_to( s, shadow);

  /* respect for fixed perturbation size for perturbation_flag smaller than 4 */
  if ( perturbation_flag <= 4 )
    k_var = k_exact;
/*    printf(" k_var = %ld\n",k_var); */
  perturbation( shadow );

  shadow_value = compute_evaluation_function( shadow );
  VERYVERBOSE ( printf("shadow_value %ld\n",shadow_value); )
/*    printf("shadow_value %ld\n",shadow_value); */
/*      print_solution( dlb ); */
/*      printf("k_var %ld",k_var); */

  local_search ( shadow );
  
/*    print_solution( dlb ); */

/*    getchar(); */

  shadow_value = compute_evaluation_function( shadow );
/*    difference = calculate_difference(s, shadow); */
/*    printf("shadow_value %ld\n",shadow_value); */
/*    print_solution( shadow ); */
/*    getchar();  */


  if ( shadow_value < best_value ) {
    time_best_found = elapsed_time( VIRTUAL ); 
    best_value = shadow_value;
    best_since_restart = shadow_value;
    iteration_best_found = iterations;
    last_improvement = iterations;
    VERBOSE ( if ( compute_evaluation_function( shadow ) != best_value) {
      printf(" Fehler im Programm, help: %ld, best_value: %ld\n",h,best_value); */
      exit(0);
    } 
    else 
      printf(" solution does verify \n"); )
/*      printf("best %ld\tcycle %ld\tk_var = %ld\ttime %.2f\n",best_value,iterations,k_var,time_best_found); */
      printf("best %ld\ttime %.2f\titerations %ld\n",best_value,time_best_found,iterations);
    fprintf(output_file,"best %ld\tcycle %ld\tsteps %ld\ttime %.2f\tk_var %ld\n",best_value,iterations,iterations,time_best_found,k_var);
    copy_from_to( shadow, best);
    copy_from_to( shadow, s);
    s_value = shadow_value;
    k_var = k_min;
    for ( h = 0 ; h < n ; h++ )
      dlb[h] = TRUE;
  }
  
  if ( shadow_value < best_since_restart ) {
    best_since_restart = shadow_value;
    last_improvement = iterations;
    k_var = k_min;
    copy_from_to( shadow, s );
    for ( h = 0 ; h < n ; h++ )
      dlb[h] = TRUE;
    VERBOSE ( printf("best since restart: %ld\tcycle %ld\tk_var = %ld\ttime = %.2f\n",best_since_restart,iterations,k_var,elapsed_time( VIRTUAL )); )
  }

  acceptance_criterion();


}

void usage()
/*    
      FUNCTION:      output usage of program
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: 
*/
{
  fprintf(stderr,"\nRobust Tabu Search, V1.0\n");
  fprintf(stderr,"usage: [- see below -] [-i dat_file] (given are long and short options)\n");
  fprintf(stderr,"--trials       -r #\t number of trials to be run on one instance\n");
  fprintf(stderr,"--time         -t #\t maximum time limit (not used for the moment)\n");
  fprintf(stderr,"--localsearch  -l #\t local search type to be run\n");
  fprintf(stderr,"                  #\t 0: first improvement, don't look bits\n");
  fprintf(stderr,"                  #\t 1: first improvement, no don't look bits\n");
  fprintf(stderr,"                  #\t 2: best improvement\n");
  fprintf(stderr,"                  #\t 3: tabu search, 1 * n\n");
  fprintf(stderr,"                  #\t 4: tabu search, 4 * n\n");
  fprintf(stderr,"                  #\t 5: tabu search, 10 * n\n");
  fprintf(stderr,"--perturbation -p #\t perturbation type to be run\n");
  fprintf(stderr,"                  #\t 0: size n / 16, fixed\n");
  fprintf(stderr,"                  #\t 1: size n / 8, fixed\n");
  fprintf(stderr,"                  #\t 2: size n / 4, fixed\n");
  fprintf(stderr,"                  #\t 3: size n / 2, fixed\n");
  fprintf(stderr,"                  #\t 4: size 3 n / 4, fixed\n");
  fprintf(stderr,"                  #\t 5: size 0.9 * n, fixed\n");
  fprintf(stderr,"                  #\t 6: size n / 16, variable\n");
  fprintf(stderr,"                  #\t 7: size n / 8, variable\n");
  fprintf(stderr,"                  #\t 8: size n / 4, variable\n");
  fprintf(stderr,"                  #\t 9: size n / 2, variable\n");
  fprintf(stderr,"                  #\t 10: size 3 n / 4, variable\n");
  fprintf(stderr,"                  #\t 11: size 0.9 * n, variable\n");
  fprintf(stderr,"--acceptance   -a #\t acceptance type to be run\n");
  fprintf(stderr,"                  #\t 0: ConstTemp, T = 0 (Better)\n");
  fprintf(stderr,"                  #\t 1: ConstTemp, T = 0.001 : low\n");
  fprintf(stderr,"                  #\t 2: ConstTemp, T = 0.005 : low-med\n");
  fprintf(stderr,"                  #\t 3: ConstTemp, T = 0.015 : med\n");
  fprintf(stderr,"                  #\t 4: ConstTemp, T = 0.05  : high\n");
  fprintf(stderr,"                  #\t 5: ConstTemp, T = 1.0   : random (kind of)\n");
  fprintf(stderr,"                  #\t 6: Restart (T = 0), quick, 1n\n");
  fprintf(stderr,"                  #\t 7: Restart (T = 0), med, 3n\n");
  fprintf(stderr,"                  #\t 8: Restart (T = 0), med, 5n\n");
  fprintf(stderr,"                  #\t 9: Restart (T = 0), high, 10n\n");
  fprintf(stderr,"--tabu         -x #\t length of one tabu search run as multiple of instance size\n");
  fprintf(stderr,"--kmax         -k #\t max kick strength for VNS\n");
  fprintf(stderr,"--Tinit        -m #\t Initial Temperature\n");
  fprintf(stderr,"--alpha        -a #\t alpha for annealing schedule\n");
  fprintf(stderr,"--optimal      -o #\t stop when hitting a solution of that quality\n");
  fprintf(stderr,"--seed         -s #\t Seed for random number generator: if absent,\n");
  fprintf(stderr,"                  #\t a seed is somehow generated...\n");
  fprintf(stderr,"--input        -i #\t input file\n");
  fprintf(stderr,"--help         -h \t\t help: prints this information\n");

}

void set_default_parameters()
/*    
      FUNCTION:      set default parameter values
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: assigns some values to global variables
*/
{
  max_trials         = 1;
  tabu_search_length = 1000;
  time_limit         = 5.;
  max_iterations     = 3;
  T_init             = 5.; /* Accept sol 5% worse than current with prob 1/E */
  T_end              = 0.001;
  T_current          = T_init;
  T_steps            = 25;
  k_min              = 3;
  k_var              = 3;
  k_max              = 30;
  alpha              = 0.9;
  optimal            = 1;
}

void checkOutOfRange(int value, int MIN, int MAX, char *optionName){
  if ((value<MIN)||(value>MAX)){
    fprintf(stderr,"Error: Option `%s' out of range\n",optionName);
    exit(1);
  }
}

void check_options(int argc, char *argv[]) {
  
  while (1){
    int option_index = 0;

    static struct option long_options[] = {
      {"time",        1, 0, 't'},
      {"trials",      1, 0, 'r'},
      {"tabu",        1, 0, 'x'},
      {"localsearch", 1, 0, 'l'},
      {"perturbation",1, 0, 'p'},
      {"acceptance",  1, 0, 'a'},
      {"input",       1, 0, 'i'},
      {"kmax",        1, 0, 'k'},
      {"Tinit",       1, 0, 'm'},
      {"dlb",         1, 0, 'd'},
      {"optimal",     1, 0, 'o'},
      {"seed",        1, 0, 's'},
      {"help",        0, 0, 'h'},
      {0,             0, 0,  0 }
    };
    
    opt = getopt_long(argc, argv, "t:r:l:p:a:i:h:k:m:os:",
		       long_options, &option_index);
    if (opt == -1)
      break;
    
    switch (opt){
    case 0:
      fprintf(stderr,"Error: Confused by input on long option...\n"), exit(1);
      break;
      
    case 't':
      time_limit = atof(optarg);
      break;
      
    case 'r':
      max_trials = atol(optarg);
      checkOutOfRange(max_trials,1,1000,"r (trials)");
      break;
      
    case 'l':
      local_search_flag = atol(optarg);
      checkOutOfRange(local_search_flag,0,5,"l (local_search_flag)");
      break;

    case 'p':
      perturbation_flag = atol(optarg);
      checkOutOfRange(perturbation_flag,0,11,"l (perturbation_flag)");
      break;

    case 'a':
      acceptance_flag = atol(optarg);
      checkOutOfRange(acceptance_flag,0,9,"l (acceptance_flag)");
      break;

    case 'x':
      tabu_search_length = atol(optarg);
      checkOutOfRange(tabu_search_length,0,100000,"l (tabu_search_length)");
      break;

    case 'k':
      k_max = atol(optarg);
      checkOutOfRange(k_max,10,500,"k (k_max)");
      break;

    case 'm':
      T_init = atof(optarg);
      checkOutOfRange(T_init,0.0,1000000.,"m (Tinit)");
      break;

    case 'o':
      optimal = atol(optarg);
/*        checkOutOfRange(optimal,0,LONG_MAX/2,"o (optimal)"); */
      break;

    case 'd':
      dlb_flag = atol(optarg);
      checkOutOfRange(max_trials,0,1,"d (dlb_flag)");
      break;

    case 's':
      seed = atol(optarg);
      /* MAO: I'm not checking if in range... I trust you :-) */
      break;

    case 'i':
      strncpy(name_buf,optarg,strlen(optarg));
      break;

    case 'h':
      usage();
      fprintf(stderr,"No more help yet, sorry.\n"), exit(1);
      break;
      
    default:
      fprintf(stderr,"Error: Confused by input... %ld\n",opt), exit(1);
      
    }
  }
}  

void init_program( long int argc, char *argv[] )
/*    
      FUNCTION:      all the stuff for starting the program: input, parameters, etc.
      INPUT:         argument list from program call
      OUTPUT:        none
      (SIDE)EFFECTS: 
*/
{
  FILE *qap_file;
  long int i;
  char buf[LINE_BUF_LEN], *c;    /* input file name */

  setbuf(stdout,NULL);
/*    printf(PROG_ID_STR); */

  /* MAO: change if negative... */
  if (seed<0) 
    seed = (long int) time(NULL); /* initialize random number generator */

  /* MAO: I want to know the seed!!! */
  printf("seed %ld\n", seed);

  if( (qap_file = fopen(name_buf, "r")) == NULL) {
    fprintf(stderr, "error opening input file %s\n",name_buf);
    exit(1);	
  }
  
  read_problem_size( qap_file );
  
  read_best_known_value( qap_file );
  
  d = read_matrix( qap_file, n);
  d_symmetric_flag = check_symmetry ( d, n );
  null_diagonal_flag = check_null_diagonal ( d, n );   
  VERYVERBOSE ( print_matrix( d ); )
    
  f = read_matrix( qap_file, n);
  f_symmetric_flag = check_symmetry ( f, n );
  if ( null_diagonal_flag )
    ; /* if one matrix has already null diagonal we need not check the other */
  else
    null_diagonal_flag = check_null_diagonal ( f, n );
  VERYVERBOSE ( print_matrix( f ); )
    
  VERBOSE ( printf("d_symmetric_flag %ld, f_symmetric_flag %ld, null_diagonal_flag %ld\n",d_symmetric_flag, f_symmetric_flag, null_diagonal_flag); )

  make_symmetric_flag = XOR(d_symmetric_flag,f_symmetric_flag);
  if ( make_symmetric_flag && null_diagonal_flag ) {
    if ( !d_symmetric_flag )
      make_matrix_symmetric ( d, n );
    else if ( !f_symmetric_flag )
      make_matrix_symmetric ( f, n );
    else {
      fprintf(stderr,"One matrix should have been symmetric\n");
      exit(1);
    }
  }
     
  if ( ( c = strrchr(name_buf,'/') ) != NULL ) {
    c++;
    i = 0;
    while ((buf[i] = *c) != '\0') {
      i++;
      c++;
    }
/*      printf("%s\n",buf); */
    strcpy(name_buf, buf);
  }
  
  sprintf(buf,"%s.LSMC.v.1.0.t%.2f.rep",name_buf,time_limit);
  
  if( (output_file = fopen(buf, "w")) == NULL) {
    fprintf(stderr, "error opening output file %s\n","qap.rep");
    exit(1);	
  }
  
  sprintf(buf,"%s.LSMC.MH-ls.v.1.0.t%.2f.sum",name_buf,time_limit);
  
  if( (output_summary_file = fopen(buf, "w")) == NULL) {
    fprintf(stderr, "error opening output file %s\n","qap.sum");
    exit(1);	
  }
 
  /* finally set these parameters to the right value, they were defined as multiples of instance size */

  tabu_search_length *= n;
  iterations_before_aspiration = (long int)(3. * n * n);
  best = malloc( n * sizeof(long int) );
  shadow = malloc( n * sizeof(long int) );
  best_in_try = malloc( max_trials * sizeof(long int) );
  best_found_at = malloc( max_trials * sizeof(long int) );
  t_best_found = malloc( max_trials * sizeof(double) );
/*    k_max = MIN( 2. * n / 3., 50); */
}

void update_temp( )
{

  T_current = alpha * T_current;
  if ( T_current < T_end ) {
    if ( iterations - last_improvement < 2 * n && no_accepted > 1 ) {
      T_current = alpha * T_current;
      if ( no_accepted > 3 ) {
	T_end = T_end / 2.;
/*  	printf("New end temp %f\n",T_end); */
      }
    }
    else {
      T_current = T_init;
/*        printf("INCREASE TEMP AGAIN, Iterations %ld, current %ld\n",iterations,s_value); */
      re_heat_start_it = iterations;
    }
  }
/*    printf("New Temp %f\n",T_current); */
}

void start_trial( )
{
  long int i;

  VERYVERBOSE ( printf("start trial .. \n"); )
  start_timers();
  time_best_found = 0.0;
  iteration_best_found = 0;
  last_improvement = 0;
  iterations = 0;
  fprintf(output_file,"begin try %li\n",trials);
  VERBOSE ( printf("begin try %li\n",trials); )
  s = generate_random_vector( n );
  dlb = malloc( n * sizeof(long int) );
  for ( i = 0 ; i < n ; i++ ) {
    dlb[i] = FALSE;
  }
/*    k_var = 3; */
  no_accepted = 0;
  re_heat_start_it = 0;
  VERYVERBOSE ( printf(".. done\n"); )  
}

void end_trial( )
{

/*    printf("END TRY %ld, BEST FOUND %ld, ITERATIONS %ld\n",trials,best_value,iterations);  */
  fprintf(output_file,"end try %ld\t total time %f\n",trials, elapsed_time( VIRTUAL ));
  fprintf(output_summary_file,"best %ld\titeration %ld\ttime to best %.2f\ttotal time %.2f\n",best_value, iteration_best_found, time_best_found,elapsed_time( VIRTUAL ));
  best_in_try[trials-1] = best_value;
  best_found_at[trials-1] = iteration_best_found;
  t_best_found[trials-1] = time_best_found;

  VERBOSE ( printf("end try %li\n",trials); )
  free ( s );
}

void end_program( )
{
  fprintf(output_file,"end problem %s\n",name_buf);
  VERBOSE ( printf("end problem %s\n",name_buf); )
}

void determine_restart() {

  restart_iterations = (long int)( (log(T_end) - log(T_init)) / log (alpha)) * T_steps * 2;
/*    printf("restart_iterations %ld\n",restart_iterations); */


}


void test_temperatures ()
{
  long int temp_fixed_flag, temp_init_flag, temp_end_flag;
  long int i, j;

/*    printf("test_temps\n"); */
  j = 0;
  temp_fixed_flag = FALSE;
  temp_init_flag = FALSE;
  temp_end_flag = FALSE;
  T_init = 0.02 * s_value;
  T_end =  0.000005 * s_value;
  
  while ( !temp_fixed_flag ) {
    j++;
/*      printf("round %ld\n",j); */
    if ( !temp_init_flag ) {
      no_accepted = 0;
      T_current = T_init;
      for ( i = 1; i <= 100; i++ ) {
	ils();
	iterations++;
	if ( termination_condition() )
	  return;
      }

/*        printf("acceptance ratio %f\n",(double)no_accepted / 100.); */

      if ( 20 <= no_accepted && no_accepted <= 30 )
	temp_init_flag = TRUE;
      else if ( no_accepted < 20 ) {
	temp_init_flag = FALSE;
	T_init *= 2;
/*  	printf("T_init %f\n",T_init); */
      }
      else if ( no_accepted > 30 ) {
	temp_init_flag = FALSE;
	T_init /= 2;
/*  	printf("T_init %f\n",T_init); */
      }
    }
    if ( !temp_end_flag ) {
   
      no_accepted = 0;
      T_current = T_end;
      for ( i = 1; i <= 100; i++ ) {
	ils();
	iterations++;
      if ( termination_condition() )
	return;
      }

/*        printf("acceptance ratio %f\n",(double)no_accepted / 100.); */

      if ( 2 <= no_accepted && no_accepted <= 5 )
	temp_end_flag = TRUE;
      else if ( no_accepted < 2 ) {
	temp_end_flag = FALSE;
	T_end *= 2.;
/*  	printf("T_end %f\n",T_end); */
      }
      else if ( no_accepted > 10 ) {
	temp_end_flag = FALSE;
	T_end /= 2.;
/*  	printf("T_end %f\n",T_end); */
      }
    }

    temp_fixed_flag = temp_end_flag * temp_init_flag;
/*      getchar(); */
  }
/*    printf("\nT_INIT %f, T_END %f\n\n",T_init, T_end); */
}


void exit_program( void ) 
{
  long int best_tour_length, worst_tour_length;
  double   t_avgbest, t_stdbest;
  double   avg_sol_quality = 0., avg_cyc_to_bst = 0., stddev_best, stddev_iterations;

  fprintf(output_summary_file,"\n\nMax tries %ld",max_trials);

  best_tour_length = best_of_vector( best_in_try , max_trials );
  worst_tour_length = worst_of_vector( best_in_try , max_trials );
  avg_cyc_to_bst = mean( best_found_at , max_trials );
  stddev_iterations = std_deviation( best_found_at, max_trials, avg_cyc_to_bst );
  avg_sol_quality = mean( best_in_try ,max_trials );
  stddev_best = std_deviation( best_in_try, max_trials, avg_sol_quality);
  t_avgbest = meanr( t_best_found, max_trials );
/*    printf(" t_avgbest = %f\n", t_avgbest ); */
  t_stdbest = std_deviationr( t_best_found, max_trials, t_avgbest);

  printf("%ld\n",best_tour_length);
  fprintf(output_summary_file,"\nAverage-Best: %.2f\t Average-Iterations: %.2f", avg_sol_quality, avg_cyc_to_bst);
  fprintf(output_summary_file,"\nStddev-Best: %.2f \t Stddev Iterations: %.2f", stddev_best, stddev_iterations);
  fprintf(output_summary_file,"\nBest try: %ld \n", best_tour_length);
  fprintf(output_summary_file,"\nAvg.time-best: %.2f stddev.time-best: %.2f\n", t_avgbest, t_stdbest); 

  if ( best_known > 0 ) {
    fprintf(output_summary_file," excess best = %f, excess average = %f, excess worst = %f\n",(double)(best_tour_length - best_known) / (double)best_known,(double)(avg_sol_quality - best_known) / (double)best_known,(double)(worst_tour_length - best_known) / (double)best_known);
  }
}

void restart( void ) {

  long int i;

/*    printf("Restart %ld\n",iterations); */
  free ( s );
  s = generate_random_vector( n );
  for ( i = 0 ; i < n ; i++ ) {
    dlb[i] = FALSE;
  }
  local_search ( s );
  s_value = compute_evaluation_function( s );
  best_since_restart = s_value;
  last_improvement = iterations;
  T_current = T_init;
  re_heat_start_it = iterations;
}
/* --- main program ------------------------------------------------------ */

int main(int argc, char *argv[]) {
  

  /* variables for main programs are declared above */
  

  start_timers();             /* start timing routines */

  set_default_parameters();
  
  check_options(argc, argv);
  
  init_program(argc, argv);   /* initialize all important data */

  if ( local_search_flag < 6 && local_search_flag >= 0 )
    ;
  else {
    fprintf(stderr,"local_search_flag outside limits\n Abort\n"), exit(1);
  }

  if ( local_search_flag == 0 )
    dlb_flag = TRUE;
  else
    dlb_flag = FALSE;

  if ( local_search_flag == 3 )
    tabu_search_length = n;
  else if ( local_search_flag == 4 )
    tabu_search_length = 4 * n;
  else if ( local_search_flag == 5 ) 
    tabu_search_length = 10 * n;
  
/*    printf("perturbation_flag %ld\n",perturbation_flag); */

  switch ( perturbation_flag ) {
    case 0:
      k_exact = n / 16;
      k_min = k_exact; k_max = k_exact; k_var = k_exact;
      break;
    case 1:
      k_exact = n / 8;
      k_min = k_exact; k_max = k_exact; k_var = k_exact;
      break;
    case 2:
      k_exact = n / 4;
      k_min = k_exact; k_max = k_exact; k_var = k_exact;
      break;
    case 3:
      k_exact = n / 2;
      k_min = k_exact; k_max = k_exact; k_var = k_exact;
      break;
    case 4:
      k_exact = 3 * n / 4;
      k_min = k_exact; k_max = k_exact; k_var = k_exact;
      break;
    case 5:
      k_exact = 0.9 * n ;
      k_min = k_exact; k_max = k_exact; k_var = k_exact;
      break;
    case 6:
      k_exact = n / 16;
      k_min = 3;
      k_max = k_exact;
      k_var = k_min;
      break;
    case 7:
      k_exact = n / 8;
      k_min = 3;
      k_max = k_exact;
      k_var = k_min;
      break;
    case 8:
      k_exact = n / 4;
      k_min = 3;
      k_max = k_exact;
      k_var = k_min;
      break;
    case 9:
      k_exact = n / 2;
      k_min = 3;
      k_max = k_exact;
      k_var = k_min;
      break;
    case 10:
      k_exact = 3 * n / 4;
      k_min = 3;
      k_max = k_exact;
      k_var = k_min;
      break;
    case 11:
      k_exact = 0.9 * n ;
      k_min = 3;
      k_max = k_exact;
      k_var = k_min;
      break;
    default:
      fprintf(stderr,"You have to specify a perturbation (option -p, see --help for more info)!!\n Abort\n"), exit(1);
  }      

  switch ( acceptance_flag ) {
    case 0:
      T_constant = 0.00000000000000001;
      T_current = T_constant;
      restart_iterations = INT_MAX;
      break;
    case 1:
      T_constant = 0.001;
      T_current = T_constant;
      restart_iterations = INT_MAX;
      break;
    case 2:
      T_constant = 0.005;
      T_current = T_constant;
      restart_iterations = INT_MAX;
      break;
    case 3:
      T_constant = 0.015;
      T_current = T_constant;
      restart_iterations = INT_MAX;
      break;
    case 4:
      T_constant = 0.05;
      T_current = T_constant;
      restart_iterations = INT_MAX;
      break;
    case 5:
      T_constant = 1.0;
      T_current = T_constant;
      restart_iterations = INT_MAX;
      break;
    case 6:
      T_constant = 0.00000000000000001;
      restart_iterations = n;
      break;
    case 7:
      T_constant = 0.00000000000000001;
      restart_iterations = 3 * n;
      break;
    case 8:
      T_constant = 0.00000000000000001;
      restart_iterations = 5 * n;
      break;
    case 9:
      T_constant = 0.00000000000000001;
      restart_iterations = 10 * n;
      break;
    default:
      fprintf(stderr,"You have to specify an acceptance criterion (option -a, see --help for more info)!!\n Abort\n"), exit(1);
  }      

  print_header ();
/*    printf("before starting trials\n"); */

  assert ( seed < IM );

  for ( trials = 1 ; trials <= 1 ; trials ++ ) {

    start_trial();

    print_solution( s );

    local_search ( s );
    s_value = compute_evaluation_function( s );
/*      printf("s_value %ld\n",s_value); */
    best_value = s_value;
    copy_from_to(s, best);
    copy_from_to(s, shadow);
    time_best_found = elapsed_time( VIRTUAL ); 
    iteration_best_found = iterations;
    last_improvement = iterations;
    printf("best %ld\ttime %.2f\n",best_value,time_best_found);
    fprintf(output_file,"best %ld\tcycle %ld\tsteps %ld\ttime %.2f\tk_var %ld\n",best_value,iterations,iterations,time_best_found,k_var);

    /* T_init = 0.025 * s_value; */
    /* T_current = T_init; */
    /* T_end =  MAX( 2, 0.000005 * s_value); */
    /* T_steps = MIN ( n / 8, 10); */
/*      printf("\nT_INIT %f, T_END %f\n\n",T_init, T_end); */

/*      test_temperatures(); */

/*    if ( acceptance_flag > 6 )
      determine_restart();
*/
    
    while ( !termination_condition() ) {
      
      ils();
      iterations++;

      if ( iterations - last_improvement > restart_iterations ) {
	restart();
      }
    
      /*      if ( !(iterations % T_steps) && acceptance_flag > 6 ) {
  	k_min = 3;
	update_temp();
      }
*/
      /*      if ( !(iterations % 50) && acceptance_flag > 6 ) {
	VERBOSE ( printf("acceptance ratio %f\n",(double)no_accepted / 50.); )
  	printf("acceptance %ld\n",no_accepted);
	if ( no_accepted <= 1 ) {
	  T_current = T_init;
 	  printf("INCREASE TEMP AGAIN, Iterations %ld, current %ld\n",iterations,s_value);
	  re_heat_start_it = iterations;
	  T_current = T_init;
	  printf("INCREASE TEMP AGAIN, Iterations %ld, current %ld\n",iterations,s_value);
	  re_heat_start_it = iterations;
  	  restart();
	  no_accepted = 0;
	}
	no_accepted = 0;
      }
      */
    }

    end_trial();
    fflush ( output_file );
    fflush ( output_summary_file );

  }

  fprintf(output_file,"end problem %s\n",name_buf);

  exit_program();
	
  fflush ( output_file );
  fflush ( output_summary_file );
  
  free ( d );  /* first matrix */
  free ( f );  /* second matrix */
  free ( p );  /* solution vector */
  /* free ( move_values ); */
  /* free ( tabu_values ); */
  return(0);
}











