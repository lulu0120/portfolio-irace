
/*****************************************************************************/
/*                                                                           */
/*      Version:  1.00   Date: 22/04/96   File: statistics.c                 */
/* Last Version:                          File:                              */
/* Changes:                                                                  */
/* 22/04/96 Created                                                          */
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
#include <stdlib.h>
#include <math.h>
#include "stat.h"

double mean( long int *values, const long int max ) {
  long int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += (double)values[j];
  }
  m = m / (double)max;
  return m;
}

double meanr( double *values, const long int max ) {
  long int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += values[j];
  }
  m = m / (double)max;
  return m;
}

double std_deviation( long int *values, long int i, double mean ) {

  long int j;
  double   help, dev;

  if (i <= 1)
    return 0.;
  dev = 0.;
  for ( j = 0 ; j < i; j++ ) {
    help = ((double)values[j] - mean);
    dev += help * help;
  }
  return sqrt(dev/(double)(i - 1));
}

double std_deviationr( double *values, long int i, double mean ) {

  long int j;
  double   help, dev;

  if (i <= 1)
    return 0.;
  dev = 0.;
  for ( j = 0 ; j < i ; j++ ) {
    help = ((double)values[j] - mean);
    dev += help * help;
  }
  return sqrt(dev/(double)(i - 1));
}

long int best_of_vector( long int *values, long int i ) {

  long int min, k;

  k = 0;
  min = values[k];
  for( k = 1 ; k < i ; k++ ) {
    if( values[k] < min ) {
      min = values[k];
    }
  }
  return min;
}

long int worst_of_vector( long int *values, long int i ) {

  long int max, k;

  k = 0;
  max = values[k];
  for( k = 1 ; k < i ; k++ ) {
    if( values[k] > max ){
      max = values[k];
    }
  }
  return max;
}
