// extern "C" {
#include <R.h>
// }
#include <vector>
#include <algorithm>

using namespace std;

double half_range_mode( double *start, double *end, double beta, int diag ) {

  // The end pointer is one step beyond the data...

  double w, w_prime;
  double *last, *new_start, *new_end;
  vector<int> counts, J;
  vector<double> w_range;
  int i, s, e;
  int N, N_prime, N_double_prime;
  double lo, hi;

  last = end - 1;
  N = end - start;

  // How many elements are in the set? Terminate recursion appropriately...

  switch ( N ) {

  case 1:
    return *start;

  case 2:
    return .5 * ( *start + *last );

  // Main recursive code begins here

  default:
    
    w = beta * ( *last - *start ); 

    // If all values are identical, return immediately...

    if ( w == 0 ) return *start;

    // If we're at the end of the data, counts can only get worse, so there's no point in continuing...
    e = 0;
    for( s = 0; s < N && e < N; s++ ) {
      while ( e < N && start[ e ] <= start[ s ] + w ) { e++; }
      counts.push_back( e - s );
    }

    // Maximum count, and its multiplicity

    N_prime = *( max_element( counts.begin(), counts.end() ) );

    for ( i = 0; i < (int) counts.size(); i++ ) if ( counts[i] == N_prime ) J.push_back( i );
    
    // Do we have more than one maximal interval?

    if ( J.size() == 1 ) { 
      // No... the interval's unique.
      new_start = start + J[0];
      new_end = start + J[0] + N_prime;
    }
    
    else {

      // Yes.. What's the smallest range?
      for ( i = 0; i < (int) J.size(); i++ ) w_range.push_back( start[ J[i] + N_prime - 1 ] - start[ J[i] ] );
      w_prime = *( min_element( w_range.begin(), w_range.end() ) );

      // Set new start and end. We skip the more cumbersome V.min and V.max of the Bickel algorithm

      i = 0;
      while( w_range[ i ] > w_prime ) i++;
      new_start = start + J[i];
      new_end = start + J[i] + N_prime;
      
      // If there are any more maximal-count, minimal-range intervals, adjust
      // new_end accordingly.
      for ( i++; i < (int) J.size(); i++ ) if ( w_range[ i ] == w_prime ) new_end = start + J[i] + N_prime; 
      
    }
    
    // Adjustments in rare cases where the interval hasn't shrunk. Trim one end,
    // the other, or both if lo == hi. Originally, this was inside the else
    // block above. With discrete data with a small number of levels, it is
    // possible, however for |J| = 1 AND N_double_prime = N, leading to an
    // infinite recursion.
    
    N_double_prime = new_end - new_start;
    
    if (N_double_prime == N ) {
      lo = new_start[1] - new_start[0];
      hi = new_start[ N - 1 ] - new_start[ N - 2 ];
      if ( lo <= hi ) { new_end--; }
      if ( lo >= hi ) { new_start++; }
    }

    // Diagnostic output if requested

    if (diag) Rprintf( "N = %i, N'' = %i, w = %.4f, |J| = %i\n", N, N_double_prime, w, J.size() );

    // Clean up and then go in recursively

    counts.clear(); J.clear(); w_range.clear();

    return half_range_mode( new_start, new_end, beta, diag ); 

  }
  
}




extern "C" {

  void half_range_mode( double *data, int *n, double *beta, int *diag, double *M ) {

    // We assume that that data is already sorted for us...

    *M = half_range_mode( data, data + *n, *beta, *diag );

  }
    
}
