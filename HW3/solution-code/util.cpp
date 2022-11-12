#include <stdexcept>
#include <math.h>
#include "util.h"

using namespace std;

double Combination(int n, int k)
{
  if(n < k || k < 0){
    throw invalid_argument("Comination: the number of elements should be bigger than selection numbers AND two numbers should be positive\n");
    return EXIT_FAILURE;
  }
  double result = 1e308;
  if (pow(result, 1.0/k) < n) {
    throw invalid_argument("The Combination number may be bigger than the maxium of double precision\n");
    return EXIT_FAILURE;
  }
  double n_d = (double)n;
  result = 1.0;
  int k_small = min(k, n - k);
  for(double j = (double)k_small; j > 0; j--){
    result /= j;
    result *= n_d;
    n_d --;
  }
  return result;
}

double DoubleFactorial(int n){
  if(n< -1)
    throw invalid_argument("DoubleFactorial: Input should be 1 at least\n");
  if(n ==0 || n== -1)
    return 1;
  double result = 1;
  while(n > 1){
    result *= n;
    n -= 2;
  }
  return result;
}
