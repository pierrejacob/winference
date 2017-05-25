#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List swapsweep(IntegerVector permutation, const NumericMatrix Cp, double totalcost){
int n = Cp.cols();
double currentcost, proposedcost;
int perm_i, perm_j;
for (int i = 0; i < n - 1; i ++){
  int perm_i = permutation[i];
  for (int j = i+1; j < n; j++){
    perm_j = permutation[j];
    currentcost = Cp(i, perm_i) + Cp(j, perm_j);
    proposedcost = Cp(i, perm_j) + Cp(j, perm_i);
    if (proposedcost < currentcost){
      permutation[i] = perm_j;
      permutation[j] = perm_i;
      perm_i = perm_j;
      totalcost = totalcost - currentcost + proposedcost;
    }
  }
}
return List::create(Named("totalcost")=totalcost,
                    Named("permutation") = permutation);
}
