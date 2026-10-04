#include <Rcpp.h>
using namespace Rcpp;
using std::vector;
#include "ARIbrain.h"


// Calculates the size of the concentration set at a fixed alpha
// [[Rcpp::export]]
int findConcentration(Rcpp::NumericVector& allp,         // vector of all p-values (unsorted!)
                      Rcpp::IntegerVector& ORD,          // sorted orders for non-decreasing p-values
                      double               simesfactor,  // simesfactor at h(alpha)
                      int                  h,            // h(alpha)
                      double               alpha,        // alpha itself
                      int                  m)            // size of the problem
  
{
  // from m-h we increase z until we fulfil the condition
  int z = m - h;
  if (z > 0)  // h=m implies z=0
  {
    while ((z < m) & (simesfactor * allp[ORD[z-1] - 1] > (z - m + h + 1) * alpha))
    {
      z++;
    }
  }
  return z;
}


// Find function for disjoint set data structure
// (1) Old recursive version
// int Find(int              x,
//          std::vector<int> &parent)
// {
//   if (parent[x] != x)
//   {
//     parent[x] = Find(parent[x], parent);
//   }
//   
//   return parent[x];
// }
// (2) iterative version (more stable)
int Find(int               x,
         std::vector<int>& parent)
{
  while (parent[x] != x)
  {
    parent[x] = parent[parent[x]];
    x         = parent[x];
  }
  
  return x;
}



// Union function for disjoint set data structure
// Extra: we keep track of the lowest entry of each disjoint set
// That way we can find the lower set to merge with
void Union(int x,
           int y,
           std::vector<int>& parent,
           std::vector<int>& lowest,
           std::vector<int>& rank)
{
  int xRoot = Find(x, parent);
  int yRoot = Find(y, parent);
  
  // if x and y are already in the same set (i.e., have the same root or representative)
  if (xRoot == yRoot) return; // Note: this never happens in our case
  
  // x and y are not in same set, so we merge
  if (rank[xRoot] < rank[yRoot])
  {
    parent[xRoot] = yRoot;
    lowest[yRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
  else if (rank[xRoot] > rank[yRoot])
  {
    parent[yRoot] = xRoot;
    lowest[xRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
  else
  {
    parent[yRoot] = xRoot;
    rank[xRoot]++;
    lowest[xRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
}

// Calculate the category for each p-value
int getCategory(double p,               // p-value for which we need the category
                double simesfactor,     // simesfactor at h(alpha)
                double alpha,           // alpha itself
                int    m)               // size of the problem
{
  if (p==0 || simesfactor==0)
    return 1;
  else
    if (alpha == 0)
      return m+1;
    else
    {
      double cat = (simesfactor / alpha) * p;
      return static_cast<int> (std::ceil(cat));
    }
}

// Calculates the lower bound to the number of false hypotheses
// Implements the algorithm based on the disjoint set structure
// [[Rcpp::export]]
Rcpp::IntegerVector findDiscoveries(Rcpp::IntegerVector& idx,          // indices in set I (from 1)
                                    Rcpp::NumericVector& allp,         // all p-values
                                    double               simesfactor,  // simesfactor at h(alpha)
                                    int                  h,            // h(alpha)
                                    double               alpha,        // alpha
                                    int                  k,            // size of I
                                    int                  z,            // size of concentration set
                                    int                  m)            // size of the problem
{
  // calculate categories for the p-values
  std::vector<int> cats(k);
  for (int i=0; i<k; i++)
  {
    cats[i] = getCategory(allp[idx[i]-1], simesfactor, alpha, m);
  }
  
  // find the maximum category needed
  int maxcat = std::min(z-m+h+1, k);
  int maxcatI = 0;
  for (int i=k-1; i >= 0; i--)
  {
    if (cats[i] > maxcatI)
    {
      maxcatI = cats[i];
      if (maxcatI >= maxcat) break; 
    }
  }
  maxcat = std::min(maxcat, maxcatI);
  
  // prepare disjoint set data structure
  std::vector<int> parent(maxcat+1);
  std::vector<int> lowest(maxcat+1);
  std::vector<int> rank(maxcat+1, 0);
  for (int i=0; i <= maxcat; i++)
  {
    parent[i] = i;
    lowest[i] = i;
  }

  // The algorithm proper. See pseudocode in paper
  Rcpp::IntegerVector discoveries(k+1,0);
  int lowestInPi;
  for (int i=0; i < k; i++)
  {
    if (cats[i] <= maxcat)
    {
      lowestInPi = lowest[Find(cats[i], parent)];
      if (lowestInPi == 1)
      {
        discoveries[i+1] = discoveries[i]+1;
      }
      else
      {
        discoveries[i+1] = discoveries[i];
        Union(lowestInPi-1, Find(cats[i], parent), parent, lowest, rank);
      }
    }
    else
    {
        discoveries[i+1] = discoveries[i];
    }
  }
  
  return discoveries;
}

