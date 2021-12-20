#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
#include "ARIbrain.h"

// Macros:
// 1) Convert coordinates to index
#define xyz2index(x, y, z, DIMS) ( (z-1)*DIMS[0]*DIMS[1] + (y-1)*DIMS[0] + (x-1) )
// 2) Compute size of 3D image
#define ndims(DIMS) ( DIMS[0]*DIMS[1]*DIMS[2] )


//------------------------- (1) FIND ALL STCS -------------------------//

// Convert voxel index to [x y z] coordinates
std::vector<int> index2xyz(int                  index,
                           Rcpp::IntegerVector& DIMS)
{
    std::vector<int> XYZ;
    XYZ.reserve(3);
    XYZ.push_back( index % DIMS[0] + 1 );
    XYZ.push_back( ((index-(XYZ[0]-1))/DIMS[0]) % DIMS[1] + 1 );
    XYZ.push_back( (index-(XYZ[0]-1)-(XYZ[1]-1)*DIMS[0]) / (DIMS[0]*DIMS[1]) + 1 );
    
    return XYZ;
}

// Check if a voxel is in the mask
bool xyz_check(int                  x,
               int                  y,
               int                  z,
               int                  index,
               Rcpp::IntegerVector& DIMS,
               Rcpp::IntegerVector& MASK)
{
    return (x >= 1 && x <= DIMS[0] && \
            y >= 1 && y <= DIMS[1] && \
            z >= 1 && z <= DIMS[2] && \
            MASK[index] > 0);
}

// Find valid neighbours of a voxel
std::vector<int> findNeighbours(Rcpp::IntegerVector& MASK,   // 3D mask
                                Rcpp::IntegerVector& DIMS,   // image dimensions
                                int                  index,  // voxel index of interest
                                int                  conn)   // connectivity criterion
{
    // compute [x y z] coordinates based on voxel index
    std::vector<int> XYZ = index2xyz(index, DIMS);
    
    // xyz coordinate adjustment vectors
    int DX[26] = {1,-1,0, 0,0, 0,  1,-1, 1,-1,1,-1, 1,-1,0, 0, 0, 0,  1,-1, 1,-1, 1,-1, 1,-1};
    int DY[26] = {0, 0,1,-1,0, 0,  1, 1,-1,-1,0, 0, 0, 0,1,-1, 1,-1,  1, 1,-1,-1, 1, 1,-1,-1};
    int DZ[26] = {0, 0,0, 0,1,-1,  0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1,  1, 1, 1, 1,-1,-1,-1,-1};
    
    // find all valid neighbours of a voxel
    std::vector<int> IDS;
    IDS.reserve(conn);
    for (int i = 0; i < conn; i++)
    {
        int      id = xyz2index(XYZ[0]+DX[i], XYZ[1]+DY[i], XYZ[2]+DZ[i], DIMS);
        bool inmask = xyz_check(XYZ[0]+DX[i], XYZ[1]+DY[i], XYZ[2]+DZ[i], id, DIMS, MASK);
        if (inmask)
        {
            IDS.push_back(MASK[id]-1);
        }
    }
    
    return IDS;
}

// Find the adjacency list for all in-mask voxels
std::vector< std::vector<int> > findAdjList(Rcpp::IntegerVector& MASK,    // 3D mask
                                            Rcpp::IntegerVector& INDEXP,  // indices of sorted p-values
                                            Rcpp::IntegerVector& DIMS,    // image dimensions
                                            int                  m,       // number of in-mask voxels
                                            int                  conn)    // connectivity criterion
{
    std::vector< std::vector<int> > ADJ(m);
    for (int i = 0; i < m; i++)
    {
        std::vector<int> IDS = findNeighbours(MASK, DIMS, INDEXP[i], conn);
        ADJ[i] = IDS;
    }
    
    return ADJ;
}



//// Find function for disjoint set data structure
//// (1) recursive version (path compression) - may result in stack overflow when recursion is too deep
//int Find(int               x,
//         std::vector<int>& PARENT)
//{
//    if (PARENT[x] != x)
//    {
//        PARENT[x] = Find(PARENT[x], PARENT);
//        return PARENT[x];
//    }
//    else
//    {
//        return x;
//    }
//}
//// (2) iterative version (path halving) - reduce the amount of heap memory
////int Find(int               x,
////         std::vector<int>& PARENT)
////{
////    while (PARENT[x] != x)
////    {
////        PARENT[x] = PARENT[PARENT[x]];
////        x         = PARENT[x];
////    }
////
////    return x;
////}
//
//// Union function for disjoint set data structure (union by rank)
//// Extra: we keep track of the lowest entry of each disjoint set
//// That way we can find the lower set to merge with
//void Union(int               x,
//           int               y,
//           std::vector<int>& PARENT,
//           std::vector<int>& LOWEST,
//           std::vector<int>& RANK)
//{
//    int xRoot = Find(x, PARENT);
//    int yRoot = Find(y, PARENT);
//
//    // if x and y are already in the same set (i.e., have the same root or representative)
//    if (xRoot == yRoot) return; // Note: this never happens in our case
//
//    // x and y are not in same set, so we merge
//    if (RANK[xRoot] < RANK[yRoot])
//    {
//        PARENT[xRoot] = yRoot;
//        LOWEST[yRoot] = std::min(LOWEST[xRoot], LOWEST[yRoot]);
//    }
//    else if (RANK[xRoot] > RANK[yRoot])
//    {
//        PARENT[yRoot] = xRoot;
//        LOWEST[xRoot] = std::min(LOWEST[xRoot], LOWEST[yRoot]);
//    }
//    else
//    {
//        PARENT[yRoot] = xRoot;
//        RANK[xRoot]++;
//        LOWEST[xRoot] = std::min(LOWEST[xRoot], LOWEST[yRoot]);
//    }
//}
//
//// Calculate the category for each p-value
//int getCategory(double p,            // p-value for which we need the category
//                double simesfactor,  // simesfactor at h(alpha)
//                double alpha,        // alpha itself
//                int    m)            // size of the problem
//{
//    if ( p == 0 || simesfactor == 0 )
//    {
//        return 1;
//    }
//    else
//    {
//        if (alpha == 0)
//        {
//            return (m+1);
//        }
//        else
//        {
//            double cat = (simesfactor / alpha) * p;
//            return static_cast<int>(std::ceil(cat));
//        }
//    }
//}
//
//// Calculates the size of the concentration set at a fixed alpha
//int findConcentration(Rcpp::NumericVector& P,            // vector of p-values (sorted!)
//                      double               simesfactor,  // simesfactor at h(alpha)
//                      int                  h,            // h(alpha)
//                      double               alpha,        // alpha itself
//                      int                  m)            // size of the problem
//
//{
//  // from m-h we increase z until we fulfil the condition
//  int z = m-h;
//  if (z > 0)  // h=m implies z=0
//  {
//    while (z < m && simesfactor * p[z-1] > (z-m+h+1) * alpha)
//    {
//      z++;
//    }
//  }
//
//  return z;  // concentration set: {1, ..., z}
//}
//
//// Calculates the lower bound to the number of false hypotheses
//// Implements the algorithm based on the disjoint set structure
//std::vector<int> findDiscoveries(Rcpp::IntegerVector& IDS,         // indices in set I (not sorted)
//                                 Rcpp::NumericVector& P,           // all p-values (sorted!)
//                                 double               simesfactor, // simesfactor at h(alpha)
//                                 int                  h,           // h(alpha)
//                                 double               alpha,       // alpha
//                                 int                  k,           // size of I
//                                 int                  m)           // size of the problem
//{
//    // calculate categories for the p-values
//    std::vector<int> CATS;
//    CATS.reserve(k);
//    for (int i = 0; i < k; i++)
//    {
//        CATS.push_back(getCategory(P[IDS[i]-1], simesfactor, alpha, m));
//    }
//
//    // find the maximum category needed
//    int z = findConcentration(P, simesfactor, h, alpha, m);
//    int maxcat  = std::min(z-m+h+1, k);
//    int maxcatI = 0;
//    for (int i = k-1; i >= 0; i--)
//    {
//        if (CATS[i] > maxcatI)
//        {
//            maxcatI = CATS[i];
//            if (maxcatI >= maxcat) break;
//        }
//    }
//    maxcat = std::min(maxcat, maxcatI);
//
//    // prepare disjoint set data structure
//    std::vector<int> PARENT;
//    std::vector<int> LOWEST;
//    std::vector<int> RANK;
//    PARENT.reserve(maxcat+1);
//    LOWEST.reserve(maxcat+1);
//    RANK.reserve(maxcat+1);
//    for (int i = 0; i <= maxcat; i++)
//    {
//        PARENT.push_back(i);
//        LOWEST.push_back(i);
//        RANK.push_back(0);
//    }
//
//    // The algorithm proper. See pseudocode in paper "Hommel's procedure in linear time".
//    Rcpp::IntegerVector DISC(k+1, 0);
//    int lowestInPi;
//    for (int i = 0; i < k; i++)
//    {
//        if (CATS[i] <= maxcat)
//        {
//            lowestInPi = LOWEST[Find(CATS[i], PARENT)];
//            if (lowestInPi == 1)
//            {
//                DISC[i+1] = DISC[i] + 1;
//            }
//            else
//            {
//                DISC[i+1] = DISC[i];
//                Union(lowestInPi-1, Find(CATS[i], PARENT), PARENT, LOWEST, RANK);
//            }
//        }
//        else
//        {
//            DISC[i+1] = DISC[i];
//        }
//    }
//
//    return DISC;
//}





// Union by size
void UnionBySize(int                  x,
                 int                  y,
                 std::vector<int>&    PARENT,
                 Rcpp::IntegerVector& SIZE)
{
    int xRoot = Find(x, PARENT);
    int yRoot = Find(y, PARENT);
    
    // if x and y are already in the same set
    if (xRoot == yRoot) return;
    
    // x and y are not in same set, so we merge
    if (SIZE[xRoot] < SIZE[yRoot])
    {
        PARENT[xRoot] = yRoot;
    }
    else
    {
        PARENT[yRoot] = xRoot;
    }
}

// Compute all supra-threshold clusters (STCs)
// [[Rcpp::export]]
Rcpp::List findClusters(Rcpp::IntegerVector& MASK,    // 3D mask
                        Rcpp::IntegerVector& INDEXP,  // index vector of sorted p-values
                        Rcpp::IntegerVector& DIMS,    // image dimensions
                        Rcpp::IntegerVector& SIZE,    // size of subtrees
                        Rcpp::IntegerVector& ROOT,    // forest roots (>=0 for roots & -1 for not)
                        int                  m,       // number of in-mask voxels
                        int                  conn)    // connectivity criterion
{
    // find neighbours for all voxels
    std::vector< std::vector<int> > ADJ = findAdjList(MASK, INDEXP, DIMS, m, conn);
    
    // prepare disjoint set data structure
    std::vector<int> PARENT, FORESTROOT;
    PARENT.reserve(m);
    FORESTROOT.reserve(m);
    for (int i = 0; i < m; i++)
    {
        PARENT.push_back(i);
        FORESTROOT.push_back(i);
    }
    
    // initialize output: a list of children for all voxels
    std::vector< std::vector<int> > CHILD;
    CHILD.reserve(m);
    // initialize child vector for a single voxel
    std::vector<int> CHD;
    CHD.reserve(conn);
    
    // loop through all voxels
    for (int i = 0; i < m; i++)
    {
        // get neighbour list of i
        std::vector<int> IDS = ADJ[i];
        
        // loop through all neighbours of i
        for (int j = 0; j < IDS.size(); j++)
        {
            if (IDS[j] < i)  // check if IDS[j] has a smaller rank, i.e., a smaller p-value
            {
                int jRoot = Find(IDS[j], PARENT);  // representative of the tree containing IDS[j]
                int     k = FORESTROOT[jRoot];     // forest root of the tree containing jRoot&IDS[j]
                
                if (k != i)
                {
                    UnionBySize(i, jRoot, PARENT, SIZE);
                    int iRoot = Find(jRoot, PARENT);  // new disjoint-set root of the merged set
                    FORESTROOT[iRoot] = i;            // store forest root of the just merged set
                    SIZE[i] += SIZE[k];               // add an edge from i to k to the forest
                    
                    // put a heavy child in front
                    if (CHD.size() == 0 || SIZE[CHD[0]] >= SIZE[k])
                    {
                        CHD.push_back(k);
                    }
                    else
                    {
                        CHD.push_back(CHD[0]);
                        CHD[0] = k;
                    }
                }
            }
        }
        
        // update child list
        CHILD.push_back(CHD);
        CHD.resize(0);
    }
    
    // find forest roots
    int j = 0;
    for (int i = 0; i < m; i++)
    {
        if (PARENT[i] == i)
        {
            ROOT[j] = FORESTROOT[i];
            j++;
        }
    }
    
    return Rcpp::wrap(CHILD);
}


//------------------------- (2) COMPUTE TDPS FOR ALL STCS -------------------------//

// 1) Iterative pre-order traversal to find descendants of v (including v)
//std::vector<int> descendants(int                  v,
//                             Rcpp::IntegerVector& SIZE,   // subtree size for all nodes
//                             Rcpp::List&          CHILD)  // list of children for all voxels
//{
//    std::vector<int> DESC;
//    DESC.reserve(SIZE[v]);
//    DESC.push_back(v);
//
//    int i = 0;
//    while (i < DESC.size())
//    {
//        Rcpp::IntegerVector CHD = Rcpp::as<Rcpp::IntegerVector>(CHILD[DESC[i]]);
//        for (int j = 0; j < CHD.size(); j++)
//        {
//            DESC.push_back(CHD[j]);
//        }
//        i++;
//    }
//
//    return DESC;
//}

// 2) Iterative post-order traversal to find descendants of v (including v).
// Note: When we pop a vertex from the stack, we push that vertex again as a value and
// then all its children in reverse order on the stack. If we pop a value, it means that
// all its children have been fully explored and added to the descendants, so we append
// the current value to the descendants too.
Rcpp::IntegerVector descendants(int                  v,
                                Rcpp::IntegerVector& SIZE,   // subtree size for all nodes
                                Rcpp::List&          CHILD)  // a list of children for all voxels
{
    Rcpp::IntegerVector DESC(SIZE[v], 0);
    
    int len = 0;          // track the number of found descendants
    int top = SIZE[v]-1;  // track the top of the stack
    DESC[top] = v;        // push v on the stack (stack grows from right to left)
    
    while (top < DESC.size())  // while stack is non-empty
    {
        v = DESC[top];       // pop the top of the stack
        top++;
        if (v < 0)           // if ~v is a value, append ~v to the descendants
        {
            DESC[len] = ~v;  // bitwise negation is used as minus doesn't work for zero
            len++;
        }
        else
        {
            top--;
            DESC[top] = ~v;  // push v as a value
            
            // push all children in reverse order
            Rcpp::IntegerVector CHD = Rcpp::as<Rcpp::IntegerVector>(CHILD[v]);  // convert the component of Rcpp::List input to its C++ equivalent with Rcpp::as(), which can be achieved through an implicit call to Rcpp::as(), see below.
            //Rcpp::IntegerVector CHD = CHILD[v];
            int j = CHD.size() - 1;
            while (j >= 0)
            {
                top--;
                DESC[top] = CHD[j];
                j--;
            }
        }
    }
    
    return DESC;
}

// Compute the TDP bounds of the heavy path starting at v
void heavyPathTDP(int                  v,       // start of the heavy path
                  int                  par,     // parent of v (-1 indicates no parent)
                  int                  m,
                  int                  h,       // h(alpha)
                  double               alpha,   // alpha
                  double               simesh,  // simesfactor at h(alpha)
                  Rcpp::NumericVector& P,       // all p-values (sorted!)
                  Rcpp::IntegerVector& SIZE,    // subtree size for all nodes
                  Rcpp::List&          CHILD,   // a list of children for all nodes
                  Rcpp::NumericVector& TDP)     // TDP bounds
{
    Rcpp::IntegerVector HP  = descendants(v, SIZE, CHILD);
    for (int i = 0; i < HP.size(); i++)
    {
        HP[i] = HP[i] + 1;
    }
    Rcpp::IntegerVector NUM = findDiscoveries(HP, P, simesh, h, alpha, HP.size(), m);
    
    while (true)  // walk down the heavy path
    {
        // check if v represents an STC
        if (par == -1 || P[v] != P[par])
        {
            TDP[v] = ((double) NUM[SIZE[v]]) / ((double) SIZE[v]);
        }
        else
        {
            TDP[v] = -1;  // invalid STCs get TDP of -1
        }

        // check if v is a leaf
        Rcpp::IntegerVector CHD = CHILD[v];
        if (CHD.size() == 0) break;
        
        // update v & its parent
        par = v;
        v = CHD[0];
    }
}

// Find the start of every heavy path & compute the TDPs of that heavy path
// start of heavy path: 1) root of F;
//                      2) non-root node that is not the 1st heavy child
// [[Rcpp::export]]
Rcpp::NumericVector forestTDP(int                  m,       // number of vertices
                              int                  h,       // h(alpha)
                              double               alpha,   // alpha
                              double               simesh,  // simesfactor at h(alpha)
                              Rcpp::NumericVector& P,       // all p-values (sorted!)
                              Rcpp::IntegerVector& SIZE,    // subtree size for all nodes
                              Rcpp::IntegerVector& ROOT,    // all roots of the forest
                              Rcpp::List&          CHILD)   // a child list for all nodes
{
    Rcpp::NumericVector TDP(m);
    for (int i = 0; i < m; i++)
    {
        if (ROOT[i] >= 0)  // check if ROOT[i] is a root
        {
            heavyPathTDP(ROOT[i], -1, m, h, alpha, simesh, P, SIZE, CHILD, TDP);
        }

        Rcpp::IntegerVector CHD = CHILD[i];
        int j = 1;
        while (j < CHD.size())  // for each child of i except the 1st
        {
            heavyPathTDP(CHD[j], i, m, h, alpha, simesh, P, SIZE, CHILD, TDP);
            j++;
        }
    }

    return TDP;
}


//------------------------- (3) PREPARE ADMISSIBLE STCS -------------------------//

// Set up L: a list of representative of admissible STCs
// [[Rcpp::export]]
Rcpp::IntegerVector queryPreparation(int                  m,      // number of vertices
                                     Rcpp::IntegerVector& ROOT,   // all roots of the forest
                                     Rcpp::NumericVector& TDP,    // all TDP bounds
                                     Rcpp::List&          CHILD)  // a children list for all vertices
{
    std::vector<int> L;  // an index vector of representatives of admissible STCs
    L.reserve(m);
    std::vector<double> STACK;
    STACK.reserve(m*2);
    
    // loop through all roots
    int i = 0;
    while (ROOT[i] >= 0)  // check if ROOT[i] is a forest root
    {
        STACK.push_back(ROOT[i]);  // walk down the forest from ROOT[i]
        STACK.push_back(-1);       // maximum seen TDP on the path to ROOT[i]: non-existent
        while (STACK.size() > 0)
        {
            double q = STACK.back();         // maximum seen TDP on the path to v
            STACK.pop_back();
            int    v = int(STACK.back());
            STACK.pop_back();
            
            // check if v has higher TDP than its ancestors
            if (TDP[v] > q) L.push_back(v);  // note: q>=-1 & invalid STCs have TDP=-1
            
            Rcpp::IntegerVector CHD = CHILD[v];
            for (int j = 0; j < CHD.size(); j++)
            {
                STACK.push_back(CHD[j]);
                STACK.push_back(std::max(TDP[v], q));
            }
        }
        i++;
    }
    
    // sort L in ascending order of TDP
    std::vector< std::pair<double,int> > PAIRED;  // (tdp, index)
    for (int i = 0; i < L.size(); i++)
    {
        PAIRED.push_back(std::make_pair(TDP[L[i]], L[i]));
    }
    std::sort(PAIRED.begin(), PAIRED.end());
    
    Rcpp::IntegerVector RCPPL(L.size());
    for (int i = 0; i < L.size(); i++)
    {
        RCPPL[i] = PAIRED[i].second;
    }

    return RCPPL;
}


//------------------------- (4) FORM CLUSTERS USING g(amma) -------------------------//

// Find leftmost index i in L such that TDP[L[i]] >= g
// return size(L) if no such index exists;
// run linear search & binary search in parallel
int findLeft(double               g,    // a TDP threshold (gamma)
             Rcpp::IntegerVector& L,    // an index list of all admissible vertices (sorted on TDP)
             Rcpp::NumericVector& TDP)  // all TDP bounds
{
    int right = L.size();
    int   low = 0;
    int  high = right;
    while (low < high)
    {
        int mid = (low+high)/2;  // (1) binary search part (using integer division)
        if (TDP[L[mid]] >= g)
        {
            high = mid;
        }
        else
        {
            low = mid + 1;
        }
        
        right--;                     // (2) linear search part
        // no need to guard against right<0 as right>=0 will always be true
        if (TDP[L[right]] < g) return (right+1);
    }
    
    return low;
}

// Answer the query, i.e., find maximal STCs under the TDP condition
// [[Rcpp::export]]
Rcpp::IntegerVector answerQuery(int                  m,
                                double               g,
                                Rcpp::IntegerVector& L,
                                Rcpp::IntegerVector& SIZE,
                                Rcpp::NumericVector& TDP,
                                Rcpp::List&          CHILD)
{
    if (g < 0) g = 0;  // constrain TDP threshold gamma to be non-negative
    
    int clus = 1;
    int left = findLeft(g, L, TDP);
    
    Rcpp::IntegerVector MARK(m, 0);
    for (int i = left; i < L.size(); i++)
    {
        if (MARK[L[i]] == 0)
        {
            Rcpp::IntegerVector DESC = descendants(L[i], SIZE, CHILD);
            for (int j = 0; j < DESC.size(); j++)
            {
                MARK[DESC[j]] = clus;
            }
            clus++;
        }
    }
    
    return MARK;
}
