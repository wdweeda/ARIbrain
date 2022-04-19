#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
#include "ARIBrain.h"

//------------------------- (0) PREPARATIONS FOR 3D INPUTS -------------------------//

// Macros:
// 1) Convert xyz coordinates to index
#define xyz2index(x, y, z, DIMS) ( (z-1)*DIMS[0]*DIMS[1] + (y-1)*DIMS[0] + (x-1) )
// 2) Compute size of 3D image
#define ndims(DIMS) ( DIMS[0]*DIMS[1]*DIMS[2] )

// Convert voxel index to [x y z] coordinates
std::vector<int> index2xyz(int                  index,
                           Rcpp::IntegerVector& DIMS)
{
    std::vector<int> XYZ;
    XYZ.reserve(3);
    XYZ.push_back( index % DIMS[0] + 1 );
    //XYZ.push_back( ((index-(XYZ[0]-1))/DIMS[0]) % DIMS[1] + 1 );
    XYZ.push_back( (index/DIMS[0]) % DIMS[1] + 1 );
    //XYZ.push_back( (index-(XYZ[0]-1)-(XYZ[1]-1)*DIMS[0]) / (DIMS[0]*DIMS[1]) + 1 );
    XYZ.push_back( index/(DIMS[0]*DIMS[1]) + 1 );
    
    return XYZ;
}

// Convert several voxel indices to an xyz-coordinate matrix
// [[Rcpp::export]]
Rcpp::IntegerMatrix ids2xyz(Rcpp::IntegerVector& IDS,
                            Rcpp::IntegerVector& DIMS)
{
    Rcpp::IntegerMatrix XYZS(IDS.size(), 3);
    for (int i = 0; i < IDS.size(); i++)
    {
        std::vector<int> XYZ = index2xyz(IDS[i], DIMS);
        XYZS(i,0) = XYZ[0];
        XYZS(i,1) = XYZ[1];
        XYZS(i,2) = XYZ[2];
    }
    
    return XYZS;
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
            MASK[index] != 0);
}

// Find valid neighbours of a voxel
std::vector<int> findNeighbours(Rcpp::IntegerVector& MASK,   // 3D mask of original orders (1:m)
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
// [[Rcpp::export]]
Rcpp::List findAdjList(Rcpp::IntegerVector& MASK,    // 3D mask of original orders (1:m)
                       Rcpp::IntegerVector& INDEXP,  // voxel indices of unsorted p-values
                       Rcpp::IntegerVector& DIMS,    // image dimensions
                       int                  m,       // number of in-mask voxels
                       int                  conn)    // connectivity criterion
{
    Rcpp::List ADJ(m);
    for (int i = 0; i < m; i++)
    {
        std::vector<int> IDS = findNeighbours(MASK, DIMS, INDEXP[i], conn);
        ADJ[i] = Rcpp::IntegerVector(IDS.begin(), IDS.end());
    }
    
    return ADJ;
}


//------------------------- (1) FIND ALL STCS (USING SORTING RANKS) -------------------------//

// Union function of a disjoint-set data structure based on the "union by size" technique.
// The disjoint sets will represent the components of a forest, and the data structure is
// augmented to keep track of the forest root of each component. UnionBySize(i,j) merges
// the sets S_i and S_j (where S_k denotes the set containing k), and assigns the forestroot
// of S_i to be the forestroot of the resulting union.
void UnionBySize(int                  i,
                 int                  j,
                 std::vector<int>&    PARENT,
                 std::vector<int>&    FORESTROOT,
                 Rcpp::IntegerVector& SIZE)
{
    int irep = Find(i, PARENT);
    int jrep = Find(j, PARENT);
    
    // if i and j are already in the same set
    if (irep == jrep) return;
    
    // i and j are not in same set, so we merge
    int iroot = FORESTROOT[irep];
    int jroot = FORESTROOT[jrep];
    if (SIZE[iroot] < SIZE[jroot])
    {
        PARENT[irep] = jrep;
        FORESTROOT[jrep] = iroot;
    }
    else
    {
        PARENT[jrep] = irep;
    }
    SIZE[iroot] += SIZE[jroot];
}

// Compute all supra-threshold clusters (STCs)
// [[Rcpp::export]]
Rcpp::List findClusters(int                  m,     // number of nodes
                        Rcpp::List&          ADJ,   // a list of neighbours for all nodes (unsorted!)
                        Rcpp::IntegerVector& ORD,   // sorted orders for non-decreasing p-values
                        Rcpp::IntegerVector& RANK)  // sorting ranks for all p-values
{
    // initialize output (1): a list of children for all nodes
    Rcpp::List CHILD(m);
    // initialize output (2): a vector of sizes of subtrees
    Rcpp::IntegerVector SIZE(m, 1);
    // initialize output (3): a list of forest roots
    std::list<int> ROOT;
    
    // initialize a child list for a single node
    std::list<int> CHD;
    // // using std::forward_list instead of std::list (Note: must include header <forward_list>)
    // std::forward_list<int> CHD;
    
    // prepare disjoint set data structure
    std::vector<int> PARENT, FORESTROOT;
    PARENT.reserve(m);
    FORESTROOT.reserve(m);
    for (int i = 0; i < m; i++)
    {
        PARENT.push_back(i);
        FORESTROOT.push_back(i);
    }
    
    // loop through all nodes in the ascending order of p-values
    for (int i = 0; i < m; i++)
    {
        // find neighbours for node with the ith smallest p-value
        Rcpp::IntegerVector IDS = ADJ[ORD[i]-1];
        
        // loop through all its neighbours
        for (int j = 0; j < IDS.size(); j++)
        {
            if (RANK[IDS[j]] < i+1)  // check if the neighbour has a smaller rank
            {
                int jrep = Find(RANK[IDS[j]]-1, PARENT);  // representative of the tree
                int    k = FORESTROOT[jrep];              // forest root of the tree
                
                if (k != i)
                {
                    // Merge S_i and S_{jrep}
                    UnionBySize(i, jrep, PARENT, FORESTROOT, SIZE);
                    
                    // put a heavy child in front (using std::list)
                    if (CHD.empty() || SIZE[CHD.front()] >= SIZE[k])
                    {
                        CHD.push_back(k);
                        // CHD.insert_after(CHD.begin(), k);  // for std::forward_list
                    }
                    else
                    {
                        CHD.push_front(k);
                        // CHD.push_back(CHD[0]);  // for std::vector
                        // CHD[0] = k;
                    }
                }
            }
        }
        
        // update child list
        CHILD[i] = Rcpp::IntegerVector(CHD.begin(), CHD.end());
        CHD.resize(0);
    }
    
    // find forest roots
    for (int i = 0; i < m; i++)
    {
        if (PARENT[i] == i)
        {
            ROOT.push_back(FORESTROOT[i]);
        }
    }

    return Rcpp::List::create(Rcpp::Named("CHILD") = CHILD,
                              Rcpp::Named("ROOT") = Rcpp::IntegerVector(ROOT.begin(), ROOT.end()),
                              Rcpp::Named("SIZE") = SIZE);
}


//------------------------- (2) COMPUTE TDPS FOR ALL STCS -------------------------//

// Iterative post-order traversal to find descendants of v (including v).
// Note: When we pop a vertex from the stack, we push that vertex again as a value and
// then all its children in reverse order on the stack. If we pop a value, it means that
// all its children have been fully explored and added to the descendants, so we append
// the current value to the descendants too.
Rcpp::IntegerVector descendants(int                  v,      // sorting rank (0:m-1)
                                Rcpp::IntegerVector& SIZE,   // subtree sizes for all nodes
                                Rcpp::List&          CHILD)  // a list of children for all nodes
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
            Rcpp::IntegerVector CHD = Rcpp::as<Rcpp::IntegerVector>(CHILD[v]);  // convert the component of Rcpp::List input to its C++ equivalent with Rcpp::as(), which can be achieved through an implicit call to Rcpp::as() using the code below.
            // Rcpp::IntegerVector CHD = CHILD[v];
            for (int j = CHD.size() - 1; j >= 0; j--)
            {
                top--;
                DESC[top] = CHD[j];
            }
        }
    }
    
    return DESC;
}

// Compute the TDP bounds of the heavy path starting at v
void heavyPathTDP(int                  v,       // start of the heavy path (0:m-1)
                  int                  par,     // parent of v (-1 indicates no parent)
                  int                  m,       // number of all nodes
                  int                  h,       // h(alpha)
                  double               alpha,   // alpha
                  double               simesh,  // simesfactor at h(alpha)
                  Rcpp::NumericVector& P,       // all p-values (unsorted!)
                  Rcpp::IntegerVector& ORD,     // sorted order for non-decreasing p-values
                  Rcpp::IntegerVector& SIZE,    // subtree sizes for all nodes
                  Rcpp::List&          CHILD,   // a list of children for all nodes
                  Rcpp::NumericVector& TDP)     // TDP bounds
{
    Rcpp::IntegerVector HP = descendants(v, SIZE, CHILD);
    for (int i = 0; i < HP.size(); i++)
    {
        HP[i] = ORD[HP[i]];  // convert sorting rank to original order
    }
    Rcpp::IntegerVector NUM = findDiscoveries(HP, P, simesh, h, alpha, HP.size(), m);
    
    while (true)  // walk down the heavy path
    {
        // check if v represents an STC
        if (par == -1 || P[ORD[v]-1] != P[ORD[par]-1])
        {
            TDP[v] = ((double) NUM[SIZE[v]]) / ((double) SIZE[v]);
        }
        else
        {
            TDP[v] = -1;  // invalid STCs get TDP of -1
        }

        // check if v is a leaf
        if (SIZE[v] == 1) break;
        
        // update v & its parent
        par = v;
        Rcpp::IntegerVector CHD = CHILD[v];
        v = CHD[0];
    }
}

// Find the start of every heavy path & compute the TDPs of that heavy path
// start of heavy path: 1) root of F;
//                      2) non-root node that is not the 1st heavy child
// [[Rcpp::export]]
Rcpp::NumericVector forestTDP(int                  m,       // number of all nodes
                              int                  h,       // h(alpha)
                              double               alpha,   // alpha
                              double               simesh,  // simesfactor at h(alpha)
                              Rcpp::NumericVector& P,       // all p-values (unsorted!)
                              Rcpp::IntegerVector& ORD,     // sorted order for non-decreasing p-values
                              Rcpp::IntegerVector& SIZE,    // subtree size for all nodes
                              Rcpp::IntegerVector& ROOT,    // all roots of the forest
                              Rcpp::List&          CHILD)   // a child list for all nodes
{
    Rcpp::NumericVector TDP(m);
    // loop through all roots
    for (int i = 0; i < ROOT.size(); i++)
    {
        heavyPathTDP(ROOT[i], -1, m, h, alpha, simesh, P, ORD, SIZE, CHILD, TDP);
    }
    // loop through all nodes
    for (int i = 0; i < m; i++)
    {
        Rcpp::IntegerVector CHD = CHILD[i];
        for (int j = 1; j < CHD.size(); j++)
        {
            heavyPathTDP(CHD[j], i, m, h, alpha, simesh, P, ORD, SIZE, CHILD, TDP);
        }
    }

    return TDP;
}


//------------------------- (3) PREPARE ADMISSIBLE STCS -------------------------//

// construct a comparator for the below sorting step
struct compareBy
{
    Rcpp::NumericVector& value;
    compareBy(Rcpp::NumericVector& val) : value(val) {}
    bool operator() (int i, int j) {return value[i] < value[j];}
};

// Set up ADMSTC: a list of representative of admissible STCs
// [[Rcpp::export]]
Rcpp::IntegerVector queryPreparation(int                  m,      // number of vertices
                                     Rcpp::IntegerVector& ROOT,   // all roots of the forest
                                     Rcpp::NumericVector& TDP,    // all TDP bounds
                                     Rcpp::List&          CHILD)  // a children list for all vertices
{
    std::vector<int> ADMSTC;  // a vector of representatives of admissible STCs
    ADMSTC.reserve(m);
    std::vector<double> STACK;
    STACK.reserve(m*2);
    
    // loop through all roots
    for (int i = 0; i < ROOT.size(); i++)
    {
        STACK.push_back(ROOT[i]);  // walk down the forest from ROOT[i]
        STACK.push_back(-1);       // maximum seen TDP on the path to ROOT[i]: non-existent
        while (STACK.size() > 0)
        {
            double q = STACK.back();  // maximum seen TDP on the path to v
            STACK.pop_back();
            int    v = int(STACK.back());
            STACK.pop_back();
            
            // check if v has higher TDP than its ancestors
            if (TDP[v] > q) ADMSTC.push_back(v);  // note: q>=-1 & invalid STCs have TDP=-1
            
            Rcpp::IntegerVector CHD = CHILD[v];
            for (int j = 0; j < CHD.size(); j++)
            {
                STACK.push_back(CHD[j]);
                STACK.push_back(std::max(TDP[v], q));
            }
        }
    }
    
    // // use a lambda function for the sorting step (To compile, require C++11 and use the command: g++ –std=c++11 ARIBrain.cpp; to construct R package, add the following to the DESCRIPTIONS file to enable C++11: SystemRequirements: C++11)
    // std::sort(ADMSTC.begin(), ADMSTC.end(), [&TDP](int v, int w) {return TDP[v] < TDP[w];});
    
    // sort ADMSTC in ascending order of TDP using the comparator
    std::sort(ADMSTC.begin(), ADMSTC.end(), compareBy(TDP));

    return Rcpp::IntegerVector(ADMSTC.begin(), ADMSTC.end());
}


//-------------------------- (4) FORM CLUSTERS USING gamma --------------------------//

// Find leftmost index i in ADMSTC such that TDP[ADMSTC[i]] >= g
// return size(ADMSTC) if no such index exists;
// run linear search & binary search in parallel;
// gamma>=0 is needed because inadmissible STCs have been assigned TDP -1.
int findLeft(double               gamma,    // a TDP threshold (gamma)
             Rcpp::IntegerVector& ADMSTC,   // a list of all admissible vertices (sorted on TDP)
             Rcpp::NumericVector& TDP)      // all TDP bounds
{
    int right = ADMSTC.size();
    int   low = 0;
    int  high = right;
    while (low < high)
    {
        int mid = (low+high)/2;  // (1) binary search part (using integer division)
        if (TDP[ADMSTC[mid]] >= gamma)
        {
            high = mid;
        }
        else
        {
            low = mid + 1;
        }
        
        right--;                 // (2) linear search part
        // no need to guard against right<0 as right>=0 will always be true
        if (TDP[ADMSTC[right]] < gamma) return (right+1);
    }
    
    return low;
}

// Answer the query, i.e., find maximal STCs under the TDP condition.
// gamma>=0 is needed because inadmissible STCs have been assigned TDP -1.
// [[Rcpp::export]]
Rcpp::List answerQuery(double               gamma,
                       Rcpp::IntegerVector& ADMSTC,
                       Rcpp::IntegerVector& ORD,
                       Rcpp::IntegerVector& SIZE,
                       Rcpp::IntegerVector& MARK,
                       Rcpp::NumericVector& TDP,
                       Rcpp::List&          CHILD)
{
    if (gamma < 0) gamma = 0;  // constrain TDP threshold gamma to be non-negative
    
    // initialise output: a list of sorting rank vectors for all clusters
    std::list<Rcpp::IntegerVector> ANS;
    
    int left = findLeft(gamma, ADMSTC, TDP);
    
    for (int i = left; i < ADMSTC.size(); i++)
    {
        if (MARK[ORD[ADMSTC[i]]-1] == 0)
        {
            // append a cluster to ANS
            Rcpp::IntegerVector DESC = descendants(ADMSTC[i], SIZE, CHILD);
            ANS.push_back(DESC);
            // mark the corresponding voxels
            for (int j = 0; j < DESC.size(); j++)
            {
                MARK[ORD[DESC[j]]-1] = 1;
            }
        }
    }
    
    // clear marks back to 0
    for(std::list<Rcpp::IntegerVector>::iterator it = ANS.begin(); it != ANS.end(); ++it)
    {
        for(int j = 0; j < (*it).size(); j++)
        {
            MARK[ORD[(*it)[j]]-1] = 0;
        }
    }
    
    return Rcpp::wrap(ANS);
}

// Counting sort in descending order of cluster sizes.
// [[Rcpp::export]]
Rcpp::IntegerVector counting_sort(int                  n,          // #{clusters}
                                  int                  maxid,      // max(cluster size)
                                  Rcpp::IntegerVector& CLSTRSIZE)  // unsorted cluster sizes
{
    // initialise output sorted indices for descending cluster sizes
    Rcpp::IntegerVector SORTED(n, 0);
    std::vector<int>    COUNT(maxid+1, 0);
    
    // store count of each cluster size
    for (int i = 0; i < n; i++)
    {
        COUNT[CLSTRSIZE[i]]++;
    }
    
    // find cumulative frequency
    for (int i = maxid; i > 0; i--)
    {
        COUNT[i-1] += COUNT[i];
    }
    
    for (int i = 0; i < n; i++)
    {
        SORTED[COUNT[CLSTRSIZE[i]] - 1] = i;
        COUNT[CLSTRSIZE[i]]--;
    }
    
    return SORTED;
}
