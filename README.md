---
output:
  rmarkdown::html_vignette: default
  html_document: default
---
# ARIbrain
All-Resolution Inference (ARI) for brain imaging is an R-package to estimate True Discovery Proportions (TDP) for brain imaging analysis derived clusters of functional MRI activation. The TDP gives the lower-bound (with a certain confidence, usually 95%) of the number of truly active voxels within each cluster.

# Introduction
`ARIbrain` is the R-package for All-Resolution Inference (ARI) in neuroscience. It allows researchers to estimate the True Discovery Proportion (TDP) of any cluster in a statistical map derived from a (functional) MRI experiment. Statistical maps can be derived using your favorite fMRI analysis package (e.g. SPM, FSL, AFNI). It is convenient for the output to be in NIfti format, as this can be read in by the package. Alternatively you could use an R array as input as well (as the nifti files will be converted to an array internally). 

ARIbrain can be used in two different 'modes'. Using `ARI`, we show how to compute lower bound for proportion of active voxels (or any other spatially located units) within given clusters. Alternatively, with `ARIcluster`, we show how to find maximal clusters under the given threshold of true discovery proportion (TDP).

ARIbrain requires R to run and the 'ARIbrain'-package to be installed. For non R-users the easiest way to install R is in combination with Rstudio. You can find the instructions how to install R and RStudio.

## References
The original paper introducing ARI can be found here: (https://doi.org/10.1016/j.neuroimage.2018.07.060), and the efficient algorithm used to perform ARICluster analysis can be found here: (https://doi.org/10.48550/arXiv.2206.13587).

## Installing the 'ARIbrain' package

### Step 1, Installing R and Rstudio.
Go [here](https://posit.co/download/rstudio-desktop/) to download R and Rstudio from the Posit website (you need the RStudio Desktop FREE version). First install R and then install RStudio.

### Step 2, Installing ARIbrain
You can install the stable version of ARIbrain from [CRAN](https://cran.r-project.org/web/packages/ARIbrain/), or use the *Tools > Install packages* option from Rstudio (select CRAN Repository and search for ARIbrain, leave the install dependencies option checked), or use the `install.packages('ARIbrain')` command in R/Rstudio. The development version of ARIbrain can be downloaded from this GitHub repository using the 'devtools' package. First install this package using `install.packages('devtools')`, and then install ARIbrain using the following command: `devtools::install_github('wdweeda/ARIbrain')`.

### Step 3, Running ARI in RStudio
After installing, open RStudio (if not already open), and load the ARIbrain package by typing `library(ARIbrain)`. This will load the package for usage. For easy access to the files it is convenient to change to the working directory of where your statistics maps (in Nifti format) of interest are located by typing `setwd('workdirpath')` where `workdirpath` is the path to your working directory (e.g. `'/Users/wouter/fmri'` or `'c:/Users/wouter/fmridir'`)

# ARI analysis in R
There are two main flavors of TDP estimation using ARI either providing clusters to the analysis and estimating TDPs for these clusters, or the other way around, setting a minimal TDP level and letting ARI estimate the largest clusters wit at least that TDP level.

## ARI using pre-defined clusters
ARI can caluculate TDPs for any cluster provided (with full FWER control). These can be clusters defined by, for example, a cluster-forming threshold, or clusters from an anatomical atlas. The basic input for ARI is a map of (2-sided) p-values and a map of cluster indices (0's for non-clusters and integer values for each voxel that belong to a specific clusters, e.g. 1's for cluster 1, 2's for cluster 2, etc.)

### Main ARI syntax 
If you are familiar with neuroimaging analysis is R here is the main syntax for the R function `ARI` (type `?ARIbrain::ARI` for more details). Below the syntax we will continue with an example of a 'standard' analysis.

`ARI(Pmap, clusters, mask=NULL, alpha=0.05, Statmap=function(ix) -qnorm(Pmap[ix]), summary_stat=c("max", "center-of-mass"), silent=FALSE)`

The main input parameters of `ARI()` are:   

- `Pmap`: the map of p-values, 
- `clusters`: the map of cluster indices.

Others optional maps (parameters) are:   

- `mask`: the map of logicals (not mandatory, but useful),
- `Statmap`: the map of statistics (usually z-scores or t-values).

The function accepts input map formats of character file names or 3D arrays. Therefore the minimal syntax is   
`ARI(Pmap, clusters)`

### Define clusters

The clusters can be defined *a priori*, on the basis of previous knowledges or on the basis of anatomical regions. Clusters of such a kind are usually called ROIs. There are no limitations to the number of ROIs that can be evaluated in the same analysis; the lower bounds for each ROI is valid simultaneously for all estimates (i.e. corrected for multiplicity). 

Even more interestingly, the clusters can be defined on the basis of the same data. This is valid as `ARI` allows for circular analysis, still controlling for multiplicity of inferences.

### Example: 'standard' cluster analysis
If you have an output z-map (i.e., containing z-statistics) of a contrast/analysis of interest and want to to a 'standard' cluster-extent analysis. We first need to load the statistics file into R and threshold the image at a certain Z value (e.g., 3.1) to form clusters. Make sure your input file is the _unthresholded_ map of statistics. Use the following commands to load the file and threshold the map into clusters.

```
zdat <- readNifti('zstat1.nii.gz')
clus31 <- ARIbrain::cluster_threshold(zdat>3.1)
```
The z-statistics data is now loaded into an r-object called `zdat`, the clusters that are formed based on a z-statistics value larger than 3.1 (`zdat>3.1`) which are subsequently stored in the `clus31` object.

Now we need to calculate the p-values (preferably 2-sided) from the z-values. We can do that using the following command:
```
pvals2 <- pnorm(abs(zdat), lower.tail = F)*2
```
Finally, we can estimate the TDPs for the clusters. Best practice is to also give a brain-mask for the in-brain voxels `mask = zdat!=0` (else the method will correct over all voxels including the voxels outside the brain).
```
ari_out <- ARIbrain::ARI(Pmap = pvals2, clusters = clus31, mask=zdat!=0, Statmap = zdat)
```
The `ari_out` output object contains a table with the TDPs (in the ActiveProp column) for all clusters and include locations for the maximum.
```
A hommel object for 199918 hypotheses.
Simes inequality is assumed.
Use p.adjust(), discoveries() or localtest() to access this object.

With 0.95 confidence: at least 25616 discoveries.
1588 hypotheses with adjusted p-values below 0.05.

       Size FalseNull TrueNull ActiveProp dim1 dim2 dim3     Stat
cl92   4470      2668     1802 0.59686801   24   43   58 6.964157
cl91   2644      1187     1457 0.44894100   58   64   63 5.973241
cl90   2226       991     1235 0.44519317   67   33   15 6.258505
cl89   1903       829     1074 0.43562796   29   35   21 6.346173
cl88   1019       227      792 0.22276742   23   80   53 5.484357
cl87    800       324      476 0.40500000   32   68   65 6.333398
cl86    407        59      348 0.14496314   30   76   38 5.569287
cl85    337        43      294 0.12759644   62   75   39 5.331001
cl84    297         0      297 0.00000000   63   92   46 4.748392
cl83    146         0      146 0.00000000   34   22   27 4.258435
cl82    128         0      128 0.00000000   51   18   37 4.881685
cl81     81         0       81 0.00000000   16   50   27 4.635155
cl80     38         0       38 0.00000000   46   53   42 4.503975
cl0  185422     11120   174302 0.05997131   48   19   30 4.937728
```
The output also includes `cl0` which is a 'null' cluster containing all in-mask voxels that are not in any other cluster. If this gives a TDP > 0 this means that there is some activitiy somewhere in the brain not captured in the clusters. Usually this number is relatively low (<1%) If this number is really high, you might want to redefine clusters, as you might have missed some area of activation.

The following paragraphs contain additional ways to define clusters.

## FSL command-line add-on
You can also append your FSL cluster analysis with TDPs. For this go the cope directory of interest in your multilevel gfeat analysis in the terminal using cd 'gfeatdir/copedir/', where 'gfeatdir/copedir/' is the cope directory of your multilevel analysis. Also download the [get_tdp.R](https://github.com/wdweeda/ohbm2023_edu_course/blob/main/practicals/aribrain/get_tdp.R) file (for example in the 'download' directory. The command line has the following input

```
True Discovery Proportions (TDP) using ARIbrain version 1.2.a

Usage:
Rscript get_tdp.R --zstat=<filename> --cluster=<filename> [options]

Compulsory arguments:
--zstat/tstat   filename of t/z-statistics file (nifti)
--cluster       filename of cluster-index file (nifti)
--df            if tstat is specified, df is also needed

Optional arguments:
[method]
--alpha         nominal alpha level of TDPs, default is 0.05 (%95 confidence)
[output]
--outfile       optional output file of TDP values (nifti)
--outtable      optional output table of TDP values (txt)
--intable       optional input table (txt) of clusters (TDPs will be added with _tdp)
--inhtml        optional input html cluster file (TDPs will be added with _tdp)
[options]
--verbose       print additional information on TDPs
--quiet         suppresses almost all output to console
--help          prints this message
```

The following command will append your cluster.html file (with extension _tdp.html) and save a text file with TDPs and a nifti file with voxels having the TDP value of the cluster they belong to. In the cope directory run: `Rscript /download/get_tdp.R --zstat=./stats/zstat1.nii.gz --cluster=cluster_mask_zstat1.nii.gz --alpha=0.05 --outtable=tdptable.txt --outfile=tdpclus.nii.gz --inhtml=cluster_zstat1_std.html`

Ouput will look approximately like this:
```
True Discovery Proportions (TDP) using ARIbrain version 1.2.a
Calculated assuming Simes' inequality with 95% confidence.

 > converted z-stats to 2-sided p-values
 > writing TDP table
 > writing cluster nifti with TDP values
 > adding TDPs to cluster html

<TDP output>
      Size FalseNull TrueNull ActiveProp dim1 dim2 dim3      Stat
cl5   9161      7134     2027 0.77873595   15   60   37 10.095283
cl4   8851      6531     2320 0.73788273   77   58   39 10.245214
cl3    555        40      515 0.07207207   36   65   41  5.100245
cl2    380        67      313 0.17631579   69   70   50  5.211789
cl1    195        53      142 0.27179487   72   61   62  5.915697
cl0 146612      7048   139564 0.04807246   44   87   29  3.539190

```

# ARICluster analysis

Here we show an analysis where clusters are defined by a TDP threshold. Using a sufficiently high TDP threshold leads to achieving better spatial localisation. In contrast to classical cluster inference by providing a fixed cluster-forming threshold (CFT), `ARICluster` uses flexible CFTs, each defind by the TDP threshold, and ensures all derived clusters are maximal and obtain the TDP meeting or exceeding the pre-specified TDP threshold.

## Syntax and parameters
The syntax of the function includes two steps:

1. Create an ARIBrainCluster-class object (type `?ARIbrain::ARIBrainCluster` for more details).

    `ARIBrainCluster(Pmap, mask, conn=18, alpha=0.05)`

    The main input argument of `ARIBrainCluster()` is:   

    - `Pmap`: the map of p-values.

    Other optional arguments are:   

    - `mask`: the map of numerics/logicals (not mandatory, but useful),
    - `conn`: the connectivity criterion: face (8), edge (18) and vertex (26),
    - `alpha`: the significance level.
    
    The output of `ARIBrainCluster()` is:
    
    - `ARIBrainCluster-class`: the ARIBrainCluster-class object.

2. Answer queries given a TDP threshold (type `?ARIbrain::TDPQuery` for more details).

    `TDPQuery(ARIBrainCluster-class, threshold)`

    The input arguments of `TDPQuery()` are:   

    - `ARIBrainCluster-class`: the ARIBrainCluster-class object,
    - `threshold`: the TDP threshold.
    
    The output of `TDPQuery()` is:
    
    - `TDPBrainClusters`: the TDPBrainClusters object.
    
3. Change the size of a chosen cluster by increasing or decreasing its TDP (type `?ARIbrain::TDPChange` for more details).

    `TDPChange(TDPBrainClusters, v, tdpchg)`
    
    The main input arguments of `TDPChange()` are:  

    - `TDPBrainClusters`: the TDPBrainClusters object.
    - `v`: index or XYZ coordinates of a node within the chosen cluster.
    
    Another optional argument is:
    
    - `tdpchg`: a sufficient TDP change.
    
    The output of `TDPChange()` is:
    
    - `TDPBrainClusters`: the TDPBrainClusters object.
    
The function accepts input map formats of character file names or 3D arrays. Therefore the minimal syntax is
```{r, eval=FALSE}
aricluster <- ARIBrainCluster(Pmap)
tdpclusters <- TDPQuery(aricluster, threshold)
tdpchanges <- TDPChange(tdpclusters, v=c(X,Y,Z))
```

Others methods can be used to show the resulting cluster information:
    
  - `summary(TDPBrainClusters, rest=FALSE)`  # summarize cluster information & show summary table 
  - `print(TDPBrainClusters)`                # print summary table
  - `show(TDPBrainClusters)`                 # display summary table
  - `length(TDPBrainClusters)`               # compute the number of clusters
  - `TDPBrainClusters[[i]]`                  # obtain 3D voxel indices of the ith largest cluster
  - `TDPBrainClusters[1:i]`                  # obtain 3D voxel indices of the top i largest cluster
    
and to write out the cluster image
    
  - `writeClusters(TDPBrainClusters, file, template)`

## Select proper TDP threshold

Finding clusters with non-zero TDP threshold `threshold` indicates the presence of some signal in each cluster, however, `threshold` of 40%, 70% and 90% could be characterized as weak, moderate and strong spatial localisation, respectively.

# Additional syntax
## ARI examples
### Create `cluster.nii.gz` with FSL

You simply need to run on the shell:

`cluster -z zstat1.nii.gz -t 3.2 -o cluster.nii.gz`

This will create `cluster.nii.gz` that you need.

*hint*: In case it retun an error message like  
`cluster: error while loading shared libraries: libutils.so: cannot open shared object file: No such file or directory`
type into the shell (replacing the path with your own path of the file fsl.sh):

`source /etc/fsl/5.0/fsl.sh`

and try again.

Get a complete help for FSL at  
<https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster>

### Compute thresholds and clusters on the basis of concentration set (optimal threshold)

```{r}
library(RNifti)
library(hommel)
library(ARIbrain)

Tmap <- RNifti::readNifti(system.file("extdata", "zstat.nii.gz", package="ARIbrain"))
Pmap <- RNifti::readNifti(system.file("extdata", "pvalue.nii.gz", package="ARIbrain"))
mask <- RNifti::readNifti(system.file("extdata", "mask.nii.gz", package="ARIbrain"))
mask <- mask!=0

# compute p-value threshold (thr_p) and z-score threshold (thr_z)
hom <- hommel::hommel(Pmap[mask])
thr_p <- hommel::concentration(hom)
thr_z <- -qnorm(thr_p)

# define clusters
Tmap[!mask] <- 0
clstr <- cluster_threshold(Tmap>thr_z)
table(clstr)
```


### Nifti (nii) inputs

```{r}
library(ARIbrain)

pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARIbrain")
cluster_name <- system.file("extdata", "cluster_th_3.2.nii.gz", package="ARIbrain")
zstat_name <- system.file("extdata", "zstat.nii.gz", package="ARIbrain")
mask_name <- system.file("extdata", "mask.nii.gz", package="ARIbrain")

res_ARI <- ARI(Pmap=pvalue_name, clusters=cluster_name, mask=mask_name, Statmap=zstat_name)

str(res_ARI)
```

### Array inputs

```{r}
library(RNifti)
library(ARIbrain)

Tmap <- RNifti::readNifti(system.file("extdata", "zstat.nii.gz", package="ARIbrain"))
# compute p-values from test statistic (refering to Normal distribution, right-side alternative)
Pmap <- pnorm(-Tmap)

# read the mask file
mask <- RNifti::readNifti(system.file("extdata", "mask.nii.gz", package="ARIbrain"))
# make sure that mask is a logical map
mask <- mask!=0

# create clusters using a threshold equal to 3.2
Tmap[!mask] <- 0
clstr <- cluster_threshold(Tmap>3.2)
table(clstr)

res_ARI <- ARI(Pmap, clusters=clstr, mask=mask, Statmap=Tmap)

str(res_ARI)
```

## ARICluster examples

```{r}
library(ARIbrain)

Pmap <- RNifti::readNifti(system.file("extdata", "pvalue.nii.gz", package="ARIbrain"))
mask <- RNifti::readNifti(system.file("extdata", "mask.nii.gz", package="ARIbrain"))
mask <- mask!=0

# (1) create an ARIBrainCluster-class object
aricluster <- ARIBrainCluster(Pmap = Pmap, mask = mask)

# (2) answer query: find all maximal clusters given a TDP threshold
tdpclusters <- TDPQuery(aricluster, threshold = 0.7)
summary(tdpclusters)

# (3) change query: change the size of chosen cluster based on request
tdpchanges <- TDPChange(tdpclusters, v = c(17,57,38), tdpchg =  0.01)  # increase TDP bound / decrease cluster size
summary(tdpchanges)
tdpchanges <- TDPChange(tdpchanges,  v = c(17,57,38), tdpchg =  0.01)  # further update the chosen cluster
summary(tdpchanges)

tdpchanges <- TDPChange(tdpchanges,  v = c(17,57,38), tdpchg = -0.01)  # decrease TDP bound / increase cluster size
summary(tdpchanges)
```
