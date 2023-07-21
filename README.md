# ARIbrain
All-Resolution Inference (ARI) for brain imaging

# Introduction
`ARIbrain` is the R-package for All-Resolution Inference (ARI) in neuroscience. It allows researchers to estimate the True Discovery Proportion (TDP) of any cluster in a statistical map derived from a (functional) MRI experiment. Statistical maps can be derived using your favorite fMRI analysis package (e.g. SPM, FSL, AFNI). It is convenient for the output to be in NIfti format, as this can be read in by the package. Alternatively you could use an R array as input as well (as the nifti files will be converted to an array internally). 

ARIbrain can be used in two different 'modes'. Using `ARI`, we show how to compute lower bound for proportion of active voxels (or any other spatially located units) within given clusters. Alternatively, with `ARIcluster`, we show how to find maximal clusters under the given threshold of true discovery proportion (TDP).

ARIbrain requires R to run and the 'ARIbrain'-package to be installed. For non R-users the easiest way to install R is in combination with Rstudio. You can find the instructions how to install R and RStudio.

## Installing the 'ARIbrain' package

### Step 1, Installing R and Rstudio.
Go [here](https://posit.co/download/rstudio-desktop/) to download R and Rstudio from the Posit website (you need the RStudio Desktop FREE version). First install R and then install RStudio.

### Step 2, Installing ARIbrain
You can install the stable version of ARIbrain from [CRAN](https://cran.r-project.org/web/packages/ARIbrain/), or use the *Tools > Install packages* option from Rstudio (select CRAN Repository and search for ARIbrain, leave the install dependencies option checked), or use the `install.packages('ARIbrain')` command in R/Rstudio. The development version of ARIbrain can be downloaded from this GitHub repository using the 'devtools' package. First install this package using `install.packages('devtools')`, and then install ARIbrain using the following command: `devtools::install_github('wdweeda/ARIbrain')`.

# ARI analysis

First, we show an analysis where clusters are defined by a supra-threshold-statistic rule. This is the typical case of cluster-wise analysis followed by a multiplicity correction based on Random Field Theory. Here we follow an alternative way: we provide lower bound for proportion for the estimate of active voxels.

## Syntax and parameters
The syntax of the function is (type `?ARIbrain::ARI` for more details)

`ARI(Pmap, clusters, mask=NULL, alpha=0.05, Statmap=function(ix) -qnorm(Pmap[ix]), summary_stat=c("max", "center-of-mass"), silent=FALSE)`

The main input parameters of `ARI()` are:   

- `Pmap`: the map of p-values, 
- `clusters`: the map of cluster indices.

Others optional maps (parameters) are:   

- `mask`: the map of logicals (not mandatory, but useful),
- `Statmap`: the map of statistics (usually z-scores or t-values).

The function accepts input map formats of character file names or 3D arrays. Therefore the minimal sintax is   
`ARI(Pmap, clusters)`

## Define clusters

The clusters can be defined *a priori*, on the basis of previous knowledges or on the basis of anatomical regions. Clusters of such a kind are usually called ROIs. There are no limitations to the number of ROIs that can be evaluated in the same analysis; the lower bounds for each ROI is valid simultaneously for all estimates (i.e. corrected for multiplicity). 

Even more interestingly, the clusters can be defined on the basis of the same data. This is true, since `ARI` allows for circular analysis, still controlling for multiplicity of inferences.

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

## ARI examples

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

# ARICluster analysis

Here we show an analysis where clusters are defined by a TDP threshold. Using a sufficiently high TDP threshold leads to achieving better spatial localisation. In contrast to classical cluster inference by providing a fixed cluster-forming threshold (CFT), `ARICluster` uses flexible CFTs, each defind by the TDP threshold, and ensures all derived clusters obtain the TDP meeting or exceeding the pre-specified TDP threshold.

## Syntax and parameters
The syntax of the function includes two steps:

1. Create an ARIBrainCluster object (type `?ARIbrain::ARIBrainCluster` for more details).

    `ARIBrainCluster(Pmap, mask, conn=18, alpha=0.05)`

    The main input parameter of `ARIBrainCluster()` is:   

    - `Pmap`: the map of p-values.

    Others optional maps (parameters) are:   

    - `mask`: the map of numerics/logicals (not mandatory, but useful),
    - `conn`: the connectivity criterion: face (8), edge (18) and vertex (26),
    - `alpha`: the significance level.

2. Answer queries given a TDP threshold (type `?ARIbrain::TDPquery` for more details).

    `TDPquery(ARIBrainCluster, gamma)`

    The input parameters of `TDPquery()` are:   

    - `ARIBrainCluster`: the ARIBrainCluster object,
    - `gamma`: the TDP threshold.

    Others methodss can be used to summarize the resulting cluster information:
    
    - `summaryCluster(ARIBrainCluster, TDPquery, rest=FALSE)`
    - `writeCluster(ARIBrainCluster, TDPquery, file="aribrain.nii.gz", template=NULL)`

The function accepts input map formats of character file names or 3D arrays. Therefore the minimal sintax is
```{r, eval=FALSE}
ari <- ARIBrainCluster(Pmap)
TDPquery(ari, gamma)
```

## Select proper TDP threshold

Finding clusters with non-zero TDP threshold `gamma` indicates the presence of some signal in each cluster, however, `gamma` of 40%, 70% and 90% could be characterised as weak, moderate and strong spatial localisation, respectively.

## ARIcluster examples

### Nifti (nii) inputs

```{r}
library(ARIbrain)

pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARIbrain")
mask_name <- system.file("extdata", "mask.nii.gz", package="ARIbrain")

# (1) create an ARIBrainCluster object
ari <- ARIBrainCluster(Pmap=pvalue_name, mask=mask_name)

# (2) answer queries: find all maximal clusters given a TDP threshold
res <- TDPquery(ari, gamma=0.7)
# res@clusterlist   # access cluster list
# print cluster summary table
summaryCluster(ari, res, rest=TRUE)
# write cluster image
writeCluster(ari, res, file="ari.nii.gz", template=mask_name)
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

# (1) create an ARIBrainCluster object
ari <- ARIBrainCluster(Pmap, mask=mask)

# (2) answer queries: find all maximal clusters given a TDP threshold
res <- TDPquery(ari, gamma=0.7)
# res@clusterlist   # access cluster list
# print cluster summary table
summaryCluster(ari, res, rest=TRUE)
# write cluster image
writeCluster(ari, res, file="ari.nii.gz", template=mask)
```
