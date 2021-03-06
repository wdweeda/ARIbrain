# ARIbrain
All-Resolution Inference (ARI) for brain imaging

# Introduction
`ARIbrain` is the package for All-Resolution Inference in neuroscience.

Here we show how to compute lower bound for proportion of active voxels (or any other spatially located units) within given clusters.

The clusters can be defined *a priori*, on the basis of previous knowledges or on the basis of anatomical regions. Clusters of such a kind are usually called ROIs. There are no limitations to the number of ROIs that can be evaluated in the same analysis; the lower bounds for each ROI is valid simultaneously for all estimates (i.e. corrected for multiplicity). 

Even more interestigly, the clusters can be defined on the basis of the same data. This is true, since the `ARI` allows for circular analysis, still controlling for multiplicity of inferences.

In the following we show an analysis where clusters are defined by a supra-threshold-statistic rule. This is the typical case of cluster-wise analysis followed by a multiplicity correction based on Random Field Theory. Here we follow an alternative way: we provide lower bound for proportion for the estimate of active voxels.


## Syntax and parameters
The syntax of the function is (type `?ARIbrain::ARI` for more details)

`ARI(Pmap, clusters, mask = NULL, alpha = 0.05, Statmap = function(ix)  -qnorm(Pmap[ix]), summary_stat = c("max", "center-of-mass"),  silent = FALSE)`


The main input parameters of `ARI()` are:   

- `Pmap`: the map of p-values and 
- `clusters`: the map of cluster index.

The function accepts both character file names and 3D arrays. Therefore the minimal syntax is   
`ARI(Pmap, clusters)`

Others maps (parameters) are:   

- `mask` which is a 3D array of logicals (i.e.`TRUE`/`FALSE` means in/out of the brain). Alternatively, it may be a (character) nifti file name.  If omitted, all voxels are considered.  
- `Statmap` which is a 3D array of statistics (usually t-values) on which the summaries are based. File name is also accepted.



#  <a name="nii"> Performing the analysis from nifti (nii) files </a>

In order to perfom the analysis you need:   

- a `zstat.nii.gz` containing the test statistic used in the analysis 
- a `mask.nii.gz` (not mandatory, but usefull)
- a `cluster.nii.gz` image of cluster index.

## Making the map cluster.nii.gz with FSL

You simply need to run on the shell:

`cluster -z zstat1.nii.gz -t 3.2 -o cluster.nii.gz`

This will create the `cluster.nii.gz` that you need.

*hint*: In case it retun an error message like  
`cluster: error while loading shared libraries: libutils.so: cannot open shared object file: No such file or directory`,  
type into the shell (replacing the path with your own path of the file fsl.sh):  
`source /etc/fsl/5.0/fsl.sh`  
and try again.


Get a complete help for FSL at  
<https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster>


# ARI analysis
```{r}
library(ARIbrain)
pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARIbrain")
cluster_name <- system.file("extdata", "cluster_th_3.2.nii.gz", package="ARIbrain")
zstat_name <- system.file("extdata", "zstat.nii.gz", package="ARIbrain")
mask_name <- system.file("extdata", "mask.nii.gz", package="ARIbrain")
res_ARI=ARI(Pmap = pvalue_name, clusters= cluster_name,
    mask=mask_name, Statmap = zstat_name)
str(res_ARI)
```


# other ARI examples


## using arrays


```{r}
library(RNifti)
Tmap = readNifti(system.file("extdata", "zstat.nii.gz", package="ARIbrain"))
# compute p-values from Test statistic (refering to Normal distribution, right-side alternative)
Pmap=pnorm(-Tmap)
#Read the mask file. 
mask = RNifti::readNifti(system.file("extdata", "mask.nii.gz", package="ARIbrain"))
# Make sure that it is a logical map by: ()!=0
mask=mask!=0
#Create Clusters using a threshold equal to 3.2
Tmap[!mask]=0
clstr=cluster_threshold(Tmap>3.2)
table(clstr)
res_ARI=ARI(Pmap,clusters = clstr,mask = mask,Statmap = Tmap)
```




## Define threshold and clusters on the basis of concentration set (optimal threshold)


```{r}
hom=hommel::hommel(Pmap[mask])
(thr_p=hommel::concentration(hom))
(thr_z=-qnorm(thr_p))
Tmap[!mask]=0
clstr=cluster_threshold(Tmap>thr_z)
table(clstr)
res_ARI_conc=ARI(Pmap,clusters = clstr,mask = mask,Statmap = Tmap)
```
