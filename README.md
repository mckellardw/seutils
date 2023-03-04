# seutils
Utility functions for Seurat


#### Wrapper/helper functions written to smooth out Seurat pipelines... See [this wonderful website](https://satijalab.org/seurat/index.html) from the Satija Lab to learn more about `Seurat`!

## **seutils.R** (seurat utils)
`Seurat` utility functions, to make your life easier.

### Generic Seurat utils

#### `Features()`
Quickly grab feature names from any `Assay` within a `Seurat` object - analogous to `Seurat::Cells()`

#### `grepGene()`
Quick function to look for genes in a `Seurat` object that match a certain pattern. Can also set a `filter.pattern` to exclude genes which contain the string. Built on `grep` and `Seurat`

#### `ens2gene()`
Convert a list of ensembl IDs to gene symbols. Compatible with outputs from `biomaRt`.

#### `npcs()`
Calculate how many principal components account for X% of the variance.

#### `collapseMultimappers()`
Quickly collapse multimapper genes in a `Seurat` assay (non-unique features, marked by a "." in the feature name). Note, the literal regex used is `\\.`

#### `seuPreProcess()`
Generic single-cell pipeline, based on [this vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) from the Satija Lab

#### `AddCellTypeIdents()`
Add cell type labels given a clustering output


### Spatial tools
#### `addSpatialLocation()`
Add spatial location for each bead/spot within a Seurat object. Written to work with the whitelists from 10x Genomics (for Visium) and Curio Biosciences (for Seeker/SlideSeq)

#### `removeSpatialSinglets()`
Strategy to remove singlets based on spatial position/nearest neighbors - filters out beads/spots with fewer than `K` neighbors with `D` distance units

#### `rotateClockwise90N()`
Rotate the spatial embedding of a sample clockwise, 90 degrees, `N` times (i.e. N=2 to rotate 180 degrees clockwise). Useful if you have Visium/SlideSeq samples that were sectioned in different orientations.

#### `spatialRel2Abs()`
Convert relative (row/col) positions to absolute (X/Y) positions

#### `subsetLoupeJson()`
Subset Visium Seurat object based on spots manually selected via loupe browser


### Single-cell
TODO

## **seuplots.R** (seurat plots)
Additional plots not included in Seurat including

#### `visListPlot()`
Wrapper function for neatly plotting multiple Visium objects and multiple features into a grid. Input is a `list()` of `Seurat` objects. Includes automatic color scaling and lots of options for customization. Built on `patchwork`. Combine with `coord_fixed()` to maintain true spatial coordinates.
```
visListPlot <- function(
  seu.list, # list of Seurat objects
  features=NULL, # genes/metadata entries/etc to plot
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial',
  reduction="space",
  slot="data",
  legend.position="bottom",
  pt.size=1,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  colormap="viridis", # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.85), # color for na.value (spot where gene is not detected)
  verbose=FALSE
)
```
**Usage notes on `visListPlot()`:**
- Looks for the spatial position in the `reduction` passed (Check `SEU@reductions`), which means this also works with PCA/UMAP/PHATE/etc embeddings!
- Only works with **continuous** variables for now (i.e. gene expression, nCounts, deconvolution estimates) - does NOT work with cluster IDs
- Keeps the same color scale for all samples in each feature; use `colormap.same.scale=T` to force the same color scale across all samples AND features (i.e, when working with scaled gene expression values)
- The `na.value` is the color used for cells/spots/beads where the feature value is either missing or `NA`. I recommend using the default light gray option (`gray(0.85)`)


#### `visCoMap()`
Plot co-expression of any two genes/features, where "co-expression" can be abstracted to an generic function (default is just `prod`, to see expression of Gene_A times expression of Gene_B). Built on top of visListPlot(), so it also takes a list of `Seurat` objects.
```
visCoMap <- function(
  seu.list,
  features=NULL, # pair of genes
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial', # Either a single assay, or a pair of assays (useful for spliced vs. unspliced plotting)
  reduction="space",
  slot="data", # Either a single slot, or a pair of slots
  legend.position="bottom",
  pt.size=1,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  colormap=NULL, # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.69), # color for na.value (spot where gene is not detected)
  comap.fxn = prod,
  coex.name = NULL, # plot title for computed co-expression values
  verbose=FALSE
)
```
#TODO
- Volcano plots (ggvolcano_v2)
- Silhouette plot

