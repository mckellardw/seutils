###############################
# seutils -  Seurat utility functions
# Written by: David McKellar
# version: 1.0
###############################

########################################
## Helpers for grabbing/searching feature names
########################################

# Quickly return genes/feature names from a Seurat object
Features <- function(
  SEU,
  assay=NULL
){
  if(is.null(assay)){
    assay <- SEU@active.assay
  }

  return(rownames(GetAssayData(SEU,assay=assay)))
}

# Check gene names for a pattern, using grep
grepGenes <- function(
  SEU,
  pattern="", # pattern to look for
  assay=NULL,
  filter.pattern=NULL, # pattern to remove
  sort.by = c("expression","abc"),
  verbose=T
){
  require(Seurat)
  require(dplyr)
  
  if(is.null(pattern)){
    if(verbose){message("Need a pattern to look for!")}
    return(NULL)
  }
  if(is.list(SEU)){
    if(verbose){message("Don't pass a list!")}
    return(NULL)
  }
  if(is.null(assay)){
    assay <- SEU@active.assay
  }

  genes = SEU@assays[[assay]]@counts@Dimnames[[1]]



  if(length(pattern)>1){
    if(verbose){
      message(paste0("Found ", length(genes), " features in the assay '",assay,"'..."))
      message(paste0("Looking for multiple patterns in these features..."))
    }
    out.genes <- lapply(
      pattern,
      FUN = function(PAT) genes[grep(pattern=PAT,x=genes)] #get genes
    )
    out.genes <- unlist(out.genes)
  }else{
    if(verbose){
      message(paste0("Found ", length(genes), " features in the assay '",assay,"'..."))
      message(paste0("Looking for '", pattern, "' in these features..."))
    }
    out.genes <- genes[grep(pattern=pattern,x=genes)] #get genes
  }

  if(length(out.genes)==0){
    message("Nothing found!\n")
    return(NULL)
  }

  if(!is.null(filter.pattern)){ # filter out filter.pattern
    if(verbose){message(paste0("Removing features containing '", filter.pattern,"' from output...\n"))}
    for(fp in filter.pattern){
      out.genes <- out.genes[!grepl(pattern=fp, x=out.genes)]
    }
  }

  if(is.null(sort.by[1])){
    if(verbose){message("Output genes not sorted...")}
  }else if(length(out.genes)==1){
    # Do nothing
  }else if(sort.by[1]=="expression"){
    if(verbose){message("Output genes sorted by expression...")}

    out.genes <- GetAssayData(
      SEU,
      assay=assay,
      slot="counts"
    )[out.genes,] %>%
        rowSums() %>% # sum counts across all cells
        sort(decreasing=T) %>% # sort genes according to total expression
        names() # get gene names back
  }else if(sort.by[1]=="abc"){
    if(verbose){message("Output genes sorted alphabetically...")}

    out.genes <- sort(out.genes, decreasing=T)
  }

  # Return matching gene names!
  return(
    out.genes
  )
}

# Convert ensembl IDs to gene IDs using a biomaRt reference
ens2gene <- function(
  ens=NULL, # vector of ensembl IDs to convert
  biomart.info=NULL, # biomaRt database
  ens.colname="ensembl_gene_id",
  gene.colname="mgi_symbol",
  min.cells,
  ncores=1,
  force.unique=F, # switch to make gene names unique (adds ".1", ".2", etc.)
  verbose=F
){
  require(dplyr)

  if(is.null(ens)){
    message("Need ensembl IDs to convert!")
    return(NULL)
  }
  if(is.null(biomart.info)){
    message("Need biomaRt reference for conversion!")
    return(NULL)
  }
  if(!ens.colname %in% colnames(biomart.info)){
    message("ensembl ID column not found in biomaRt reference. Check input for 'ens.colname'")
    return(NULL)
  }
  if(!gene.colname %in% colnames(biomart.info)){
    message("Gene ID column not found in biomaRt reference. Check input for 'gene.colname'")
    return(NULL)
  }

  if(ncores==1){
    out <- lapply(
      ens,
      FUN = function(ENS){
        tmp = biomart.info[biomart.info[[ens.colname]]==ENS,]
        if(nrow(tmp==1)){ # only found one entry that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])==1)){ # only found one unique gene name that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])>1)){ # only found one unique gene name that matches
          if(verbose){message(paste0("Found multiple matches for ", ENS,", returning the first one I see..."))}
          feat = tmp[1,gene.colname]
        }else{
          if(verbose){cat(paste0("Nothing found for ", ENS,", returning ensembl ID"))}
          feat = ENS
        }

        if(feat==""){
          if(verbose){message(paste0("No gene name found in '", gene.colname,"', for",ENS, ", returning ensembl ID."))}
          feat = ENS
        }

        return(feat)
      }
    ) %>%
      unlist()
  }else if(ncores>1){
    require(parallel)
    if(verbose){message(paste0("Running on ", ncores," threads..."))}

    out <- mclapply(
      ens,
      FUN = function(ENS){
        tmp = biomart.info[biomart.info[[ens.colname]]==ENS,]
        if(nrow(tmp==1)){ # only found one entry that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])==1)){ # only found one unique gene name that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])>1)){ # only found one unique gene name that matches
          if(verbose){message(paste0("Found multiple matches for ", ENS,", returning the first one I see..."))}
          feat = tmp[1,gene.colname]
        }else{
          if(verbose){cat(paste0("Nothing found for ", ENS,", returning ensembl ID"))}
          feat = ENS
        }

        if(feat==""){
          if(verbose){message(paste0("No gene name found in '", gene.colname,"', for",ENS, ", returning ensembl ID."))}
          feat = ENS
        }

        return(feat)
      },
      mc.cores = ncores
    ) %>%
      unlist()
  }

  if(force.unique){
    out <- make.unique(out)
  }

  return(
    out
  )
}


# Add a new assay from Market Matrix (.mtx) file to a Seurat object 
addNewAssayMtx <- function(
  SEU,
  assay=NULL,
  mtx=NULL,
  cells=NULL,
  features=NULL,
  current.cells=T, #whether to only add counts for current cells/barcodes in the Seurat object
  normalize.col=NULL, # metadata column used to normalize counts (useful for splitting counts across biotypes/etc)
  
  # ReadMtx() params
  cell.column = 1,
  feature.column=2,
  cell.sep = "\t",
  feature.sep = "\t",
  mtx.transpose=F,
  
  # CreateAssayObject() params
  min.cells=1, 
  
  
  verbose=F
){
  # arg checks
  if(is.null(mtx)){
    message("Need to pass mtx filepath!")
    return(SEU)
  }else if(is.null(cells)){
    message("Need to pass cells filepath!")
    return(SEU)
  }else if(is.null(features)){
    message("Need to pass features filepath!")
    return(SEU)
  }
  
  # lib requirements
  require(Seurat)
  
  # read mat
  tmp.mat <- ReadMtx(
    mtx = mtx,
    cells = cells,
    features = features,
    feature.column = feature.column,
    mtx.transpose = mtx.transpose
  ) 
  
  # get cell list for output matrix
  keep.cells <- intersect(Cells(SEU), colnames(tmp.mat))
  if(sum(Cells(SEU) %in% colnames(tmp.mat))==0){
    message("Cells not found in the new count matrix!")
    return(SEU)
  }else{
    if(verbose){
      cat(paste0(
        length(keep.cells)," out of ",length(Cells(SEU)), " cells found in new count matrix.\n"
      ))
    }
  }
  
  #TODO: comment
  if(current.cells){
    tmp.mat <- tmp.mat[,keep.cells]
  }else{
    message("Have not implemented current.cells=F yet!")
    tmp.mat <- tmp.mat[,keep.cells]
  }
  
  # Add zeroes for any cells not found in the new matrix
  if(length(Cells(SEU)) > length(keep.cells)){
    num.zeroes <- length(Cells(SEU)) - length(keep.cells)
    missing.cells <- Cells(SEU)[!Cells(SEU) %in% keep.cells]
    
    if(verbose){
      cat(paste0("Adding zeroes for ",num.zeroes," barcodes not found in the Seurat object...\n"))
    }
    
    # Add zeroes to the matrix for missing cell barcodes
    zero.mat <- matrix(
      0,
      nrow=nrow(tmp.mat),
      ncol=num.zeroes,
      dimname = list(
        rownames(tmp.mat),
        missing.cells
      ) 
    )
    tmp.mat <- cbind(
      tmp.mat,
      zero.mat
    )
    
    # Re-order matrix to match the input Seurat object
    tmp.mat <- tmp.mat[,Cells(SEU)]
  }
  
  
  SEU[[assay]] <- CreateAssayObject(
    counts = tmp.mat,
    min.cells = min.cells
  )
  
  # Custom normalization
  ##TODO- add multi-column normalization (sum of columns)
  if(!is.null(normalize.col)){
    tmp.denominator <- SEU@meta.data[,normalize.col]
    scale.factor <- 10000 #TODO: add as param?
    tmp.mat.normalized <- log1p(tmp.mat/tmp.denominator*scale.factor)
    
    print(table(rownames(tmp.mat.normalized) %in% rownames(tmp.mat)))
    
    SEU[[assay]]@data <- tmp.mat.normalized
  }
  
  return(SEU)
}

# Change feature names in an Assay in a Seurat object
changeAssayNames <- function(
  SEU,
  assay=NULL,
  new.features=NULL,
  verbose=F
){
  # lib requirements
  require(Seurat)
  
  # arg checks
  if(is.null(assay)){
    message("Need an assay to change!")
    return(SEU)
  }else if(is.null(new.features)){
    message("Need new feature names!")
    return(SEU)
  }
  
  # check feature lists to ensure lengths are the same
  old.features <- Features(
    SEU,
    assay=assay
  )
  
  if(length(new.features)!=length(old.features)){
    message("Feature list is the wrong length!")
    if(verbose){
      cat(paste0("Trying to replace ", length(old.features), " features with ", length(new.features), "!"))
    }
  }
  
  if(verbose){
    cat(paste0("Warning: removing normalized and scaled values if they were previously computed!\n"))
  }
  
  counts.mat <- GetAssayData(
    SEU, 
    assay=assay, 
    slot = "counts"
  )
  rownames(counts.mat) <- new.features
  
  SEU[[assay]] <- CreateAssayObject(
    counts = counts.mat
  )
  
  return(SEU)
}


# Split Assay into multiple Assays based on gene biotype, or other feature metadata
#TODO - finish building/rebuilding Seurat objects/assays and fix return()s
splitAssay <- function(
  SEU,
  assay=NULL, # Assay to be chopped up
  new.assay.names=NULL, # Vector of assay names
  feature.labels=NULL, # Vector of labels which are used to chop up the Assay
  return.new.seurat=FALSE, #Whether or not to return a fresh Seurat object (memory saving feature); all meta.data is transerred over
  verbose=F
){
  require(Seurat)

  # param checks
  if(is.null(feature.labels)){
    message("Need labels to split Assay by!")
  }
  if(is.null(assay)){
    assay <- SEU@active.assay
    if(verbose){message(paste0("Using ", assay, " as `assay`..."))}
  }

  #Check if features are same length as labels
  if(length(Features(SEU,assay=assay)) != length(feature.labels)){
    message("Features and feature labels are different lengths!")

    if(return.new.seurat){
      return(NULL)
    }else{
      return(SEU)
    }
  }

  # Extract matrix, and split based on unique labels in `feature.labels`
  mat <- GetAssayData(
    SEU, 
    assay=assay, 
    slot="counts"
  ) 

  split.list <- list()
  for(lab in unique(feature.labels)){
    split.list[[lab]] <- mat[feature.labels == lab, ]
  }

  # Return new assay(s)!
  if(return.new.seurat){
      return(NULL)
    }else{
      return(SEU)
    }
}


# Collapse cell/nuclei/spot counts for multimapped genes - these are features with a period ("\\.") in their name
collapseMultimappers <- function(
  SEU,
  assay=NULL,
  new.assay.name=NULL,
  verbose=F
){

  if(is.null(new.assay.name)){
    new.assay.name = paste0(assay,"_collpased")
    message("Using ",new.assay.name, " as new.assay.name...")
  }
  if(is.null(new.assay.name)){
    message("Need new.assay.name!")
    return(SEU)
  }

  SEU@active.assay <- assay

  multi.feats <- grepGenes( #Find genes with a period in them
    SEU,
    assay = assay,
    pattern="\\.",
    sort.by="abc",
    verbose=verbose
  )
  if(length(multi.feats)==0){
    message("No multimappers found!")
    return(SEU)
  }

  multi.patterns <- stringr::str_split( #extract actual gene names
    multi.feats,
    pattern = "\\.",
    n = 2
  ) %>%
    lapply(FUN=function(X) X[1]) %>%
    unlist() %>%
    unique()

  if(verbose){
    message(paste0("Found ", length(multi.patterns), " multimappers and ", length(multi.feats)," loci..."))
  }

  # Collapse counts for each gene
  mat.multi <- GetAssayData(
    SEU,
    assay=assay,
    slot="counts"
  )

  collapsed.list <- lapply(
    multi.patterns,
    FUN=function(X){
      tmp.genes = rownames(mat.multi)[grep(rownames(mat.multi),pattern=X)]
      tmp.mat = mat.multi[tmp.genes,]

      if(length(tmp.genes)==1){
        return(tmp.mat)
      }else{
        return(colSums(tmp.mat))
      }
    }
  )
  collapsed.mat <- do.call(rbind, collapsed.list) %>% as.sparse()
  rownames(collapsed.mat) <- multi.patterns

  # Add new assay with collapsed counts + the rest of the genes
  if(verbose){cat(paste0("Adding back ", nrow(collapsed.mat), " features...\n"))}

  solo.feats <- rownames(SEU)[!rownames(SEU)%in%c(multi.feats,multi.patterns)]

  out.mat <- rbind(
    GetAssayData(SEU,assay=assay, slot="counts")[solo.feats,],
    collapsed.mat
  )
  SEU[[new.assay.name]] <- CreateAssayObject(counts=out.mat)

  SEU@active.assay <- new.assay.name

  # Return Seurat object!
  return(SEU)
}


# Add spatial location for each cell/spot/bead in a Seurat object, given a whitelist of barcods/locations (X- & Y- coordinates)
#TODO: generalize whitelist formatting/reading in; add option for passing vector instead of file
#TODO: finish `as.fov` for better seurat handling
addSpatialLocation <- function(
  SEU,
  bc_map="~/txg_snake/resources/whitelists/visium-v1_coordinates.txt",
  bc_prefix=NULL, # prefix on cell barcodes; useful if cells have been merged/renamed
  loupe.json=NULL, # path to a .json file which contains the absolute spatial position (TXG Visium data only)
  reduction.name = "space",
  scale.factor = 1, # Factor to **multiply** the spatial coordinates by; useful for pixel-to-unit transformations
  assay=NULL,
  move.to.origin=F, # set lower bounds of spatial coordinates to zero
  as.fov=F, # Whether to add the coordinates as an `FOV`; needed for spatially-aware analyses w/ Seurat
  as.reduction=T, # Whether to add the spatial location as a reduction; useful for plotting
  verbose=F
){
  require(Seurat)
  require(dplyr)

  # Check assays
  if(is.null(assay) | !assay %in% Assays(SEU)){
    assay=SEU@active.assay
    if(verbose){message("Using the default assay...")}
  }

  # Read in whitelist coordinates
  if(file.exists(bc_map)){
  bc.coords <- read.table(
    bc_map,
    row.names = 1,
    col.names = c(
      "X",
      "Y"
    )
  )
  }else{
    cat("Whitelist file not found!")
    return(SEU)
  }
  
  if(scale.factor != 1){
    if(verbose)(message(paste0("Scaling coordinates by a factor of ", scale.factor)))
    bc.coords$X <- bc.coords$X * scale.factor
    bc.coords$Y <- bc.coords$Y * scale.factor
  }
  
  if(move.to.origin){
    if(verbose)(message(paste0("Moving to origin...\n")))
    bc.coords$X <- bc.coords$X - min(bc.coords$X)
    bc.coords$Y <- bc.coords$Y - min(bc.coords$Y)
  }
  
  #Add prefix to whitelist, if required
  if(!is.null(bc_prefix)){
    if(verbose)(message(paste0("Adding ", bc_prefix, " as a prefix to barcodes")))
    rownames(bc.coords) <- paste0(bc_prefix, rownames(bc.coords))
  }
  
  # Build reduction based on spot barcodes & whitelist
  bcs <- Cells(SEU)
  
  if(verbose){
    message(paste0(
      table(bcs %in% rownames(bc.coords))["TRUE"], " out of ", length(bcs), " barcodes found in whitelist"
    ))
  }

  #TODO- check/remove '-1'

  tmp.mat <- lapply(
    bcs,
    FUN=function(BC) bc.coords[BC,]
  ) %>%
    do.call(what=rbind)

  colnames(tmp.mat) <- paste0(reduction.name, 1:2)
  rownames(tmp.mat) <- bcs

  # Add spatial coordinates as a `reduction`
  if(as.reduction){
    # Add reduc to seurat obj
    SEU[[reduction.name]] <- CreateDimReducObject(
      embeddings=as.matrix(tmp.mat),
      key = paste0(reduction.name,"_"),
      assay=assay,
      global=TRUE
    )
    
    if(!is.null(loupe.json)){
      #TODO
      #read json in as a dataframe/mat
      subsetLoupeJson(
        SEU = SEU,
        reduction = reduction.name,
        json_path = loupe.json,
        verbose = verbose
      )
      # replace the entries in `tmp.mat` (relative positions) with the absolute positions
    }
  }
  
  # Add spatial coordinates as a `fov`
  if(as.fov){
    #TODO
  }
  
  return(SEU)
}


# Rotate spatial coordinates clockwise 90 degrees, N times 
rotateClockwise90N <- function(
    SEU,
    reduction = "space",
    N = 1, # number of times to rotate 90 degrees
    verbose=F
){
  require(Seurat)
  require(dplyr)
  
  for(i in 1:N){
    # grab coordinates
    tmp <- SEU[[reduction]]@cell.embeddings
    
    # map (x,y) to (-y,x)
    SEU[[reduction]]@cell.embeddings[,1] <- -1*tmp[,2]
    SEU[[reduction]]@cell.embeddings[,2] <- tmp[,1]
    
    # translocate coordinates back into 1st quadrant
    trans.factor = SEU[[reduction]]@cell.embeddings[,2]%>% na.omit()%>%range()%>%sum()%>%abs() # sum of min & max; abs() to make positive
    SEU[[reduction]]@cell.embeddings[,2] <- SEU[[reduction]]@cell.embeddings[,2] + trans.factor
  }
  
  return(SEU)
}


# Remove spots/beads that have fewer than `K` neighbors within `D` units of another spot/bead
removeSpatialSinglets <- function(
    SEU,
    reduction = "space",
    K = 1, # minimum number of nearest neighbors within distance D
    D = NULL, # minimum distance of a nearest neighbor; default is 1/100th of the x-axis range()
    verbose=F
){
  require(Seurat)
  require(dplyr)
  
  if(is.null(D)){
    D <- SEU[[reduction]]@cell.embeddings%>%na.omit()%>%range()%>%diff()/100
  }
  
  # Find euclidean distance matrix
  distmat <- as.matrix(dist(SEU[[reduction]]@cell.embeddings, method = "euclidean"))
  diag(distmat) <- NA # set diagonal to NA, to ignore self-distances
  
  # Find spots/beads without a nearest neighbor within D units
  cells.remove <- apply(
    distmat,
    MARGIN = 1,
    # FUN = function(X) min(X, na.rm = TRUE) >= D
    FUN = function(X) na.omit(X)%>%subset(subset = X<D)%>%length() < K
  )
  
  # Filter SEU
  if(length(cells.remove)>0){
    SEU <- subset(
      SEU,
      cells = Cells(SEU)[!cells.remove]
    )
  }else{
    cat("No singlets found!")
  }
  
  return(SEU)
}


# Convert relative (row/col) positions to absolute (X/Y) positions
spatialRel2Abs <- function( 
    SEU,
    reduction.rowcol.in="space",
    reduction.xy.out="space.xy",
    assay=NULL,
    json_path=NULL,
    verbose=FALSE
){
  require(jsonlite)
  
  # param checks
  if(is.null(SEU)){
    message("Need Seurat object to subset...")
  }
  
  if(is.null(json_path)){
    message("Need `json_path`...")
  }
  
  if(is.null(assay)){
    assay=SEU@active.assay
    if(verbose){message("Using the default assay...")}
  }else if(!assay %in% Assays(SEU)){
    assay=SEU@active.assay
    if(verbose){message(paste0("Assay `",assay, "` not found. Using the default assay..."))}
  }
  
  # read in json
  df <- jsonlite::fromJSON(
    txt = json_path,
    flatten = T
  )
  
  rc.coords <- df$oligo[!is.na(df$oligo$tissue),c("row", "col")]
  xy.coords <- df$oligo[!is.na(df$oligo$tissue),c("x", "y")]
  
  row_col <- SEU@reductions[[reduction.rowcol.in]]@cell.embeddings
  
  out.xy <- list()
  missing.spot.count = 0
  for(i in 1:nrow(row_col)){
    tmp.coords <- xy.coords[rc.coords[,1]==row_col[i,2] & rc.coords[,2]==row_col[i,1], ]
    if(is.null(tmp.coords)){
      out.xy[[ Cells(SEU)[i] ]] <- c(0,0)
      missing.spot.count = missing.spot.count + 1
    }else{
      out.xy[[ Cells(SEU)[i] ]] <- tmp.coords
    }
  }
  
  out.xy <- do.call(rbind,out.xy) #concatenate list into data.frame
  if(verbose){message("  ",c(ncol(SEU)-missing.spot.count)," out of ", ncol(SEU)," spots found...")}
  
  # add new reduction
  colnames(out.xy) <- paste(reduction.xy.out,c(1,2),sep="_")
  SEU[[reduction.xy.out]] <- CreateDimReducObject(
    embeddings=as.matrix(out.xy),
    assay=assay,
    key = paste0(reduction.xy.out,"_"),
    global=TRUE
  )
  
  gc()
  
  return(SEU)
}


# Subset Visium Seurat object based on spots manually selected via loupe browser
subsetLoupeJson <- function(
  SEU,
  reduction="space",
  spatial_format=c("row_col","x_y"),
  json_path=NULL,
  verbose=FALSE
){
  require(jsonlite)
  require(dplyr)
  
  # param checks
  if(is.null(SEU)){
    message("Need Seurat object to subset...")
  }
  
  if(is.null(json_path)){
    message("Need `json_path`...")
  }
  
  spatial_format <- spatial_format[1]
  
  # read in json
  df <- jsonlite::fromJSON(
    txt = json_path,
    flatten = T
  )
  
  if(spatial_format == "row_col"){
    out.coords <- df$oligo[!is.na(df$oligo$tissue),c("row", "col")]
  }else if(spatial_format == "x_y"){
    out.coords <- df$oligo[!is.na(df$oligo$tissue),c("x", "y")]
  }else{
    message("Incorrect `spatial_format` given...")
    return(SEU)
  }
  
  # paste coord pairs to make them unique for easy filtering
  out.coords <- apply(
    out.coords,
    MARGIN = 1, # each row-col or x-y pair
    FUN = function(COORDS){
      paste(COORDS[1], COORDS[2],sep = "_")
    }
  )%>%
    as.vector()
  
  # Get cell barcodes to keep 
  out.cells <- list()
  for(i in 1:nrow(SEU@reductions[['space']]@cell.embeddings)){
    tmp.cell = SEU@reductions[['space']]@cell.embeddings[i,]
    tmp.coords = paste(tmp.cell[2], tmp.cell[1],sep = "_")
    if(tmp.coords %in% out.coords){
      out.cells <- append(out.cells, rownames(SEU@reductions[['space']]@cell.embeddings)[i])
    }
  }
  
  # Subset Seurat object
  SEU <- subset(
    SEU,
    cells = unlist(out.cells)
  )
  
  gc()
  
  return(SEU)
}

########################################
## General Seurat workflow helpers
########################################
# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
npcs <- function(
  SEU,
  var.total=0.95,
  reduction="pca"
){
  if(is.null(SEU@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }

  tmp.var <- (SEU@reductions[[reduction]]@stdev)^2
  var.cut <- var.total*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }

  return(n.pcs)
}


# Preprocessing wrapper function
#   (1) NormalizeData(SEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
#   (2) FindNeighbors, RunUMAP (optional), FindClusters (optional)
seuPreProcess <- function(
  SEU=NULL,
  assay='RNA',
  n.pcs=50,
  res=0.8,
  nVarFeatures = 2000,
  run.clustering = F,
  run.umap=F,
  verbose=F
){
  if(is.null(SEU)){
    cat("Need a Seurat object to preprocess!\n")
    return(NULL)
  }
  if(!assay %in% Assays(SEU)){
    cat(paste0(assay, " not found in the seurat object! Not preprocessed.\n"))
    return(SEU)
  }

  # NormalizeData(SEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')

  SEU = NormalizeData(
    SEU
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = nVarFeatures,
    verbose = verbose
  ) %>% ScaleData(
    assay = assay
  ) 
  
  if(run.clustering){
    SEU = RunPCA(
      SEU,
      assay = assay,
      reduction.name = pca.name,
      reduction.key = pca.key,
      verbose = verbose,
      npcs = n.pcs
    )
    
    SEU = FindNeighbors(
      SEU,
      reduction = pca.name,
      dims = 1:n.pcs,
      force.recalc = TRUE,
      verbose = verbose
    )
    
    if(length(res)<1){
      message("Need resolution parameter to run clustering...\n Returning Seurat object without clusters.")
    }else{
      for(tmp.res in res){
        SEU = FindClusters(
          object = SEU,
          resolution = tmp.res,
          graph.name = paste0(assay,"_snn"),
          verbose = verbose
        )
      }
    }
  }
  
  if(run.umap){
    if(!pca.name %in% Reductions(SEU)){
      SEU = RunPCA(
        SEU,
        assay = assay,
        reduction.name = pca.name,
        reduction.key = pca.key,
        verbose = verbose,
        npcs = n.pcs
      )
    }
    
    umap.name = paste0('umap_', assay)
    
    SEU = RunUMAP(
      SEU,
      reduction = pca.name,
      dims = 1:n.pcs,
      verbose = verbose,
      reduction.name=umap.name
    )
    SEU@reductions[[umap.name]]@misc$n.pcs <- n.pcs
  }
  
  gc()
  
  return(
    tryCatch(
      SEU,
      error=function(e) NULL
    )
  )
}


# Add a new ident seurat metadata filed) based on a list of cell types
#
#     object:     seurat object
#     old.idents: name of the idents metadata you will be assigning cell types to
#     new.idents: vector of cell types, in order of cluster number
#     newName:    string of the new idents name
#
AddCellTypeIdents <- function(
    SEU=NULL, 
    old.name, 
    new.name=NULL, 
    new.idents, 
    verbose=FALSE
){
  old.idents = as.vector(names(table(SEU[[old.name]])))
  
  if(is.null(new.name)){
    cat("**Need a new.name for the new idents**\n")
  }else{
    SEU[[new.name]] <- as.vector(SEU[[old.name]])
    for(i in 1:length(old.idents)){
      if(verbose){cat("Adding ", new.idents[i],"...", sep = "")}
      SEU[[new.name]][ SEU[[new.name]]==old.idents[i] ] <- new.idents[i]
      if(verbose){cat("Done!\n", sep = "")}
    }
  }
  
  return(SEU)
}


# A 1-function wrapper for DoubletFinder
RunDoubletFinder <- function(
  SEU,
  NCORES=1,
  to.filter=F
){
  
  # Estimated Doublet Rate for each dataset
  edr <- estimateDoubletRate.DWM(seur.list = seu.list)/100 #use your own known EDR here
  
  for(i in 1:length(seu.list)){
    cat(' --------------------------------------------\n',
        '### DoubletFinder for dataset number ', i, '###\n',
        '--------------------------------------------\n')
    
    n.pcs = npcs(SEU)
    
    # SEU@reductions$umap_RNA@misc$n.pcs.used
    
    ## pK Identification (no ground-truth)
    bcmvn<- paramSweep_v3_DWM(
      seu=SEU, 
      PCs = 1:n.pcs, 
      num.cores = ncores
    ) %>% summarizeSweep(#sweep.res.list, 
      GT = FALSE
    ) %>% find.pK(#sweep.stats
    ) 
    
    # Pull out max of bcmvn
    pK <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)])) # ugly, but functional...
    
    ## Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(SEU$seurat_clusters) 
    
    nExp_poi <- round(edr*length(colnames(SEU)))  
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    gc()
  }
  
    rm(bcmvn, pK, nExp_poi, nExp_poi.adj)
    gc()
  
  
  SEU <- 
    doubletFinder_V3.DWM_v2(
      seu=SEU, 
      PCs = 1:SEU@reductions$umap_RNA@misc$n.pcs.used, 
      pN = 0.25, 
      pK= pK, 
      nExp = nExp_poi.adj,  
      reuse.pANN = F,
      classification.name='DF.individual', 
      pANN.name='DF.pANN.individual'
    )
  gc()
  
  #Filter seurat objects
  if(to.filter){
    return()
  }else{
    return(SEU)
  }
  
}


# Borrowed/adapted from the Marioni Lab, DropletUtils package (https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R)
#   (Had R version issues getting it to work as a dependency)
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
write_sparse <- function(
  path, # name of new directory
  x, # matrix to write as sparse
  barcodes=NULL, # cell IDs, colnames
  features=NULL, # gene IDs, rownames
  overwrite=F,
  verbose=F
){
  require(utils,quietly = T)
  require(Matrix,quietly = T)
  require(R.utils,quietly = T)

  if(!dir.exists(path)){
    dir.create(
      path,
      showWarnings=verbose,
      recursive = T
    )
  }

  if(is.null(barcodes)){
    barcodes=colnames(x)
  }
  if(is.null(features)){
    features=rownames(x)
  }
  # gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

  # gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
  mhandle <- file.path(path, "matrix.mtx")
  bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
  fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
  on.exit({
    close(bhandle)
    close(fhandle)
  })

  if(overwrite){
    if(verbose){message("Overwriting old files if they exist...")}

    if(file.exists(paste0(mhandle,".gz"))){
      file.remove(paste0(mhandle,".gz"))
    }
    if(file.exists(bhandle)){
      file.remove(bhandle)
    }
    if(file.exists(fhandle)){
      file.remove(fhandle)
    }
  }

  writeMM(x, file=mhandle)
  write(
    barcodes, 
    file=bhandle
  )
  write.table(
    features, 
    file=fhandle, 
    row.names=FALSE, 
    col.names=FALSE, 
    quote=FALSE, 
    sep="\t"
  )

  # Annoyingly, writeMM doesn't take connection objects.
  gzip(mhandle)

  return(NULL)
}


# TODO- haven't used this in a while, make sure it works...
# calculate entropy across groups in a seurat object
#     output: returns data.frame ("group.by" rows by 1 col)
seu_entropy <- function(
  SEU,
  group.by="sample",
  entropy.on="factorIDs",
  out.name=NULL,
  weighted=T,
  norm2one=T,
  verbose=T
){
  group.levels <- unlist(unique(SEU[[group.by]]))

  sub.meta <- SEU@meta.data[,c(group.by,entropy.on)]

  entropy.out <- list()
  for(lev in group.levels){
    if(verbose){cat("calculating entropy for ", lev, "... ",sep = "")}
    tmp <- sub.meta[sub.meta[[group.by]]==lev,] # subset metadata by sample

    entro.levels <- unlist(unique(sub.meta[[entropy.on]]))

    perc <- table(tmp[[entropy.on]])/nrow(tmp) # find proportions for each entropy level
    if(weighted){
      w <- (table(sub.meta[[entropy.on]])/nrow(sub.meta))[perc!=0] # weights for present levels
    }else{
      w <- rep(1,length(table(sub.meta[[entropy.on]])))[perc!=0]
    }

    perc <- perc[perc!=0] # remove zeroes

    entropy.out[[lev]] <- sum(-1*perc*log2(perc)*w) #calculate entropy

    if(verbose){cat("Done!\n",sep = "")}
  }

  # re-format entropy output
  if(is.null(out.name)){
    out.name <- paste0("entropy_",group.by,".by.",entropy.on)
  }
  entropy.out <- t(as.data.frame(entropy.out))
  colnames(entropy.out) <- out.name

  #replace NaN's with zeroes
  # entropy.out[is.nan(entropy.out)] <- 0

  if(norm2one){
    entropy.out <- entropy.out/log2(length(unique(sub.meta[[entropy.on]])))
  }

  return(entropy.out)
}


# TODO- haven't used this in a while, make sure it works...
# Calculate silhouette
seu_silhouette <- function(
  SEU,
  group.by,
  reduction,
  #TODO- add graph input as option
  meta.name.out=NULL,
  dims=1:10
){
  require(cluster)

  if(is.null(meta.name.out)){
    meta.name.out=paste0('sil.',reduction,'.',group.by)
  }
  if(is.null(group.by)){
    message("group.by is missing!")
    return(SEU)
  }
  if(is.null(reduction)){
    message("reduction is missing!")
    return(SEU)
  }

  # Calculate silhouette coefficient
  sil.out <- silhouette(
    x = as.numeric(x = as.factor(x = unlist(SEU[[group.by]]))),
    dist = dist(x = Embeddings(object = SEU, reduction=reduction)[,dims])
  )
  SEU[[meta.name.out]] <- sil.out[,3]

  return(SEU)
}

########################################
## biomaRt helper functions
########################################

# Basic function to convert human to mouse gene names
#   From @leonfodoulian (https://github.com/satijalab/seurat/issues/462)
#   https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
#
#   x = list of genes to be converted
#   Usage: mouse.genes <- lapply(X = human.genes, ConvertHumanGeneListToMM)
#
ConvertHumanGeneListToMM <- function(x){
  require(biomaRt)

  # Load human ensembl attributes
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Load mouse ensembl attributes
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # Link both datasets and retrieve mouse genes from the human genes
  genes.list = biomaRt::getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x ,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = TRUE
  )

  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])

  return(mouse.gene.list)
}

# Add gene biotype percentages to a seurat object, given a biomaRt object.
seu_biotypes_pct <- function(
  SEU,
  biomart=NULL, # biomaRt or data.frame containing gene biotypes
  gene.colname="GeneSymbol",
  biotype.colname="Biotype",
  add.as=c("metadata","assay"), # how percent features should be added
  assay="RNA",
  prefix="pct.",
  scale=c(1,100),
  verbose=TRUE
){

  if(is.null(biomart)){
    message("Need a list of gene biotypes! Nothing done.")
    return(SEU)
  }

  if(add.as[1]=="assay"){
    message("add.as='assay' is not yet implemented")
    new.assay.name = paste0(assay,".biotypes")
    return(SEU)
  }

  if(verbose){message(paste0("Adding gene biotype percentage values as ", add.as, " ..."))}

  biotypes = unique(biomart[[biotype.colname]])
  for(biotype in biotypes){
    tmp.mart =  biomart[biomart[[biotype.colname]]==biotype,] #subset out current gene biotype
    tmp.feat = unique( #get unique gene names which are present in the SEU object
      tmp.mart[[gene.colname]][tmp.mart[[gene.colname]]%in%Features(SEU, assay=assay)]
    )

    if(length(tmp.feat)==0){
      if(verbose){message(paste0("  No ", biotype, " genes found..."))}
    }else{
      if(add.as[1]=="metadata"){
        SEU <- PercentageFeatureSet(
          SEU,
          col.name = paste0(prefix, biotype),
          assay = assay,
          features = tmp.feat
        )
        
        if(scale[1] == 1){ #[0,1]
          SEU[[paste0(prefix, biotype)]] <- SEU[[paste0(prefix, biotype)]]/100
        }else if(scale[1] == 100){ #[0,100]
          #Do nothing...
        }else{
          message("Given scale was not found. Scaling to 100...")
        }
      }else if(add.as[1] == "assay"){
        #TODO
        message("Not implemented yet...")
      }else{
        message("`add.as` option not found... Try again.")
      }
    }
  }
  return(SEU)
}

# Compute gene biotype diversity on a seurat object, given a biomaRt object.
## Uses `vegan::diversity`
seu_biotypes_div <- function(
    SEU,
    biomart=NULL, # biomaRt or data.frame containing gene biotypes
    gene.colname="GeneSymbol",
    biotype.colname="Biotype",
    add.as=c("metadata","assay"), # how percent features should be added
    assay="RNA",
    slot="counts",
    index = c("shannon", "simpson", "invsimpson"),
    prefix=NULL,
    verbose=TRUE
){
  require(vegan)
  
  index = index[1]
  if(!index %in% c("shannon", "simpson", "invsimpson")){
    message("Given `index` is not supported by `vegan::diversity`")
    return(SEU)
  }
  
  if(is.null(biomart)){
    message("Need a list of gene biotypes! Nothing done.")
    return(SEU)
  }
  
  if(is.null(prefix)){
    prefix = paste0(index,".")
  }
  
  if(add.as[1]=="assay"){
    message("add.as='assay' is not yet implemented")
    new.assay.name = paste0(assay,".biotypes")
    return(SEU)
  }
  
  if(verbose){message(paste0("Adding gene biotype percentage values as ", add.as, " ..."))}
  
  biotypes = unique(biomart[[biotype.colname]])
  for(biotype in biotypes){
    tmp.mart =  biomart[biomart[[biotype.colname]]==biotype,] #subset out current gene biotype
    tmp.feat = unique( #get unique gene names which are present in the SEU object
      tmp.mart[[gene.colname]][tmp.mart[[gene.colname]] %in% Features(SEU, assay=assay)]
    )
    tmp.mat <- GetAssayData(
      SEU,
      assay = assay,
      slot = slot
    )[tmp.feat,]
    
    # if(slot=="counts"){ # normalize prior to computing diversity
    #   if(verbose){message("Normalizing counts...")}
    #   tmp.mat <- apply(
    #     tmp.mat,
    #     MARGIN = 2, 
    #     FUN = function(X) X/sum(X)
    #   )
    # }
    
    if(length(tmp.feat)==0){
      if(verbose){message(paste0("  No ", biotype, " genes found..."))}
    }else{
      # Calculate Shannon diversity index for each column
      tmp.div <- apply(
        tmp.mat,
        MARGIN = 2, 
        FUN = diversity, 
        index = index
      )
      
      if(add.as[1]=="metadata"){
        SEU[[paste0(prefix, biotype)]] <- tmp.div
      }else if(add.as[1] == "assay"){
        #TODO
        message("Not implemented yet...")
      }else{
        message("`add.as` option not found... Try again.")
      }
    }
  }
  return(SEU)
}
########################################
## Other helpers...
########################################

# Run PHATE on reduction, with a Seurat objects
seuPHATE <- function(
  SEU=NULL,
  reduction="pca",
  ndims=50,
  reduction.name=NULL, # output reduction name
  reduction.key=NULL, # output reduction key
  # phate() defaults
  ndim = 2,
  knn = 5,
  decay = 40,
  n.landmark = 2000,
  gamma = 1,
  t = "auto",
  mds.solver = "sgd",
  knn.dist.method = "euclidean",
  knn.max = NULL,
  init = NULL,
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max = 100,
  npca = 100,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed = NULL,
  potential.method = NULL,
  k = NULL,
  alpha = NULL,
  use.alpha = NULL
){
  require(Seurat)
  require(phateR)

  dims <- SEU@reductions[[reduction]]@stdev %>% order(decreasing = T)
  dims <- dims[1:ndims]

  # run PHATE
  tmp.phate <- phate(
    SEU@reductions[[reduction]]@cell.embeddings[,dims], #cells x reduc dims
    ndim = ndim,
    knn = knn,
    decay = decay,
    n.landmark = n.landmark,
    gamma = gamma,
    t = t,
    mds.solver = mds.solver,
    knn.dist.method = knn.dist.method,
    knn.max = knn.max,
    init = init,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed = seed,
    potential.method = potential.method,
    k = k,
    alpha = alpha,
    use.alpha = use.alpha
  )

  # Set reduction name & key for SEU
  if(is.null(reduction.name)){
    reduction.name <- paste0("phate_",reduction)
  }
  if(is.null(reduction.key)){
    reduction.key <- paste0("phate_",reduction,"_")
  }

  colnames(tmp.phate$embedding) <- paste0(reduction.key, 1:2)
  tmp.phate$params$data <- NULL

  SEU[[reduction.name]] <- CreateDimReducObject(
    embeddings = tmp.phate$embedding,
    key = reduction.key,
    assay = 'RNA',
    misc = tmp.phate$params
  )

  # Add std dev to reduction
  SEU@reductions[[reduction.name]]@stdev <-
    apply(SEU@reductions[[reduction.name]]@cell.embeddings, 2, sd) #find std dev for phate vals

  return(SEU)
}
