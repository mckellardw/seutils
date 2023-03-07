#############################################################################
# seuplots -  Seurat-compatible plotting functions
# Written by: David McKellar
#############################################################################

#############################################################################
## QC plots  ----------------------------------------------------------------

# Sorted scatter plot showing the number of unique molecules (UMIs) captured in each cell/spot/bead/etc.
kneePlot <- function(
    SEU,
    nUMI = "nCount_STAR", #nUMIs metadata column
    assay = NULL,
    alpha = 0.2,
    umi.threshold = NULL,
    umi.threshold.color = "red",
    color.by = NULL, #TODO
    verbose = FALSE
){
  #TODO
  # Warning message:
  #   In xtfrm.data.frame(x) : cannot xtfrm data frames
  
  if(is.null(assay)){
    assay <- SEU@active.assay
  }
  
  # add plottting order
  # SEU$plotting_order <- order(SEU[[nUMI]], decreasing = TRUE)
  
  out.plot <- ggplot(
    SEU@meta.data[order(SEU[[nUMI]], decreasing = TRUE), ],
    aes(
      x = 1:ncol(SEU),
      y = .data[[nUMI]]
    )
  ) +
    geom_point(
      alpha=alpha
    )+
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if(!is.null(umi.threshold) & is.numeric(umi.threshold)){
    cell.count = sort(table( SEU[[nUMI]] > umi.threshold )) # nCells in cell.count[1]
    
    out.plot +
      geom_hline(
        yintercept = umi.threshold,
        color = umi.threshold.color
      )+
      geom_vline(
        xintercept = cell.count[1],
        color = umi.threshold.color
      )
  }else{
    return()
  }
}

## Spatial/Visium plot functions --------------------------------------------


# Generate feature plots given a list of Visium/Seurat objects
seuListPlot <- function(
  seu.list,
  features=NULL,
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial',
  reduction="space",
  slot="data",
  legend.position="bottom",
  pt.size=1,
  order=FALSE,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  x.window = NULL, # list of xlim parameters for each sample; if length is one, # TODO:will use the same params for each samples
  y.window = NULL, # list of ylim parameters for each sample; if length is one, # TODO:will use the same params for each samples
  flip_row_col=F, #Default is samples in rows, features in cols. Set to `T` for samples in cols, features in rows
  colormap="viridis", # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.85), # color for na.value (spot where gene is not detected)
  min.value=10^-100, #minimum value use to label "NA" spots
  verbose=FALSE
){
  require(dplyr)
  require(Seurat)
  require(ggplot2)
  require(viridis)
  
  if(is.null(seu.list)){message("List of Seurat objects to plot is NULL!")}
  
  if(verbose){cat(paste0("Plotting data, using the assay ",assay,"\n"))}
  
  # Check colormap
  if(length(colormap)==1){ # single colormap option
    if(!colormap %in% c("magma","inferno","plasma","viridis","cividis","A","B","C","D")){
      message(paste0("Color palette '", colormap, "' not found!"))
      return(NULL)
    }else{
      colormap = rep(
        x=colormap,
        length=length(features)
      )%>%
        as.list()
    }
  }else if(length(colormap)>1 & !is.list(colormap)){ 
    colormap = rep(
      x=list(colormap),
      length=length(features)
    )
  }else if(is.list(colormap)){ # list for feature-specific colormaps
    if(length(colormap)!= length(features)){
      message(paste0("Wrong number of color palettes!"))
      return(NULL)
    }
  }
  
  # Alternate feature names
  if(is.null(alt.titles)){
    alt.titles=features
  }
  
  # Check slot, assay length(s)
  if(length(assay)==1){
    assay = rep(assay,length(features))
  }else if(length(assay) != length(features)){
    message("`assay` and `features` lengths don't match!")
  }
  
  if(length(slot)==1){
    slot = rep(slot,length(features))
  }else if(length(slot) != length(features)){
    message("`slot` and `features` lengths don't match!")
  }
  
  # Check x.window & y.window params
  if(is.null(x.window)){
    x.window <- rep(list(c(NA,NA)),length(seu.list))
  }else if(!is.list(x.window) & length(x.window)==2){
    x.window <- rep(list(x.window),length(seu.list))
  }else if(is.list(x.window)&length(x.window)==length(seu.list)){
    #Do nothing
  }else{
    message("Error with `x.window` - incorrect parameterization!")
  }
  
  if(is.null(y.window)){
    y.window <- rep(list(c(NA,NA)),length(seu.list))
  }else if(!is.list(y.window) & length(y.window)==2){
    y.window <- rep(list(y.window),length(seu.list))
  }else if(is.list(y.window)&length(y.window)==length(seu.list)){
    #Do nothing
  }else{
    message("Error with `y.window` - incorrect parameterization!")
  }
  
  # Check for spatial info in Reductions
  ##TODO: fix the assumption to look just in seu.list[[1]]
  if(is.null(reduction)){
    reductions <- Reductions(seu.list[[1]])[1]
    if(verbose){message(paste0("Using ", reduction, " as spatial location..."))}
  }else if(! reduction %in% Reductions(seu.list[[1]])){
    reductions <- Reductions(seu.list[[1]])[1]
    if(verbose){message(paste0("Reduction not found. Using ", reduction, " as spatial location..."))}
  }
  
  # Check for genes
  tmp.features = paste0("tmp.",features) # place-holder name for features; allows assay-specific feature plotting in FeaturePlot
  
  seu.list <- lapply(
    seu.list,
    FUN = function(SEU){
      for(i in 1:length(features)){
        if(features[i] %in% colnames(SEU@meta.data)){
          SEU@meta.data[,tmp.features[i]] <- SEU@meta.data[,features[i]]
          
        }else if(!features[i] %in% Features(SEU, assay = assay[i])){
          # Stick a bunch of zeroes into the metadata to act like an undetected gene
          SEU@meta.data[,tmp.features[i]] <- rep(0,nrow(SEU@meta.data))
          
        }else{
          SEU@meta.data[,tmp.features[i]] <- GetAssayData(
            SEU,
            assay=assay[i],
            slot=slot[i]
          )[features[i],]
        }
      }
      return(SEU)
    }
  )
  
  # Get expression limits for each gene, across all datasets
  gene.lims <- mapply(
    FUN = function(FEAT, ASS, SLOT){
      out.max <- lapply(
        seu.list,
        FUN = function(SEU){
          if(FEAT %in% colnames(SEU@meta.data)){
            return(
              max(
                SEU@meta.data[[FEAT]],
                na.rm = TRUE
              )
            )
            
          }else if(FEAT %in% Features(SEU, assay=ASS)){
            return(
              max(
                GetAssayData(SEU,assay=ASS,slot=SLOT)[FEAT,],
                na.rm = TRUE
              )
            )
            
          }else{
            if(verbose){message(FEAT, " not found!")}
            return(0)
          }
        }
      ) %>% 
        unlist() %>% 
        max()
      
      if(verbose){message(paste0("Using ", out.max, " as the max value for ", FEAT,"\n"))}
      return(c(min.value, out.max))
    },
    SIMPLIFY = F,
    FEAT = features,
    ASS = assay,
    SLOT = slot
  )
  
  # Set colormap scale to the global min/max of all features 
  if(colormap.same.scale){
    gene.lims <- lapply(
      gene.lims,
      FUN=function(X) c(min.value, max(unlist(gene.lims)) )
    )
  }
  
  # Get plot heights
  #TODO- build in coord_fixed() 
  if(abs.heights){
    # Select axis for setting plot heights/widths; dependent on orientation of the plot grid
    if(flip_row_col){
      column = 1
    }else{
      column = 2
    }
    
    heights <- lapply(
      seu.list,
      FUN=function(SEU) abs(
        diff(
          range(
            SEU@reductions[[reduction]]@cell.embeddings[,column]
          )
        )
      )
    ) %>% unlist()
    if(verbose){
      message(paste0("Using these plot ",c("widths","heights")[column],":"))
      print(heights)
    }
  }else{
    heights <- rep(1,length(features))
  }
  
  # Plot
  plot.list <- list()
  for(i in 1:length(features)){
    tmp <- lapply(
      seu.list,
      FUN = function(SEU){
        tmp.plot = FeaturePlot(
          SEU,
          slot = slot[i],
          order = order,
          features = tmp.features[i],
          pt.size = pt.size,
          reduction = reduction
        ) +
          theme(
            plot.margin = unit(rep(0,4), "inches"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            plot.title = element_blank(),
            legend.position = legend.position,
            legend.title = element_text(size=font.size,face="bold", hjust=0.5),
            legend.text = element_text(size=font.size,face="bold")
          )
        #TODO
        # lims(
        #   x=x.window[[i]],
        #   y=y.window[[i]]
        # )
        
        # set colormap
        if(length(colormap[[i]])==1){
          if(verbose){message(paste0("Using viridis palette ", colormap))}
          tmp.plot = tmp.plot + 
            scale_color_viridis(
              limits = unlist(gene.lims[i]),
              option = colormap[[i]],
              direction = colormap.direction,
              na.value = na.value
            )
        }else{
          if(verbose){message(paste0("Using custom color gradient"))}
          
          # Flip colormap if direction is `-1`
          if(colormap.direction==-1){
            colormap = rev(colormap[[i]])
          }
          
          tmp.plot = tmp.plot + 
            scale_color_gradientn(
              limits = unlist(gene.lims[i]),
              colors = colormap[[i]],
              na.value = na.value
            )
        }
        return(tmp.plot)
      }
    )
    
    if(!flip_row_col){
      tmp[[1]] <- tmp[[1]] +
        theme(
          plot.title = element_text(
            size=font.size,
            face="bold.italic",
            vjust=1
          )
        ) +
        labs(title=alt.titles[i])
    }
    
    plot.list[[i]] <- tmp
  }
  
  if(flip_row_col){ #samples=columns; features=rows
    
    # Add sample titles
    for(i in 1: length(plot.list[[1]])){
      plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
        theme(
          plot.title = element_text(
            size=font.size,
            face="bold",
            vjust=1
          )
        ) +
        labs(title=sample.titles[i])
    }
    
    # Feature titles on y-axis
    for(j in 1:length(plot.list)){
      plot.list[[j]][[1]] <- plot.list[[j]][[1]] +
        theme(
          axis.title.y = element_text(
            size=font.size,
            face="bold.italic",
            color="black",
            hjust=0.5,
            vjust=0.5,
            angle=axis.title.angle.y
          )
        )
      
      if(!is.null(alt.titles)){ # add feature titles
        plot.list[[j]][[1]] <- plot.list[[j]][[1]] +
          labs(y=alt.titles[j])
      }
    }
    
    # TEST
    # if(abs.heights){
    #   plot.list <- lapply(
    #     plot.list,
    #     FUN = function(PLOT) align_patches(PLOT, align)
    #   )
    # }
    
    #Wrap
      plot.list <- lapply(
        plot.list,
        FUN = function(X){
          wrap_plots(
            X,
            nrow=1,
            heights=heights,
            guides="collect"
          )&theme(
            legend.position=legend.position,
            legend.margin = margin(0,0,0,0,"inches")
          )
        }
      )
  }else{ #samples=rows; features=columns
    # Feature axis title
    for(i in 1:length(plot.list[[1]]) ){
      plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
        theme(
          axis.title.y = element_text(
            size = font.size,
            face = "bold",
            color = "black",
            hjust = 0.5,
            vjust = 0.5,
            angle = axis.title.angle.y
          )
        )
      
      if(!is.null(sample.titles)){ # add sample titles
        plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
          labs(y=sample.titles[i])
      }
    }
    
    # TEST
    # if(abs.heights){
    #   plot.list <- lapply(
    #     plot.list,
    #     FUN = function(PLOT) align_patches(PLOT, align)
    #   )
    # }
    
    # Sample axis title
    plot.list <- lapply(
      plot.list,
      FUN = function(X){
        wrap_plots(
          X,
          ncol=1,
          widths=heights,
          guides="collect"
        )&theme(
          legend.position=legend.position,
          legend.margin = margin(0,0,0,0,"inches")
        )
      }
    )
  }
  
  # # TEST
  # if(abs.heights){
  #   plot.list <- lapply(
  #     plot.list,
  #     FUN = function(PLOT) align_patches(PLOT, align)
  #   )
  # }
  
  if(combine){
    return(
      wrap_plots(
        plot.list,
        nrow=nrow,
        ncol=ncol,
        guides=NULL
      )
    )
  }else{
    return(plot.list)
  }
}


# Plot co-expression of any two genes/features
## Built on top of seuListPlot()
seuCoMap <- function(
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
  na.value=gray(0.85), # color for na.value (spot where gene is not detected)
  min.value=10^-100, #minimum value use to label "na" spots
  comap.fxn = prod,
  coex.name = NULL, # plot title for computed co-expression values
  include.scatter = F,
  scatter.group.by= 'orig.ident',
  scatter.theme= NULL, # default is theme_minimal()
  verbose=FALSE
){
  require(dplyr)
  require(Seurat)
  require(ggplot2)
  require(viridis)
  
  if(is.null(seu.list)){message("List of Seurat objects to plot is NULL!")}
  
  if(verbose){cat(paste0("Plotting Visium data, using the assay ",assay,"!\n"))}
  
  # Check colormap
  if(is.null(colormap)){
    colormap=list(
      mckolors$blue_ramp, 
      mckolors$red_ramp, 
      mckolors$purple_ramp
    )
  }
  
  # Alternate feature names
  if(is.null(alt.titles)){
    alt.titles=features
  }
  
  # Check `features` length
  if(length(features)!=2){
    message(paste0("Only supports plotting co-=expression of 2 features... for now..."))
  }
  
  # Check slot, assay length(s)
  if(length(assay)==1){
    assay = rep(assay,2)
  }
  if(length(slot)==1){
    slot = rep(slot,2)
  }
  
  # Compute co-expression values and add back to Seurat objects (as entry in meta.data)
  if(is.null(coex.name)){
    coex.name = "coexpression.tmp.values"
  }
  
  maps.out <- lapply(
    seu.list,
    FUN = function(SEU){
      tmp.coex <- list()
      
      DefaultAssay(SEU) <- assay[1]
      if(features[1] %in% Features(SEU,assay=assay[1])){
        tmp.coex[[1]] <- FetchData(
          object = SEU,
          vars = features[1],
          slot = slot[1]
        ) 
      }else{
        if(verbose){message(paste0(features[1], " is missing from `", assay[1], "`"))}
        tmp.coex[[1]] <- rep(0, ncol(SEU))
        
      }
      
      DefaultAssay(SEU) <- assay[2]
      if(features[2] %in% Features(SEU,assay=assay[2])){
        tmp.coex[[2]] <- FetchData(
          object = SEU,
          vars = features[2],
          slot = slot[2]
        )  
      }else{
        if(verbose){message(paste0(features[2], " is missing from `", assay[2], "`"))}
        tmp.coex[[2]] <- rep(0, ncol(SEU))
      }
      
      tmp.coex <- do.call(
        cbind,
        tmp.coex
      ) %>%
        apply(
          # X = tmp.coex,
          MARGIN = 1,
          FUN = prod # Function for co-expression
        ) 
      
      # SEU <- AddMetaData(
      #   object = SEU,
      #   metadata = tmp.coex,
      #   col.name = coex.name
      # )
      SEU@meta.data[,coex.name] <- tmp.coex
      return(SEU)
    }
  ) %>%
    seuListPlot( # Plot all three with seuListPlot
      # seu.list=seu.list,
      features=c(features,coex.name), # pass two features plus the co-expression values
      alt.titles=c(alt.titles, coex.name),
      sample.titles=sample.titles,
      reduction=reduction,
      assay=c(assay,assay[1]), #extra assay entry for the coex data
      slot=c(slot,slot[1]),#extra assay entry for the coex data
      legend.position=legend.position,
      pt.size=pt.size,
      font.size=font.size,
      axis.title.angle.y=axis.title.angle.y,
      combine=combine,
      abs.heights=abs.heights,
      nrow=nrow,
      ncol=ncol,
      colormap=colormap,
      colormap.direction=colormap.direction,
      colormap.same.scale=colormap.same.scale,
      na.value=na.value,
      min.value=min.value,
      verbose=verbose
    )
  
  #Scatter plot
  if(include.scatter){
    if(is.null(scatter.theme)){
      scatter.theme <- theme_minimal()
    }
    
    # Get expression limits for each gene, across all datasets
    gene.lims <- mapply(
      FUN = function(FEAT, ASS, SLOT){
        out.max <- lapply(
          seu.list,
          FUN = function(SEU){
            if(FEAT %in% Features(SEU, assay=ASS)){
              return(max(GetAssayData(SEU,assay=ASS,slot=SLOT)[FEAT,]))
            }else if(FEAT %in% colnames(SEU@meta.data)){
              return(max(SEU@meta.data[,FEAT]))
            }else{
              if(verbose){message(FEAT, " not found!")}
              return(0)
            }
          }
        ) %>% 
          unlist() %>% 
          max()
        
        return(c(0, out.max))
      },
      SIMPLIFY = F,
      FEAT = features,
      ASS = assay,
      SLOT = slot
    )
    
    scatter.out <- lapply(
      seu.list,
      FUN=function(VIS){
        
        # tmp.nonzeros <- GetAssayData(
        #   VIS,
        #   assay=assay,
        # )[c(features[1],features[2]),]%>%
        #   apply(
        #     MARGIN=2,
        #     FUN=function(X) data.frame(
        #       FEAT1 = X[1]>0 & X[2]==0,
        #       FEAT2 = X[1]==0 & X[2]>0,
        #       DP = X[1]>0 & X[2]>0
        #     )
        #   )%>%
        #   do.call(what=rbind)
        
        # tmp.label <- paste0(
        #   "Notch1+: ",table(tmp.nonzeros$NOTCH)["TRUE"],"/",sum(tmp.nonzeros),"\n",
        #   "Gm13568+:",table(tmp.nonzeros$GM)["TRUE"],"/",sum(tmp.nonzeros),"\n",
        #   "Double Pos: ",table(tmp.nonzeros$DP)["TRUE"],"/",sum(tmp.nonzeros)
        # )
        
        # FeatureScatter(
        #   VIS,
        #   shuffle = T,
        #   jitter = F,
        #   feature1 = features[1],
        #   feature2 = features[2],
        #   plot.cor = F,
        #   pt.size = pt.size,
        #   cols=c("#c42a2e", "#d8a837", "#2e6198"),#TODO: parameterize
        #   group.by = scatter.group.by
        # )+
        # geom_label(
        #   x=1,
        #   y=5.2,
        #   size = small.font/ggplot2::.pt,
        #   # size = small.font/4,
        #   label=tmp.label
        # )+
        # 
        # tmp.df <- cbind(
        #   VIS@meta.data,
        #   t(
        #     GetAssayData(
        #       VIS,
        #       assay=assay[1],
        #       slot=slot[1]
        #     )[features[1],] 
        #   ),
        #   t(
        #     GetAssayData(
        #       VIS,
        #       assay=assay[2],
        #       slot=slot[2]
        #     )[features[2],] 
        #   )
        # )
        # 
        # #Plot
        # ggplot(
        #   tmp.df[sample(rownames(tmp.df)),],
        #   aes_string(
        #     x=features[1],
        #     y=features[2]
        #   )
        # )+
        #   geom_point(
        #     alpha=0.7
        #   )+
        
        #For undetected features...
        for(i in 1:length(features)){
          if(!features[i] %in% Features(VIS,assay=assay[i])){
            features[i] <- paste0("tmp_",features[i])# in case gene starts with a number...
            VIS@meta.data[[features[i]]]<- rep(0, length(Cells(VIS)))
          }
        }
        
        tmp.plot <- FeatureScatter(
          VIS,
          shuffle = T,
          jitter = F,
          feature1 = features[1],
          feature2 = features[2],
          plot.cor = T,
          pt.size = pt.size,
          cols = "black",
          # cols=c("#c42a2e", "#d8a837", "#2e6198"),#TODO: parameterize
          group.by = scatter.group.by
        )+
          geom_smooth(
            method="lm",
            formula = "y ~ x",
            color="black"
          )+
          xlim(gene.lims[[features[1]]])+ 
          ylim(gene.lims[[features[2]]])+
          labs(
            x=stringr::str_remove(features[1],pattern="tmp_"),
            y=stringr::str_remove(features[2],pattern="tmp_")
          )+
          scale_color_manual(
            values = c("#c42a2e", "#d8a837", "#2e6198"),#TODO: parameterize
          )+
          guides(
            color = guide_legend(override.aes = list(size=2))
          )+
          scatter.theme+
          theme(
            axis.title.y = element_text(
              face="bold.italic",
              hjust=0.5,
              size=font.size
            ),
            axis.title.x = element_blank(),
            legend.position = "right",
            plot.margin = unit(rep(0,4),"cm")
          )+
          coord_fixed(
            ratio=max(gene.lims[[features[1]]])/max(gene.lims[[features[2]]])
          )
        
        tmp.plot$layers[[1]]$aes_params$alpha <- 0.5
        
        return(tmp.plot)
      }
    )
    
    scatter.out[[length(scatter.out)]] <- scatter.out[[length(scatter.out)]]+
      theme(
        axis.title.x = element_text(
          face="bold.italic",
          hjust=0.5,
          size=font.size
        )
      )
    
    scatter.out <- wrap_plots(
      scatter.out,
      ncol=1,
      heights=rep(1,length(scatter.out)),
      guides="collect"
    )&theme(
      legend.position="none"#TODO
    )
    
    return(
      # scatter.out
      wrap_plots(
        maps.out,
        scatter.out,
        nrow=1,
        widths=c(3,1.5)
      )
    )
  }else{
    return(maps.out)
  }
}

# Generate split violin plots for Visium data to look at co-occurence of cell types & gene expression
#   Co-occurence is based on BayesPrism theta values
visCoocVlnPlot<-function(
  SEU,
  assay.ge='Spatial',# gene expression assay
  assay.ct, #cell type assay
  slot.ge="data",
  features.lr, # ligand-receptor(in that order)  pair
  features.ct=NULL, # 2 celltypes to compare
  bp.thresh=0.1, # Lower threshold for bayesprism theta values to determine presence of a cell type
  scale.data=T,
  scale.min=-2.5,
  scale.max=2.5,
  width=0.9,
  #TODO
  legend.position="bottom",
  pt.size=0,
  font.size=8,
  verbose=T
){
  require(ggplot2)
  require(dplyr)

  #Check inputs
  if(is.null(features.lr)){
    message("Need genes to plot (features.lr)...")
    return(NULL)
  }
  if(is.null(features.ct)|length(features.ct)<2){
    message("Need 2 cell types to compare (features.ct)...")
    return(NULL)
  }

  # build df
  df <- data.frame(
    cell=Cells(SEU),
    GetAssayData(SEU, assay=assay.ge,slot=slot.ge)[features.lr[1:2],]%>%t(),
    GetAssayData(SEU, assay=assay.ct)[features.ct[1:2],]%>%t()
  )
  colnames(df)<-c("cell",features.lr[1:2],features.ct[1:2])

  df$cooc <- rep("L-/R-",nrow(df))

  df$cooc[df[features.ct[1]]>=bp.thresh & df[features.ct[2]]>=bp.thresh] <- "L+/R+" #both present
  df$cooc[df[features.ct[1]]>=bp.thresh & df[features.ct[2]]<bp.thresh] <- "L+/R-" #ligand cell type only
  df$cooc[df[features.ct[1]]<bp.thresh & df[features.ct[2]]>=bp.thresh] <- "L-/R+" # receptor cell type only

  if(verbose){
    print(table(df$cooc))
  }

  #TODO- scale gene expression data
  # if(scale.data){
  #   df[features.lr[1]] <- scale(x = df[features.lr[1]])
  #   df[features.lr[2]] <- scale(x = df[features.lr[2]])
  #
  #   df[features.lr[1]] <- MinMax(data = df[features.lr[1]], min = scale.min, max = scale.max)
  #   df[features.lr[2]] <- MinMax(data = df[features.lr[2]], min = scale.min, max = scale.max)
  # }

  # melt df for ggplot
  df <- reshape2::melt(
    df,
    measure.vars=features.ct,
    variable.name="cell_type",
    value.name = "theta"
  ) %>% reshape2::melt(
    measure.vars=features.lr[1:2],
    variable.name="gene",
    value.name = "log_norm_expr"
  )

  # plot!
  out.plot <-ggplot(
    df,
    aes(
      x=cooc,
      y=log_norm_expr,
      fill=gene
      # group=cooc
    )
  )+
    .geom_split_violin(
      width=width,
      show.legend = T
    )+
    labs(title = paste0(features.ct[1]," -> ", features.ct[2]))+
    vln.theme+
    theme(
      legend.position = "right",
      plot.title = element_text(color="black", face="bold")
    )

    return(out.plot)
}

#############################################################################
## Differential Gene Expression Plots ---------------------------------------

# Build volcano plot after running FindMarkers()
#   *Note- requires Seurat v4 or greater (b/c colnames in the output file changed for some reason...)
ggVolcano <- function(
  markers=NULL,
  expression=NULL,
  seu=NULL,
  logFC_filter = 1,
  neg.log.pval.Thresh=50,
  pct.thresh = 0.4,
  plotTitle='Volcano Plot',
  pseudocount=10^-300,
  fill.cols = c("#000000","#cf4c59","#2c30a8"),
  xlim=NULL, ylim=NULL,
  genes=NULL,
  gene.text.size=6, 
  repel=T, 
  nudge_x=-0.1,
  dot.scale=1,
  line.width=0.5, 
  segment.color="gray",
  pt.size=1,
  pt.alpha=1
){
  # TODO: add split.by - generalize?

  require(ggrepel)

  if(is.null(rownames(markers))){
    stop("Differential gene expression matrix needs gene names as rownames()...")
  }
  if(plotTitle=='Volcano Plot'){
    cat('Warning: no plot title given...\n')
  }

  # Build data.frame ####
  df <- data.frame(
    p_val_adj = -log10(markers$p_val_adj+pseudocount),
    avg_log2FC = markers$avg_log2FC,
    genes = rownames(markers),
    p_val_filter = -log10(markers$p_val_adj+pseudocount)>neg.log.pval.Thresh, # Values that pass p value threshold
    FC_filter.low = markers$avg_log2FC < -1*logFC_filter,
    FC_filter.high= markers$avg_log2FC > logFC_filter
  )

  if(is.null(seu)){
    df$pct =
      apply(markers,1,
            function(X){
              if(as.numeric(X["avg_log2FC"])>0){
                (return(as.numeric(X["pct.1"])))
              }else{
                return(as.numeric(X["pct.2"]))
              }
            }
      )

    df$expr.filter = df$pct > pct.thresh # which genes to label
  }else{
    #TODO: add in gene expression data
  }

  # Add colors ####
  df$cols <- rep(x='none', times=length(df$genes))
  df$cols[df$p_val_filter & df$FC_filter.high] <- 'high'
  df$cols[df$p_val_filter & df$FC_filter.low] <- 'low'
  df$cols <- factor(df$cols,levels=c("none", "high", "low"))

  message(table(df$colors!="none" & df$expr.filter)[2], "genes drawn \n")
  # Build ggplot object ####
  out.gg <- ggplot(
    data = df,
    aes(x=avg_log2FC, y=p_val_adj)
  ) +
    geom_hline(# pval dotted line
      yintercept = neg.log.pval.Thresh,
      linetype='dashed',
      color='lightgray', 
      size = line.width
    )  +
    geom_vline(# logFC dotted line
      xintercept = c(logFC_filter, -logFC_filter),
      linetype='dashed',
      color='lightgray', 
      size = line.width
    ) +
    geom_vline(
      xintercept = 0, 
      color='black',
      size = line.width
    ) +
    geom_point(
      pch=21, size=pt.size, alpha=pt.alpha,
      aes(
        fill = cols, col=cols
        # size = pct # scale dot size according to relative pct expression
      )
    )
  if(repel){
    out.gg <- out.gg +
      geom_text_repel(
        data=df[df$cols!="none" & df$expr.filter,],
        size=gene.text.size*(5/14), #convert to same scale as normal ggplot
        segment.size=0.25,
        segment.alpha = 0.8,
        segment.color = segment.color,
        point.padding = 0.4,
        aes(
          label=genes,
          col=cols
        )
      )
  }else{
    out.gg <- out.gg + 
      geom_text(
        data=df[df$cols=="both" & df$expr.filter,],
        size=gene.text.size*(5/14), #convert to same scale as normal ggplot
        # position=position_dodge(width = 1),
        nudge_x = nudge_x,
        aes(
          label=genes,
          col=cols
        )
    )
  }
  # Finish plot ####
  out.gg <- out.gg +
    labs(
      x="log2_FC",
      y="-log10(adj._p_val)",
      size = "Pct. Expr."
    ) +
    ggtitle(plotTitle) +
    guides(color=FALSE, fill=FALSE) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),#element_line(colour = "black", size = line.width),
      axis.ticks=element_line(colour = "black", size = line.width),
      panel.border = element_rect(color='black', size=line.width, fill=NA),
      plot.title = element_text(face="bold", hjust=0.5)
    ) +
    scale_fill_manual(values=fill.cols,aesthetics = c("col","fill"))+
    scale_x_continuous(expand=c(0.01,0.01))+
    scale_y_continuous(expand=c(0.01,0.01))

  return(out.gg)
}
