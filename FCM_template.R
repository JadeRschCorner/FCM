---
# title: "Analysis of Flow Cytometry"
# Description: Use Phenoflow package adapted from http://rprops.github.io/PhenoFlow/ to analyze microbial FCM data
# author: "Natchaya Luangphairin"
# date last revised: "5/18/2023"
# output: R Script
---

# install and load pacman if not already installed
install.packages("pacman")
library(pacman)
# This will install the libraries if they are not already installed, and then load them into your R session.
p_load(devtools,Phenoflow,flowviz,ggplot2,flowAI) 
install_github("rprops/Phenoflow_package") 

# set a fix seed to ensure reproducible analysis
set.seed(777)

  ###########################################################
  ################### 1. Load data ###################
  ###########################################################
    data(flowData)
    #path = "test_data"
    #flowData <- read.flowSet(path = path, pattern=".fcs")

    # quick look at the dimensions of our flowData
    attributes(flowData)

  ########################################################
  ###################  2. Denoise data ###################
  ########################################################
    # Select phenotypic features of interest and transform parameters 
    # In this case we chose two fluorescent parameters and two scatter parameters in their height format (-H)
    flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                       `SSC-H`=asinh(`SSC-H`), 
                                       `FL3-H`=asinh(`FL3-H`), 
                                       `FSC-H`=asinh(`FSC-H`))
    param=c("FL1-H", "FL3-H","SSC-H","FSC-H")


    ### Create a PolygonGate for denoising the dataset
    ### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
    sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3),ncol=2, nrow=4)
    colnames(sqrcut1) <- c("FL1-H","FL3-H")
    polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

    ###  Gating quality check
    xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
           scales=list(y=list(limits=c(0,14)),
                       x=list(limits=c(6,16))),
           axis = axis.default, nbin=125, 
           par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)


    ### Isolate only the cellular information based on the polyGate1
    flowData_transformed <- Subset(flowData_transformed, polyGate1)

    ### Extract metadata from sample names
    metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData),"_"), rbind)))
    colnames(metadata) <- c("Cycle_nr", "Location", "day", "timepoint", "Staining", "Reactor_phase", "replicate")


    # the phenotypic intensity values of each cell are normalized to the [0,1] range based on the maximum FL1-H intensity value over the data set
    # Note: depending on the max value of the primary fluorescence parameter (FL1-H in this case), large differences in observed alpha/beta diversity may occur.
    #       If you observe these variations, please try increasing the grid resolution of the fingerprints from `128` to `256`. 
    summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
    maxval <- max(summary[,"FL1-H"]) #Replace with the column representing the green fluorescence channel (e.g. "FITC-H")
    mytrans <- function(x) x/maxval
    flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                      `FL3-H`=mytrans(`FL3-H`), 
                                      `SSC-H`=mytrans(`SSC-H`),
                                      `FSC-H`=mytrans(`FSC-H`))

  ###########################################################
  ###################  23. Fingerprinting ###################
  ###########################################################
    ### Randomly resample to the lowest sample size
    #flowData_transformed <- FCS_resample(flowData_transformed, replace=TRUE)
    ### Calculate fingerprint with bw = 0.01
    fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                       bw=0.01,normalize=function(x) x)

  #####################################################################
  ###################  4. Calculate alpha diversity ###################
  #####################################################################
    # R is the number of bootstraps. d is a rounding factor which is used to remove spurious density values from the dataset.
    ### Calculate Diversity from normalized fingerprint 
    Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)  

    # Add the argument plot=TRUE in case a quick plot of the D2 diversity values with their errors is desired. 
    # A bit more tidied up figure can be easily achieved as such:
    p1 <- ggplot(data = Diversity.fbasis, aes(x = as.numeric(as.character(metadata$day)), y = D2, color = metadata$Reactor_phase))+
      geom_point(size = 8, alpha = 0.7)+
      geom_line()+
      scale_color_manual(values = c("#a65628", "red", 
                                    "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
      theme_bw()+
      labs(color = "Reactor phase", y = "Phenotypic diversity (D2)", x = "Days")+
      geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
    print(p1)


    ### Important: Diversity_rf function ###
      # This function allows more accurate diversity and error estimation by bootstrapping from the initial FCS files, as well as denoising
      # Diversity assessment with cleaning
      Diversity.clean <- Diversity_rf(flowData_transformed, param = param, R = 3, R.b = 3,
      cleanFCS = TRUE)

      # Note that we only took 3 bootstraps here. It is recommended to increase R and R.b to 100.

      p2 <- ggplot(data = Diversity.clean, aes(x = as.numeric(as.character(metadata$day)), y = D2, color = metadata$Reactor_phase))+
        geom_point(size = 8, alpha = 0.7)+
        geom_line()+
        scale_color_manual(values = c("#a65628", "red", 
                                      "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
        theme_bw()+
        labs(color = "Reactor phase", y = "Phenotypic diversity (D2)", x = "Days")+
        geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
      print(p2)

      ### Export diversity estimates to .csv file in the chosen directory
      write.csv2(file="results.metrics.csv", Diversity.clean)


  ###################################################################
  ###################  5. Beta diversity analysis ###################
  ###################################################################
    # Beta-diversity assessment of fingerprint
    beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

    # Plot ordination
    plot_beta_fcm(beta.div, color = metadata$Reactor_phase, labels="Reactor phase") + 
      theme_bw() +
      geom_point(size = 8, alpha = 0.5)  


  ###############################################################
  ###################  6. Extract cell counts ###################
  ###############################################################
    ### Creating a rectangle gate for counting HNA and LNA cells
    rGate_HNA <- rectangleGate("FL1-H"=c(asinh(20000), 20)/maxval,"FL3-H"=c(0,20)/maxval, 
                               filterId = "HNA bacteria")
    ### Normalize total cell gate
    sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3)/maxval,ncol=2, nrow=4)
    colnames(sqrcut1) <- c("FL1-H","FL3-H")
    polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

    ### Check if rectangle gate is correct, if not, adjust rGate_HNA
    xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
           scales=list(y=list(limits=c(0,1)),
                       x=list(limits=c(0.4,1))),
           axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                              cex=2), smooth=FALSE)
    ### Extract the cell counts
    a <- flowCore::filter(flowData_transformed, rGate_HNA) 
    HNACount <- summary(a);HNACount <- toTable(HNACount)
    s <- flowCore::filter(flowData_transformed, polyGate1)
    TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

    ### Extract the volume
    vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$`$VOL`))/1000

    ### Store the data
    results_counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                                 Total.cells = TotalCount$true/vol, HNA.cells = HNACount$true/vol)



    ### Exporting cell counts to .csv file to working directory
    write.csv2(file="results.counts.csv", results_counts)

    ### Plot cell density
    ggplot(data = results_counts, aes(x = as.numeric(as.character(metadata$day)), y = Total.cells, color = metadata$Reactor_phase))+
      geom_point(size = 8, alpha = 0.9)+
      scale_color_manual(values = c("#a65628", "red", 
                                    "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
      theme_bw()+
      labs(color = "Reactor phase", y = "Total cell density (cells/µL)", x = "Days")  