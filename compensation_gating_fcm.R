---
# title: "Analyze flow cytometry data in R"
# Description: "Script to set compensation and live/dead gating"
# author: Natchaya Luangphairin
# date last revised: "8/25/2024"
# output: R Script
---


# Overview ----------------------------------------------------------------
#1) Load .fcs file
#2) Compensation (unstained help set background and autofluorescence baseline, while single-stain controls help compensate for spillover)
#3) Cleaning (Before running the analysis, rename all files to add "_" to spaces); Quality Control
#4) Transformation 
#5) Gating (Important to identify the cell population: singlets, doublets, live/dead cells analysis)
#6) Visualize the results

#Samples: 
##Unstained: use as negative control for background and autofluorescence signal, 
##Single-stained (SYBR and PI): use for compensation
##Co-stained (SYBR/PI): Use for live/dead gating to differentiate live cells (SYBR+, PI-)

# Increase memory by using cache for larger files
knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE,
                      echo = TRUE)



# Load Required Packages --------------------------------------------------
if (!require("pacman")) install.packages(pacman)
devtools::install_github("DillonHammill/CytoExploreR")
pacman::p_load("tidyverse","knitr","ggplot2","remotes","BiocManager","devtools","kableExtra","dplyr")
#BiocManager must be downloaded first before the following can be downloaded:
pacman::p_load("flowCore","flowWorkspace","openCyto","flowAI","ggcyto","CytoML","CytoExploreR")
citation("CytoExploreR") # for citation

# Save FCS Files ----------------------------------------------------------
#Unstained and FMO Controls: Unstained, FITC channel (SGI), and PerCP channel (PI)
#Here we did not use co-stain for compensation because co-stain won't allow us to isolate spillover effects of each fluorophore individually. Co-stains will be used later to gate live/dead cells via activation samples
myfiles <- "gating_samples/compensation"
fs <- read.flowSet(path = myfiles, pattern = ".fcs", ignore.case = TRUE, truncate_max_range.names = TRUE)

fs #flowSet for compensation/gating controls 
sampleNames(fs) #find out sample names
#Create a vector with the new names
new_names <- c("Compensation-Unstained", "Compensation-FITC", "Compensation-PerCP")  # Modify this according to data. 
#Assign the new names to the flowSet
sampleNames(fs) <- new_names

spillover(fs[["Compensation-FITC"]]) #need to choose one, so get spillover info from first file
fs_comp <- compensate(fs, spillover(fs[["Compensation-FITC"]])$`$SPILLOVER`) #notice we compensate with compensation matrix from first fcs file. If want to compensate each file comparing to their own must use fsApply. 

# automatically go through all parameter and looks for cells that are outside the dynamic range and remove events with anomalies (bad events)
fs_comp_clean <- flow_auto_qc(fs_comp) # takes a while; remove bad anomalies "out-of-range" events

# Compensation of FCS Files
cyto_save(fs_comp_clean,
          save_as = "Compensation-Samples")

#gs_linear <- GatingSet(fs_comp_clean)

# Compensation of Fluorescent Spillover -----------------------------------
##Prepare compensation controls (includes: Unstained, Single Stain, Full Range Stain)
##Setup compensation controls and mark FL1 (FITC) channel as SYBR Green I, FL3 (PerCP) channel as Propidium Iodide
gs_linear <- cyto_setup("Compensation-Samples",
                 gatingTemplate = "GACCompensation-gatingTemplate.csv",
                 markers = "GACExperiment-Marker.csv",
                 details = "GACExperiment-Details.csv") 

# Data transformations
##It is recommended that the compensation controls be appropriately transformed before proceeding.
##Transform: SYBR Green I (FL1-A), Propidium Iodide (FL2-A), and any other fluorescence channels with a wide dynamic range or compensation applied.
##Do not transform: Forward scatter (FSC) and side scatter (SSC) unless there is a compelling reason to do so.
all_columns <- colnames(gs_linear) #Get the column names of the flowFrame
exclude_columns <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "Width", "Time") #Identify columns to exclude by name
include_columns <- setdiff(all_columns, exclude_columns) #Create a vector of column names to include (exclude the identified columns)

# Visualize each transformation -------------------------------------------
#https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html
gs_log <- cyto_transform(cyto_copy(gs_linear), 
                         channels = "SYBR Green I-A",
                         type = "log",
                         copy = TRUE)
gs_arc <- cyto_transform(cyto_copy(gs_linear), 
                         channels = "SYBR Green I-A",
                         type = "arcsinh",
                         copy = TRUE)
gs_biex <- cyto_transform(cyto_copy(gs_linear), 
                          channels = "SYBR Green I-A",
                          type = "biex",
                          copy = TRUE)
gs_logicle <- cyto_transform(cyto_copy(gs_linear), 
                             channels = "SYBR Green I-A",
                             type = "logicle",
                             copy = TRUE)

# Save the plot with the correct file path and extension
cyto_plot_save("Gating-plots/Transformations-4.png", 
               width = 10, 
               height = 10)

# Set up the grid for a 2x2 layout
cyto_plot_custom(c(2, 2))

# LOG
p_log <- cyto_plot(gs_log[["Compensation-FITC"]],
                   parent = "root",
                   channels = "FL1-A",
                   title = "Log Transformation",
                   density_fill = "deepskyblue")
# ARCSINH
p_arc <- cyto_plot(gs_arc[["Compensation-FITC"]],
                   parent = "root",
                   channels = "FL1-A",
                   title = "Arcsinh Transformation",
                   density_fill = "deeppink")
# BIEX
p_biex <- cyto_plot(gs_biex[["Compensation-FITC"]],
                    parent = "root",
                    channels = "FL1-A",
                    title = "Biexponential Transformation",
                    density_fill = "green")
# LOGICLE
p_logicle <- cyto_plot(gs_logicle[["Compensation-FITC"]],
                       parent = "root",
                       channels = "FL1-A",
                       title = "Logicle Transformation",
                       density_fill = "yellow")

# Finalize plot layout
cyto_plot_complete()


# Apply transformation ----------------------------------------------------
# Good transformation will have least distortion near zero, symmetrical, and can handle dynamic ranges of samples
# Good transformation allows for clear distinction of microbial clusters
## for GAC work, consider either logicle or arcsinh
cyto_plot_save("Gating-plots/Transformations-Final.png", 
               width = 10, 
               height = 10)
cyto_plot_profile(gs_arc[["Compensation-FITC"]],
                  parent = "root",
                  channels = cyto_fluor_channels(gs_arc),
                  header = "Arcsinh Transformers",
                  density_fill = "deeppink")
cyto_plot_complete()

# can see that arcsinh works best for water samples, agreed with literature
# Using cyto_transform automatically applies transformation to the gating set
gs_arc <- cyto_transform(cyto_copy(gs_linear), 
                     channels=include_columns, 
                     type = "arcsinh") 

gs_logicle <- cyto_transform(cyto_copy(gs_linear), 
                            channels=include_columns, 
                            type = "logicle") 

time_check <- autoplot(fs_comp_clean, x = "Time", y = "SSC-A", bin = 456) +
  coord_cartesian(ylim = c(0, 5000000)) +
  labs(x = "Time", y = "SSC-A") +
  theme_bw() +
  theme(legend.position = "none") 

ggsave("Gating-plots/compensation_time_qc_plot.png", 
       plot = time_check, 
       width = 12, 
       height = 4, 
       dpi = 300)


single_cell_check <- autoplot(fs_comp_clean, x = "FSC-A", y = "FSC-H", bin = 456) +
  coord_cartesian(xlim = c(0, 5000000) , ylim = c(0, 6000000)) +
  labs(x = "FSC-A", y = "FSC-H") +
  theme_bw() +
  theme(legend.position = "none") 

ggsave("Gating-plots/compensation_single_cell_check_plot.png", 
       plot = single_cell_check, 
       width = 12, 
       height = 4, 
       dpi = 300)


# Gate Singlet Cells (Single Cells Gate 98.90%)
cyto_gate_draw(gs_logicle,
               parent = "root",
               alias = "Single Cells",
               channels = c("FSC-A", "FSC-H"),
               axes_limits = "auto")


# Automated computation of Spillover Matrix
spill <- cyto_spillover_compute(gs_logicle,
                                parent = "Single Cells",
                                spillover = "Spillover-Matrix.csv",
                                axes_limits = "auto") #select final peak

spill

# Interactively Edit Spillover Matrices (
# For SGI, primary fluroscence for SGI is FL1 and spillover fluroescence is FL2 and FL3. For PI, primary is FL3 and spillover is FL2 and FL1.
# Open CytoExploreR spillover matrix editor
spill <- cyto_spillover_edit(gs_logicle,
                             parent = "Single Cells",
                             spillover = "Spillover-Matrix.csv")


# Visualize compensation
cyto_plot_compensation(gs_logicle,
                       parent = "Single Cells")

cyto_plot_compensation(gs_logicle,
                       parent = "Single Cells",
                       spillover = "Spillover-Matrix.csv",
                       compensate = TRUE)


### DONE WITH COMPENSATION CONTROLS ###


# Live/Dead Gate ----------------------------------------------------------
# Load multiple .fcs file 
myfiles <- "gating_samples/live_dead/rep1"
#if packages already loaded and used library() no need to type in flowCore:: and command every single time, can just type in command
fs <- read.flowSet(path = myfiles, pattern = ".fcs", ignore.case = TRUE)

fs # flowSet for compensation/gating controls 
sampleNames(fs) # find out sample names
# Create a vector with the new names
new_names <- c("Gating-GAC-Unstained", "Gating-GAC-SGI", "Gating-GAC-PI", "Gating-GAC-SGIPI")  # Modify this according to data. 
# Assign the new names to the flowSet
sampleNames(fs) <- new_names

spillover(fs[["Gating-GAC-SGI"]]) # need to choose one, so get spillover info from first file
fs_comp <- compensate(fs, spillover(fs[["Gating-GAC-SGI"]])$`$SPILLOVER`) # notice we compensate with compensation matrix from first fcs file. If want to compensate each file comparing to their own must use fsApply. 

# QC to check three properties of flow cytometry: 1) flow rate, 2) signal acquisition, 3) dynamic range.
fs_comp_clean <- flow_auto_qc(fs_comp) # takes a while; remove bad anomalies "out-of-range" events from data

# Compensation of FCS Files
cyto_save(fs_comp_clean,
          save_as = "GACLiveDead-Samples")

# Create GatingSet to add LiveDead gate
gs <- cyto_setup("GACLiveDead-Samples",
                 gatingTemplate = "GACLiveDead-gatingTemplate.csv",
                 markers = "GACLiveDead-Marker.csv",
                 details = "GACLiveDead-Details.csv")

# Apply compensation spillover
# Use cyto_gatingTemplate_apply() if want to re-using gates from before.  
# cyto_gatingTemplate_apply() useful when set multiple gate; must be done sequentially in order: unstained control, main cell, single cell, live/dead within single cell, buffer control)
gs_comp_clean <- cyto_compensate(gs,
                      spillover = "Spillover-Matrix.csv") 

# Data transformations
# Get the column names of the flowFrame
all_columns <- colnames(gs_comp_clean)
# Identify columns to exclude by name
exclude_columns <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "Width", "Time")
# Create a vector of column names to include (exclude the identified columns)
include_columns <- setdiff(all_columns, exclude_columns)

# Apply arcsinh transformation to fluorescent channels
#The arcsinh transformation is generally the best choice for flow cytometry data from water samples due to its ability to handle a wide range of fluorescence intensities, accommodate low signals, and deal with background noise and autofluorescence effectively.
gs_comp_clean_trans <- cyto_transform(cyto_copy(gs_comp_clean), channels=include_columns, type = "arcsinh")

cyto_gatingTemplate_active() #to check which gate file is being active

# Gate Live Cells
# set gate manually using method to draw and Negated gates
# Single Cells
cyto_gate_draw(gs_comp_clean_trans,
               parent = "root",
               alias = "Single Cells",
               channels = c("FSC-A", "FSC-H"),
               type = "polygon",
               axes_limits = "auto")

# Extract unstained and stain control 
NIL <- cyto_extract(gs_comp_clean_trans, "Single Cells")[["Gating-GAC-Unstained"]] # for overlay
SGI_only <- cyto_extract(gs_comp_clean_trans, "Single Cells")[["Gating-GAC-SGI"]] # for overlay
PI_only <- cyto_extract(gs_comp_clean_trans, "Single Cells")[["Gating-GAC-PI"]] # for overlay
SGI_PI_costain <- cyto_extract(gs_comp_clean_trans, "Single Cells")[["Gating-GAC-SGIPI"]] # for overlay

# Subset the GatingSet to only include the "Gating-GAC-SGI" sample
gs_sgi <- gs_comp_clean_trans["Gating-GAC-SGI"]
# Plot SYBR Green I channel for visualization
cyto_plot(gs_sgi, parent = "Single Cells", channels = "SYBR Green I-A")
# Draw a quadrant gate for SYBR Green I (FL1-A) vs PI (FL3-A). Use cyto_gate_remove() to remove gate
cyto_gate_draw(gs_sgi,
               parent = "Single Cells",
               alias = c("Intact Cells", "Damaged or Dead Cells"), #SGI+/PI- Intact Cells
               channels = c("SYBR Green I-A", "Propidium Iodide-A"),
               type = "polygon",
               negate = "TRUE")

# Explore Gated Samples
gs_pop_get_stats(gs_sgi)
cyto_plot_gating_tree(gs_sgi,
                      stat = "freq")


### DONE WITH LIVE/DEAD GATING ###

