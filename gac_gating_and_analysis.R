---
# title: "Analyze flow cytometry in R"
# Description: "Script to set analyze FCM data"
# author: Natchaya Luangphairin
# date last revised: "8/19/2024"
# output: R Script
---

# Analyze Actual Data 
#1) Load .fcs file
#2) Compensation (unstained help set background and autofluorescence baseline, while single-stain controls help compensate for spillover)
#3) Cleaning (Before running the analysis, rename all files to add "_" to spaces); Quality Control
#4) Transformation 
#5) Gating (Important to identify the cell population: singlets, doublets, live/dead cells analysis)
#6) Visualize the results

#https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Visualisations.html

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


# Load FCS files ----------------------------------------------------------
myfiles <- "raw_gating_samples/live_dead/gas_pbs_5min_sample3"
#if packages already loaded and used library() no need to type in flowCore:: and command every single time, can just type in command
fs <- read.flowSet(path = myfiles, pattern = ".fcs", ignore.case = TRUE)

fs # flowSet for compensation/gating controls 
sampleNames(fs) # find out sample names
# Create a vector with the new names
#new_names <- c("Gating-GAC-Unstained", "Gating-GAC-SGI", "Gating-GAC-PI", "Gating-GAC-SGIPI")  # Modify this according to data. 
# Assign the new names to the flowSet
#sampleNames(fs) <- new_names

# Compensate Spillover ----------------------------------------------------
# Initial spillover from instrument
spillover(fs[[1]]) # need to choose one, so get spillover info from first file
fs_comp <- compensate(fs, spillover(fs[[1]])$`$SPILLOVER`) # notice we compensate with compensation matrix from first fcs file. If want to compensate each file comparing to their own must use fsApply. 

# Quality Control ---------------------------------------------------------
# QC to check three properties of flow cytometry: 1) flow rate, 2) signal acquisition, 3) dynamic range.
fs_comp_clean <- flow_auto_qc(fs_comp)

# Compensation of FCS Files
cyto_save(fs_comp_clean,
          save_as = "GACStudy-Samples")

# Samples as GatingSet
gs <- cyto_setup("GACStudy-Samples") #IF want to create new gating template add gatingTemplate = "GACLiveDead-gatingTemplate.csv"

gs_comp <- cyto_compensate(gs,
                           spillover = "Spillover-Matrix.csv") #reuse from spillover set by compensation_gating_fcm.R

# Clean and Extract metadata ----------------------------------------------
# Naming scheme: YYYYMMDD_WaterMatrix_SonicationTime_Replicate_Dye.fcs
# Pattern to extract the data, accounting for potential naming schemes and put in separate columns

# Extract sample names once for efficiency
sample_names <- sampleNames(gs)

# Create a column for the date of the experiment
pData(gs)$date <- sub(".*([0-9]{8}).*", "\\1", sample_names)

# Create a column for reactor or unique sample names to group later
pData(gs)$reactor <- sub(".*\\d{8}_([a-zA-Z0-9]+).*", "\\1", sample_names)
pData(gs)$reactor[!grepl("\\d{8}_[a-zA-Z0-9]+", sample_names)] <- NA

# Create a column for matrix media (DI, MQ, 1xPBS, AutoclavedFeed)
pData(gs)$matrix <- NA
pData(gs)$matrix[grepl("di", sample_names, ignore.case = TRUE)] <- "DI"
pData(gs)$matrix[grepl("mq", sample_names, ignore.case = TRUE)] <- "MQ"
pData(gs)$matrix[grepl("pbs", sample_names, ignore.case = TRUE)] <- "1xPBS"
pData(gs)$matrix[grepl("feed", sample_names, ignore.case = TRUE)] <- "AutoclavedFeed"

# Create a column for replicate type (duplicate, triplicate)
pData(gs)$replicate <- "sample"
pData(gs)$replicate[grepl("_\\d+d|_d_|mq2|di2|pbs2|feed2", sample_names, ignore.case = TRUE)] <- "duplicate"
pData(gs)$replicate[grepl("_\\d+t|_t_|mq3|di3|pbs3|feed3", sample_names, ignore.case = TRUE)] <- "triplicate"

# Extract sonication times from the sample names
pData(gs)$sonication_time_min <- as.numeric(sub(".*_(\\d+).*", "\\1", sample_names))
pData(gs)$sonication_time_min[grepl("di|mq|pbs|feed", sample_names, ignore.case = TRUE)] <- 0
pData(gs)$sonication_time_min[grepl("gac\\d+", sample_names, ignore.case = TRUE)] <- 5

# Create a column for dye (SYBR, PI, both)
pData(gs)$dye <- NA
pData(gs)$dye[grepl("sybr", sample_names, ignore.case = TRUE)] <- "sybr"
pData(gs)$dye[grepl("pi", sample_names, ignore.case = TRUE)] <- "pi"
pData(gs)$dye[grepl("both", sample_names, ignore.case = TRUE)] <- "both"

# Extract any other string pattern that could be of interest
pData(gs)$other <- sub(".*?([a-zA-Z]{2,}).*", "\\1", sample_names)
pData(gs)$other[!grepl("[a-zA-Z]{2,}", sample_names)] <- NA

# Create a description column based on conditions
pData(gs)$description <- NA
pData(gs)$description[grepl("di|mq|pbs|feed", sample_names, ignore.case = TRUE)] <- "blank"
pData(gs)$description[grepl("gac\\d+", sample_names, ignore.case = TRUE)] <- "matrix_study"
pData(gs)$description[grepl(".*_(\\d+).*", sample_names, ignore.case = TRUE)] <- "sonication_study"

# View the updated pData, channels, markers
pData(gs)
colnames(gs[[1]]) #renaming (optional): colnames(gs)[colnames(gs)=="FL1-A"] <- "SYBR Green I-A"
cyto_channels(gs)
cyto_markers(gs)
cyto_fluor_channels(gs)

gs_comp_clean <- gs_comp

# Transformation ----------------------------------------------------------
##Do not transform: Forward scatter (FSC) and side scatter (SSC) unless there is a compelling reason to do so.
all_columns <- colnames(gs_comp_clean) #Get the column names of the flowFrame
exclude_columns <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "Width", "Time") #Identify columns to exclude by name
include_columns <- setdiff(all_columns, exclude_columns) #Create a vector of column names to include (exclude the identified columns)
gs_comp_clean_trans <- cyto_transform(cyto_copy(gs_comp_clean), channels=include_columns, type = "arcsinh")

#Apply live/dead gating template set by compensation_gating_fcm.R and edit if needed using cyto_gate_edit()
gs_gated <- cyto_gatingTemplate_apply(gs_comp_clean_trans, gatingTemplate = "GACLiveDead-gatingTemplate.csv")
cyto_gatingTemplate_active()#use cyto_gate_draw() to draw new gate

pData(gs_gated)


# Filter Data and Apply Live/Dead Gate ------------------------------------
# Extract the pData from the GatingSet
pdata <- pData(gs_gated)
colnames(pdata)
# Subset pData where dye is "both" and sample contains "replicate"
subset_indices <- which(pdata$reactor %in% c("s1") & pdata$replicate %in% c("sample"))
# Subset the GatingSet using the identified indices
gs_subset <- gs_gated[subset_indices]

# #Time Check
# cyto_plot(gs_subset,
#           parent = "root",
#           channels = c("Time", "SSC-A"),
#           title = "2-D Scatter Time Check",
#           axes_limits = "auto")
# 
# #FSC-A vs. SSC-A Loose Cell Gate to gate clusters and select targeted cells (optional)
# cyto_plot(gs_subset,
#           parent = "root",
#           channels = c("FSC-A", "FSC-H"),
#           title = "2-D Scatter Loose Cells",
#           axes_limits = "auto")
# 
# #Single cells
# cyto_plot(gs_subset,
#           parent = "root",
#           channels = c("FSC-A", "FSC-H"),
#           title = "2-D Scatter Single Cells",
#           axes_limits = "auto")



# Live/Dead ---------------------------------------------------------------
#If unstained control shows intact cells, means there's autofluorescence
#autofluoirescence can be caused by debris, cell metabolic cofactors, EPS, photoreophic organisms

# Save png image of gating scheme after plotting
plot_name <- paste0("Gating-plots/",pData(gs_subset)$reactor[1], "_", pData(gs_subset)$replicate[1], "_Gating-Scheme.png")
cyto_plot_save(plot_name,
               width = 10, #1349x1228 = 1349/300, 1228/300 = width=4.50, height=4.09
               height = 10,
               res = 300 #high resolution for publication
               )

  # 2D scatter plot
  #when save set at 1000 and maintain aspect ratio
NIL <- cyto_extract(gs_subset, "root")[[1]] # for overlay
SGI_only <- cyto_extract(gs_subset, "root")[[2]] # for overlay
PI_only <- cyto_extract(gs_subset, "root")[[3]] # for overlay
SGI_PI_costain <- cyto_extract(gs_subset, "root")[[4]] # for overlay


  cyto_plot(gs_subset,
            parent = "root", #select "root", "Intact Cells", or "Damaged or Dead Cells"
            alias = c("Intact Cells", "Damaged or Dead Cells"),
            #parent = "Intact Cells",
            #parent = "Damaged or Dead Cells",
            channels = c("SYBR Green I-A", "Propidium Iodide-A"), ##For 1D density distribution define only one channel
            #group_by = c("sonication_time", "dye"),
            negate = TRUE, # statistic for events outside gates
            #overlay = NIL,
            label_text_size = 0.8,
            xlim = c(0,10e6),
            ylim = c(0,10e6),
            gate_line_type = 1,
            gate_line_width = 2,
            gate_line_col = "red",
            gate_fill = "white",
            gate_fill_alpha = 0,
            contour_lines = 15,
            title = c("Unstained Control (Background/Autofluorescence)\nGAC_1xPBS_5minSonication",
                      "SYBR Green I Single-Stained (Total Cells)\nGAC_1xPBS_5minSonication",
                      "Propidium Iodide Single-Stained (Damaged or Dead Cells)\nGAC_1xPBS_5minSonication",
                      "SGI/PI Co-Stained (Intact/Damaged or Dead Cells)\nGAC_1xPBS_5minSonication"),
            #axes_limits = "auto",
            point_size = 1,
            point_col_alpha = 1)

# Finalize plot layout
cyto_plot_complete()


# In case need to edit gate
cyto_gatingTemplate_active() #to check which gate file is being active
cyto_gate_edit(gs_gated, 
               parent = "root",
               alias = "Intact Cells",
               gatingTemplate = "GACLiveDead-gatingTemplate.csv")


#Cluster
cyto_plot(gs_subset[[4]],
          parent = "Intact Cells",
          channels = c("FSC-A", "SSC-A"),
          title = "2-D Scatter Cluster",
          xlim = c(0,1e5),
          ylim = c(0,1e5),
          display = 0.05,
          point_size = 2,
          point_col_alpha = 1,
          contour_lines = 15)


# Explore Gated Samples
gs_pop_get_stats(gs_subset)
cyto_plot_gating_tree(gs_subset[[4]],
                      stat = "freq")



# get population stats for downstream analysis
ps <- gs_pop_get_count_with_meta(gs_gated)
ps <- ps %>% mutate(percent_of_parent = Count/ParentCount)

styled_table <- ps %>%
  select(sampleName, Population, Count, ParentCount, percent_of_parent) %>%
  head() %>%
  mutate(across(everything(), ~cell_spec(.x, background = "white", color = "black"))) %>%
  kable("html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F)

styled_table








# Visualize the results using ggcyto --------------------------------------
# Compare before and after transformation
# Filter the data frame to include only rows where matrix == "1xPBS", replicate == "sample", description == "sonication_study"
pData(gs)
filtered_fs <- pData(gs_gated) %>% 
  dplyr::filter(matrix == "1xPBS" & replicate == "sample" & description == "sonication_study")

fs_comp_clean_trans <- filtered_fs

autoplot(gs_comp_clean[[1]]) #before transformation
autoplot(gs_comp_clean_trans[[4]]) #after transformation
# plot FSC-A vs. SSC-A
autoplot(gs_comp_clean_trans, x="FSC-A", y="SSC-a", bins = 256) + facet_wrap(~sonication_time_min)
# plot FSC-A vs. FL1-H or other markers
autoplot(fs_comp_clean_trans, x="FSC-A", y="FL1-H", bins = 256) + facet_wrap(~sonication_time_min)
# plot Time vs. FSC-A to check for doublets or aggregates
autoplot(fs_comp_clean_trans, x="Time", y="FSC-A", bins = 256) + facet_wrap(~sonication_time_min)









# Code Dump ---------------------------------------------------------------
gs_sgipi <- gs_comp_clean_trans["Gating-GAC-SGIPI"]
sample_costained <- cyto_gatingTemplate_apply(gs_sgipi, gatingTemplate = "GACLiveDead-gatingTemplate.csv")

# Gating scheme
cyto_plot_gating_scheme(sample_costained,
                        back_gate = TRUE,
                        gate_track = TRUE)

gs_pop_get_stats(sample_costained)
cyto_plot_gating_tree(sample_costained,
                      stat = "freq")

GAC_controls_trans <- cyto_transform(cyto_copy(gs_comp_clean), channels=include_columns, type = "arcsinh")
GAC_controls <- cyto_gatingTemplate_apply(GAC_controls_trans, gatingTemplate = "GACLiveDead-gatingTemplate.csv")
gs_pop_get_stats(GAC_controls)
cyto_plot_gating_tree(GAC_controls,
                      stat = "freq")
cyto_plot_gating_scheme(GAC_controls,
                        back_gate = TRUE,
                        gate_track = TRUE)


# get population stats for downstream analysis
ps <- gs_pop_get_count_with_meta(GAC_controls)
ps <- ps %>% mutate(percent_of_parent = Count/ParentCount)

styled_table <- ps %>%
  select(sampleName, Population, Count, ParentCount, percent_of_parent) %>%
  head() %>%
  mutate(across(everything(), ~cell_spec(.x, background = "black", color = "white"))) %>%
  kable("html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F)

styled_table






# Gating set must be done after data is transformed
gs <- gs_gated

# Auto Gating -------------------------------------------------------------
## Cell gate ---------------------------------------------------------------
# FSC vs SSC removes debris
fs_data <- gs_pop_get_data(gs)
noneDebris_gate <- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels = c("FSC-A","SSC-A")))
gs_pop_add(gs, noneDebris_gate, parent = "root", name = "noneDebris_gate")
recompute(gs)
autoplot(gs, x = "FSC-A", y= "SSC-A", "noneDebris_gate", bins = 256)

# plot each set
for (i in 1:length(gs)) {
  # generate the plot for each sample
  p <- autoplot(gs[[i]], x = "FSC-A", y= "SSC-A", "noneDebris_gate", bins = 256)
  
  # print the plot
  print(p)
  
  # Optionally, you can save each plot to a file
  # ggsave(filename=paste("plot_",i,".png,sep=""),plot=p)
}


## Setup single gate -------------------------------------------------------
#cell gate
# FSC-A vs. FSC-H Singlet gate removes doublets and aggregates
fs_data <- gs_pop_get_data(gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels = c("FSC-A","FSC-H")))
gs_pop_add(gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(gs)
autoplot(gs, x = "FSC-A", y = "FSC-H", "singlets", bins = 256)

# plot each set
for (i in 1:length(gs)) {
  # generate the plot for each sample
  p <- autoplot(gs[[i]], x = "FSC-A", y= "FSC-", "singlets", bins = 256)
  
  # print the plot
  print(p)
  
  # Optionally, you can save each plot to a file
  # ggsave(filename=paste("plot_",i,".png,sep=""),plot=p)
}


