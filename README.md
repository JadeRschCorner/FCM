# Flow Cytometry for Water Quality Analysis
Script to analyze flow cytometry results for sonication, degradation, and other water quality study applied to drinking water analysis.

Note: Other methods for microbial counts include ATP kits and heterotrophic plate counts (HPC). ATP test does not directly measure bacteria, but higher ATP concentration indicates higher number of living cells. HPC can be labor intensive and prone to underestimation, since not all water bacteria grow on conventional agar plates.

## Flow Cytometry Functional Data Analysis
R packages and scripts:
- For cytometry analysis: CytoExploreR: Interactive Analysis of Cytometry Data [Dillon Hammill (2021)](https://dillonhammill.github.io/CytoExploreR/index.html) <br/>
- For analysing flow cytometry experiments with model based clustering: functional principal component analysis adapted from http://rprops.github.io/PhenoFlow/

Templates:
- To help water quality researchers standardize their analyses, the [Eawag template](https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/unpublished-pdfs/Accuri-WP-Assessing-Water-Quality.pdf) for the BD Accuri C6 is available at bdbiosciences.com

Methods used for WQ analysis adapted from:
- See the [BD Accuri White Paper (Gatza et al., 2013)](https://www.umces.edu/sites/default/files/accuri-wp-assessing-water-quality.pdf), Assessing Water Quality with the BD Accuri™ C6 Flow Cytometer

## Abbreviations commonly used in flow cytometry:
1.  **FSC**: Forward Scatter is related to the size of the cell. The larger the cell, the more light is scattered in the forward direction.
    -   **FSC.A**: Forward Scatter Area
    -   **FSC.H**: Forward Scatter Height

2.  **SSC**: Side Scatter reflects the complexity or granularity of the cell. For example, cells with internal structures, like granules, will scatter more light to the side.
    -   **SSC-A**: Side Scatter Area
    -   **SSC-H**: Side Scatter Height
      
3.  **FITC**: Fluorescein Isothiocyanate is a fluorescent dye used to label cells or components within cells. It's not explicitly listed in your abbreviations, but it's a common fluorochrome in flow cytometry.

4.  **APC**: Allophycocyanin is another fluorescent dye used in flow cytometry.

5.  **PE**: Phycoerythrin is a bright fluorescent dye commonly used in flow cytometry.

6.  **PerCP**: Peridinin Chlorophyll Protein Complex is a water soluble carotenoid pigment found in photosynthetic dinoflagellates.

7.  **Time**: This usually refers to the time parameter, indicating the duration of the cell's passage through the laser beam, which can be useful for identifying doublets or aggregates.
Each of these parameters provides different information about the cells being analyzed, allowing for detailed characterization of cell populations based on size, granularity, and the presence of specific markers or dyes.

8. Other abbreviations are likely referring to specific fluorescent markers or dyes used in the analysis. Each of these markers is designed to bind to a specific cell component or molecule and will emit light at a specific wavelength when excited by a laser. 

## Packages
```{r}
# Install and load packages and libraries
#install.packages("pacman")
library(pacman)
p_load("tidyverse","knitr","ggplot2","remotes","BiocManager","devtools","kableExtra","dplyr")
# BiocManager must be downloaded first before the following can be downloaded:
p_load("flowCore","flowWorkspace","openCyto","flowAI","ggcyto","CytoML","CytoExploreR","CytoExploreRData")
devtools::install_github("DillonHammill/CytoExploreR")
devtools::install_github("DillonHammill/CytoExploreRData")

# Increase memory by using cache for larger files
knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
```

## Load FlowSet for multiple .fcs files
```{r}
# set your working directory
setwd("C:/Users/nluan/Box Sync/FCM_data")
# load flowset
myfiles <- "./fcm_results"
fs <- read.flowSet(path = myfiles, pattern = ".fcs", truncate_max_range.names = T)
fs

# extract metada
pData(fs)[1:5,]
```

### Explore data and attributes
```{r}
# Data
pData(fs)

# Attributes
names(fs[[1]])
colnames(fs[[1]])
exprs(fs[[1]])
keyword(fs[[1]])

# Check parameters data
fs[[1]]@parameters@data

# Show us our channels and markers
cyto_channels(fs)
cyto_markers(fs)
cyto_fluor_channels(fs)
```

## Step 1: Compensation
In flow cytometry, target-specific signal is detected in the form of Fluorescence. Compensation must be performed to correct for spectral overlaps between channels (e.g. FITC emmission peak overlaps with PE emission peak). 
  
- Proper controls: *Unstained Controls* to take into account intrinsic 'autofluorescence of a given type of cell; *Single Stained Controls* one for each color (e.g. SYBR only, PI only)
  - Some common sources of Fluorescence for flow include: Fluorochromes (e.g. FITC), Fluorescent proteins (to report gene expression), Fluroescent dyes (such as DNA or RNA binding dyes e.g. SYBR Green I (DNA), SYBR Green II (RNA), Propidium Iodide, etc.
    
  - **SYBR® Green I (SYBR)** is an asymmetrical cyanine dye that preferentially stains doublestranded DNA to form a DNA:dye complex that, when excited (λmax = 497 nm) by the 488-nm blue laser of the BD Accuri C6, emits green (λmax = 520 nm) and red light that can be measured in the FL1 and FL3 detectors.
    - SYBR stains all cells and allows for analysis of "total bacteria" (live and dead).   
    - Visualized using 2-D FL1-A (emission filter 533/30) vs FL3-A (emission filter 670 LP) log-scale density plot.
    - Bacteria are gated on the FL1-A vs FL3-A plots using a lower limit of 2,000 on FL1 based on real samples[BD Accuri White Paper (Gatza et al., 2013)](https://www.umces.edu/sites/default/files/accuri-wp-assessing-water-quality.pdf)
      
  - Adding **Propidium iodide (PI)** to SYBR® Green I allows for analysis of intact (viable) cells. PI only stains damaged "dead" cells, thus bacterial cells co-stained with SYBRPI will shift out of bacterial cell gate.
    
  - Together, **SYBR® Green I and PI (SYBR/PI)** can optimally discriminate bacteria with disrupted vs intact membranes. In the presence of PI, the same gate that was used to determine the total bacterial cell concentration will now include only viable, intact bacteria.
  
```{r}
# Renaming markers (optional)
colnames(fs)[colnames(fs)=="FL1-A"] <- "SYBR Green I-A"
colnames(fs)[colnames(fs)=="PerCP-A"] <- "PI-H"

# Compensation
spillover(fs[[1]]) #use spilllover from the first .fcs file
fs_comp <- compensate(fs, spillover(fs[[1]])$`$SPILLOVER`)
# if want to compensate each file comparing to their own must use fsApply
```

## Step 2: Quality Control
Detect and remove anomalies by checking 1) flow rate, 2) signal acquisition, 3) dynamic range, then continue same process with cleaning and transformation
- Takes a while; remove bad anomalies "out-of-range" events
```{r}
# Quality Check of samples
flow_auto_qc(fs_comp)
```

## Step 3: Transformation
Transform data to improve visualization and separate negative and positive events into discrete populations.
- Optimal transformation for each parameter is data-dependent.
- Methods include:
  - Log transformation
  - Arcsinh transformation
  - Biexponential transformation
  - Logicle transformation
  - See [Data Transformation in CytoExploreR](https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html)
 ```{r}
# Transform data
# Get the column names of the flowFrame
all_columns <- colnames(fs_comp_clean)
# Identify columns to exclude by name
exclude_columns <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "Width", "Time")
# Create a vector of column names to include (exclude the identified columns)
include_columns <- setdiff(all_columns, exclude_columns)

fs_comp_clean_trans <- transform(fs_comp_clean, estimateLogicle(fs_comp_clean[[1]], include_columns))
```
   
### Visualize transformed data
Compare before and after transformation

```{r}
# Compare between before and after transformation
fs_comp_clean_trans[[1]]
autoplot(fs_comp_clean[[1]]) #before
autoplot(fs_comp_clean_trans[[1]]) #after
```
  
## Step 4: Gating
Gating set should only be done after data is transformed
- Can use autogate or manual gate

Full Gating Strategy:
1. Unstained Control: Analyze an unstained sample to determine the baseline autofluorescence and background noise.
2. FSC vs. SSC: Gate the **Main Cell** Population: Identify and gate the main cell population to exclude debris and large aggregates.
3. FSC-A vs. FSC-H or FSC-W: Gate **Single Cells**: Ensure analysis is focused on single cells by gating out doublets and aggregates.
4. SYBR Green I vs. PI: Distinguish **Live and Dead** Cells within the Single-Cell Gate: Use the gated single cells to set up quadrants for live and dead cells based on specific markers and     5. Buffer Control: Use a 1xPBS buffer **control** sample (or DI, MilliQ, Autoclaved feed, depending on what water matrix used) to identify and subtract any background cell populations present in the buffer.

### Auto Gating
```{r}
# Auto gating
gs <- GatingSet(fs_comp_clean_trans)

# Main cell gate: FSC vs SSC removes debris
fs_data <- gs_pop_get_data(gs)
noneDebris_gate <- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels = c("FSC-A","SSC-A")))
gs_pop_add(gs, noneDebris_gate, parent = "root", name = "noneDebris_gate")
recompute(gs)
autoplot(gs, x = "FSC-A", y= "SSC-A", "noneDebris_gate", bins = 256)

# Single cell gate: FSC-A vs. FSC-H Singlet gate removes doublets and aggregates
fs_data <- gs_pop_get_data(gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels = c("FSC-A","FSC-H")))
gs_pop_add(gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(gs)
autoplot(gs, x = "FSC-A", y = "FSC-H", "singlets", bins = 256)

# Removing gates; can only do one at a time; example below:
#gs_pop_remove(gs, "singlets")
```
### Manual Gating
```{r}
# Create empty gatingTemplate
cyto_gatingTemplate_create("Gating_Template.csv")
# Select new gatingTemplate as active
cyto_gatingTemplate_select("Gating_Template.csv")

cyto_gatingTemplate_active() #to check which gate file is being active

# Set gate manually using method to draw and Negated gates
cyto_gate_draw(gs,
               alias = c("Background Cells","Targeted Cells"),
               channels = c("FL1-H","FL3-H"),
               type = "polygon",
               negate = TRUE,
               gatingTemplate = "Gating_Template.csv")

#Adding newly constructed gate(s) to gatingTemplate.csv

# This next line is necessary if you are re-using gates from before; useful when set multiple gate; must be done sequentially in order: unstained control, main cell, single cell, live/dead within single cell, buffer control)
cyto_gatingTemplate_apply(gs, "Gating_Template.csv")
```

## Step 5: Data Analysis and Visualizing FCM results
```{r}
# Applying gate: when read new set and want to automatically apply the fixed predetermined gate
gs <- cyto_setup(path = "./Gating",
                 gatingTemplate = "Gating_Template.csv",
                 details = "Gating-Details.csv",
                 restrict = TRUE)

# Analysis and Statistics
gs_pop_get_stats(gs)
gs_pop_get_stats(gs, "my_gate", "percent")

# Get population stats for downstream analysis
ps <- gs_pop_get_count_with_meta(gs)
ps <- ps %>% mutate(percent_of_parent = Count/ParentCount)
```

For scripts used for specific study, see folder (In progress)
