#Check if necessary packages are installed and install them if not:
CRANpackages <- c("remotes", "BiocManager", "purrr", "dplyr", "stringr", "data.table", "lubridate", "tidyr", "tidyselect", "tibble")
BiocPackages <- c("flowCore", "flowWorkspace", "openCyto")
install.packages(setdiff(CRANpackages, rownames(installed.packages())))
BiocManager::install(setdiff(BiocPackages, rownames(installed.packages())), lib=.libPaths()[1])
if(!"Luciernaga" %in% rownames(installed.packages())){remotes::install.github("DavidRach/Luciernaga")}
if(!"ggcyto" %in% rownames(installed.packages())){remotes::install.github("RGLab/ggcyto")}
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(purrr)
library(dplyr)
library(Luciernaga)

#Load helper functions for pulling metadata from the cytoSet object for the purpose of adding additional metadata to pData later
functionsList <- list.files("functions", pattern=".R", full.names=TRUE)
walk(.x=functionsList, .f=source)

#Create list of FCS file paths
storageLocation <- file.path("data")
FCS_list <- list.files(storageLocation, pattern=".fcs", full.names=TRUE)

#Combine FCS files into a gating set object
cytoSet <- load_cytoset_from_fcs(FCS_list, truncate_max_range=FALSE)
gatingSet <- GatingSet(cytoSet)

#Arcsinh transformation of the CyTOF parameters
FCS_parameters <- colnames(gatingSet)
ms_parameters <- FCS_parameters[!stringr::str_detect(FCS_parameters, "Time|Event_length")]
arcsinh <- flowjo_fasinh_trans(m=4, t=20000, a=0.7, length=512)
myArcsinhTransform <- transformerList(ms_parameters, arcsinh)
transform(gatingSet, myArcsinhTransform)

#Apply automated gates using min density
autogates_CSV <- list.files(storageLocation, pattern="Autogates.csv", full.names=TRUE)
autogates_data.frame <- data.table::fread(autogates_CSV)
gateTemplate <- gatingTemplate(autogates_data.frame)
gt_gating(gateTemplate, gatingSet)

#Create new metadata columns and add them to the gating set's pData
originalPdata <- pData(gatingSet)
dateColumn <- map(.x=c(1:length(FCS_list)), .f=getDate, data=cytoSet) |> as.character() |> lubridate::dmy()
beadsColumn <- map(.x=c(1:length(FCS_list)), .f=getBeads, data=cytoSet) |> as.character()
newPdata <- originalPdata |> mutate(Date=dateColumn) |> mutate(BeadsType=beadsColumn)
pData(gatingSet) <- newPdata

#Create stacked histogram plots by date
outputLocation <- file.path("reports")
if(!dir.exists(outputLocation)){dir.create(outputLocation, recursive = FALSE)}
BeadMFIHistograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di+", TheX="140Ce_EQ Bead", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="EQ_beads")
CenterHistograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di-", TheX="Center", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="Center")
WidthHistograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di-", TheX="Width", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="Width")
ResidualHistograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di-", TheX="Residual", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=storageLocation, filename="Residual")
OffsetHistograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di-", TheX="Offset", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="Offset")
AmplitudeHistograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di-", TheX="Amplitude", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="Amplitude")
CD45Histograms <- Utility_RidgePlots(gs=gatingSet, subset="Ce140Di-", TheX="89Y_CD45_EQ6", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="CD45")
CD19Histograms <- Utility_RidgePlots(gs=gatingSet, subset="142Nd_CD19+", TheX="142Nd_CD19", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="CD19")
CD3Histograms <- Utility_RidgePlots(gs=gatingSet, subset="142Nd_CD19-", TheX="154Sm_CD3", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="CD3")
CD4Histograms <- Utility_RidgePlots(gs=gatingSet, subset="154Sm_CD3+", TheX="174Yb_CD4", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="CD4")
CD8Histograms <- Utility_RidgePlots(gs=gatingSet, subset="154Sm_CD3+", TheX="162Dy_CD8", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="CD8")
HLADRHistograms <- Utility_RidgePlots(gs=gatingSet, subset="170Er_HLADR+", TheX="170Er_HLADR", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="HLADR")
CD16Histograms <- Utility_RidgePlots(gs=gatingSet, subset="154Sm_CD3-", TheX="145Nd_CD16", TheY="Date", TheFill="BeadsType", returntype="pdf", outpath=outputLocation, filename="CD16")

#Create table containing counts from each gate
freqData <- gs_pop_get_count_with_meta(gatingSet, "count")
freqData$Population <- basename(freqData$Population)
columnsToRemove <- c("BeadsType", "sampleName", "Parent", "ParentCount")
freqData <- freqData |> select(!tidyselect::any_of(columnsToRemove)) |> relocate("Date", .before="Population")
freqData <- freqData |> tidyr::pivot_wider(names_from=Population, values_from=Count)

#Create percent of parent columns and change column names to not include problematic symbols
freqData <- freqData |> mutate(`percentCD45`=100*`CD45+`/`Ce140Di-`) |> mutate(`percentCD3`=100*`154Sm_CD3+`/`CD45+`) |> mutate(`percentCD19`=100*`142Nd_CD19+`/`CD45+`) |> mutate(`percentCD3CD4`=100*`174Yb_CD4+`/`154Sm_CD3+`) |> mutate(`percentCD3CD8`=100*`162Dy_CD8+`/`154Sm_CD3+`) |> mutate(`percentCD7CD3pos`=100*`CD7+CD3+`/`CD45+`) |> mutate(`percentCD7CD3neg`=100*`CD7+CD3-`/`CD45+`) |> mutate(`percentCD16HLADR`=100*`HLADR+CD16+`/`154Sm_CD3-`) |> mutate(`percentCD11cHLADR`=100*`HLADR+CD11c+`/`154Sm_CD3-`) |> mutate(numCells=`Ce140Di-`) |> mutate(numBeads=`Ce140Di+`)

#Create longitudinal plots
numCellsPlot <- ggplot(freqData, aes(x=as.Date(Date), y=numCells)) + geom_point() + geom_line() + coord_cartesian(ylim=c(0,100000)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("# of cells")
numBeadsPlot <- ggplot(freqData, aes(x=as.Date(Date), y=numBeads)) + geom_point() + geom_line() + coord_cartesian(ylim=c(0,60000)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("# of beads")
CD45plot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD45)) + geom_point() + geom_line() + coord_cartesian(ylim=c(40,100)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD45+")
CD3plot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD3)) + geom_point() + geom_line() + coord_cartesian(ylim=c(60,80)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3+")
CD19plot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD19)) + geom_point() + geom_line() + coord_cartesian(ylim=c(0,16)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD19+")
CD4plot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD3CD4)) + geom_point() + geom_line() + coord_cartesian(ylim=c(60,80)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3+CD4+")
CD8plot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD3CD8)) + geom_point() + geom_line() + coord_cartesian(ylim=c(15,30)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3+CD8+")
CD7CD3posplot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD7CD3pos)) + geom_point() + geom_line() + coord_cartesian(ylim=c(50,80)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3+CD7+")
CD7CD3negplot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD7CD3neg)) + geom_point() + geom_line() + coord_cartesian(ylim=c(0,20)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3-CD7+")
CD16plot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD16HLADR)) + geom_point() + geom_line() + coord_cartesian(ylim=c(0,16)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3-CD19-CD16+HLADR+")
CD11cplot <- ggplot(freqData, aes(x=as.Date(Date), y=percentCD11cHLADR)) + geom_point() + geom_line() + coord_cartesian(ylim=c(0,60)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("%CD3-CD19-CD11c+HLADR+")

#Retrieve MFI data and generate longitudinal plots
beadsMFIcolumn <- map(.x=c(1:length(FCS_list)), .f=retrieveMFI, gs=gatingSet, subset="Ce140Di+", parameter="Ce140Di") |> as.numeric()
centerMFIcolumn <- map(.x=c(1:length(FCS_list)), .f=retrieveMFI, gs=gatingSet, subset="Ce140Di-", parameter="Center") |> as.numeric()
Ir191MFIcolumn <- map(.x=c(1:length(FCS_list)), .f=retrieveMFI, gs=gatingSet, subset="Ce140Di-", parameter="Ir191Di") |> as.numeric()
Ir193MFIcolumn <- map(.x=c(1:length(FCS_list)), .f=retrieveMFI, gs=gatingSet, subset="Ce140Di-", parameter="Ir193Di") |> as.numeric()
newPdata <- pData(gatingSet) |> mutate(beadsMFI=beadsMFIcolumn)|> mutate(centerMFI=centerMFIcolumn) |> mutate(Ir191MFI=Ir191MFIcolumn)|> mutate(Ir193MFI=Ir193MFIcolumn)
pData(gatingSet) <- newPdata
beadMFIplot <- ggplot(newPdata, aes(x=as.Date(Date), y=beadsMFI)) + geom_point() + geom_line() + coord_cartesian(ylim=c(2000,2400)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("Beads MFI")
centerMFIplot <- ggplot(newPdata, aes(x=as.Date(Date), y=centerMFI)) + geom_point() + geom_line() + coord_cartesian(ylim=c(420,480)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("Center MFI")
Ir191MFIplot <- ggplot(newPdata, aes(x=as.Date(Date), y=Ir191MFI)) + geom_point() + geom_line() + coord_cartesian(ylim=c(500,850)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("191Ir MFI")
Ir193MFIplot <- ggplot(newPdata, aes(x=as.Date(Date), y=Ir193MFI)) + geom_point() + geom_line() + coord_cartesian(ylim=c(1000,1500)) + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5)) + labs(x=NULL, y=NULL) + ggtitle("193Ir MFI")

#Create reports
LJplots <- (numCellsPlot | numBeadsPlot | beadMFIplot)/(centerMFIplot | Ir191MFIplot | Ir193MFIplot)/(CD45plot | CD7CD3posplot | CD7CD3negplot)/(CD3plot | CD4plot | CD8plot)/(CD19plot | CD16plot | CD11cplot)
ggsave(filename=file.path(outputLocation, "LJplots.pdf"), plot=LJplots, width=10, height=9, units="in")