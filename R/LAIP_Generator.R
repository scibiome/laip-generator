#############################################################
# Data generator to simulate LAIP-positive events.
# If desired the simulated events are integrated into an existing flowFrame.
#
# This script requires the blast data and data of the other used population as
# separate csv- or txt-files. These files can be created using the export
# function of gating software like FlowJo or Infinicyt.
#
# The blast data can be filtered using the LAIPx_constraints parameter.
#
# The simulated LAIP-positive events can be saved as csv-file or integrated in
# a given basic population. This basic population can be imported as csv-, txt-
# or fcs-file.
# If desired the generated data contain an additional parameter "LAIP", labeling
# the LAIP-positive events.
#
# Hint: the enriched data are reshuffled before exporting into FCS-file to avoid
# anomalies in later visualization with cen-se.
#
# LAIP1 has to be defined. LAIP2, LAIP3 and LAIP4 are optional.
#
# The used CD-marker must be present in the imported data!
# The name of the LAIP-marker has to match the corresponding column name used in
# the csv-/txt-file.
#
# Examples of known LAIPs:
#  CD34+ CD22+ CD20-
#  CD34+ CD22+ CD20+
#  CD34+ CD22- CD20+
#  CD34+ CD22+ CD20+/-
#  CD7+
#  CD33-
#  CD56+
#  CD45-
#  HLA_DR-
#  CD1+/CD5+
#
# Author: Sebastian Mueller
#
##############################################################

#' simulate_LAIP_events
#'
#' Simulates blast-events with given LAIP-expression
#'
#' @param LAIP matrix
#' @param blast_pop data.frame
#' @param pop_colnum number of columns to use in given blast_pop-data.frame
#' @param number_of_LAIPevents number of events to be simulated
#' @param LAIP_constraints additional constraints for given blasts. blast-events that do not fullfill
#' @param addNoise boolean
#' @param printScatterplot boolean
#'
#' @return data.frame
#' @importFrom Spectre make.colour.plot
#' @export
#'
#' @examples
simulate_LAIP_events <- function(LAIP, blast_pop, pop_colnum, number_of_LAIPevents, LAIP_constraints = NULL, addNoise=TRUE, printScatterplot=FALSE){

  if(is.null(LAIP_constraints) == FALSE){
    #Removing blast events that violate given LAIP_constraints
    del_events <- c()
    z <- 0
    i <- 1
    for(i in 1:nrow(LAIP_constraints)){
      #check lower boundary
      if(is.na(LAIP_constraints[i,2]) == FALSE){
        j <- 1
        for(j in 1:nrow(blast_pop)){
          if(blast_pop[j,LAIP_constraints[i,1]] < as.numeric(LAIP_constraints[i,2])){
            del_events <- append(del_events, j)
          }
        }
      }
      #check upper boundary
      if(is.na(LAIP_constraints[i,3]) == FALSE){
        j <- 1
        for(j in 1:nrow(blast_pop)){
          if(blast_pop[j,LAIP_constraints[i,1]] > as.numeric(LAIP_constraints[i,3])){
            del_events <- append(del_events, j)
          }
        }
      }

      print(paste("LAIP-constraint ", i, ": ", length(del_events) - z, " violations."))
      z <- length(del_events)
    }
    #remove redundant values
    del_events <- unique(del_events)
    #print(quantile(blast_pop[,"PE-Cy7-A"],c(0.01, 0.1, 0.4, 0.6, 0.8, 0.9, 0.99)))
    print(paste(length(del_events), "of", nrow(blast_pop),"blast events violate LAIP_constraints"))

    #remove all LAIP_constraints-violating blast-events from blast_pop
    blast_pop <- blast_pop[-del_events,]
  }

  #Sampling blast events
  if(nrow(blast_pop) < number_of_LAIPevents && addNoise == FALSE){
    stop("Cannot simulate LAIP-population larger than blast population without adding noise.")
  }

  if(addNoise){
    LAIP_data <- blast_pop[sample(1:nrow(blast_pop), number_of_LAIPevents, replace = TRUE), 1:pop_colnum]
    num_vector <- runif(n=number_of_LAIPevents, min=1-noiseIntensity, max=1+noiseIntensity)
    for(i in 1:number_of_LAIPevents){
      LAIP_data[i,] <- round(LAIP_data[i,] * num_vector[i], digits = 2)
    }
  }else{
    LAIP_data <- blast_pop[sample(1:nrow(blast_pop), number_of_LAIPevents), 1:pop_colnum]
  }

  for(i in 1:nrow(LAIP)){
    #LAIP_parameter_data <- fread(LAIP[i,2], select = LAIP[i,1])
    ref_pop <- fread(LAIP[i,2])

    #exported population data contain sometimes additional attributes like "simulated time"
    #these columns are truncated here
    ref_pop <- as.matrix(ref_pop[,1:pop_colnum])
    colnames(ref_pop) <- colnames(data)[1:pop_colnum]
    # if(addLAIPlabel){
    #   j <- colnum_data-1
    #   ref_pop <- as.matrix(ref_pop[,1:j])
    #   colnames(ref_pop) <- colnames(data)[1:j]
    # }else{
    #   ref_pop <- as.matrix(ref_pop[,1:colnum_data])
    #   colnames(ref_pop) <- colnames(data)[1:colnum_data]
    # }

    LAIP_parameter_data <- ref_pop[,LAIP[i,1]]
    density_function <- density(LAIP_parameter_data)

    #Replacing LAIP-specific parameter with simulated data
    LAIP_data[,LAIP[i,1]] <- approx(cumsum(density_function$y)/sum(density_function$y),
                                    density_function$x,
                                    runif(number_of_LAIPevents))$y

    if(printScatterplot){
      #Share of simulated LAIPpos-data within plot should be at least 5%.
      if(nrow(blast_pop) > number_of_LAIPevents*20){
        blast_limit <- number_of_LAIPevents*19
      }else{
        blast_limit <- nrow(blast_pop)
      }

      plotDir <- getDirFromFilepath(output_FileName)

      #Extract the SSC-H and the LAIP-relevant parameter from the blast and the simulated data
      plotting_data <- rbind(blast_pop[1:blast_limit,c("SSC-H",LAIP[i,1])],LAIP_data[,c("SSC-H",LAIP[i,1])])

      laipLabel <- seq(0, 0, length.out = nrow(plotting_data)-number_of_LAIPevents)
      laipLabel <- append(laipLabel,seq(1,1, length.out = number_of_LAIPevents))

      plotting_data <- cbind(plotting_data, laipLabel)

      #Perform log-transformation on LAIP-relevant parameter and renaming the column
      plotting_data[,LAIP[i,1]] <- sapply(plotting_data[,LAIP[i,1]],FUN=log_transform)
      laipColumn <- paste("log(",LAIP[i,1],")",sep="")
      colnames(plotting_data) <- c("SSC-H",laipColumn,"laipLabel")

      make.colour.plot(dat = as.data.table(plotting_data),
                       x.axis = "SSC-H",
                       y.axis = laipColumn,
                       #y.axis = "BV786-A",
                       col.type = "factor",
                       col.axis = "laipLabel",
                       title = paste("Blasts with abnormal", LAIP[i,1], "expression."),
                       plot.width = 20,
                       plot.height = 20,
                       legend.loc = "right",
                       save.to.disk = TRUE,
                       path = plotDir)
      }
  }

  return(LAIP_data)
}


#' generateFlowFrame
#'
#' Generates a flowFrame based on the given data.frame, parameters and description
#' Parameters and description are optional. Without given description the first
#' 7 columns must be non-color-parameter (time, SSC-H, etc.) data
#'
#' @importFrom flowCore flowFrame
#' @importFrom BiocGenerics combine
#' @importFrom Biobase AnnotatedDataFrame
#' @param data = data.frame containing the measured intensities. Rows correspond to cells, columns to the different measurement channels.
#' @param parameters = data.frame containing information about each column of the flowFrame. columns of parameters: name, description.
#' @param description = a list containing the meta data included in the FCS file.
#'
#' @return flowFrame-object
#' @export
#'
#' @examples generateFlowFrame(bm_data)
generateFlowFrame <- function(data, parameters=NULL, description=NULL){

  #flowFrame-objects contain 3 slots: exprs, parameters and description
  #exprs = Object of class matrix, rows = events, columns = cytometer parameter (Time, FSC-H, SSC-H, ...)
  #parameters = Object of class AnnotatedDataFrame (package "Biobase") with 4 slots
  #     varMetadata (data.frame): 1 column "labelDescription" with 5 rows
  #description = Object of class list

  #Build parameters slot

  varMetadata <- data.frame(labelDescription = c("Name of Parameter",
                                                 "Description of Parameter",
                                                 "Range of Parameter",
                                                 "Minimum Parameter Value after Transforamtion",
                                                 "Maximum Parameter Value after Transformation"))
  rownames(varMetadata) <- c("name", "desc", "range", "minRange", "maxRange")

  metadata <- NULL
  if(is.null(parameters)){
    metadata <- data.frame(name=colnames(data),desc=paste('column',colnames(data),'from dataset'))
  }else{
    for(i in 1:nrow(parameters)){
      metadata <- data.frame(name=colnames(data),desc=paste('column',colnames(data),'from dataset'))
    }
  }


  ## Create FCS file metadata - ranges, min, and max settings
  metadata$range <- apply(apply(data,2,range),2,diff)
  metadata$minRange <- apply(data,2,min)
  metadata$maxRange <- apply(data,2,max)

  # in order to create a flow frame, data for exprs needs to be delivered as matrix
  exprs <- as.matrix(data)
  parameters <- AnnotatedDataFrame(metadata, varMetadata)


  gen_FlowFrame <- new("flowFrame",exprs=exprs, parameters=parameters)

  return(gen_FlowFrame)
}

#' addIdParameter
#' Adds an event-unique Id-parameter to a given flowFrame
#' @param ff = flowFrame-object
#'
#' @return flowFrame-object
#' @importFrom BiocGenerics combine
#' @export
#'
#' @examples addIdParameter(ff_fcs)
addIdParameter <- function(ff){

  # Kopieren eines bestehenden Parameters des FlowFrames als Vorlage
  new_p <- parameters(ff)[1,]

  # Umbenennen des kopierten Parameters von $P1 in $Px mit x = Parameteranzahl + 1
  new_p_number <- as.integer(dim(ff)[2]+1)
  rownames(new_p) <- c(paste0("$P", new_p_number))

  # Merge der Parameter des FlowFrames mit dem neuen Id-Parameter
  allPars <- combine(parameters(ff), new_p)

  # Setzen des Parameternamens und -beschreibung des neuen Id-Parameters
  new_p_name <- "event_id"
  #new_p_desc <- "Unique event ID"
  new_p_desc <- NA  #Durch Setzen der Parameterbeschreibung als NA wird dieser als Nicht-Farbkanal betrachtet.
  allPars@data$name[new_p_number] <- new_p_name
  allPars@data$desc[new_p_number] <- new_p_desc

  # Erweitern der Datenmatrix des FlowFrames um die event_id-Spalte
  orig_col_names <- dimnames(ff@exprs)[[2]]
  num_events <- as.integer(dim(ff)[1])
  event_ids <- as.matrix(1:num_events, ncol=1)
  new_exprs <- cbind(ff@exprs, event_ids)
  new_par_col_name <- setNames(new_p_name,
                               paste0("$P",as.character(new_p_number),"N"))
  dimnames(new_exprs)[[2]] <- c(orig_col_names, new_par_col_name)

  # Erweitern der Keywords des FlowFrames um neuen event_id-Parameter
  new_kw <- ff@description
  new_kw["$PAR"] <- as.character(new_p_number)
  new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
  new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
  new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
  new_kw[paste0("$P",as.character(new_p_number),"R")] <- new_kw["$P1R"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- 1
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- num_events
  #new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- new_kw["flowCore_$P1Rmin"]
  #new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- new_kw["flowCore_$P1Rmax"]

  # Erstellen des neuen FlowFrame-Objektes
  new_ff <- new("flowFrame", exprs=new_exprs, parameters=allPars, description=new_kw)

  return(new_ff)
}

#' addLaipParameter
#' Adds the parameter isLAIP to the given flowFrame.
#' @param ff = flowFrame-object
#' @importFrom BiocGenerics combine
#'
#' @return
#' @export
#'
#' @examples
addLaipParameter <- function(ff){

  labelColumn <- rep(0, as.integer(dim(ff)[1]))
  # Cloning an existing parameter of the flow frame as template
  new_p <- parameters(ff)[1,]

  # Renaming the copied parameter from $P1 in $Px with x = parameter count + 1
  new_p_number <- as.integer(dim(ff)[2]+1)
  rownames(new_p) <- c(paste0("$P", new_p_number))

  # Merging the parameter of the flowFrame and the new isLAIP-parameter
  allPars <- combine(parameters(ff), new_p)

  # Set name and description of the new isLAIP-parameter
  new_p_name <- "isLAIP"
  #new_p_desc <- "Unique event ID"
  new_p_desc <- NA  #setting description-field to NA ==> isLAIP will be handled as non-color-parameter.
  allPars@data$name[new_p_number] <- new_p_name
  allPars@data$desc[new_p_number] <- new_p_desc

  # Erweitern der Datenmatrix des FlowFrames um die isLAIP-Spalte
  orig_col_names <- dimnames(ff@exprs)[[2]]
  #num_LAIP_pops <- as.integer(dim(ff)[1])
  num_LAIP_pops <- max(labelColumn)
  LAIP_pop_ids <- as.matrix(labelColumn, ncol=1)
  new_exprs <- cbind(ff@exprs, LAIP_pop_ids)
  new_par_col_name <- setNames(new_p_name,
                               paste0("$P",as.character(new_p_number),"N"))
  dimnames(new_exprs)[[2]] <- c(orig_col_names, new_par_col_name)

  # Adding keywords in FlowFrame for new isLAIP-parameter
  new_kw <- ff@description
  new_kw["$PAR"] <- as.character(new_p_number)
  new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
  new_kw[paste0("$P",as.character(new_p_number),"V")] <- "0"
  new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
  new_kw[paste0("$P",as.character(new_p_number),"R")] <- new_kw["$P1R"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- 1
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- 4 * laipLabel
  #new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- new_kw["flowCore_$P1Rmin"]
  #new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- new_kw["flowCore_$P1Rmax"]

  #Erstellen des neuen FlowFrame-Objektes
  new_ff <- new("flowFrame", exprs=new_exprs, parameters=allPars, description=new_kw)

  return(new_ff)
}


#' print_scatterplot
#' print_scatterplot is a function for debugging purpose
#' @param data
#' @param max_rows
#' @param x = x-axis
#' @param y = y-axis
#' @param log_transform_y
#' @param laipColumnName
#'
#' @return
#' @importFrom Spectre make.colour.plot
#' @export
#'
#' @examples
print_scatterplot <- function(data, max_rows=10000, x="SSC-H",y,log_transform_y = TRUE,laipColumnName=NULL){

  if(nrow(data) > max_rows){
    #Sampling events to be plotted
    data_sample <- data[sample(1:nrow(data), max_rows),]
  }else{
    data_sample <- data
    max_rows <- nrow(data)
  }

  plotDir <- getDirFromFilepath(output_FileName)

  #Extract the x and the y parameter from data_sample and the LAIP-label if given
  if(is.null(laipColumnName)){
    plotting_data <- data_sample[1:max_rows,c(x,y)]
  }else{
    plotting_data <- data_sample[1:max_rows,c(x,y,laipColumnName)]
  }

  yColumn <- y

  if(log_transform_y){
    #Perform log-transformation on y parameter and renaming the column
    #plotting_data[,y] <- sapply(plotting_data[,y],FUN=log_transform)
    plotting_data[,y] <- sapply(plotting_data[,y],FUN=log_transform_no_neg)
    yColumn <- paste("log(",y,")",sep="")

    if(is.null(laipColumnName)){
      colnames(plotting_data) <- c(x,yColumn)
    }else{
      colnames(plotting_data) <- c(x,yColumn,laipColumnName)
    }
  }

  if(is.null(laipColumnName)){
    make.colour.plot(dat = as.data.table(plotting_data),
                     x.axis = x,
                     y.axis = yColumn,
                     title = paste(x, "-", y, "-scatterplot."),
                     plot.width = 15,
                     plot.height = 12,
                     square = FALSE,
                     legend.loc = "right",
                     save.to.disk = TRUE,
                     path = plotDir)
  }else{
    make.colour.plot(dat = as.data.table(plotting_data),
                     x.axis = x,
                     y.axis = yColumn,
                     col.type = "factor",
                     col.axis = laipColumnName,
                     title = paste(x, "-", y, "-scatterplot with LAIPs."),
                     plot.width = 15,
                     plot.height = 12,
                     square = FALSE,
                     legend.loc = "right",
                     save.to.disk = TRUE,
                     path = plotDir)
  }
}

#' writeToFCS
#' writeToFCS replaces the data matrix of a given flowFrame with the given data matrix
#' and saves the flowFrame as FCS-file.
#' @param ff = flowFrame-object
#' @param data = matrix
#' @param filepath
#' @param reshuffle = boolean, if TRUE the rows of data are reshuffled
#'
#' @return
#' @importFrom flowCore write.FCS
#' @export
#'
#' @examples
writeToFCS <- function(ff, data, filepath, reshuffle = TRUE){

  #reshuffling data to avoid anomalies in cen-se visualization
  if(reshuffle){
    data <- data[sample(nrow(data)),]
  }

  ff@exprs <- data

  print(paste("Saving simulated data in ", filepath, sep=""))
  write.FCS(ff_fcs, filepath, what="numeric", delimiter = "|", endian="big")
}
