#' Timeseries
#'
#' This function allows you plot relative abundance reads over time in a timeseries plot.
#' @param data (required) A phyloseq object.
#' @param time (required) The name of the column containing the time variables, e.g. "Date".
#' @param group A variable from the associated sample data to group samples by (default: "Samples").
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: "Phylum")
#' @param tax.add Additional taxonomic levels to display for each entry, e.g. "Phylum" (default: "none").
#' @param tax.class Converts a specific phyla to class level instead, e.g. "p__Proteobacteria" (default: "none) .
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param plot.x.label Label on x-axis (default: "Time").
#' @param plot.y.label Label on y-axis (default: "Abundance").
#' @param plot.legend.title Label on ggplot legend. (default: "Sample").
#' @param plot.layout Either "ggplot" or "dygraph" (default: "dygraph").
#' @param plot.theme Chose different standard layouts choose from "normal" or "clean" (default: "normal").
#' @keywords timeseries
#' 
#'
#' @import data.table
#' @import xts
#' @import dygraphs
#' @import broom
#' @import ggplot2
#' @export 
#' @author Julie Klessner Thun Pedersen \email{julieklessnerthun@@gmail.com}

timeseries <- function(data, time, tax.aggregate="Phylum", tax.add=NULL, tax.class=NULL, tax.empty="best", plot.x.label="Time", plot.y.label="Abundance", plot.legend.title="Sample", plot.theme="normal", plot.layout="dygraph"){
  
  
  # Extract data from phyloseq object ------------------------------------------------
  
  data0 <- list(abund = as.data.frame(otu_table(data)@.Data),
                tax = data.frame(tax_table(data)@.Data, OTU =   rownames(tax_table(data))),
                sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  
  # Clean and rename taxonomy ---------------------------------------------------------
  
  data <- amp_rename(data = data0,
                     tax.class=tax.class, 
                     tax.empty=tax.empty, 
                     tax.level = tax.aggregate)
  
  # Divide data to seperate data frames ----------------------------------------------
  
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["sample"]]
  
  #Display multiple levels using tax.add argument ------------------------------------
  
  suppressWarnings(
    if (!is.null(tax.add)){
      
      if (tax.add != tax.aggregate){
        tax <- data.frame(tax, 
                          Display = apply(tax[,c(tax.add,tax.aggregate)], 1, 
                                          paste, collapse="; "))
      }
      
    } else {
      tax <- data.frame(tax, Display = tax[, tax.aggregate])
    }
  )  
  
  # Aggregate to a specific taxonomic level using tax.aggregate argument--------------
  
  abund3 <- cbind(Display = tax[,"Display"], abund) %>%
    melt(id.var = "Display", 
         value.name = "Abundance", 
         variable.name = "Sample")
  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  # Reshaping data for xts ---------------------------------------------------------
  
  keep <- c("Display","Sample", "sum")
  abund4 <- abund3[keep]
  
  abund5 <- reshape(abund4, timevar = "Sample", idvar = "Display",
                    direction = "wide")
  
  rownames(abund5) <- abund5[,1]
  abund5[,1] <- NULL
  
  # Add time variables to abundance reads ----------------------------------------
  
  time1 <- sample[, time] %>% as.Date()
  xts1 <- xts(x=t(abund5), order.by = time1, unique = T)
  
  xts2 <- apply.daily(xts1, FUN=mean)
  
  
  if(plot.layout == "dygraph"){
    q <- dygraph(xts2) %>%
      dyOptions(drawPoint = T, pointSize=2) %>%
      dyAxis("y", label = plot.y.label) %>%
      dyAxis("x", label = plot.x.label) %>%
      dyLegend(labelsSeparateLines = TRUE)
  }
  
  show(q)
  
  if(plot.layout == "ggplot"){
    
    # Reshaping for ggplot  --------------------------------------------
    
    ggxts <- tidy(xts2)  
    
    # Base plot, user may choose labels ------------------------------------------------
    
    p <- ggplot(ggxts, aes_string(x="index", y="value", col = "series"))+
      geom_line()+
      geom_point()+
      labs(x=plot.x.label, y=plot.y.label)+
      scale_colour_discrete(name=plot.legend.title)+
      theme(axis.text.x = element_text(size = 10, hjust = 1)) + 
      theme(axis.text.y = element_text(size = 12))
    
    # Setting layout if plot.theme = clean----------------------------------------------
    
    if(plot.theme == "clean"){
      
      p <- p + theme(legend.position = "none",
                     axis.text.y = element_text(size = 8, color = "black"),
                     axis.text.x = element_text(size = 8, color = "black", vjust = 0.5),
                     axis.title = element_blank(),
                     text = element_text(size = 8, color = "black"),
                     axis.ticks.length = unit(1, "mm"),
                     plot.margin = unit(c(0,0,0,0), "mm"),
                     title = element_text(size = 8))
    }
    # Output----------------------------------------------------------------------------
    p
    
  }

}
  