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
#' @keywords timeseries
#' @export

timeseries <- function(data, time, group="Sample", tax.aggregate="Phylum", tax.add=NULL, tax.class=NULL, tax.empty="best"){
  
  ## Pulling data from phyloseq 
  data <- list(abund = as.data.frame(otu_table(data)@.Data),
               tax = data.frame(tax_table(data)@.Data, OTU = rownames(tax_table(data))),
               sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  
  ##Cleaning and renaming taxonomy
  data <- amp_rename(data = data, tax.class=tax.class, tax.empty=tax.empty, tax.level = tax.aggregate)
  
  ##Dividing data to seperate data frames
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["sample"]]
  
  #Making it possible to display multiple levels using tax.add
  ##Display: contains value names
  
  suppressWarnings(
    if (!is.null(tax.add)){
      if (tax.add != tax.aggregate) {
        tax <- data.frame(tax, Display = apply(tax[,c(tax.add,tax.aggregate)], 1, paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[, tax.aggregate])
    }
  )  
  
  #Adding tax.aggregates to abundance reads (as value names)
  #cols <- tax[, "Display"]
  #colnames(abund) <- cols
  
  
  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display = tax[,"Display"], abund) %>%
    melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  #Adding time variables to abundance reads
  time1 <- sample[, time] %>% as.Date()
  abund1 <- cbind(Time=time1, abund3) ### Brug merge
  
  ##Adding group information to abundance reads
  suppressWarnings(
    if (group != "Sample"){
      if (length(group) > 1){
        abund3 <- data.frame(abund1, Group = apply(sample[,group], 1, paste, collapse = " ")) 
      } else{
        abund3 <- data.frame(abund1, Group = sample[,group]) 
      }
    } else{
      abund3 <- data.frame(abund1, Group = abund1$Sample )
    }
  )
  
  ##Plot  
  ggplot(abund3, aes_string(x="Time", y="sum", col = "Display"))+
    geom_line()+
    theme(axis.text.x = element_text(size = 10, hjust = 1)) + 
    theme(axis.text.y = element_text(size = 12)) 
}