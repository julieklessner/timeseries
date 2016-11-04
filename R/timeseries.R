#' Timeseries
#'
#' This function allows you plot relative abundance reads over time in a timeseries plot.
#' @param data A phyloseq object.
#' @param time The name of the column containing the time variables (e.g. "Date")
#' @keywords timeseries
#' @export

timeseries <- function (data, time){
   
  ##Pulling data from phyloseq
  abund <- otu_table(data)@.Data %>% t() %>% as.data.frame()
  tax <- data.frame(tax_table(data)@.Data, OTU = rownames(tax_table(data)))
  sample <- suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  
  ##Adding time variables to abundance values and melting
  timeVar <- time
  time1 <- sample[, timeVar] %>% as.Date() 
  abund1 <- cbind(Time=time1, abund)
  abund2 <- melt(abund1, id.var="Time", value.name="Abundance", variable.name = "Sample")
  
  ##Plot  
  ggplot(abund2, aes_string(x="Time", y="Abundance", col="Sample"))+
    geom_line() 
}