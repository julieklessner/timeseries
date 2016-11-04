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
  
  ##Adding time to subset and melting
  time <- get_variable(sample, time) %>% as.Date() %>% as.data.frame()
  abund3 <- cbind(time, abund)
  names(abund3)[1] <-"Time"
  abund5 <- melt(abund3, id.var="Time", value.name="Abundance", variable.name = "Sample")
  
  ##Plot  
  ggplot(abund5, aes_string(x="Time", y="Abundance", col="Sample"))+
    geom_line()  
}