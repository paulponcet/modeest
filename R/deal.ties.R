 
# #' @author 
# #' Adapted from the function written by Wolfgang Huber and Ligia Pedroso Bras 
# #' in package \pkg{genefilter}. 
# #' 
# #' @importFrom stats median
.deal.ties <-
function(ny,         # length of the data
         i,          # index
         tie.action, # action to be taken
         tie.limit,  # limit
         warn = FALSE)
{
  ## Deal with ties
  maxi <- max(i)
  mini <- min(i)
  if (maxi-mini > tie.limit * ny) {
    warning(paste("encountered a tie, and the difference between minimal and 
                   maximal value is > length('x') * 'tie.limit'",
                  "the distribution could be multimodal", sep="\n"))
  }
  
  ## Take the action specified in "tie.action"
  switch(tie.action,
         mean = mean(i),
         median = stats::median(i),
         max = maxi,
         min = mini,
         stop(sprintf("invalid value '%s' for argument 'tie.action'", tie.action)))
}
