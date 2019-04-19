#' Random Sample Consensus
#'
#' @param formula the specification for the fit
#' @param data the data to fit the model with
#' @param min_to_fit the minimum number of datum to fit the initial linear model
#' @param max_iteration the maximum number of iterations 
#' @param fit_threshold the theshold to decide if a point is fit well with the
#'                      initial model - the method of well fit detection is MSE
#' @param data_close the number of data that fit well required before 
#'                   refitting a new model
#'
#' @return a robust linear regression model
#' @export
#' 
RANSAC <- function(formula, data, min_to_fit = 10, max_iteration = 10,
                   fit_threshold = 0.5, data_close = 10) {
  
  # initial model is NULL
  best_fit <- NULL
  
  # the initial best error is massive
  # this is the error to beat!
  best_error <- 1e6
  
  # extract dependent variable from the formula
  dependent_variable <- paste(strsplit(as.character(formula), split = "~"))[[2]]
  
  for (iteration in 1:max_iteration) {
    
    # sample from the data and fit an initial model
    possible_inliers <- sample(1:nrow(data), size = min_to_fit, replace = FALSE)
    initial_model <- lm(formula, data = data[possible_inliers, ])
    
    # data which will also be considered inliers if fit well with initial model
    also_inliers <- NULL
    
    for (data_point in (1:nrow(data))[-possible_inliers]) {
      
      # generate a prediction for each other observation and see if the fit
      # is any good with MSE
      # keep track of those which do fit well
      p <- predict(initial_model, newdata = data[data_point, ])
      if ((data[data_point, dependent_variable] - p) ^ 2 < fit_threshold) {
        
        also_inliers <- c(also_inliers, data_point)
        
      }
      
      # if we've saved enough data, fit a new model with these new saves
      # and the data we used to build the initial model
      if (length(also_inliers) > data_close) {
        
        better_model <- lm(formula, data[c(possible_inliers, also_inliers), ])
        
        # extract the residual standard error from the better model
        this_error <- summary(better_model)$sigma
        
        # check if the model improved
        if (this_error < best_error) {
          
          best_fit <- better_model
          best_error <- this_error
          
        }
        
      }
      
    }
    
  }
  
  # return the best fit to the user
  best_fit
  
}
       