.pc_compute_capability_metrics <- function(object, LSL, USL, target, sigma = 3, force_normal = FALSE) {
  # This function is a wrapper around .bpc_compute_capability_metrics
  # to handle the different class structure
  
  
  # Get the class of the object
  obj_class <- class(object)[2]  # pc_normal or pc_t
  
  # Create fit structure expected by .bpc_compute_capability_metrics
  fit <- list(stanfit = NULL)
  
  if (obj_class == "pc_normal") {
    # Handle normal distribution case
    mu_sigma <- object$fit$fit
    
    # Create a temporary stanfit-like object for .bpc_compute_capability_metrics
    fit$stanfit <- list()
    class(fit$stanfit) <- "stanfit"
    
    # Set up the extract method for the stanfit object
    fit$stanfit$extract <- function(pars) {
      if (identical(pars, c("mu", "sigma"))) {
        return(list(
          mu = mu_sigma$mu,
          sigma = mu_sigma$sigma
        ))
      }
    }
    
    # Set proper class for dispatching
    class(fit) <- c("list", "bpc_normal")
  } else if (obj_class == "pc_t") {
    # Handle t distribution case
    t_params <- object$fit$fit
    
    # Create a temporary stanfit-like object for .bpc_compute_capability_metrics
    fit$stanfit <- list()
    class(fit$stanfit) <- "stanfit"
    
    # Set up the extract method for the stanfit object
    fit$stanfit$extract <- function(pars) {
      if (identical(pars, c("mu", "scale", "nu"))) {
        return(list(
          mu = t_params$mu,
          scale = t_params$sigma,  # sigma in t_params maps to scale
          nu = t_params$nu
        ))
      } else if (identical(pars, c("mu", "sigma"))) {
        # Also handle mu, sigma extraction
        return(list(
          mu = t_params$mu,
          sigma = t_params$sigma * sqrt(t_params$nu / (t_params$nu - 2.0))
        ))
      }
    }
    
    # Set proper class for dispatching
    class(fit) <- c("list", "bpc_t")
  }
  
  # Delegate to .bpc_compute_capability_metrics
  metrics <- .bpc_compute_capability_metrics(fit, LSL = LSL, USL = USL, target = target, 
                                           sigma = sigma, force_normal = force_normal)
  
  return(metrics)
}