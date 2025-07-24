filter_plots <- function(plot_names, condition, microstate) {

  keep_names <- plot_names %>%
    keep(~ str_detect(.x, condition)) %>%
    keep(~ {
      if (str_detect(.x, "transition_probability_peaks_")) {
        str_detect(.x, paste0("transition_probability_peaks_", microstate, "_"))
      } else if (str_detect(.x, "transition_probability_")) {
        str_detect(.x, paste0("transition_probability_", microstate, "_"))
      } else {
        str_detect(.x, paste0("_", microstate, "$"))
      }
    })
  
  return(keep_names)
}