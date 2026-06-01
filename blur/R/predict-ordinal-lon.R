predict_ordinal_lon <- function(object, predX, predPsi, counts, timepoints) {
  if(attr(object, "has_prev_covar")) {
    predict_with_prev_covar(...)
  } else {
    for(tt in 1:n_tp) {
      id_t <- which(timepoints == tt)
      preds_t <- predict_ordinal_cs(object,
                                    predX[id_t,],
                                    predPsi[id_t,],
                                    counts[id_t])
    }      
  }  
}
