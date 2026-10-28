predict.ppmlasso = function(object, ..., newdata, interactions = NA)
{
  if (any(lapply(newdata, class) == "factor"))
  {
    unpacknewdata = catConvert(newdata)
    newdata       = unpacknewdata$X
    cat.names     = setdiff(unique(unpacknewdata$cat.names), NA)
    use.form = as.character(object$formula)[2]
    for (i in 1:length(cat.names))
    {
      use.form = gsub(cat.names[i], paste(names(newdata)[which(unpacknewdata$cat.names == cat.names[i])], collapse = " + "), use.form)
    }
    object$formula = as.formula(paste("~", use.form))
  }
  var.0 = which(apply(newdata, 2, var) == 0)
  #if (length(var.0) > 0)
  #{
  #newdata[,var.0] = newdata[,var.0] + rnorm(dim(newdata)[1], 0, 1.e-8)
  #}
  mf    = model.frame(object$formula, data = newdata)
  mt    = attr(mf, "terms")
  X.des = if (!is.empty.model(mt)) 
    model.matrix(mt, mf)
  else matrix(, length(object$mu), 0L)
  X.var = X.des[,-1]
  if (is.null(object$s.means) == FALSE)
  {
    X.var = scale(X.var, center = object$s.means, scale = object$s.sds)
    X.des = cbind(1, X.var)
  }
  if (object$family == "area.inter")
  {
    if (any(is.na(interactions)) == TRUE)
    {
      if (is.null(object$s.means) == FALSE)
      {
        X.des = cbind(X.des, min(scale(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE)
      {
        X.des = cbind(X.des, 0)
      }
    }
    if (any(is.na(interactions)) == FALSE)
    {
      if (is.null(object$s.means) == FALSE)
      {
        X.des = cbind(X.des, scale(interactions, center = mean(object$pt.interactions), scale = sd(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE)
      {
        X.des = cbind(X.des, interactions)
      }
    }
  }
  pred.int = exp(as.matrix(X.des) %*% object$beta)
  return(pred.int)
}

setMethod("predict", "ppmlasso", predict.ppmlasso)

