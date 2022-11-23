#' @title Lasso-based variable selection
#' @import glmnet
#' @importFrom stats cor
#' @param data data of clustered results, do not contain cluster column
#' @param col scalar of which row number selected
#' @param reduction use n/log(n) dimension reduction
#' @return list contain relationship of each row
#' @export
get_interaction <- function(data, col, reduction = FALSE ){

  if (nrow(data)==2) {
    return_obj = list(ind.name = rownames(data)[col],
                      dep.name = rownames(data)[-col],
                      coefficient = cor(t(data))[1,2])

  } else{
    data <- t(data); name <- colnames(data)
    y = as.matrix(data[,col])
    x = as.matrix(data[,-col])

    n <- ncol(x)
    if (reduction == TRUE) {
      vec <- abs(apply(x, 2, cor, y))
      if (all(is.na(vec))) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
      } else{
        x = x[,order(vec, decreasing = T)[1:(n/log(n))]]
      }
    }

    if ( all(y==0) |  all(y==1) ) {
      return_obj = list(ind.name = name[col],
                        dep.name = NA,
                        coefficient = 0)
    } else{
      ridge_cv <- try(cv.glmnet(x = x, y = y,alpha = 0))
      if ('try-error' %in% class(ridge_cv)) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)

      } else{
        ridge_cv <- cv.glmnet(x = x, y = y, type.measure = "mse", nfolds = 10, alpha = 0)
        best_ridge_coef <- abs(as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1])

        fit <- cv.glmnet(x = x, y = y, alpha = 1, family = "gaussian", type.measure = "mse",
                         penalty.factor = 1/best_ridge_coef,
                         nfolds = 10, keep = TRUE, thresh=1e-10, maxit=1e6)
        lasso_coef <- coef(fit, s = fit$lambda.1se)
        return_obj = list(ind.name = name[col],
                          dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1][-1],
                          coefficient = lasso_coef@x[-1])
        if ( length(return_obj$dep.name)==0 ) {
          tmp = cor(x,y)
          return_obj$dep.name = rownames(tmp)[which.max(abs(tmp))]
          return_obj$coefficient = tmp[which.max(abs(tmp))]*1/3
        }

      }

    }


  }


  return(return_obj)
}
