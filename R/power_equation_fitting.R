#' @title remove observation with too many 0 values
#' @param data dataframe of imported dataset, must have first column as ID
#' @param x scales indicate how many 0 to remove
#' @return a dataframe without too many 0 observations
#' @examples
#' data_cleaning(matrix(c(c(0,1,1,0,0,1,1), c(2,1,0,3,5,2,2), c(1,1,3,2,4,5,1)), 3, 7), 2)
#' @importFrom stats aggregate
#' @export
data_cleaning <- function(data, x = round(ncol(data)*0.3)){
  data = aggregate(data[,2:ncol(data)], by=list(data[,1]), FUN = 'sum')
  rownames(data) = data[,1]
  data = data[,-1]
  tmp = apply(data, 1, function(c) which( as.numeric(c) != 0) )
  keep_no = which(sapply(tmp, length) >= x)
  data2 = data[keep_no,]
  return(data2)
}


#' @title match power_equation fit result for bi-variate model
#' @param result1 list object from power_equation fit
#' @param result2 list object from power_equation fit
#' @return a id match list for input dataset
#' @export
data_match <- function(result1, result2){
  matchname = intersect(rownames(result1$original_data),rownames(result2$original_data))

  new_result1 = list(original_data = result1$original_data[matchname,],
                     trans_data = result1$trans_data[matchname,],
                     power_par = result1$power_par[matchname,],
                     power_fit = result1$power_fit[matchname,],
                     Time = result1$Time)

  new_result2 = list(original_data = result2$original_data[matchname,],
                     trans_data = result2$trans_data[matchname,],
                     power_par = result2$power_par[matchname,],
                     power_fit = result2$power_fit[matchname,],
                     Time = result2$Time)
  result = list(dataset1 = new_result1, dataset2 = new_result2)
  return(result)
}

#' @title use power equation parameters to generate y values
#' @param x vector for x values
#' @param power_par matrix contain parameters for power equation
#' @return y values for given power equation parameters
#' @examples
#' power_equation(c(1,2,3,5,7), matrix(c(2,1,1,2),2,2))
#' @export
power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}


#' @title use power equation to fit observed values
#' @param x vector for x values
#' @param y vector for y valyes
#' @return nls model
#' @examples
#' power_equation_base(c(1,2,3,5,7), c(5,10,15,17,20))
#' @importFrom stats nls nls.control runif lm
#' @export
power_equation_base <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_value = min(y[y!=0])

  lmFit <- lm( log( y + runif(1, min = 0, max = min_value))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]

  model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                   control = nls.control(maxiter = 1e3, minFactor = 1e-200)))
  if( 'try-error' %in% class(model)) {
    result = NULL
  }
  else{
    result = model
  }
  return(result)
}


#' @title use power equation to fit observed values
#' @param x vector for x values
#' @param y vector for y values
#' @param maxit numeric value for maximum initial pars try
#' @return nls model
#' @examples
#' power_equation_all(c(1,2,3,5,7), c(5,10,15,17,20))
#' @export
power_equation_all <- function(x,y, maxit=1e2){
  result <- power_equation_base(x,y)
  iter <- 1
  while( is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x,y))
  }
  return(result)
}

#' @title use power equation to fit given dataset
#' @param data cleaned dataframe
#' @param n scales for how many interpolation needed
#' @param trans indicate log/log2/log10 transform dataset
#' @param thread scales for how many thread used
#' @return list contain power equation parameters and fitted data
#' @importFrom stats coef predict
#' @import parallel
#' @export
power_equation_fit <- function(data, n=30, trans = log10, thread = 2) {
  data = data[,order(colSums(data))]
  if ( is.null(trans)) {
    X = colSums(data)
    trans_data = data
  } else{
    X = trans(colSums(data+1))
    trans_data = trans(data+1)
  }
  colnames(trans_data) = X

  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("power_equation_all", "power_equation_base", "trans_data", "X"), envir = environment())
  all_model = parLapply(cl = cl, 1:nrow(data), function(c) power_equation_all(X, trans_data[c,]))
  stopCluster(cl)


  names(all_model) = rownames(data)
  no = which(sapply(all_model, length)>=1)
  all_model2 = all_model[no]
  data2 = data[no,]
  trans_data2 = trans_data[no,]

  new_x = seq(min(X), max(X), length = n)
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                       FUN.VALUE = numeric(n), USE.NAMES = TRUE))

  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2,
                power_par = power_par, power_fit = power_fit,
                Time = X)
  return(result)
}



#' @title plot power equation fitting results
#' @import ggplot2 scales
#' @importFrom reshape2 melt
#' @param result list object returned from power_equation_fit
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param n scales for how many subplots needed
#' @return plot show power curve fitting result
#' @export
power_equation_plot <- function(result, label = 10, n = 9){
  data1 = result[[2]]
  data2 = result[[4]]

  no = sample(1:nrow(data1),n)

  df_original =  reshape2::melt(as.matrix(data1[no,]))
  df_fit = reshape2::melt(as.matrix(data2[no,]))


  p <- ggplot() +
    geom_point(df_original, mapping = aes_string(x = "Var2", y = "value",colour = "Var1"),
               show.legend = F, alpha = 0.5, shape = 1) +
    geom_line(df_fit, mapping = aes_string(x = "Var2", y = "value",colour = "Var1"), size = 1.25, show.legend = F)  +
    facet_wrap(~Var1) +
    xlab("Habitat Index") + ylab("Niche Index") + theme(axis.title=element_text(size=18)) +
    theme_bw() + geom_text(df_fit, mapping = aes_string(label = "Var1"), show.legend = FALSE,
                           x = mean(df_fit$Var2), y = max(df_original$value)*0.9, check_overlap = TRUE, size = 4) +
    theme_bw() +
    theme(axis.title=element_text(size=15),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10,hjust = 0),
          panel.spacing = unit(0.0, "lines"),
          plot.margin = unit(c(1,1,1,1), "lines"),
          strip.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_blank())


  if (is.null(label)) {
    p = p
  } else {

    xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
    ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks

    xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
    if (ylabel[1] == 0) {
      ylabel2 = parse(text=c(0,paste(label,"^", ylabel[2:length(ylabel)], sep="")))
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p = p + scale_x_continuous(labels = xlabel2) + scale_y_continuous(labels = ylabel2)
  }
  return(p)
}


#' @title plot power equation fitting results for bi-variate model
#' @import ggplot2 scales
#' @importFrom reshape2 melt
#' @param result list object returned from data_match
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param n scales for how many subplots needed
#' @param show.legend show legend or not
#' @param color1 Hex Color Codes for first data
#' @param color2 Hex Color Codes for second data
#' @return plot show power curve fitting result
#' @export
bipower_equation_plot <- function(result, label = 10, n = 9, show.legend = FALSE,
                                  color1 = "#38E54D", color2 = "#FF8787"){
  data1 = result$dataset1[[2]]
  data2 = result$dataset1[[4]]
  no = sample(1:nrow(data1),n)
  df_original1 =  reshape2::melt(as.matrix(data1[no,]))
  df_original1$type = "L"
  df_fit1 = reshape2::melt(as.matrix(data2[no,]))
  df_fit1$type = "L"

  data1 = result$dataset2[[2]]
  data2 = result$dataset2[[4]]

  df_original2 =  reshape2::melt(as.matrix(data1[no,]))
  df_original2$type = "R"
  df_fit2 = reshape2::melt(as.matrix(data2[no,]))
  df_fit2$type = "R"

  df_fit = rbind(df_fit1,df_fit2)
  df_original = rbind(df_original1,df_original2)
  df_fit$Var1 = as.character(df_fit$Var1)
  df_original$Var1 = as.character(df_original$Var1)

  p <- ggplot() +
    geom_point(df_original, mapping = aes_string(x = "Var2", y = "value",
                                                 colour = factor(df_original$type)),
               show.legend = show.legend, alpha = 0.5, shape = 1) +
    geom_line(df_fit, mapping = aes_string(x = "Var2", y = "value",
                                           colour = factor(df_fit$type)),
              size = 1.25, show.legend = show.legend, shape = 1)  +
    facet_wrap(~Var1) + labs(color = "Type")  +
    xlab("Habitat Index") + ylab("Niche Index") + theme(axis.title=element_text(size=18)) +
    theme_bw() +
    geom_text(df_fit, mapping = aes_string(label = "Var1"), show.legend = show.legend,
              x = mean(df_fit$Var2), y = max(df_original$value)*0.9, check_overlap = TRUE, size = 4) +
    theme_bw() + scale_color_manual(values = c("L" = color1, "R" = color2)) +
    theme(axis.title=element_text(size=15),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10,hjust = 0),
          panel.spacing = unit(0.0, "lines"),
          plot.margin = unit(c(1,1,1,1), "lines"),
          strip.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_blank())

  if (is.null(label)) {
    p = p
  } else {
    xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
    ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks

    xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )

    if (ylabel[1] == 0) {
      ylabel2 = parse(text=c(0,paste(label,"^", ylabel[2:length(ylabel)], sep="")))
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p = p + scale_x_continuous(labels = xlabel2) + scale_y_continuous(labels = ylabel2)
  }
  return(p)
}
