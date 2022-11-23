#' @title Q-function to replace log-likelihood function
#' @importFrom mvtnorm dmvnorm
#' @param par numeric vector for parameters need to be estimated
#' @param prob_log mixture component weights(log)
#' @param omega_log latent variables(log)
#' @param X matrix for cluster
#' @param k vector for the cluster number
#' @param n1 scalar for number of column contain first trait/location etc
#' @param n2 scalar for number of column contain second trait/location etc
#' @param times1 vector for the x values or time points
#' @param times2 vector for the x values or time points
#' @importFrom mvtnorm dmvnorm
#' @return the Loglikelihood value
biQ_function <- function(par, prob_log, omega_log, X, k, n1, n2, times1, times2){
  n = dim(X)[1]; d = dim(X)[2];
  X1 = X[,1:n1]; X2 = X[,(n1+1):(n1+n2)]
  par.cov <- par[1:4]
  par.mu <- array(par[-c(1:4)], dim = c(k,2,2))

  cov1 = get_SAD1_covmatrix(par.cov[1:2], n1)
  cov2 = get_SAD1_covmatrix(par.cov[3:4], n2)
  mu1 <- power_equation(times1, par.mu[,,1][1:k,])
  mu2 <- power_equation(times2, par.mu[,,2][1:k,])
  mvn_log1 <- sapply(1:k, function(c) dmvnorm(X1, mu1[c,], cov1, log = T))
  mvn_log2 <- sapply(1:k, function(c) dmvnorm(X2, mu2[c,], cov2, log = T))
  mvn.log = mvn_log1 + mvn_log2
  tmp = sweep(mvn.log, 2, FUN = "+", STATS = prob_log) - omega_log
  Q = -sum(tmp*exp(omega_log))
  return(Q)
}

#' @title acquire initial parameters for functional clustering
#' @importFrom stats kmeans
#' @param X matrix for cluster
#' @param k vector for the cluster number
#' @param times1 vector for the x values or time points
#' @param times2 vector for the x values or time points
#' @param n1 scalar for number of column contain first trait/location etc
#' @param n2 scalar for number of column contain second trait/location etc
#' @importFrom stats cov kmeans IQR
#' @return the initial parameters for functional clustering
#' @export
biget_par_int <- function(X, k, times1, times2, n1, n2){
  n = dim(X)[1]; d = dim(X)[2];
  X1 = X[,1:n1]; X2 = X[,(n1+1):(n1+n2)]

  cov.int = c(0.9,IQR(diag(cov(X1))), 0.9,IQR(diag(cov(X2))))

  init.cluster <- kmeans(X,centers = k)
  prob <- table(init.cluster$cluster)/n

  fit1 <- lapply(1:k,function(c) power_equation_all(times1,init.cluster$centers[c,1:n1]))
  fit2 <- lapply(1:k,function(c) power_equation_all(times2,init.cluster$centers[c,(n1+1):(n1+n2)]))
  mu.par.int1 <- t(sapply(fit1, coef))
  mu.par.int2 <- t(sapply(fit2, coef))
  return_obj <- list(initial_cov_params = cov.int,
                     initial_mu_params = c(mu.par.int1,mu.par.int2),
                     initial_probibality = prob)
  return(return_obj)
}


#' @title main function for bifunctional clustering
#' @param data1 matrix or data for cluster
#' @param data2 matrix or data for cluster
#' @param k vector for the cluster number
#' @param Time1 vector for the time point
#' @param Time2 vector for the time point
#' @param trans indicate log/log2/log10 transform dataset
#' @param inv.cov matrix for directly solve cov matrix, default not given(currently not available)
#' @param initial.pars vector for manual give initial parameters, default not given
#' @param iter.max scales control iteration for EM algorithm
#' @param parscale scales control parameters scales for cov pars
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats optim
#' @return the initial parameters for functional clustering
#' @export
bifun_clu <- function(data1, data2, k, Time1 = NULL, Time2 = NULL, trans = log10, inv.cov = NULL,
                    initial.pars = NULL, iter.max = 1e2, parscale = 1e-3){
  #sort from low to high
  data1 = data1[,order(colSums(data1))]
  data2 = data2[,order(colSums(data2))]
  #attribute
  data = as.matrix(cbind(data1,data2)); n = dim(data)[1]; d = dim(data)[2]
  n1 = dim(data1)[2]; n2 = dim(data2)[2]
  eplison = 1; iter = 0;
  if (is.null(trans)) {
    times1 = as.numeric(colSums(data1))
    times2 = as.numeric(colSums(data2))
    X = data
  } else{
    times1 = as.numeric(trans(colSums(data1)+1));
    times2 = as.numeric(trans(colSums(data2)+1));
    X = trans(data+1)
  }

  if (is.null(Time1)|is.null(Time2) ) {
  } else{
    times1 = as.numeric(Time1)
    times2 = as.numeric(Time2)
  }

  # initial pars
  if (is.null(initial.pars)) {
    initial.pars = biget_par_int(X, k, times1, times2, n1, n2)
  }

  X1 = X[,1:n1]; X2 = X[,(n1+1):(n1+n2)]
  par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
  prob_log <- log(initial.pars$initial_probibality)

  #parscale = par.int[2]*par.int[4]*parscale

  while( abs(eplison) > 1e-3 && iter <= iter.max ){
    #E step
    par.mu <- array(par.int[-c(1:4)], dim = c(k,2,2))
    par.cov = par.int[1:4]
    cov1 = get_SAD1_covmatrix(par.cov[1:2], n1)
    cov2 = get_SAD1_covmatrix(par.cov[3:4], n2)
    mu1 <- power_equation(times1, par.mu[,,1][1:k,])
    mu2 <- power_equation(times2, par.mu[,,2][1:k,])
    mvn_log1 <- sapply(1:k, function(c) dmvnorm(X1, mu1[c,], cov1, log = T))
    mvn_log2 <- sapply(1:k, function(c) dmvnorm(X2, mu2[c,], cov2, log = T))

    mvn_log = mvn_log1 + mvn_log2
    mvn = sweep(mvn_log, 2, FUN = '+', STATS =  prob_log )

    omega_log = t(sapply(1:n, function(c) mvn[c,] - logsumexp(mvn[c,]) ))
    omega = exp(omega_log)

    LL.mem <- biQ_function(par = par.int, prob_log = prob_log, omega_log, X, k, n1, n2, times1, times2)

    #M step
    prob_exp = apply(omega_log, 2, logsumexp)
    prob_log = prob_exp - log(n)

    Q.maximization <- try(optim(par = par.int, biQ_function,
                                prob_log = prob_log,
                                omega_log = omega_log,
                                X = X,
                                k = k,
                                n1 = n1,
                                n2 = n2,
                                times1 = times1,
                                times2 = times2,
                                method = "BFGS",
                                lower = c(-10,-10,-10,-10,rep(-Inf,4*k)),
                                upper = c(10,10,10,10,rep(Inf,4*k)),
                                control = list(trace = TRUE,
                                               parscale = c(rep(parscale,4),rep(1,4*k)),
                                               maxit = 1e2
                                )))
    if ('try-error' %in% class(Q.maximization))
      break
    par.hat <- Q.maximization$par
    par.int = par.hat
    LL.next <- biQ_function(par = par.int, prob_log = prob_log, omega_log, X, k, n1, n2, times1, times2)
    eplison <-  LL.next - LL.mem
    LL.mem <- LL.next
    iter = iter + 1

    cat("\n", "iter =", iter, "\n", "Log-Likelihood = ", LL.next, "\n")
  }
  AIC = 2*(LL.next) + 2*(length(par.hat)+k-1)
  BIC = 2*(LL.next) + log(n)*(length(par.hat)+k-1)

  omega = exp(omega_log)
  X.clustered <- data.frame(X, apply(omega,1,which.max),check.names = F)

  return_obj <- list(cluster_number = k,
                     Log_likelihodd = LL.mem,
                     AIC = AIC,
                     BIC = BIC,
                     cov_par = par.cov,
                     mu_par = par.mu,
                     probibality = exp(prob_log),
                     omega = omega,
                     cluster = X.clustered,
                     cluster2 = data.frame(data, apply(omega,1,which.max), check.names = F),
                     Time1 = times1,
                     Time2 = times2,
                     original_data = data)
  return(return_obj)
}


#' @title parallel version for functional clustering
#' @import parallel
#' @param data1 data for cluster
#' @param data2 data for cluster
#' @param Time1 vector for the time point
#' @param Time2 vector for the time point
#' @param trans indicate log/log2/log10 transform dataset
#' @param start vector for the minimum cluster number
#' @param end vector for the maximum cluster number
#' @param iter.max scales control iteration for EM algorithm
#' @param thread scales for how many thread used
#' @return the initial parameters for functional clustering
#' @export
bifun_clu_parallel <- function(data1, data2, Time1 = NULL, Time2 = NULL, trans = log10, start, end, iter.max = 100, thread = 2){
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(mvtnorm)})
  clusterExport(cl,c(c("bifun_clu","biget_par_int","biQ_function","logsumexp","power_equation_base",
                     "power_equation_all","power_equation","get_biSAD1","get_SAD1_covmatrix"),ls()),
                envir=environment())
  result <- parLapply(cl=cl, start:end, function(c)bifun_clu(data1 = data1,
                                                             data2 = data2,
                                                             k = c,
                                                             Time1 = Time1,
                                                             Time2 = Time2,
                                                             trans = trans,
                                                             iter.max = iter.max))
  stopCluster(cl)

  return(result)
}

#' @title convert result of bifunctional clustering result
#' @param result list directly from bifun_clu_parallel function
#' @param best.k scale of BIC-determined cluster number
#' @return list contain module data and fitted data
#' @export
bifun_clu_convert <- function(result, best.k){
  cluster.result = result[[which(sapply( result , "[[" , 'cluster_number' )==best.k)]]

  times = cluster.result$Time1
  times_new = seq(min(times),max(times),length = 30)

  n1 = length(times);n2 = length(cluster.result$Time2)
  par.mu = cluster.result$mu_par[,,1]
  colnames(par.mu) = c("a","b")
  rownames(par.mu) = paste0("M",1:best.k)

  k = cluster.result$cluster_number
  mu.fit = power_equation(times_new, par.mu[1:k,])
  colnames(mu.fit) = times_new
  rownames(mu.fit) = paste0("M",1:best.k)

  df = cluster.result$cluster2[,c(1:n1,(n1+n2+1))]
  tmp = split(df,df$apply.omega..1..which.max.)
  tmp2 = lapply(tmp, function(x) { x["apply.omega..1..which.max."] <- NULL; x })

  a = list(original_data = mu.fit,
           trans_data = mu.fit,
           power_par = par.mu,
           power_fit = mu.fit,
           Module.all = tmp2)

  times = cluster.result$Time2
  times_new = seq(min(times),max(times),length = 30)

  par.mu = cluster.result$mu_par[,,2]
  colnames(par.mu) = c("a","b")
  rownames(par.mu) = paste0("M",1:best.k)

  k = cluster.result$cluster_number
  mu.fit = power_equation(times_new, par.mu[1:k,])
  colnames(mu.fit) = times_new
  rownames(mu.fit) = paste0("M",1:best.k)

  df = cluster.result$cluster2[,c((n1+1):(n1+n2+1))]
  tmp = split(df,df$apply.omega..1..which.max.)
  tmp2 = lapply(tmp, function(x) { x["apply.omega..1..which.max."] <- NULL; x })

  b = list(original_data = mu.fit,
           trans_data = mu.fit,
           power_par = par.mu,
           power_fit = mu.fit,
           Module.all = tmp2)
  return_obj = list(a = a, b = b)
  return(return_obj)
}

#' @title bifunctional clustering plot
#' @param result list directly from bifun_clu_parallel function
#' @param best.k scale of BIC-determined cluster number
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param degree scalar control transparency degree
#' @param show.legend show legend or not
#' @param color1 Hex Color Codes for first data
#' @param color2 Hex Color Codes for second data
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import scales
#' @return functional clustering plot
#' @export
bifun_clu_plot <- function(result, best.k, label = 10, degree = 1/4, show.legend = FALSE,
                           color1 = "#38E54D", color2 = "#FF8787"){
  cluster.result = result[[which(sapply( result , "[[" , 'cluster_number' )==best.k)]]
  kk = length(table(cluster.result$cluster$apply.omega..1..which.max.))
  if ( kk!= best.k) stop("Please use a smaller k or rerun functional clustering")

  times1 = cluster.result$Time1
  times2 = cluster.result$Time2
  times1_new = seq(min(times1),max(times1),length = 30)
  times2_new = seq(min(times2),max(times2),length = 30)

  n1 = length(times1);n2 = length(times2)
  timesall = c(min(c(times1,times2)),max(c(times1,times2)))

  par.mu = cluster.result$mu_par
  k = cluster.result$cluster_number
  alpha = as.numeric(table(cluster.result$cluster$apply.omega..1..which.max.))

  mu.fit1 = power_equation(times1_new, par.mu[,,1][1:k,])
  mu.fit2 = power_equation(times2_new, par.mu[,,2][1:k,])


  colnames(mu.fit1) = times1_new
  mu.fit1 = melt(as.matrix(mu.fit1))
  colnames(mu.fit1) = c("cluster","x","y")
  mu.fit1$x = as.numeric(as.character(mu.fit1$x))

  colnames(mu.fit2) = times2_new
  mu.fit2 = melt(as.matrix(mu.fit2))
  colnames(mu.fit2) = c("cluster","x","y")
  mu.fit2$x = as.numeric(as.character(mu.fit2$x))

  mu.fit1$type = "L"
  mu.fit2$type = "R"
  mu.fit = rbind(mu.fit1,mu.fit2)

  plot.df1 = cluster.result$cluster[,c(1:n1,ncol(cluster.result$cluster))]
  X1 = plot.df1[,-ncol(plot.df1)]
  plot.df1$name = rownames(plot.df1)
  colnames(plot.df1) = c(times1,"cluster","name")

  plot.df1 = melt(plot.df1, id.vars = c('cluster',"name"))
  colnames(plot.df1) = c("cluster","name", "x","y")
  plot.df1$x = as.numeric(as.character(plot.df1$x))
  plot.df1$name = paste0("a_",plot.df1$name)

  plot.df2 = cluster.result$cluster[,((n1+1):ncol(cluster.result$cluster))]
  X2 = plot.df2[,-ncol(plot.df2)]
  plot.df2$name = rownames(plot.df2)
  colnames(plot.df2) = c(times2,"cluster","name")

  plot.df2 = melt(plot.df2, id.vars = c('cluster',"name"))
  colnames(plot.df2) = c("cluster","name", "x","y")
  plot.df2$x = as.numeric(as.character(plot.df2$x))
  plot.df2$name = paste0("b_",plot.df2$name)

  plot.df1$type = "L"
  plot.df2$type = "R"

  plot.df = rbind(plot.df1,plot.df2)


  df.alpha = data.frame(cluster = 1:k, alpha =normalization(alpha)*degree)
  plot.df = merge(plot.df, df.alpha, by = "cluster")

  name.df = data.frame(label = paste0("M",1:best.k,"(",alpha ,")"),
                       x = mean(timesall), y = max(X1,X2)*0.9, cluster = 1:best.k)

  p = ggplot() + geom_point(plot.df,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(plot.df$type),
                                                 group = "name",
                                                 alpha = factor(plot.df$cluster)),
                            show.legend = show.legend, shape = 1) +
    geom_line(mu.fit, mapping = aes_string(x = "x", y = "y", colour = factor(mu.fit$type)),
              size=1.25, show.legend = show.legend)+
    facet_wrap(~cluster) + scale_alpha_manual(values = df.alpha$alpha) +
    xlab("Habitat Index") + ylab("Niche Index") + theme(axis.title=element_text(size=18)) +
    geom_text(name.df, mapping = aes_string(x = "x", y = "y",label = "label"),
              check_overlap = TRUE, size = 4) +
    theme_bw() + scale_color_manual(name = "Type",values = c("L" = color1,
                                                             "R" = color2))+
    guides(alpha = "none") +
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
