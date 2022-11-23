#' @title Q-function to replace log-likelihood function
#' @importFrom mvtnorm dmvnorm
#' @param par numeric vector for parameters need to be estimated
#' @param prob_log mixture component weights(log)
#' @param omega_log latent variables(log)
#' @param X matrix for cluster
#' @param k vector for the cluster number
#' @param times vector for the x values or time points
#' @importFrom mvtnorm dmvnorm
#' @return the Loglikelihood value
Q_function <- function(par, prob_log, omega_log, X, k, times){
  n = dim(X)[1]; d = dim(X)[2];

  par.cov <- par[1:2]
  par.mu <- matrix(par[-c(1:2)],nrow = k,ncol = 2 )

  cov = get_SAD1_covmatrix(par.cov[1:2], n = d)

  mu <- power_equation(times, par.mu[1:k,1:2])

  mvn.log <- sapply(1:k, function(c) dmvnorm(X, mu[c,], cov, log = T))

  tmp = sweep(mvn.log, 2, FUN = "+", STATS = prob_log) - omega_log
  Q = -sum(tmp*exp(omega_log))
  return(Q)
}

#' @title acquire initial parameters for functional clustering
#' @importFrom stats kmeans cov
#' @param X matrix for cluster
#' @param k vector for the cluster number
#' @param times vector for the x values or time points
#' @return the initial parameters for functional clustering
#' @export
get_par_int <- function(X, k, times){
  n = dim(X)[1]; d = dim(X)[2];
  cov.int = c(0.5,mean(diag(cov(X))))

  init.cluster <- kmeans(X,centers = k)
  prob <- table(init.cluster$cluster)/n

  fit <- lapply(1:k,function(c) power_equation_all(times,init.cluster$centers[c,]))
  mu.par.int <- t(sapply(fit, coef))

  return_obj <- list(initial_cov_params = cov.int,
                     initial_mu_params = mu.par.int,
                     initial_probibality = prob)
  return(return_obj)
}


#' @title main function for functional clustering
#' @param data matrix or data for cluster
#' @param k vector for the cluster number
#' @param Time vector for the time point
#' @param trans indicate log/log2/log10 transform dataset
#' @param inv.cov matrix for directly solve cov matrix, default not given(currently not available)
#' @param initial.pars vector for manual give initial parameters, default not given
#' @param iter.max scales control iteration for EM algorithm
#' @param parscale scales control parameters scales for cov pars
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats optim
#' @return the initial parameters for functional clustering
#' @export
fun_clu <- function(data, k, Time = NULL, trans = log10, inv.cov = NULL,
                    initial.pars = NULL, iter.max = 1e2, parscale = 1e-1){
  data = data[,order(colSums(data))]
  #attribute
  data = as.matrix(data); n = dim(data)[1]; d = dim(data)[2]
  eplison = 1; iter = 0;
  if (is.null(trans)) {
    times = as.numeric(colSums(data))
    X = data
  } else{
    times = as.numeric(trans(colSums(data)+1));
    X = trans(data+1)
  }

if (is.null(Time)) {
  } else{
    times = Time
  }
  # initial pars
  if (is.null(initial.pars)) {
    initial.pars = get_par_int(X, k, times)
  }
  par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
  prob_log <- log(initial.pars$initial_probibality)

  parscale = par.int[2]*parscale

  while( abs(eplison) > 1e-3 && iter <= iter.max ){
    #E step
    par.mu <- matrix(par.int[-c(1:2)],nrow = k,ncol = 2 )
    par.cov = par.int[1:2]
    cov = get_SAD1_covmatrix(par.cov[1:2], n = d)
    mu <- power_equation(times, par.mu[1:k,])

    mvn_log <- sapply(1:k, function(c) dmvnorm(X, mu[c,], cov, log = T))

    mvn = sweep(mvn_log, 2, FUN = '+', STATS =  prob_log )

    omega_log = t(sapply(1:n, function(c) mvn[c,] - logsumexp(mvn[c,]) ))
    omega = exp(omega_log)

    LL.mem <- Q_function(par = par.int, prob_log = prob_log, omega_log, X, k, times)

    #M step
    prob_exp = apply(omega_log, 2, logsumexp)
    prob_log = prob_exp - log(n)

    Q.maximization <- try(optim(par = par.int, Q_function,
                                prob_log = prob_log,
                                omega_log = omega_log,
                                X = X,
                                k = k,
                                times = times,
                                method = "BFGS",
                                #lower = c(0,0,rep(-Inf, 2*k)),
                                #upper = Inf,
                                control = list(trace = TRUE,maxit = 1e2,
                                               parscale = c(rep(parscale,2),rep(1,2*k))
                                )))
    if ('try-error' %in% class(Q.maximization))
      break
    par.hat <- Q.maximization$par
    par.int = par.hat
    LL.next <- Q_function(par = par.hat, prob_log = prob_log, omega_log, X, k, times)
    eplison <-  LL.next - LL.mem
    LL.mem <- LL.next
    iter = iter + 1

    cat("\n", "iter =", iter, "\n", "Log-Likelihood = ", LL.next, "\n")
  }
  AIC = 2*(LL.next) + 2*(length(par.hat)+k-1)
  BIC = 2*(LL.next) + log(n)*(length(par.hat)+k-1)

  omega = exp(omega_log)
  X.clustered <- data.frame(X, apply(omega,1,which.max), check.names = F)

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
                     Time = times,
                     original_data = data)
  return(return_obj)
}


#' @title parallel version for functional clustering
#' @import parallel
#' @param data data for cluster
#' @param Time vector for the time point
#' @param trans indicate log/log2/log10 transform dataset
#' @param start vector for the minimum cluster number
#' @param end vector for the maximum cluster number
#' @param iter.max scales control iteration for EM algorithm
#' @param thread scales for how many threads used
#' @return the initial parameters for functional clustering
#' @export
fun_clu_parallel <- function(data, Time = NULL, trans = log10, start, end, iter.max = 100, thread = 2){
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(mvtnorm)})
  clusterExport(cl,c(c("fun_clu","get_par_int","Q_function","logsumexp","power_equation_base",
                     "power_equation_all","power_equation","get_SAD1_covmatrix"),ls()),
                envir=environment())
  result <- parLapply(cl=cl, start:end, function(c) fun_clu(data = data,
                                                            Time = Time,
                                                            trans = trans,
                                                            k = c,
                                                            iter.max = iter.max))
  stopCluster(cl)

  return(result)
}



#' @title plot BIC results for functional clustering
#' @param result list directly from fun_clu_parallel function
#' @param crit either BIC or AIC for module selection
#' @param title title for the plot
#' @import ggplot2
#' @return the BIC plot
#' @export
fun_clu_BIC <- function(result, crit = "BIC", title = NULL){
  BIC.df = data.frame( k = sapply( result , "[[" , 'cluster_number' ),
                       BIC = sapply( result , "[[" , 'BIC' ) )
  AIC.df = data.frame( k = sapply( result , "[[" , 'cluster_number' ),
                       AIC = sapply( result , "[[" , 'AIC' ) )
  min.k = BIC.df$k[which.min(BIC.df$BIC)]
  min.k2 = AIC.df$k[which.min(AIC.df$AIC)]

  p = ggplot(BIC.df)+ theme_bw() + geom_line(mapping = aes_string(x = "k" ,y = "BIC"), color = 'grey10') +
    geom_vline(xintercept = min.k, linetype="dashed", color = "red", size=1.5) +
    xlab("K") + ylab("BIC") + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  p1 = ggplot(AIC.df)+ theme_bw() + geom_line(mapping = aes_string(x = "k" ,y = "AIC"), color = 'grey10') +
    geom_vline(xintercept = min.k2, linetype="dashed", color = "red", size=1.5) +
    xlab("K") + ylab("AIC") + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  if (crit != "BIC") {
    return(p1)
  } else{
    return(p)
  }
}


#' @title convert result of functional clustering result
#' @param result list directly from fun_clu_parallel function
#' @param best.k scale of BIC-determined cluster number
#' @return list contain module data and fitted data
#' @export
fun_clu_convert <- function(result, best.k){
  cluster.result = result[[which(sapply( result , "[[" , 'cluster_number' )==best.k)]]

  times = cluster.result$Time
  times_new = seq(min(times),max(times),length = 30)

  par.mu = cluster.result$mu_par
  colnames(par.mu) = c("a","b")
  rownames(par.mu) = paste0("M",1:best.k)

  k = cluster.result$cluster_number
  mu.fit = power_equation(times_new, par.mu[1:k,])
  colnames(mu.fit) = times_new
  rownames(mu.fit) = paste0("M",1:best.k)

  df = cluster.result$cluster2
  tmp = split(df,df$apply.omega..1..which.max.)
  tmp2 = lapply(tmp, function(x) { x["apply.omega..1..which.max."] <- NULL; x })

  return_obj <- list(original_data = mu.fit,
                     trans_data = mu.fit,
                     power_par = par.mu,
                     power_fit = mu.fit,
                     Module.all = tmp2)
  return(return_obj )
}

#' @title select result of functional clustering result
#' @param result_fit list directly from power_equation_fit
#' @param result_funclu list from fun_clu_convert
#' @param i scale of which cluster selected
#' @return list contain microbe data and fitted data
#' @export
fun_clu_select <- function(result_fit, result_funclu, i){


  cluster.name = rownames(result_funclu$Module.all[[i]])
  original_data = result_fit$original_data[cluster.name,]

  times = result_fit$Time
  times_new = seq(min(times),max(times),length = 30)

  par.mu = result_fit$power_par[cluster.name,]
  mu.fit = result_fit$power_fit[cluster.name,]

  return_obj <- list(original_data = original_data,
                     trans_data = mu.fit,
                     power_par = par.mu,
                     power_fit = mu.fit)
  return(return_obj )
}

#' @title functional clustering plot
#' @param result list directly from fun_clu_parallel function
#' @param best.k scalar of BIC-determined cluster number
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param degree scalar control transparency degree
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import scales
#' @return functional clustering plot
#' @export
fun_clu_plot <- function(result, best.k, label = 10, degree = 1){
  cluster.result = result[[which(sapply( result , "[[" , 'cluster_number' )==best.k)]]
  kk = length(table(cluster.result$cluster$apply.omega..1..which.max.))
  if ( kk!= best.k) stop("Please use a smaller k or rerun functional clustering")

  times = cluster.result$Time
  times_new = seq(min(times),max(times),length = 30)

  par.mu = cluster.result$mu_par
  k = cluster.result$cluster_number
  alpha = as.numeric(table(cluster.result$cluster$apply.omega..1..which.max.))

  mu.fit = power_equation(times_new, par.mu[1:k,])
  colnames(mu.fit) = times_new
  mu.fit = melt(as.matrix(mu.fit))
  colnames(mu.fit) = c("cluster","x","y")
  mu.fit$x = as.numeric(as.character(mu.fit$x))

  plot.df = cluster.result$cluster
  X = plot.df[,-ncol(plot.df)]
  plot.df$name = rownames(plot.df)
  colnames(plot.df) = c(times,"cluster","name")

  plot.df = melt(plot.df, id.vars = c('cluster',"name"))
  colnames(plot.df) = c("cluster","name", "x","y")
  plot.df$x = as.numeric(as.character(plot.df$x))

  df.alpha = data.frame(cluster = 1:k, alpha = min(alpha)/alpha*degree)
  plot.df = merge(plot.df, df.alpha, by = "cluster")

  name.df = data.frame(label = paste0("M",1:best.k,"(",alpha ,")"),
                       x = mean(range(times)), y = max(X)*0.9, cluster = 1:best.k)


  p = ggplot() + geom_point(plot.df,
                           mapping = aes_string(x = "x", y = "y",
                                                colour = factor(plot.df$cluster),
                                                alpha = factor(plot.df$cluster),
                                                group = "name"),
                           show.legend = F, shape = 1) +
    geom_line(mu.fit, mapping = aes_string(x = "x", y = "y", colour = factor(mu.fit$cluster)),
              size=1.25, show.legend = F)+
    facet_wrap(~cluster) + scale_alpha_manual(values = df.alpha$alpha) +
    xlab("Habitat Index") + ylab("Niche Index") + theme(axis.title=element_text(size=18)) +
    geom_text(name.df, mapping = aes_string(x = "x", y = "y",label = "label"),
              check_overlap = TRUE, size = 4) +
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
