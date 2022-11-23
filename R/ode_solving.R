#' @title quasi-dynamic lotka volterra model
#' @param Time vector of time point
#' @param State vector of ODE initial state
#' @param Pars vector for unknown ODE parameters
#' @param power_par matrix of power equation parameters for dependent effect
#' @return list used in ode function
qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)

  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}

#' @title least-square fit for qdODE model
#' @importFrom deSolve ode
#' @param pars vector for unknown ODE parameters
#' @param data data contain independent effect as first row and dependent effect
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @return mean-square error
qdODE_ls <- function(pars, data, Time, power_par){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(data[1,])
  fit = as.numeric(out[,2])
  ind = as.numeric(out[,(n+2)])
  sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  return(sse)
}

#' @title legendre polynomials fit to qdODE model
#' @importFrom deSolve ode
#' @param pars vector of qdODE parameters
#' @param data dataframe of observed data
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @param LOP_order scalar of LOP order
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @return list contain legendre polynomials parameters, qdODE values and LOP fitted values
qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = data[1,], y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL

  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))

  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}



#' @title wrapper for qdODE model
#' @importFrom stats optim
#' @param result result from power_equation_fit
#' @param relationship list contain variable selection results
#' @param i scalar for which id used for qdODE solving, must <= nrow
#' @param init_pars scalar for initial parameters
#' @param LOP_order scalar of LOP order
#' @param method scalar of qdODE solving methodm, cuurent only support least square
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = 6, method = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e3){
  Time = as.numeric(colnames(result$power_fit))
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit[variable,]

  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = result$power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (method == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))

      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         #lower = c(0, rep(-10,(length(pars_int))-1)),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))

      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}


#' @title wrapper for qdODE_all in parallel version
#' @param result result from power_equation_fit
#' @param reduction use n/log(n) dimension reduction
#' @param thread scales for how many threads used
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_parallel <- function(result, reduction = FALSE, thread = 2, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit"), envir=environment())
  result = parLapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                        relationship = relationship,
                                                        i = c,
                                                        maxit = maxit
  ), cl = cl)
  stopCluster(cl)
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}

#' @title convert qdODE results to plot data
#' @importFrom reshape2 melt
#' @param result list of qdODE all
#' @importFrom stats na.omit
qdODEplot_convert <- function(result){
  data = result$predict
  n = ncol(data)
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))

  plot.df = melt(data, id.vars = c("x"))

  name = levels(plot.df[,2])

  ind.name = name[grep("ind", name)]
  ind.name2 = strsplit(ind.name,split = "\\.")[[1]][2]
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  ind.df$type = "ind"
  ind.df$variable = ind.name2

  depname = levels(plot.df[,2])[grep("dep",name )]
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  dep.df$type = "dep"
  dep.df$variable = sapply(strsplit(as.character(dep.df$variable),"\\."),"[",2)


  original.df = subset(plot.df, plot.df[,2] == "y")
  original.df$type = "original"

  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  fit.df$type = "fit"

  plot.df2 = rbind(ind.df, dep.df,fit.df)

  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  name.df = name.df[-nrow(name.df),]
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2

  name.df = name.df[-which(name.df[,4] == "fit"),]

  name.df[,1] = name.df[,1]*1.002
  return_obj = list(plot.df2 = plot.df2,
                    name.df = name.df,
                    ind.name2 = ind.name2)
  return(return_obj)
}


#' @title plot single decompose plot
#' @import ggplot2
#' @param result list of qdODE all
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @export
qdODE_plot_base <- function(result,label = 10, show.legend = TRUE){
  result2 = qdODEplot_convert(result)
  plot.df2 = result2$plot.df2
  name.df = result2$name.df
  ind.name2 = result2$ind.name2

  p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df, mapping = aes_string(label = "variable", colour = "type",
                                            x = "x", y = "value"), show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) +
    xlab("Habitat Index") + ylab("Niche Index") +
    ggtitle(ind.name2) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y = element_text(hjust = 0))

  if (is.null(label)) {
    p = p + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))
  } else {
    xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
    ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks

    xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )

    if (0 %in% ylabel) {
      pos0 = which(ylabel==0)
      text=paste(label,"^", ylabel, sep="")
      text[pos0] = "0"
      if (any(na.omit(ylabel)<0)) {
        pos1 = which(ylabel<0)
        text[pos1] = paste0("-",label,"^",abs(ylabel[pos1]))
        ylabel2 = parse(text=text)
      } else{ylabel2 = parse(text=text)}
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p = p + scale_x_continuous(labels = xlabel2, limits = c(min(plot.df2$x), max(plot.df2$x)*1.005)) +
      scale_y_continuous(labels = ylabel2)+ theme(axis.text.y = element_text(hjust = 0))
  }
  return(p)
}

#' @title plot all decompose plot
#' @import ggplot2 patchwork
#' @param result list of qdODE parallel
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @param nrow scalar for subplot row number
#' @param ncol scalar for subplot column number
#' @return all effect curve decompose plot
#' @export
qdODE_plot_all <- function(result,label = 10, show.legend = TRUE, nrow = NULL, ncol = NULL){
  p = lapply(result$ode_result, qdODE_plot_base, label = label, show.legend = show.legend)
  #p = lapply(1:length(result$ode_result), function(c) qdODE_plot_base(result$ode_result[[c]], label = label))
  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))

  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Niche Index", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    plot_annotation(caption = "Habitat Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.05, 1),guides = 'collect')

  return(pp)
}

#' @title plot single decompose plot for two data
#' @import ggplot2
#' @param result1 list of qdODE all for first data
#' @param result2 list of qdODE all for second data
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @param remove.label to remove x and y label
#' @export
biqdODE_plot_base <- function(result1, result2, label = 10, show.legend = FALSE, remove.label = FALSE){
  resulta = qdODEplot_convert(result1)
  resultb = qdODEplot_convert(result2)
  plot.df1 = resulta$plot.df2
  name.df1 = resulta$name.df
  ind.name2 = resulta$ind.name2
  name.df1$x = name.df1$x*0.99

  plot.df2 = resultb$plot.df2
  name.df2 = resultb$name.df
  name.df2$x = name.df2$x*0.99

  lower1 = min(plot.df1[,3])
  upper1 = max(plot.df1[,3])
  lower2 = min(plot.df2[,3])
  upper2 = max(plot.df2[,3])
  y_min = round(min(lower1, lower2),1)-0.05
  y_max = round(max(upper1, upper2),1)+0.05

  p1 = ggplot(plot.df1, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df1, mapping = aes_string(label = "variable", colour = "type",
                                             x = "x", y = "value"), size = 3,show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) + theme_bw() +
    xlab("Habitat Index") + ylab("Niche Index") +
    theme(axis.text.y = element_text(hjust = 0)) + xlab(NULL)+
    ggtitle(ind.name2) + theme(plot.title = element_text(hjust = 1))


  s1 = p1 + scale_y_continuous(limits = c(y_min, y_max))

  p2 = ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df2, mapping = aes_string(label = "variable", colour = "type",
                                             x = "x", y = "value"), size = 3,show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) + theme_bw() +
    theme(axis.text.y = element_text(hjust = 0)) + xlab(NULL)

  s2 = p2 + scale_y_continuous(limits = c(y_min, y_max))

  if (is.null(label)) {
    p1 = p1 + scale_x_continuous(limits = c(min(plot.df1$x), max(plot.df1$x)*1.005))+
      theme(axis.title.x = element_text(vjust=1))
    p2 = p2 + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))
  } else {
    xlabel1 = ggplot_build(s1)$layout$panel_params[[1]]$x.sec$breaks
    xlabel2 = ggplot_build(s2)$layout$panel_params[[1]]$x.sec$breaks

    ylabel = ggplot_build(s1)$layout$panel_params[[1]]$y.sec$breaks

    xlabel1.2 = parse(text= paste(label,"^", xlabel1, sep="") )
    xlabel2.2 = parse(text= paste(label,"^", xlabel2, sep="") )
    if (0 %in% ylabel) {
      pos0 = which(ylabel==0)
      text=paste(label,"^", ylabel, sep="")
      text[pos0] = "0"
      if (any(na.omit(ylabel)<0)) {
        pos1 = which(ylabel<0)
        text[pos1] = paste0("-",label,"^",abs(ylabel[pos1]))
        ylabel2 = parse(text=text)
      } else{ylabel2 = parse(text=text)}
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p1 = p1 + scale_x_continuous(labels = xlabel1.2, limits = c(min(xlabel1)-0.05, max(xlabel1)+0.05)) +
      scale_y_continuous(limits = c(y_min,y_max),labels = ylabel2) +
      theme(axis.text.y = element_text(hjust = 0))+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))+
      theme(axis.title.x = element_text(vjust=1))

    p2 = p2 + scale_x_continuous(labels = xlabel2.2, limits = c(min(xlabel2)-0.05, max(xlabel2)+0.05)) +
      scale_y_continuous(limits = c(y_min,y_max),labels = ylabel2) +
      theme(axis.text.y = element_text(hjust = 0)) +
      theme(axis.text.y = element_blank(), axis.ticks.length.y = unit(-0.1,"cm")) +
      ylab(NULL)+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))

  }

  if (remove.label == TRUE) {
    p1 = p1 + theme(axis.title=element_blank())
    p2 = p2 + theme(axis.title=element_blank())
  }

  pp = p1+p2
  return(pp)
}


#' @title plot all decompose plot for two data
#' @import ggplot2
#' @param result1 list of qdODE all for first data
#' @param result2 list of qdODE all for second data
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @param remove.label to remove x and y label
#' @param nrow scalar for subplot row number
#' @param ncol scalar for subplot column number
#' @export
biqdODE_plot_all <- function(result1, result2, label = 10, show.legend = FALSE,
                             remove.label = TRUE, nrow = NULL, ncol = NULL){
  n = length(result1$ode_result)
  p = lapply(1:n, function(c)
    biqdODE_plot_base(result1$ode_result[[c]], result2$ode_result[[c]],
                      label = label, show.legend = show.legend, remove.label = remove.label))

  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))

  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Niche Index", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    plot_annotation(caption = "Habitat Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.03, 1),guides = 'collect')

  return(pp)
}
