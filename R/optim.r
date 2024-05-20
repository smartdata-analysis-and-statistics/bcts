#' @title uniroot.integer.mod
#' @description  Optimizer Uniroot integer modified from the ssanv function package https://github.com/cran/ssanv/blob/master/R/uniroot.integer.R
#'
#' @param f  function for which a root is needed
#' @param power target power value
#' @param lower minimum allowable root
#' @param upper maximum allowable root
#' @param step.power initial step size is 2^step.power
#' @param step.up if TRUE steps up from 'lower', if FALSE steps down from 'upper'
#' @param pos.side if TRUE finds integer, i, closest to the root such that f(i) >0
#' @param maxiter maximum number of iterations (default: 100)
#' @param ... additional arguments to 'f'.
#'
#' @return a list with the following elements:
#' root= the integer on the correct side of the root
#' f.root	= value of f at root  (output= main output with final power, output.test= output table with test at endpoint level)
#' iter	= number of times f was evaluated
#' table.iter = data.frame with the estimated N and power at each iteration
#' table.test = data.frame with the test at endpoint level for a given N and each simulation draw
#'
#' @keywords internal

uniroot <- function(f,
                    power,
                    lower = lower,
                    upper = upper,
                    step.power=step.power, step.up=step.up, pos.side=pos.side,
                    maxiter = 100,
                    ...) {
  iter <- 0
  table.test <- data.frame()
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
    stop("lower < upper  is not fulfilled")
  if (lower == -Inf && step.up == TRUE) stop("lower cannot be -Inf when step.up=TRUE")
  if (upper == Inf && step.up == FALSE) stop("upper cannot be Inf when step.up=FALSE")
  if (step.up) {
    f.old <- f(lower,...)
    iter <- iter + 1
    sign <- 1
    xold <- lower }
  else{
    f.old<-f(upper,...)
    iter<-iter+1
    sign<- -1
    xold<-upper
  }

  ever.switched<-FALSE
  tried.extreme<-FALSE
  while (step.power>-1){

    if ((power-f.old$power)==0) break()
    if (iter >= maxiter) stop("reached maxiter without a solution")
    xnew <- xold + sign*(2^step.power)
    if ((step.up & xnew< upper) || (!step.up & xnew> lower) ){
      f.new <- f(xnew,...)
      iter <- iter + 1
      if (!xold %in% c(table.test$n_control)) {
        table.test <- rbind(table.test, f.old$n)}
    }
    else{

      xnew <- xold
      f.new <- f.old
      step.power <- step.power - 1

      if (!tried.extreme){
        if (step.up) {
          f.extreme <- f(upper,...)
          iter <- iter + 1
          x.extreme <- upper
        } else{
          f.extreme <- f(lower,...)
          iter <- iter+1
          x.extreme <- lower
        }
        tried.extreme <- TRUE
        xswitch <- x.extreme
        f.switch <- f.extreme
        if ((power-f.extreme$power)==0){
          xold<-x.extreme
          f.old<-f.extreme
          break()
        }

        if (((power-f.old$power)*(power-f.extreme$power))>=0){
          warning("f() at extremes not of opposite sign, try to set up upper level to a higher number")
          return(list(iter=iter,f.root=f(upper,...),root=upper,table.test=table.test))
        }
      }
    }

    if ( ((power-f.old$power)*(power-f.new$power))<0){
      sign<- sign*(-1)
      ever.switched<-TRUE
      xswitch<-xold
      f.switch<-f.old
    }
    if (ever.switched){
      step.power<-step.power-1
    }

    xold<- xnew
    f.old<-f.new
    if(step.power<0){
      if(!xold%in%c(table.test$n_control)){
        table.test<-rbind(table.test,f.old$n)}
    }
  }

  if ((power-f.old$power)==0){
    root<-xold
    f.root<-f.old
  } else if ((power-f.new$power)==0){
    root<-xnew
    f.root<-f.new

  } else if ((power-f.switch$power)==0){
    root <- xswitch
    f.root <- f.switch
  } else if (pos.side){
    root <- if((power-f.new$power)>0) xnew else xswitch
    f.root<-if((power-f.new$power)>0) f.new else f.switch
  } else {
    root<-if((power-f.new$power)<0) xnew else xswitch
    f.root<-if((power-f.new$power)<0) f.new else f.switch
  }

  if(!root%in%c(table.test$n_control)){
    table.test<-rbind(table.test,f.old$n)}

  power <- c(root, f.root$power)
  names(power) <- c("n_control","power")

  return(list(power = power,
              table.test = table.test))

}
