require("genalg")
require("DEoptim")


ga.tltsb <- function(formula, slp, h=NULL, iters=2, popsize=50, lower, upper, 
					csteps=2, method="ga", verbose=FALSE){
    X <- model.matrix(formula)
    y <- model.frame(formula)[,1]
    n <- length(y)
    p <- dim(X)[2]
    ind <- rep(0,n)
    if(is.null(h)){
        h <- floor(n/2) + floor((p + 1) / 2)
    }
    if(is.null(method) ||  (method != "de" && method != "ga")){
        cat("Please select a method: 'de' for differential evolution or 'ga' for genetic algorithm\n")
        return(NULL)
    }

   cstep <- function(candidates, csteps){
        cmybetas <- candidates
        indices <- order(abs((y - X%*% cmybetas)/sqrt(slp^2+1)))[1:p]
         
         
        for (i in 1:csteps){
            # ols <- lm(y[indices] ~ x[indices,] - 1)
            # mybetas <- ols$coefficients
         
            beta0 <- mean(y[indices])-slp*mean(x[indices])
            res <- (y - x*slp-beta0)/sqrt(1+slp^2)
            res2 <- abs(res)
            o <- order(res2)
            indices <- sort(o[1:h])
            mybetas <- c(beta0,slp)
        }
        return(mybetas)
    }

    cost <- function(candidates){
        newbetas <- cstep(candidates, csteps)
        res <- (y - X %*% newbetas)/sqrt(1+slp^2)
        fitn <- sum(sort(res^2)[1:h])
        return(fitn)
    }

    best <- rep(0, p)
    if(method == "ga"){   
        ga <- rbga(stringMin = rep(lower, p), stringMax = rep(upper, p), evalFunc = cost, 
					iters = iters, popSize = popsize, verbose = verbose)
        best <- ga$population[1,]
    }else if(method == "de"){
        de <- DEoptim(fn = cost, lower = rep(lower, p), upper = rep(upper, p), 
					control = DEoptim.control(itermax = iters, NP = popsize, trace = verbose))
        best <- de$optim$bestmem       
	 }
	newbetas <- cstep(best, 10)
    res <- (y - X %*% newbetas)/sqrt(1+slp^2)
    crit <- sum(sort(res^2)[1:h])
    result <- list(coefficients = as.double(newbetas), crit = crit, method = method)
    return(result)
}





