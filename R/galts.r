require("genalg");

ga.lts<-function(formula, h=NULL, iters=2, popsize=50, lower , upper , csteps=2){

	ols<-lm(formula,x=TRUE,y=TRUE);
	x<-ols$x;
	y<-ols$y;
	n<-length(y);
	p<-dim(x)[2];
	ind<-rep(0,n);
	if(is.null(h)){
		h<-floor(n/2)+floor((p+1)/2);
	}

	cstep<-function(candidates){
		cmybetas<-candidates;
		res<-y-x%*%cmybetas;
		indices<-order(abs(res))[1:p];
		for (i in 1:csteps){
			ind<<-rep(0,n);
			ind[indices]<<-1;
			ols<-lm(formula=y~x-1, weights=ind);
			mybetas<-ols$coefficients;
			res<-y-x%*%mybetas;
			res2<-abs(res);
			o<-order(res2);
			indices<-sort(o[1:h]);
		}
		return(mybetas);
	}


	cost<-function(candidates){
		newbetas<-cstep(candidates);
		res<-y-x%*%newbetas;
		fitn<-sum(sort(res^2)[1:h]);
		return(fitn);	
	}

	
	ga<-rbga(stringMin=rep(lower,p), stringMax=rep(upper,p), evalFunc=cost, iters=iters, popSize=popsize);
	best<-ga$population[1,];
	csteps<<-10;
	newbetas<-cstep(best);	
	res<-y-x%*%newbetas;
	crit<-sum(sort(res^2)[1:h]);
	result<-list(coefficients=newbetas, crit=crit);
	return(result);
}

