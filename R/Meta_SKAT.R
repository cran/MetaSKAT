
Beta.Weights<-function(MAF,weights.beta, Cutoff=1){

	n<-length(MAF)
	weights<-rep(0,n)	
	IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
	if(length(IDX_0) == n){
		#stop("No polymorphic SNPs")
		weights<-rep(0,n)
	} else if( length(IDX_0) == 0){
		weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
	} else {
		weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
	}
	
	#print(length(IDX_0))
	#print(weights[-IDX_0])
	return(weights)
	
}



Meta_SKAT.Work<-function(re, n.g, combined.weight=TRUE, n1=NULL, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate=FALSE, Group_Idx=NULL, MAF.cutoff=1){


	# optimal =optimal.mod
	if(method=="optimal"){
		method="optimal.mod"
	}


	# Combined MAF
	p<-length(re[[1]]$MAF)
	
	MAF.Combine=0
	MAF.Groups<-list()
	Map.Groups<-rep(0,n.g)

	MAF.list<-list()
	for(i in 1:n.g){
	
		MAF.list[[i]]<-re[[i]]$MAF

	}

	for(i in 1:n.g){
	
		MAF.Combine = MAF.Combine + MAF.list[[i]] * n1[i] / sum(n1)
	}

	# Get MAF.Groups when Group_Idx != NULL
	ID.Groups = unique(Group_Idx)	
	for(j in 1:length(ID.Groups)){
		MAF.Groups[[j]] = 0;
		temp<-which(Group_Idx == ID.Groups[j])
		Map.Groups[temp]<-j
		for(i in temp){
			MAF.Groups[[j]] = MAF.Groups[[j]] + MAF.list[[i]] * n1[i] / sum(n1[temp])
		}
	}


	for(i in 1:n.g){
		if(combined.weight == TRUE){
			weight1<-Beta.Weights(MAF.Combine,weights.beta, MAF.cutoff)
		} else {
			j<-Map.Groups[i]
			weight1<-Beta.Weights(MAF.Groups[[j]],weights.beta, MAF.cutoff)
		} 

		re[[i]]$Score =  re[[i]]$Score * weight1
		re[[i]]$SMat.Summary =  t(t(re[[i]]$SMat.Summary * weight1) * weight1)

		if(!is.null(re[[i]]$Score.Resampling)){
			re[[i]]$Score.Resampling =  re[[i]]$Score.Resampling * weight1
		}
	}

	re.method<-SKAT:::SKAT_Check_Method(method,r.corr)

	if(!is.separate){
		re.score<-Meta_SKAT.Work.OneUnit(re, n.g)
	} else {
		re.score<-Meta_SKAT.Work.Groups(re, n.g, ID.Groups, Group_Idx)
	}


	re<-Met_SKAT_Get_Pvalue(re.score$Score, re.score$SMat.Summary, re.method$r.corr, re.method$method, re.score$Score.Resampling)
	return(re)
	
}


####################################################
#
# Assume SMat and Setinfo are already aligned

MetaSKAT_withlist<-function(SMat.list, Info.list, n.cohort, n.each, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, MAF.cutoff=1){

	re<-list()
	p<-length(Info.list[[1]]$MAF)
	for(i in 1:n.cohort){
		
		idx_miss<-which(is.na(Info.list[[i]]$MAF))
		if(length(idx_miss) > 0){
			Info.list[[i]]$Score[idx_miss] = 0
			Info.list[[i]]$MAF[idx_miss] = 0
		}

		re1<-list( Score=Info.list[[i]]$Score, SMat.Summary = SMat.list[[i]], MAF=Info.list[[i]]$MAF )  
		re[[i]]<-re1	
		
	}

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.cohort
	}

	re = Meta_SKAT.Work(re, n.cohort, combined.weight, n1=n.each, weights.beta=weights.beta, method=method, r.corr=r.corr, is.separate=is.separate, Group_Idx=Group_Idx, 
	MAF.cutoff=MAF.cutoff)

	return(re)

}

MetaSKAT_MSSD_OneSet<-function(Cohort.Info, SetID, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, MAF.cutoff=1){

	n.cohort = Cohort.Info$n.cohort
	n.each=rep(0,n.cohort)
	for(i in 1:n.cohort){
		n.each[i]=Cohort.Info$EachInfo[[i]]$header$N.Sample
	}
	
	
	temp<-Get_META_Data_OneSet(Cohort.Info, SetID)
	temp1<-Get_META_Data_OneSet_Align(temp$SMat.list, temp$Info.list, temp$IsExistSNV, n.cohort)

	
	re<-MetaSKAT_withlist(temp1$SMat.list, temp1$Info.list, n.cohort, n.each, combined.weight=combined.weight, 
	weights.beta=weights.beta, method=method, r.corr= r.corr, is.separate = is.separate, Group_Idx=Group_Idx, 
	MAF.cutoff=MAF.cutoff)
	
	return(re)

}


MetaSKAT_MSSD_ALL<-function(Cohort.Info, ...){

	n.cohort = Cohort.Info$n.cohort
	n.set<-length(Cohort.Info$Set_unique)
	
	pval<-rep(NA,n.set)
	for(i in 1:n.set){
		SetID=Cohort.Info$Set_unique[i]
		out<-try(MetaSKAT_MSSD_OneSet(Cohort.Info,SetID, ... ), silent = TRUE)
		if(class(out)!= "try-error"){
			pval[i]<-MetaSKAT_MSSD_OneSet(Cohort.Info,SetID, ... )$p.value
		}
	}
	
	re<-data.frame(SetID=Cohort.Info$Set_unique, p.value=pval)
	return(re)

}


##################################################################
#
#	Genotype matrix Z should be matched with y and X
#
MetaSKAT_wZ<-function(Z, obj, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL){

	if(is.matrix(Z)!= TRUE){
		stop("ERROR: Z is not a matrix!")
	}

	if(class(obj)!= "META_NULL_Model"){
		stop("ERROR: obj class is not META_NULL_Model!")
	}

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	if(length(IDX_MISS) > 0){

		msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )

		warning(msg,call.=FALSE)
		Z<-SKAT:::Impute(Z,impute.method="fixed")
	} 
	

	n.g<-obj$n.g
	re1<-list()
	for(i in 1:n.g){

		ID<-obj$ID[[i]]
		Z1<-as.matrix(Z[ID,])
		res<-obj$out[[i]]$res
		res.out<-obj$out[[i]]$res.out
		X1<-obj$out[[i]]$X1

		if(obj$out_type=="C"){
			s2<-obj$out[[i]]$s2
			re1[[i]]<-Meta_SKAT_SaveData_Linear(res,Z1 , X1, s2, res.out)
		} else if (obj$out_type=="D"){
			pi_1<-obj$out[[i]]$pi_1
			re1[[i]]<-Meta_SKAT_SaveData_Logistic(res,Z1 , X1, pi_1, res.out)
		} else {
			stop("ERROR: out_type is wrong!")
		}
		
	}

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.g
	}

	re = Meta_SKAT.Work(re1, n.g, combined.weight, n1=obj$n.each, weights.beta=weights.beta, method=method, r.corr=r.corr,is.separate=is.separate, Group_Idx=Group_Idx)

	return(re)

}


