MetaSSDR.env <- new.env()
assign("META_SSD_FILE_OPEN_Read.isInit", 0, envir=MetaSSDR.env)

# file MetaInfo has following columns
# setid snpid score MAF missing rate, allele1, allele2, QC, Startpos
# return full path of File.MetaInfo
Open_MSSD_File_2Read<-function(File.MSSD.vec, File.MInfo.vec){

	MetaSKAT_Is_IsLittleEndian()
	
	err_code<-0
	n.cohort<-length(File.MSSD.vec)
	
	if(length(File.MSSD.vec) != length(File.MInfo.vec)){
		stop("Different numbers of MSSD and Meta Info files!")
	}
	
	cat("Number of cohorts = ", n.cohort, "\n")
	
	# Check the existence of files 
	for(i in 1:n.cohort){
		File.MSSD.vec[i]<-normalizePath(File.MSSD.vec[i] ,mustWork =FALSE)
		File.MInfo.vec[i]<-normalizePath(File.MInfo.vec[i] ,mustWork =FALSE)	
	
		SKAT:::Check_File_Exists(File.MSSD.vec[i])
		SKAT:::Check_File_Exists(File.MInfo.vec[i])
	}
	
	# Read files
	temp<-.C("META_MSSD_Read_Open_Init", as.integer(n.cohort), as.integer(err_code))
	err_code<-temp[[2]]
	Print_Error_CODE(err_code)	
	
	re<-list()
	for(i in 1:n.cohort){
		file.idx<-i-1
		File.MSSD<-normalizePath(File.MSSD.vec[i] ,mustWork =FALSE)
		File.MetaInfo<-normalizePath(File.MInfo.vec[i] ,mustWork =FALSE)
	
		temp<-.C("META_MSSD_Read_Open", as.integer(file.idx), as.character(File.MSSD), as.integer(err_code))	
		data.info<-Read_Info_File(File.MetaInfo)
		
		err_code<-temp[[3]]
		Print_Error_CODE(err_code)
	
		re[[i]]<-data.info
	}
	
	# Get unique sets
	Set_unique<-NULL
	for(i in 1:n.cohort){
		Set_unique<-union(Set_unique, re[[i]]$set_unique)
	}	
	
	assign("META_MSSD_FILE_OPEN_Read.isInit", 1, envir=MetaSSD.env);
	info<-list(n.cohort=n.cohort, Set_unique = Set_unique, EachInfo=re)
	
	return(info)
}

Close_MSSD_Files<-function(){

	if(get("META_MSSD_FILE_OPEN_Read.isInit", envir=MetaSSDR.env) == 1){
		temp<-.C("META_MSSD_Read_Close")
		Msg<-sprintf("Close the opened MSSD files \n")
		cat(Msg)
		assign("META_MSSD_FILE_OPEN_Read.isInit", 0, envir=MetaSSDR.env);
	} else{
		Msg<-sprintf("No opened MSSD file!\n");
		cat(Msg)		
	}

}

Read_Read_MSSD_File<-function(cohort_idx, start, nmarker){


	size=nmarker * nmarker

	err_code<-0
	temp<-.C("META_MSSD_GetData", as.integer(cohort_idx-1), double(size), as.integer(start), as.integer(nmarker), as.integer(err_code))

	err_code<-temp[[5]]
	Print_Error_CODE(err_code)	
	
	SMat.out<-matrix(temp[[2]],byrow=TRUE, nrow=nmarker)
	return(SMat.out)
}


Get_META_Data_OneSet<-function(Cohort.Info, SetID){

	SMat.list<-list()
	Info.list<-list()
	IsExistSNV<-rep(0,Cohort.Info$n.cohort)
	
	for(i in 1:Cohort.Info$n.cohort){
		idx<-Cohort.Info$EachInfo[[i]]$hash_set[[SetID]]
		
		if(is.null(idx)){
			IsExistSNV[i]<-0
			
		} else {
			Info.list[[i]]<-Cohort.Info$EachInfo[[i]]$Info[idx,]
		
			start<-Info.list[[i]]$StartPOS[1]
			nmarker<-dim(Info.list[[i]])[1]
		
			cohort_idx<-i
			SMat.list[[i]]<-Read_Read_MSSD_File(cohort_idx, start, nmarker)
			
			IsExistSNV[i]<-1
		}
		
	}
	
	return(list(SMat.list=SMat.list, Info.list=Info.list, IsExistSNV=IsExistSNV))
	
}



Get_META_Data_OneSet_Align<-function(SMat.list, Info.list, IsExistSNV,  n.cohort){

	SnpID.all<-NULL
	for(i in 1:n.cohort){
		if(IsExistSNV[i] == 1){
			SnpID.all<-union(SnpID.all, Info.list[[i]]$SNPID)
		}
	}
	SnpID.all<-unique(SnpID.all)
	n.all<-length(SnpID.all)
	data.master<-data.frame(SNPID = SnpID.all, IDX=1:length(SnpID.all))
	
	SMat.list.new<-list()
	Info.list.new<-list()
	
	for(i in 1:n.cohort){
	
		if(IsExistSNV[i] == 1){	
			data1<-Info.list[[i]]
			data1$IDX1=1:length(data1$SNPID)
			data2<-merge(data.master, data1, by.x="SNPID", by.y="SNPID", all.x=TRUE)
	
			IDX<-which(!is.na(data2$IDX1))
			IDX1<-data2$IDX1[IDX]
	
			SMat.list.new[[i]]<-matrix(rep(0,n.all* n.all), ncol=n.all)
			SMat.list.new[[i]][IDX,IDX]<-SMat.list[[i]][IDX1,IDX1]
		
			Info.list.new[[i]]<-data2
		} else {
			SMat.list.new[[i]]<-matrix(rep(0,n.all* n.all), ncol=n.all)
			data2<-data.master
			data2$Score= rep(0, n.all)
			data2$MAF = rep(0, n.all)			
			Info.list.new[[i]]<-data2
			
		}
	}
	
	return(list(SMat.list = SMat.list.new, Info.list=Info.list.new))
	
}

Read_Info_File_Header<-function(File.MetaInfo){



	N.Sample=-1
	N.Sets=-1
	N.SNPs=-1	
	N.SNPs.unique=-1
	for(i in 1:10){
		temp<-read.delim(File.MetaInfo, nrows=1, skip=i-1, header=FALSE, stringsAsFactors=FALSE)[1,1]
		if(substr(temp, start=1, stop=1) != "#"){
			break
		}

		temp1<-strsplit(temp, "=")
		if(temp1[[1]][1] =="#N"){
			N.Sample = as.numeric(temp1[[1]][2])
		} else if(temp1[[1]][1] =="#nSets"){
			N.Sets = as.numeric(temp1[[1]][2])
		} else if(temp1[[1]][1] =="#nSNPs"){
			N.SNPs = as.numeric(temp1[[1]][2])
		} else if(temp1[[1]][1] =="#nSNPs.unique"){
			N.SNPs.unique = as.numeric(temp1[[1]][2])
		}
		
	}

	info.header<-list(N.Sample=N.Sample, N.Sets=N.Sets, N.SNPs=N.SNPs, N.SNPs.unique=N.SNPs.unique)
	msg<-sprintf("%d samples, %d sets, %d SNPs and %d unique SNPs\n", 
	N.Sample, N.Sets, N.SNPs, N.SNPs.unique)

	cat(msg)
	return(info.header)

}

Read_Info_File<-function(File.MetaInfo){

	header=Read_Info_File_Header(File.MetaInfo)
	
	info<-read.table(File.MetaInfo, header=TRUE, stringsAsFactors=FALSE)
	info$SetID<-as.character(info$SetID)
	x<-unique(info$SetID_numeric)
	x.n<-length(x)
	
	vec<-as.vector(info$SetID_numeric)
	idx.int.end<-findInterval(x,vec)
	idx.int.start<-c(1,idx.int.end[-x.n]+1)

	hash_set<-new.env()
	set_unique<-rep(" ",x.n)
	for(i in 1:x.n){
	
		idx.start<-idx.int.start[i]
		idx.end<-idx.int.end[i]
		
		SetID1<-info$SetID[idx.start]
		SetID2<-info$SetID[idx.end]
		
		if(SetID1 != SetID2){
			msg<-sprintf("Error: Read_Set_Info, SetIDs don't match [%s] [%s]\n", SetID1, SetID2)
			stop(msg)
		}
		set_unique[i]<-SetID1
		val<-idx.start:idx.end
		hash_set[[SetID1]]<-val
	}
	
	re<-list(Info=info, hash_set=hash_set, set_unique=set_unique, header=header)
	return(re)
}

#	Hash example
#
# hash1<-new.env()
#for(i in 1:10000){
#	ID<-sprintf("Set%06d",i)
#	val<-i:rbinom(1,20,0.5)	
#	hash1[[ID]]<-val
#}
