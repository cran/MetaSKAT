MetaSSD.env <- new.env()
Print_Error_CODE<-function(err_code){

	if(err_code > 0){
		stop("Error code = ", err_code, "\n")
	}

}

##################################################################
#
#	Open and Close Files

assign("META_BED_FILE_OPEN.isOpen", 0, envir=MetaSSD.env)
assign("META_BED_FILE_OPEN.FileName", "", envir=MetaSSD.env)
assign("META_MSSD_FILE_OPEN_Write.isOpen", 0, envir=MetaSSD.env)
assign("META_MSSD_FILE_OPEN_Write.FileName", "", envir=MetaSSD.env)

#META_BED_FILE_OPEN.isOpen = 0
#META_BED_FILE_OPEN.FileName = ""

#META_MSSD_FILE_OPEN_Write.isOpen = 0
#META_MSSD_FILE_OPEN_Write.FileName = ""


Open_BED_File<-function(File.Bed, File.Bim, N.Sample){

	err_code<-0
	File.Bed<-normalizePath(File.Bed ,mustWork =FALSE)
	File.Bim<-normalizePath(File.Bim ,mustWork =FALSE)

	
	SKAT:::Check_File_Exists(File.Bed)
	SKAT:::Check_File_Exists(File.Bim)

	
	# read bim file
	cat("Read Bim file\n")
	BimInfo<-try(read.table(File.Bim, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(BimInfo)=="try-error"){
		stop("Error in Bim file!") 
	}
	nMarker<-dim(BimInfo)[1]
	colnames(BimInfo)<-c("Chr","SnpID", "cM","base","a1","a2")
	BimInfo$idx<-1:nMarker
	
	cat("Bim file has", nMarker , "markers\n")
	
	
	# open bed file
	if(get("META_BED_FILE_OPEN.isOpen", envir=MetaSSD.env) == 1){
		Close_BED_File()
	}
	temp<-.C("META_BED_Init", as.character(File.Bed), as.integer(N.Sample), as.integer(nMarker), as.integer(err_code))

	error_code<-temp[[4]]
	Print_Error_CODE(error_code)
	
	assign("META_BED_FILE_OPEN.isOpen", 1, envir=MetaSSD.env);
	assign("META_BED_FILE_OPEN.FileName",File.Bed, envir=MetaSSD.env)


	return(BimInfo)

}

Close_BED_File<-function(){

	err_code<-0
	if(get("META_BED_FILE_OPEN.isOpen", envir=MetaSSD.env) == 1){
		temp<-.C("META_BED_Close", as.integer(err_code))
		Msg<-sprintf("Close the opened Bed file: %s\n"
		,get("META_BED_FILE_OPEN.FileName", envir=MetaSSD.env));
		
		cat(Msg)
		assign("META_BED_FILE_OPEN.isOpen", 0, envir=MetaSSD.env);
	} else{
		Msg<-sprintf("No opened SSD file!\n");
		cat(Msg)		
	}

}

Read_From_BED_File<-function(IDX, N.Sample){

	if(get("META_BED_FILE_OPEN.isOpen", envir=MetaSSD.env) == 0){
		stop("BED file hasn't been opened")
	}

	
	err_code<-0
	N.SNP<-length(IDX)
	size<-N.SNP * N.Sample

	Z<-rep(11,size)

	temp<-.C("META_BED_Read",as.integer(IDX),as.raw(Z),as.integer(N.SNP), as.integer(err_code))

	error_code<-temp[[4]]
	Print_Error_CODE(error_code)
	
	Z.out<-t(matrix(as.numeric(temp[[2]]),byrow=TRUE, nrow=N.SNP))
	
	IDX_MISS<-which(Z.out == 9)
	if(length(IDX_MISS) > 0){
		Z.out[IDX_MISS]<-NA
	}
	MAF<-colMeans(Z.out, na.rm=TRUE)/2
	
	return(list(Z=Z.out, MAF=MAF))

}

# file MetaInfo has following columns
# setid snpid score MAF missing rate, allele1, allele2, QC
# return full path of File.MetaInfo
Open_Write_MSSD_File<-function(File.MSSD, File.MetaInfo){

	err_code<-0
	File.MSSD<-normalizePath(File.MSSD ,mustWork =FALSE)
	File.MetaInfo<-normalizePath(File.MetaInfo ,mustWork =FALSE)
	

	
	# open bed file
	if(get("META_MSSD_FILE_OPEN_Write.isOpen", envir=MetaSSD.env) == 1){
		Close_Write_MSSD_File()
	}
	
	temp<-.C("META_MSSD_Write_Init", as.integer(err_code))
	err_code<-temp[[1]]
	Print_Error_CODE(err_code)
	
	temp<-.C("META_MSSD_Write_Open", as.character(File.MSSD), as.integer(err_code))

	err_code<-temp[[2]]
	Print_Error_CODE(err_code)
	
	assign("META_MSSD_FILE_OPEN_Write.isOpen", 1, envir=MetaSSD.env);
	assign("META_MSSD_FILE_OPEN_Write.FileName",File.MSSD, envir=MetaSSD.env)

	return(File.MetaInfo)

}

Close_Write_MSSD_File<-function(){

	err_code<-0
	if(get("META_MSSD_FILE_OPEN_Write.isOpen", envir=MetaSSD.env) == 1){
		temp<-.C("META_MSSD_Write_Close", as.integer(err_code))
		Msg<-sprintf("Close the opened MSSD file: %s\n"
		,get("META_MSSD_FILE_OPEN_Write.FileName", envir=MetaSSD.env));
		
		cat(Msg)
		assign("META_MSSD_FILE_OPEN_Write.isOpen", 0, envir=MetaSSD.env);
	} else{
		Msg<-sprintf("No opened MSSD files!\n");
		cat(Msg)		
	}

}


Write_To_MSSD_File<-function(x){

	err_code<-0
	if(length(which(is.na(x))) > 0){
		stop("Error: Summary stat has NA values")
	}
	x1<-x[lower.tri(x, diag = TRUE)]
	size<-length(x1)
	
	temp<-.C("META_MSSD_Write", as.double(x1), as.integer(size),  as.integer(err_code))
	
	err_code<-temp[[3]]
	Print_Error_CODE(err_code)
	
	byte_used<-size* 4 + 4
	
	return(byte_used)
	
}

Check_Saved<-function(nsets, StartPos){


	err_code<-0
	nsets1<-0
	
	############################
	# check # of sets
	
	temp<-.C("META_MSSD_Num_Sets", as.integer(nsets), as.integer(err_code))
	err_code<-temp[[2]]
	Print_Error_CODE(err_code)
	nsets1<-temp[[1]]
	if(nsets1 != nsets){
		msg<-sprintf("Error: MSSD file save- number of sets [%d] [%d] \n", nsets, nsets1)
		stop(msg)
	}
	
	############################
	# check Startpos
	pos1<-rep(0,nsets)
	size1<-rep(0,nsets)
	
	temp<-.C("META_MSSD_GetStart_Pos", as.integer(pos1), as.integer(size1), as.integer(err_code) )
	err_code<-temp[[3]]
	Print_Error_CODE(err_code)
	
	pos1<-temp[[1]]
	if( sum((pos1-StartPos)^2) != 0){
		msg<-sprintf("Error: MSSD file save- startpos \n")
		cat("[POS1]", pos1, "\n")
		cat("[POS2]", StartPos, "\n")
		stop(msg)
	}
	
	############################
	# CRC check	

	temp<-.C("META_MSSD_Check_Saved", as.integer(err_code) )
	err_code<-temp[[1]]
	Print_Error_CODE(err_code)
		
	cat("Save was done successfully!\n")
}


Write_To_MetaInfo_Header<-function(File.MetaInfo, n.all, n, nSets, nSNPs, nSNPs.unique){

	header<-rbind(sprintf("#N.ALL=%d",n.all),
	sprintf("#N=%d",n),
	sprintf("#nSets=%d",nSets),
	sprintf("#nSNPs=%d",nSNPs),
	sprintf("#nSNPs.unique=%d",nSNPs.unique),
	sprintf("SetID SetID_numeric SNPID Score MAF MissingRate Allele1 Allele2 MinorAllele PASS StartPOS"))
	
	write.table(header, File.MetaInfo, col.names=FALSE, row.names=FALSE, quote=FALSE)

}

Write_To_MetaInfo<-function(File.MetaInfo, SetID, SetID_num, SNPID, Score, MAF, missing_rate, a1, a2, Minor.Allele, PASS, StartPos){


	nSNP<-length(a1)
	data.w<-data.frame(SetID=rep(SetID,nSNP), SetID_num = rep(SetID_num,nSNP), SNPID=SNPID, Score=Score
	, MAF=MAF, MissingRate=missing_rate, Allele1=a1, Allele2=a2, MinorAllele=Minor.Allele, PASS, StartPos=rep(StartPos,nSNP))
	write.table(data.w, File.MetaInfo, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep=" ")

}



Generate_Meta_Files<-function(obj, File.Bed, File.Bim, File.SetID, File.MSSD, File.MInfo, N.Sample, data=NULL){

	MetaSKAT_Is_IsLittleEndian()
	
	Magic_Byte<-1
	# Read MAP file and SetID file
	File.SetID<-normalizePath(File.SetID ,mustWork =FALSE)
	SKAT:::Check_File_Exists(File.SetID)


	##########################################
	# 	Read files 
	if(N.Sample == -1){
		N.Sample<-length(obj$res)
	}


	# read setid
	cat("Read SetID file\n")
	SetInfo<-try(read.table(File.SetID, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(SetInfo)=="try-error"){
		stop("Error in SetID file!") 
	}
	colnames(SetInfo)<-c("SetID","SnpID")
	SetIDs<-unique(SetInfo[,1])
	nSets<-length(SetIDs)
	nMarker.Sets<-dim(SetInfo)[1]

	cat("SetID file has", nSets, "sets\n")

	# open bim
	BimInfo<-Open_BED_File(File.Bed, File.Bim, N.Sample)

	# open MSSD
	File.MInfo = Open_Write_MSSD_File(File.MSSD, File.MInfo)

	###########################################
	# 	Merge two files 

	cat("Merge datasets and get set info\n")

	SetInfo$SetID.Num<-as.numeric(as.factor(SetInfo$SetID))
	SetInfo$IDX<-1:length(SetInfo$SetID)
	
	Data.all<-merge(SetInfo, BimInfo, by.x="SnpID", by.y="SnpID", all.x=FALSE, all.y=FALSE)
	ID.nofind = which(is.na(Data.all$idx))
	
	ord<-order(Data.all$SetID.Num, Data.all$IDX)
	Data.all1<-Data.all[ord,]

	###################################################
	# 	Generate Intervals
	
	SetID.Num.unique<-unique(Data.all1$SetID.Num)

	x<-1:length(SetID.Num.unique)
	x.n<-length(x)
	vec<-as.vector(Data.all1$SetID.Num)
	idx.int.end<-findInterval(x,vec)
	idx.int.start<-c(1,idx.int.end[-x.n]+1)
	
	SNPID.Num<-length(Data.all1$SnpID)
	SNPID.Num.unique<-length(unique(Data.all1$SnpID))


	Write_To_MetaInfo_Header(File.MInfo, N.Sample, N.Sample, x.n, SNPID.Num, SNPID.Num.unique)
	byte_used<-Magic_Byte
	StartPos<-rep(0, length(x))
	for(i in x){
		idx.start<-idx.int.start[i]
		idx.end<-idx.int.end[i]
		SetID<-Data.all1$SetID[idx.start]
		
		IDX<-idx.start:idx.end
		
		SNPID<-Data.all1$SnpID[IDX]
		SNPIndex<-Data.all1$idx[IDX]
		a1.org<-Data.all1$a1[IDX]
		a2.org<-Data.all1$a2[IDX]
		
		
		out<-Read_From_BED_File(SNPIndex, N.Sample)
		Z<-out$Z
		MAF<-out$MAF

		Allele1<-a1.org
		Allele2<-a2.org
		Minor.Allele<-a2.org		
		IDX1<-which(MAF > 0.5)
		
		if(length(IDX1) > 0){
			Minor.Allele[IDX1]<-a1.org[IDX1]
			MAF[IDX1]<-1-MAF[IDX1]
		}
		
		#
		#	change allele 1 and 2 by alphabetical order 
		#	currently allel1 is major

		IDX1<-which(Allele1 > Allele2)
		if(length(IDX1) > 0){
			Allele1[IDX1]<-a2.org[IDX1]
			Allele2[IDX1]<-a1.org[IDX1]
			Z[,IDX1] = 2-Z[,IDX1]
		}
		
		
		
		nSNP<-length(Allele1)
		PASS1<-rep("PASS",nSNP)
		
		re<-Meta_SKAT_SaveData(Z, obj, SetID=SetID, impute.method = "fixed")
		Write_To_MetaInfo(File.MInfo, SetID, i, SNPID, re$Score, MAF, re$missing_rate, Allele1, Allele2, Minor.Allele, PASS1 ,byte_used )
		StartPos[i]<-byte_used
		
		byte1 = Write_To_MSSD_File(re$SMat.Summary)
		byte_used = byte_used + byte1
		
	}
	
	Check_Saved(length(x), StartPos)

	
}

