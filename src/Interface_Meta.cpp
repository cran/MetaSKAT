#include <R.h>
#include "Meta.h"

/************************************************
 
 
 ************************************************/

extern "C" {

/*  1: little endian, 0: big endian (only little endian is valid)*/
void IsLittleEndian(int * re){
        
    *re = endian(); 
        
}

/*************************************************************
	MAT file for read
    "MAT" and "MSSD" are the same
*************************************************************/
void	META_MSSD_Read_Open_Init(int * pNumPop, int * err){
	
	*err = Init_MatFile_Read(pNumPop[0]);
}

void	META_MSSD_Read_Open(int * pidx, char ** pFileName, int * err){

	*err = Open_NewMatFile_Read(pidx[0], pFileName[0]);
}

void	META_MSSD_GetData(int * pidx, double * mat, int * pstart, int * pnmarker, int *err){

	*err =Mat_GetData(pidx[0], mat, pstart[0], pnmarker[0]);
}

void	META_MSSD_Read_Close(int * pidx, int * err){
	*err = Close_MatFile(pidx[0]);

}

/*************************************************************
	MAT file for save

*************************************************************/

    
void META_MSSD_Write_Init(int * err){
        
        *err = Mat_Init_Save();
}
    
void	META_MSSD_Write_Open(char ** pFileName, int *err){

	*err =Open_NewMatFile_Save(pFileName[0]);
}

void	META_MSSD_Write(double * mat, int * psize,  int *err){

	*err = Mat_PutData(mat, psize[0]);
}


void	META_MSSD_Check_Saved(int *err){

	*err = Mat_Check_Saved();
}

void	META_MSSD_Num_Sets(int * nsets, int *err){

	*nsets = Mat_Num_Sets();
}

void	META_MSSD_GetStart_Pos(int * pos, int * size, int *err){

	*err = Mat_GetStart_Pos(pos, size);

}
    
void    META_MSSD_Write_Close(int * err){
    *err = Close_Write_MatFile();
}
    

/*************************************************************
	Bed file for read

*************************************************************/
void	META_BED_Init(char** pfilename, int* pNSample, int* pNSnp, int *err){

	*err = Bed_Init(pfilename[0], pNSample[0], pNSnp[0]);
	
}
void	META_BED_Read(int * pIdxs, unsigned char * Genotype, int * pnum, int *err){

	*err = Bed_ReadData(pIdxs, pnum[0], Genotype);
	
}
    
void META_BED_Close(int *err){
    
    *err = Bed_Close();
        
}


} // extern "C"

