

#include "Read_BED.h"
#include "MatFile.h"
#include "error_messages.h"

MatFile * 	g_pMatFile_Read;
MatFile 	g_MatFile_Save;
int		g_nPop;		/* Number of Population */
BedFile		g_BedFile;

/****************************************************
 endian is important. The program is only valid
 under little endian.
 *****************************************************/
#define LITTLE_ENDIAN_DEF 1
#define BIG_ENDIAN_DEF    0


int endian() {
    int i = 1;
    char *p = (char *)&i;
    
    if (p[0] == 1)
        return LITTLE_ENDIAN_DEF;
    else
        return BIG_ENDIAN_DEF;
}



/*************************************************************
	MAT file for read

*************************************************************/
int	Init_MatFile_Read(int NumPop){

	g_nPop = NumPop;
	g_pMatFile_Read = new MatFile[NumPop];
	if(!g_pMatFile_Read){
		return ERORR_MAT_INIT;
	}

	return NO_ERRORS;

}

int	Open_NewMatFile_Read(int idx, const char * FileName){

	int re;
	if(idx >= g_nPop){
		return ERORR_MAT_IDX_OUTOFBOUND;
	}
	re = g_pMatFile_Read[idx].OpenToRead(FileName);
	return re;

}

int	Mat_GetData(int idx, double * mat, int start, int nmarker){

	int re;
	if(idx >= g_nPop){
		return ERORR_MAT_IDX_OUTOFBOUND;
	}
	re = g_pMatFile_Read[idx].GetData(mat, start, nmarker);
	return re;

}

int Close_MatFile(int idx){
    
	int re;
	if(idx >= g_nPop){
		return ERORR_MAT_IDX_OUTOFBOUND;
	}
	re = g_pMatFile_Read[idx].Close();
	return re;
    
}


/*************************************************************
	MAT file for Save

*************************************************************/

int Mat_Init_Save(){
    
    int re = g_MatFile_Save.Init(); 
    return re;
}

int	Open_NewMatFile_Save(const char * FileName){

	int re = g_MatFile_Save.OpenToSave(FileName);
	return re;
}

int	Mat_PutData(double * mat, int size){

	int re = g_MatFile_Save.PutData(mat, size);
	return re;
}


int	Mat_Check_Saved(){

	int re = g_MatFile_Save.CheckSavedData();
	return re;
}

int	Mat_Num_Sets(){

	int re = g_MatFile_Save.GetNum_Sets();
	return re;
}

int	Mat_GetStart_Pos(int * pos, int * size){

	int re = g_MatFile_Save.GetStart_Pos(pos, size);
	return re;
}

int	Close_Write_MatFile(){
    
	int re = g_MatFile_Save.Close();
	return re;
}

/*************************************************************
	Bed file for read

*************************************************************/
int	Bed_Init(char* filename, int NSample, int NSnp){

	int re = g_BedFile.Init(filename, NSample, NSnp);
	return re;
}
int	Bed_ReadData(int * pIdxs, int num, unsigned char * Genotype){

	int re =g_BedFile.ReadData(pIdxs, num, Genotype);
	return re;
}

int	Bed_Close(){
    
	int re =g_BedFile.Close();
	return re;
}

