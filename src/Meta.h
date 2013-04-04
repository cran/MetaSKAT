
#ifndef _META_H_
#define _META_H_


/*************************************************************
	MAT file for read


*************************************************************/
int	Init_MatFile_Read(int NumPop);
int	Open_NewMatFile_Read(int idx, const char * FileName);
int	Mat_GetData(int idx, double * mat, int start, int nmarker);
int Close_MatFile(int idx);

/*************************************************************
	MAT file for save

*************************************************************/

int Mat_Init_Save();
int	Open_NewMatFile_Save(const char * FileName);
int	Mat_PutData(double * mat, int size);
int	Mat_Check_Saved();
int	Mat_Num_Sets();
int	Mat_GetStart_Pos(int * pos, int * size);
int	Close_Write_MatFile();

/*************************************************************
	Bed file for read

*************************************************************/
int	Bed_Init(char* filename, int NSample, int NSnp);
int	Bed_ReadData(int * pIdxs, int num, unsigned char * Genotype);
int Bed_Close();

/***************************
 etc
 **************************/

int endian();

#endif


