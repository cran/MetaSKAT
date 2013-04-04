#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>


#include "MatFile.h"
#include "error_messages.h"


uint32_t xcrc32 (const unsigned char *buf, int len);



MatFile::MatFile(){

}

MatFile::~MatFile(){

	Close();
}

int MatFile::Init(){
    
    Close();
    m_save_startpos.clear();
    m_file_read.clear();
    m_file_save.clear();
    
    return 0;
    
}


int MatFile::OpenToSave(const char * FileName){

	char MagicNum[2];
	MagicNum[0] = MAGIC_NUM_MAT;
	MagicNum[1] = NULL;

    /* Close the opened file first*/
    Close();
    
	m_filename_save = FileName;
	m_file_save.open(m_filename_save.c_str(), std::ios::binary);

	if (!this->m_file_save)
	{
		return ERORR_MAT_OPEN_FILE4SAVE;
	}

	m_file_save. write ( MagicNum , 1 );
	m_save_startpos.push_back(1);
	return 0;

}


int MatFile::OpenToRead(const char * FileName){

	int re;
	char MagicNum[2];
	MagicNum[0] = MAGIC_NUM_MAT;
	MagicNum[1] = NULL;

    /* Close the opened file first*/
    Close();
    
	m_filename_read = FileName;
	m_file_read.open(m_filename_read.c_str(), std::ios::binary);

	if (!this->m_file_read)
	{
		return ERORR_MAT_OPEN_FILE4READ;
	}
	re = Check();
	return re;


}

int 	MatFile::Close(){

	if(m_file_read.is_open()){
		m_file_read.close();
	}
	if(m_file_save.is_open()){
		m_file_save.close();
	}

	return 0;
}


int 	MatFile::Check(){

	// 1 bit magic number
	m_file_read.read(m_magic_number,1 * sizeof(char)); // first byte

	if(m_magic_number[0] != MAGIC_NUM_MAT){
		return ERORR_MAT_MAGIC_NUM;
	}

	return NO_ERRORS;

}

/******************************************************
	PutData to the file
		one crc + matrix ( size=n *(n+1)/2)
******************************************************/
int 	MatFile::PutData(double * mat, int size){

	int len,i;
	uint32_t crc;
	char crc_c[4];
	if(size > MAX_SIZE_MAT){
		return ERORR_MAT_MAX_SIZE;
	}
	if(!m_file_save.is_open()){
		return ERORR_MAT_FILEOPEN;
	}
    
    /* convert double to float to save space */
	len=sizeof(float) * size;
    float * temp = (float *)m_buffer_crc;
    for(i=0;i<size;i++){
        temp[i] = (float)mat[i];
    }


	crc = xcrc32 ( ( unsigned char*)m_buffer_crc, len);
	memcpy(crc_c, &crc, 4);
    //printf("CRC1 :%d\n",crc);
    //printf("len:%d\n",len);

	m_file_save.write(crc_c, sizeof(uint32_t));
	m_file_save.write((char *)m_buffer_crc, len);

	/* add the starting position of the last vector */
	int last = m_save_startpos.back();
	m_save_startpos.push_back(last+len +sizeof(uint32_t));
	return NO_ERRORS;

}



/******************************************************
 Get to the file
 one crc + matrix ( size=n *(n+1)/2)
 mat is a vector of n*n
 n: number of marker
 ******************************************************/

int 	MatFile::GetData(double * mat, int start, int nmarker){

	int len, size;
	char buff1[10];
	uint32_t crc1, crc2;

    size =   (nmarker * (nmarker +1))/2; 
    len=size * sizeof(float);
    //Rprintf("len [%d], size[%d], nstart[%d], mat[%p]\n", len,size, start, mat);

	if(len > MAX_SIZE_MAT){
		return ERORR_MAT_MAX_SIZE;
	}
	if(!m_file_read.is_open()){
		return ERORR_MAT_FILEOPEN;
	}

	
	m_file_read.seekg(start ,std::ios::beg);
	m_file_read.read(buff1,sizeof(unsigned int)); // first 4 byte

	memcpy(&crc1, buff1, sizeof(unsigned int));

	m_file_read.read(m_buffer_crc, len);
	crc2 = xcrc32 (( unsigned char*) m_buffer_crc, len);

    //Rprintf("Get-CRC1 :%d\n",crc1);
    //Rprintf("Get-CRC2 :%d\n",crc2);
    //Rprintf("Get-len:%d\n",len);
	/* Check CRC */
	if(crc1 != crc2){
		
        Rprintf("CRC1 :%d\n",crc1);
        Rprintf("CRC2 :%d\n",crc2);
        return ERORR_MAT_CRC;
        
	}
 
    //Rprintf("Start Loop\n");


    int id1, id2, k, i,j;
    k=0;
    float *temp = (float *) m_buffer_crc;
   
    for(i=0;i< nmarker; i++){
        for(j=i;j<nmarker;j++){
            id1=i * nmarker +j;
            id2=j * nmarker +i;
            mat[id1] = mat[id2] = (double)temp[k];
            k++;
            
        }
    }

    //Rprintf("End Loop [%d][%d][%d][%f][%f][%p]\n",k,id1,id2, mat[id1], mat[id2], mat);
	return NO_ERRORS;

}


int 	MatFile::CheckCRC(int start, int size){
    
	int len;
	char buff1[10];
	uint32_t crc1, crc2;
    
    len=size * sizeof(float);

	if(len > MAX_SIZE_MAT){
		return ERORR_MAT_MAX_SIZE;
	}
	if(!m_file_read.is_open()){
		return ERORR_MAT_FILEOPEN;
	}
    
	
	m_file_read.seekg(start ,std::ios::beg);
	m_file_read.read(buff1,sizeof(unsigned int)); // first 4 byte
    
	memcpy(&crc1, buff1, sizeof(unsigned int));
    
	m_file_read.read(m_buffer_crc, len);
	crc2 = xcrc32 (( unsigned char*) m_buffer_crc, len);
    
    //Rprintf("Get-CRC1 :%d\n",crc1);
    //Rprintf("Get-CRC2 :%d\n",crc2);
    //Rprintf("Get-len:%d\n",len);
	/* Check CRC */
	if(crc1 != crc2){
		return ERORR_MAT_CRC;
	}
	return NO_ERRORS;
    
}





/*************************************************
	Use this file after close all files
**************************************************/
int 	MatFile::CheckSavedData(){


	Close();

	OpenToRead(m_filename_save.c_str());


	// Run
	int i, *ppos, *psize, re, nmarker;
	int NumSet = GetNum_Sets();
	ppos = new int[NumSet];
	psize = new int[NumSet];
	GetStart_Pos(ppos, psize);

	for(i=0;i<NumSet;i++){
        
        re = CheckCRC(ppos[i], psize[i]);
		if(re != NO_ERRORS){
			delete [] ppos ;
			delete [] psize ;
			Close();
			return re;
		}
	}

    delete [] ppos ;
    delete [] psize ;
	Close();
	return NO_ERRORS;
}


int 	MatFile::GetNum_Sets(){

	int len = m_save_startpos.size() -1;
	if(len < 0){
		len=0;
	}
	return len;

}

int 	MatFile::GetStart_Pos(int * pos, int * size){

	int i;
	int len = m_save_startpos.size();
	for(i=0;i<len-1;i++){
		pos[i] = m_save_startpos[i];
		size[i] = (m_save_startpos[i+1] - m_save_startpos[i])/4-1;
	}
	return 0;
}

