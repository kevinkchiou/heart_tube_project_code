#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parfile_tools.h"

int locate_whitespace(char *str,int *idx){
	int i=0,count=0;
	for(i=0;i<strlen(str);i++){if(str[i]==' '){count++;}}
	idx=(int*)malloc(count*sizeof(int));
	i=0;count=0;
	for(i=0;i<strlen(str);i++){if(str[i]==' '){idx[count]=i;count++;}}
	return count;
}

int trim_whitespace(char *out,char *in){
	int i=0,j=0;
	//not computationally efficient but is code efficient
	while(i<strlen(in)){
		if(in[i]==' '||in[i]=='\t'){out[j]=' ';i++;j++;;}
		else{out[j]=in[i];i++;j++;}
		if((out[j-1]==' ')&&(in[i]==' '||in[i]=='\t')){j--;}
		//want to return error but can't assess size inside function
		//if(j>=sizeof(out)/sizeof(out[0])){puts("Output array too short!\n");return 1;}
	}
	out[j]='\0';
	return 0; //success
}

char **split_string_whitespace(int *sz,char *in){
	int *idx=NULL,szidx=0,i=0;
	int j,count=0;
	char **out,temp[150];

	trim_whitespace(temp,in);
	szidx=locate_whitespace(temp,idx)+1;//one more word than whitespace
	out=(char**)malloc((szidx+1)*sizeof(char*));
	for(i=0;i<szidx;i++){out[i]=(char*)malloc(20*sizeof(char));}
	for(i=0;i<szidx;i++){
		j=0;
		while(1){
			out[i][j]=temp[count];count++;j++;
			if(temp[count]==' '||temp[count]=='\0'){out[i][j]='\0';count++;break;}
		}
	}
	if(idx)free(idx);
	*sz=szidx;
	return out;
}

void free_split_string(char **a,const unsigned int size){
	int i;
	for(i=0;i<size;i++){free(a[i]);}
	free(a);
}

