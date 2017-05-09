#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "definesim.h"

Myonode **allocmyoweb(int num,int size){
	Myonode **a;
	int i;
	a=(Myonode**)malloc(num*sizeof(Myonode*));
	if(!a){printf("Not enough memory!\n");fflush(stdout);return NULL;}
	for(i=0;i<num;i++){a[i]=allocmyonode(size);}
	return a;
}

Myonode *allocmyonode(int size){
	Myonode *a;
	a=(Myonode*)malloc(sizeof(Myonode));
	if(!a){printf("Not enough memory!\n");fflush(stdout);return NULL;}
	a->state=0;a->tau=10000.0;a->t=1.0;
	a->size=size;
	return a;
}

void freemyonode(Myonode *a){free(a);a=NULL;}
void freemyoweb(Myonode **a,int num){
	int i;
	for(i=0;i<num;i++){freemyonode(a[i]);}
	free(a);
	a=NULL;
}

Grid *allocgrid(int num){
	double *temp;
	int *itemp,i;
	Grid *a;
	a=(Grid*)malloc(sizeof(Grid));
	if(!a){printf("Not enough memory!\n");fflush(stdout);return NULL;}
	a->num=num;
	temp=(double*)malloc(num*sizeof(double));
	if(!temp){printf("Not enough memory!\n");fflush(stdout);return NULL;}
	for(i=0;i<num;i++){temp[i]=0.;}
	a->u=temp;
	itemp=(int*)malloc(num*sizeof(int));
	if(!itemp){printf("Not enough memory!\n");fflush(stdout);return NULL;}
	a->n=itemp;
	return a;
}
void freegrid(Grid *pg){
	free(pg->u);
	free(pg->n);
	free(pg);
}

void parameter_input_file(){
	FILE *fp;
	char temp[100];
	char *p;
	int loc,run_num;

	fp=fopen("parameters.dat","r");
	if(fp==NULL){
		printf("cannot open parameters.dat! using defaults and creating file...\n");

	}
	else{ //read in parameters from ofp
		while(1){
			if(feof(fp)){break;}
			fgets(temp,100,fp); //gets the file line by line instead of space-separated
			//what to do with each string entry
			if(temp[0]=='#'){//ignore as comment
				puts("the following line is a comment and unused!");
				puts(temp);
			}
			else{
				if(strcmp(temp,"E")){//set stiffness in next fscanf
				}
				if(strcmp(temp,"E0")){

				}
				if(strcmp(temp,"num_grid_pts")){

				}
			}
		}
	}
}

/*
int extract_bounding_values(char str[],double *a,double *b){
	char *p,cat1[10],cat2[10];
	int i=0;

	p=(char*)strchr(str,'-');
	if(p!=NULL){
		strncpy(cat1,str,p-str+1);cat1[p-str+1]='\0';
		i=p-str+2; //for copying second half
		while(str[i]!='\0'){cat2[i-p+str-2]=str[i];i++;}
		cat2[i]='\0';
		*a=atof(cat1);*b=atof(cat2);
		if(*b < *a){printf("error in ordering! we go from %lf to %lf!\n",*a,*a);}
		return 0;
	}
	else{return 1;}
}
*/
