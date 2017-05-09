#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "definesim.h"
#include "measure.h"

double strain1d(const double *gridvec,const int gridpt){
	double yp,yn;
	if(gridpt>num_grid_pts || gridpt<0){printf("strain(): indexing error!\n");return 0;}
	
	apply_boundary_conditions(gridvec,gridpt,bctype,&yp,&yn);
	//midpt definition, derivative normalized to cell size
	return (fabs(yn-gridvec[gridpt])+fabs(gridvec[gridpt]-yp))*(myoweb[0]->size)/2.0;
	//right definition - probably better with left delta-fcns
	//return (fabs(yn-gridvec[gridpt]));
}

double stress1d(const double *gridvec,const int gridpt){
	double E,lin,str;

	E=params->E;lin=params->lin;
	str=strain1d(gridvec,gridpt);
	//if(str>0.0001){printf("strain = %lf\n",strain1d(gridvec,gridpt));}
	return (E-gridpt*lin)*str;
}

int apply_boundary_conditions(const double *gridvec,const int gridpt,const int type,double *yp,double *yn){

	switch(type){
		case 0: //free boundary
			if(gridpt-1<0){*yp=gridvec[gridpt];}
			else{*yp=gridvec[gridpt-1];}
			if(gridpt+1>num_grid_pts-1){*yn=gridvec[gridpt];}
			else{*yn=gridvec[gridpt+1];}
			break;

		case 1: //fixed boundary
			if(gridpt-1<0){*yp=0.;}
			else{*yp=gridvec[gridpt-1];}
			if(gridpt+1>num_grid_pts-1){*yn=0.;}
			else{*yn=gridvec[gridpt+1];}
			break;

		default:
			printf("unrecognized case - reverting to default (free boundary)!\n");
			if(gridpt-1<0){*yp=gridvec[gridpt];}
			else{*yp=gridvec[gridpt-1];}
			if(gridpt+1>num_grid_pts-1){*yn=gridvec[gridpt];}
			else{*yn=gridvec[gridpt+1];}
			break;
	}
	return 0; //return success
}

int local_stiffness(const int gridpt,const Paras *p,double *E){

	double prefac=2.0/1.7; //for p->E = avg E with linear dependence
	*E = prefac*p->E + ((double)gridpt)*p->lin;
	*E = p->E; //for now before we start varying things
	return 0;
}

int calc_avgvel(const double newvel,double *avgvel,int *velct){
	double av;
	int vc;

	av=*avgvel;vc=*velct;
	*avgvel = (av*((double)vc) + newvel)/((double)(vc+1));
	*velct=vc+1;
	return 0;
}

double calc_velocity(int *prevloc,const double dt){
	int i,sz;
	double *y,v;
	double myosize=0.01; //size in millimeters (10 um)

	sz=myoweb[0]->size;
	y=(double*)malloc(num_grid_pts*sizeof(double));
	for(i=0;i<num_grid_pts;i++){y[i]=gridweb->u[i];}
	i=find_max_loc(y);
	v=myosize*((double)(i - *prevloc) / (((double)sz)*dt)); //normalized to myocyte size
	*prevloc=i;
	free(y);
	return v;
}

int find_max_loc(const double *gridvec){
	int i,j=0;
	double ms=0.,cs=-1.;
	double (*compFcn)(const double*,const int); //comparison function - stress or strain
	int choice=1;

	if(choice==0){compFcn=&strain1d;}
	else{compFcn=&stress1d;}
	for(i=1;i<num_grid_pts-1;i++){
		cs=(*compFcn)(gridvec,i);
		if(cs>ms){ms=cs;j=i;}
	}
	return j;
}

double calc_vel_activation(int *p_myo,const double dt){
	double v;
	int prev,curr;
	double myosize=0.01; //size in millimeters (10um)

	prev=*p_myo;curr=find_active_location("right");
	v=myosize*((double)(curr-prev))/dt; //normalized to myocyte size

	*p_myo=curr;
	return v;
}

int find_active_location(char *dir){
	int i,prev,next,flag=0;
	int retval=-1;

	if(strcmp(dir,"right")!=0 && strcmp(dir,"left")!=0){puts("Error - default to right");strcpy(dir,"right");}
	for(i=0;i<num_myo-1;i++){
		if(strcmp(dir,"right")==0){prev=i;next=i+1;}
		else{prev=i+1;next=i;}

		if(myoweb[prev]->state!=0 && myoweb[next]->state==0){flag++;retval=prev;}
	}
	if(flag>1){puts("find_active_location(): more than one location!");}
	return retval;
}
