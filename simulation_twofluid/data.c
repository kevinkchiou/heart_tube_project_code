#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "definesim.h"
#include "threshold_wave_simulation.h"
#include "data.h"
#include "measure.h"

int print_grid_state(const char *FILENAME,const int idnum){
	int i;
	FILE *ofp;
	char FNAME[100];

	snprintf(FNAME,100,"%s_%.04d.output",FILENAME,idnum);
	ofp=fopen(FNAME,"w");

	fprintf(ofp,"%d %d\n",num_grid_pts,num_myo);
	for(i=0;i<num_grid_pts;i++){fprintf(ofp,"%lf ",gridweb->u[i]);}
	fprintf(ofp,"\n");
	for(i=0;i<num_myo;i++){fprintf(ofp,"%d ",myoweb[i]->gridloc);}
	fprintf(ofp,"\n");
	for(i=0;i<num_myo;i++){fprintf(ofp,"%d ",myoweb[i]->state);}

	fclose(ofp);
	return 0;
}

int import_grid_state(const char *FILENAME){
	int i,n1,n2;
	FILE *ifp;
	Paras *p;


	ifp=fopen(FILENAME,"r");
	fscanf(ifp,"%d %d",&n1,&n2);
	init_sys(n2,(int)(n1/n2),p,200.,0.,1.,20.);
	for(i=0;i<num_grid_pts;i++){fscanf(ifp,"%lf",&(gridweb->u[i]));gridweb->n[i]=-1;}
	for(i=0;i<num_myo;i++){fscanf(ifp,"%d",&(myoweb[i]->gridloc));gridweb->n[myoweb[i]->gridloc]=i;}
	for(i=0;i<num_myo;i++){fscanf(ifp,"%d",&(myoweb[i]->state));}

	fclose(ifp);
	return 0;
}

int print_diagnostic_information(const char *FILENAME,const int idnum){

	int i;
	double eps,sig;
	double *y;
	FILE *ofp;
	char FNAME[100];

	snprintf(FNAME,100,"%s_%.04d.diag",FILENAME,idnum);
	ofp=fopen(FNAME,"w");

	y=(double*)malloc(num_grid_pts*sizeof(double));
	for(i=0;i<num_grid_pts;i++){y[i]=gridweb->u[i];}

	fprintf(ofp,"%d\n",num_grid_pts);
	for(i=0;i<num_grid_pts;i++){fprintf(ofp,"%lf\n",strain1d(y,i));}

	free(y);
	fclose(ofp);
	return 0;
}
