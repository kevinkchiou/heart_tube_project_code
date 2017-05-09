#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "definesim.h"
#include "threshold_wave_simulation.h"
#include "data.h"
#include "measure.h"

Myonode **myoweb;
Grid *gridweb;
double E0;
int num_grid_pts;
int num_myo;
int bctype=0; //0 - free boundary. 1 - fixed boundary
int modfcntype=3; //0 - no modification. 1 - sigmoid. 2 - thetafcn, 3 - double sigmoid
Paras *params;

int main(int argc, char *argv[]){
	double dE=5.,E,Eminfrac,Emaxfrac;
	char FILENAME[100],FNAME[100];
	char PFNAME[100]; //parameter filename, default is "parameters.dat"
	FILE *ofp;

	if(argc>1){strcpy(PFNAME,argv[1]);}
	else{strcpy(PFNAME,"parameters.dat");}

	if(argc<2){sprintf(FNAME,"output/default");}
	else{sprintf(FNAME,"output/%s",argv[1]);}
	if(argc<3){Eminfrac=0.1;}
	else{Eminfrac=atof(argv[2]);}
	if(argc<4){Emaxfrac=3.1;}
	else{Emaxfrac=atof(argv[3]);}
	if(argc<5){E0=100.;}
	else{E0=atof(argv[4]);}

	if(dE<4.9){dE=0.1*E0;}
	E=Eminfrac*E0;
	printf("E0=%f,Eminfrac=%f,Emaxfrac=%f\n",E0,Eminfrac,Emaxfrac);
	if(modfcntype==0){ofp=fopen("avgvel_nomod.dat","w");}
	else{ofp=fopen("avgvel_mod3.dat","w");}
	while(E<Emaxfrac*E0){
		sprintf(FILENAME,"%s_%.04f",FNAME,E);
		threshold_wave_sim(FILENAME,300,9,E,0.1,2.);
		//threshold_wave_sim(FILENAME,nummyo,gridspace,E,thresh,gamma);
		E+=dE;
	}
	fclose(ofp);

	return 0;
}

int threshold_wave_sim(char *FILENAME,int nummyo,int gridspace,double E,double thresh,double gamma){

	double r,vel;
	int i,j=0;
	int prevloc=0; //used for velocity calculation
	int length;
	Paras *p;
	FILE *ofp;

	//initialize
	p=(Paras*)malloc(sizeof(Paras));
	//global version just in case
	params=p; //global pointer to same allocated structure
	//num_myo = number of myocytes, grid is # of grid pts between myocytes
	init_sys(nummyo,gridspace,p,E,0.,thresh,gamma); //lin=0.

	//run stuff
	double t=0.,tf,dt;
	int count=0;

	//avg velocity stuff
	double avgvel=0.;
	int velct=0;
	int velflag=0;

	dt=0.01; //evolution time step
	tf=10.; //final time
	//print_grid_state(FILENAME,count);
	//print_diagnostic_information(FILENAME,count);
	count++;
	while(t<tf){
		evolgsl_wave(dt,p); //gsl integration
		//evoleuler_wave(dt,p); //stupid euler integrator

		//if(vel>0.01 || vel<0.0){printf("newvel = %lf\n",vel);fflush(stdout);}

		//uncomment one of the below to pick velocity method (grid vs myocyte)
		//difference seems to be <1%
		vel=calc_velocity(&prevloc,dt);length=num_grid_pts;
		//vel=calc_vel_activation(&prevloc,dt);length=num_myo;

		if(prevloc<0.95*length && prevloc>0.05*length && velflag==0){calc_avgvel(vel,&avgvel,&velct);}
		else if(prevloc>0.95*length){velflag=1;}//this kills velocity tracking

		//print_grid_state(FILENAME,count);
		//print_diagnostic_information(FILENAME,count);
		//output between calls to evol();
		t+=dt;count++;
	}
	printf("average velocity for E=%f: %f\n",E,avgvel);
	if(modfcntype==0){ofp=fopen("avgvel_nomod.dat","a");}
	else{ofp=fopen("avgvel_mod3.dat","a");}
	fprintf(ofp,"%f %f\n",E/E0,avgvel);
	fclose(ofp);

	//free up the data structures
	freemyoweb(myoweb,num_myo);
	freegrid(gridweb);
	free(p);

	return 0;
}

int evolfunc_wave(double t,const double y[],double f[],void *par){
	double yp,yn;
	double r,dr,test;
	int i,idx,sz;
	Paras *p;

	p=(Paras*)par;

	dr=p->lin;
	sz=myoweb[0]->size; //scale laplacian to cell size
	//need to use par structure to extract stiffness
	for(i=0;i<num_grid_pts;i++){
		r=p->E - i*(p->lin);
		apply_boundary_conditions(y,i,bctype,&yp,&yn);//gives yp and yn with boundary conditions
		f[i]=r*sz*sz*(yp-2.*y[i]+yn); //laplacian
		f[i]+=dr*sz*sz*(yn-yp)/2.; //drift
	}
	return 0;
}

//jacobian term - not necessary for explicit integration
int evoljac_wave(double t,const double y[],double *dfdy,double dfdt[],void *par){
}

int evoleuler_wave(const double Delta_t,const Paras *p){
	double *y,*f,dt,t,df;
	int i;

	y=(double*)malloc(num_grid_pts*sizeof(double));
	f=(double*)malloc(num_grid_pts*sizeof(double));
	for(i=0;i<num_grid_pts;i++){y[i]=gridweb->u[i];}
	t=0;dt=1e-4;
	while(t<Delta_t){
		if(evolfunc_wave(t,y,f,(void*)p)!=0){printf("Failure: evolfunc_wave()!\n");}//passive f[i]
		//this next line does it via delta-function modification. Thus not integrated against dt
		//for(i=0;i<num_myo;i++){if(myoweb[i]->state==1){initiate_excitation(y,i,p->gamma);}} //active f[i]
		check_myo_state(y,p,dt);//checks state of myocytes
		for(i=0;i<num_grid_pts;i++){
			//integration stability
			df=f[i]*dt;
			if(df>0.05*y[i] && y[i]>0.001){printf("df=%lf, y[%d]*dt=%lf\n",df,i,y[i]);}
			y[i]+=df;
		}
		t+=dt;
	}
	for(i=0;i<num_grid_pts;i++){gridweb->u[i]=y[i];}

	free(y);free(f);
	return 0;
}

void evolgsl_wave(const double Delta_t,const Paras *p){
	int i;
	double *y;

	y=(double*)calloc(num_grid_pts,sizeof(double));

	//initialize p and y here

	const gsl_odeiv2_step_type *type_bsimp=gsl_odeiv2_step_bsimp;
	const gsl_odeiv2_step_type *type_rkf45=gsl_odeiv2_step_rkf45;
	const gsl_odeiv2_step_type *type_rk8pd=gsl_odeiv2_step_rk8pd;

	gsl_odeiv2_step *s=gsl_odeiv2_step_alloc(type_rk8pd,num_grid_pts);
	gsl_odeiv2_control *c=gsl_odeiv2_control_yp_new(1e-4,1e-4);
	gsl_odeiv2_evolve *e=gsl_odeiv2_evolve_alloc(num_grid_pts);

	gsl_odeiv2_system sys= {evolfunc_wave,evoljac_wave,num_grid_pts,(void*)p};

	double t=0.;
	double h=1e-4;

	int status,count=0,t0;
	while(t<Delta_t){
		t0=t; //initialize before applying evolution (cell clock tracking)
		for(i=0;i<num_grid_pts;i++){y[i]=gridweb->u[i];} //initialize vector
		check_myo_state(y,p,t-t0);
		status=gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,Delta_t,&h,y);
		for(i=0;i<num_grid_pts;i++){gridweb->u[i]=y[i];}
		//consider using regular eulerian integrator here.
		//for(i=0;i<num_myo;i++){if(myoweb[i]->state==1){initiate_excitation(y,i,p->gamma);}}
		if(status!=GSL_SUCCESS){printf("failure: GSL_SUCCESS==false!\n");fflush(stdout);break;}
		count++;
	}

	for(i=0;i<num_grid_pts;i++){gridweb->u[i]=y[i];}//output to global structure

	gsl_odeiv2_evolve_free(e);
	gsl_odeiv2_control_free(c);
	gsl_odeiv2_step_free(s);
	free(y);
}

int check_myo_state(double *gridvec,const Paras *p,const double dt){
	int idx,state,i;
	double str,meas;
	Myonode *pm;

	for(i=0;i<num_myo;i++){
		pm=myoweb[i];
		pm->t-=dt; //update clock
		state=pm->state;
		idx=pm->gridloc;
		if(state==0){
			//stress and strain thresholds
			//str=modification_function(modfcntype,idx,p)*stress1d(gridvec,idx,p);
			str=modification_function(modfcntype,idx,p)*strain1d(gridvec,idx);

			meas=str;
			if(meas > p->thresh){pm->state=1;}
		}
		//else if(state==1 && pm->t<0.){pm->t=pm->tau;pm->state=2;}
		if(state==1){initiate_excitation(gridvec,i,p->gamma);pm->t=pm->tau;pm->state=2;}
		if(state==2 && pm->t<0.){pm->state=0;pm->t=0.1;}
	}
	return 0;
}

int initiate_excitation(double *y,const int myonum,const double gamma){
	Myonode *pm;
	int size,i,loc,type;
	double a,*z;

	type = 2;

	pm=myoweb[myonum];
	//printf("myocyte %d is firing!\n",myonum);
	size=((int)(0.5*pm->size));loc=pm->gridloc;
	z=(double*)malloc(pm->size*sizeof(double));
	//technically is Q/Gamma due to time derivative rescaling
	if(falloff_function(z,size,gamma,type)==1){printf("error in falloff_function()!\n");}
	else{for(i=-size;i<=size;i++){y[loc+i]+=modification_function(modfcntype,pm->gridloc,params)*z[i+size];}}
	free(z);
	return 0;
}

int falloff_function(double *y,const int size,const double gamma,const int type){
	int i,x;
	double norm,a,fx;
	
	norm=1;
	switch(type){
		case 0: //gaussian
			for(i=0;i<size;i++){norm+=2*exp(-0.5*(i+1)*(i+1));}
			a=gamma/norm;
			y[size]=a;
			for(i=0;i<size;i++){fx=exp(0.5*(-i-1)*(i+1));y[size+(i+1)]=a*fx;y[size-(i+1)]=a*fx;}
			break;

		case 1: //exponential
			for(i=0;i<size;i++){norm+=2*exp(-i-1);}
			a=gamma/norm;
			y[size]=a;
			for(i=0;i<size;i++){fx=exp(-i-1);y[size-(i+1)]=a*fx;y[size+(i+1)]=a*fx;}
			break;

		case 2: //delta
			y[size]=gamma*(2*size+1); //proper scaling for gamma
			break;

		case 3: //delta-cell-size
			a=gamma/((double)size+1.);
			y[size]=a;
			for(i=0;i<size;i++){y[size-(i+1)]=a;y[size+(i+1)]=a;}
			break;

		default: //gaussian
			printf("unrecognized case - using default: gaussian");
			for(i=0;i<size;i++){norm+=2*exp(-0.5*(i+1)*(i+1));}
			a=gamma/norm;
			y[size]+=a;
			for(i=0;i<size;i++){fx=exp((-i-1)*(i+1)/2.);y[size-(i+1)]=a*fx;y[size+(i+1)]=a*fx;}
			break;
	}
	return 0;
}

double modification_function(const int type,const int gridpt,const Paras *p){

	double a=0.75,b=0.25,delta=1.8; //max,min,shift location
	double k=1.5; //falloff factor for sigmoidal case
	double E;
	double E2,d2=0.3,k2=2.,f,g;

	local_stiffness(gridpt,p,&E);
	E/=E0; //in relation to native stiffness
	E2=E;

	switch(type){
		case 0: //no modification
			return 1.0;

		case 1: //sigmoid
			E-=delta;
			//divide by things closer to one
			if(E<0){return (a+b*exp(2.*k*E))/(1+exp(2.*k*E));}
			else{return (a*exp(-2.*k*E)+b)/(exp(-2.*k*E)+1);}
			
		case 2: //theta-fcn
			if(E>delta){return a;}
			else{return b;}

		case 3: //double sigmoid
			E2-=d2;
			E-=delta;
			if(E<0){f=(a+b*exp(2.*k*E))/(1.+exp(2.*k*E));}
			else{f=(a*exp(-2.*k*E)+b)/(exp(-2.*k*E)+1.);}
			if(E2<0){g=exp(2.*k2*E2)/(1.+exp(2.*k2*E2));}
			else{g=1./(exp(-2.*k2*E2)+1.);}
			return f*g;
	}
}

int init_sys(const int nummyo,int gridspace,Paras *p,double E,double lin,double thresh,double gamma){
	int i,j=0;

	num_myo=nummyo;
	if(gridspace%2==0){printf("gridspace should be odd, changing to %d!\n",gridspace+1);gridspace++;}
	num_grid_pts=num_myo*gridspace+2;//extra grid point on each side of end of lattice
	//printf("num_myo = %d, num_grid_pts = %d, gridspace = %d!\n",num_myo,num_grid_pts,gridspace);
	gridweb=allocgrid(num_grid_pts);
	myoweb=allocmyoweb(num_myo,gridspace);
	//printf("our modulo shift is %d...\n",((int)(gridspace/2)));fflush(stdout);
	for(i=1;i<num_grid_pts-1;i++){ //initialize
		if(j>num_myo){printf("some sort of error in counting!\n");fflush(stdout);}
		if((i-((int)(gridspace/2))-1)%gridspace==0){gridweb->n[i]=j;myoweb[j]->gridloc=i;j++;}
		else{gridweb->n[i]=-1;}
	}
	p->E = E; //avg stiffness for spatially linear stiffness dependence
	p->lin = lin;
	//p->E = 1.3;
	//p->lin = -2.0/1.7*(p->E)*0.3/num_grid_pts;
	p->thresh = thresh;
	p->gamma = gamma;

	//pacemaker(s)
	//myoweb[((int)(num_myo/2))]->state=1;
	myoweb[0]->state=1;
	//myoweb[num_myo-1]->state=1;
	//for(i=0;i<((int)(num_myo/20));i++){myoweb[i]->state=1;} //square wave
	//instead of doing a "physical" pacemaker we just initialize a gaussian curve
	/*
	int end=gridspace*5,sig=gridspace,mid;
	mid=((int)(end/2));
	for(i=0;i<end;i++){gridweb->u[i]=gridspace*p->gamma*exp((-i+mid)*(i-mid)/(2*sig*sig));}
	//put the ones under the curve into refractory (may be unnecessary)
	//for(i=0;i<end;i++){if(myoweb[i]->gridloc<=end){myoweb[i]->state=2;myoweb[i]->t=1.0;}}
	*/
	return 0;
}
