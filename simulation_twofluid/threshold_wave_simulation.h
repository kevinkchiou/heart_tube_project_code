int threshold_wave_sim(char *FILENAME,int nummyo,int gridspace,double E,double thresh,double gamma);
int evolfunc_wave(double t,const double y[], double f[],void *par);
int evoljac_wave(double t,const double y[],double *dfdy,double dfdt[],void *par);
int evoleuler_wave(const double dt,const Paras *parameters);
void evolgsl_wave(const double dt,const Paras *parameters);
int check_myo_state(double *gridvec,const Paras *p,const double dt);

int initiate_excitation(double *gridvec,const int myonum,const double excite_strength);
int falloff_function(double *gridvec,const int size,const double strength,const int type);
double modification_function(const int type,const int gridpt,const Paras *p);

int init_sys(const int nummyo,int gridspace,Paras *p,double E,double lin,double thresh,double gamma);
