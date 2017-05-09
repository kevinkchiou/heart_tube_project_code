double strain1d(const double *gridvec,const int gridpt);
double stress1d(const double *gridvec,const int gridpt);
int apply_boundary_conditions(const double *gridvec,const int gridpt,const int type,double *yp,double *yn);
int stiffness_params(const int gridpoint,const Paras *p,double *E,double *dE);
int calc_avgvel(const double newvel,double *avgvel,int *velct);
double calc_velocity(int *prev_gridpt,const double dt);
int find_max_loc(const double *gridvec);

double calc_vel_activation(int *prev_myonum,const double dt);
int find_active_location(char *dir);
