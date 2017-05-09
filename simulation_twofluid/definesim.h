typedef struct myonode
{
	int size;
	int gridloc;
	int state;//0 for primed, 1 for contracting, 2 for refractory
	double tau; //refractory time
	double t; //current time

} Myonode;

typedef struct grid
{
	int num;
	double *u;
	int *n;
} Grid;

typedef struct paras
{
	double E;
	double lin;
	double thresh;
	double gamma;
} Paras;

Myonode **allocmyoweb(int num,int size);
Myonode *allocmyonode(int size);
void freemyonode(Myonode *a);
void freemyoweb(Myonode **a,int num);
Grid *allocgrid(int num);
void freegrid(Grid *gridpointer);

void parameter_input_file();
//int extract_bounding_values(char str[],double *a,double *b);

extern Myonode **myoweb;
extern Grid *gridweb;

extern double E0;
extern int num_grid_pts;
extern int num_myo;

extern int bctype;
extern int modfcntype;
extern Paras *params;
