#define MAXNAME 50
#define MAXPROTLEN  10000
#define AA     20
#define MAXSTR 1000
# define DELTA         4 /* chem dist excluding form contact counting */
# define RLOW   1.0
# define RHIGH  6.4
#define MAXATOMS 150000
#define MAXSEQS 5000
#define MAXCHAINS 20
#define MAXCONTACTS 30
/*#define DEBUG
#define DEBUG_SEQ
#define DEBUG_PDB
#define DEBUG_ENERGY*/

static int R_atoms[]={1,7,4,4,2,5,5,0,6,4,4,5,4,7,3,2,3,10,8,3};
static char aa_oneletter[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V' };

typedef int       ivectpmty[AA];
typedef double    dvectpmty[AA];
typedef dvectpmty dmattpmty[AA];
typedef short boolean;

void Memory(void);
void Open_files(void);
int ReadAlignment(FILE *ali_input);
int set_eijmat(float *eijmat);
int Readpdb(FILE *pdbin);
int Seq_to_number(char *code,char *res,int *resnum);
int Readq (FILE *qmutname, double *qmut);
int Calculate_energy(int *seqnum , double *energy);
int MutateSeq(int *seq,double *qmatrix,int *mut);
int MakeMut(double *qmatrix);
double Calculate_score(int *seqnumT, int *seqnumA);
int Accept_trial(double score,float beta);
int Number_to_seq(int *seqnum,char *seq);
int Setrate(double *qmatrix, double *freq);
int Calculate_freq(double *pmatrix,double *freq);
double Calculate_rate(double *qmatrix,double *freq);
int Setp_cumulative(double *pmut);
char run(int *A,int *T,int *nsyn);
char Evaluate(void);
int diag(double *qmatrix,double *pmatrix,double time);
int Init_vector(int *vector,int vecsize);
int Accumul_sust(int *seqA,int *seqT,int *Nsus_pos);
int Count(char ch,double *nonaccepted, double *accepted, double *nothing,int *mut,int *nsyn, double *N_nsyn, double *N_syn);
int Print_out(FILE *outS,FILE *outPS,FILE *outPSw,FILE *outSE,FILE *outSEw,FILE *outPP);
int read_cmat(int *cmat);
int matrix(int *Nsus_pos);


/*diagonalization routines*/
void elmhes(double  a[AA][AA], int v[AA] , int n);
void eltran(double a[AA][AA], double zz[AA][AA], int v[AA], int n);
void hqr2(int n, int low, int hgh, int *err, double h[AA][AA], double zz[AA][AA], double wr[AA], double wi[AA]);
static void cdiv(double ar,double ai,double br,double bi,double *cr,double *ci);
void luinverse(dmattpmty omtrx, dmattpmty imtrx, int size);
void instantRate(dmattpmty a,double *f);
void mtrev(dmattpmty r, double *f);
void SetMatrix(double *matrix, double len);

