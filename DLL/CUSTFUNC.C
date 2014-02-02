#define INCL_WIN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <os2.h>
#include <float.h>
// the standard number of data columns
#define STDCOL 5
#define XSLOT 0
#define YSLOT 1
#define ASLOT 2
#define CALCSLOT -3
#define BESTSLOT -2
#define FLAGSLOT -1

/* static data and parameters */
static double *data;                    /* start of complete data array  */

/* The data is all in a contigous block starting with datx then daty followed by
 * zero or more additional data columns as read in from the data file. The last
 * three slots are: datt , the current best results and the flag array. 
 * Each point not in an exclusion zone has its corresponding flag TRUE. 
 */


static double *datx;                    /* x data read in from data file */
static double *daty;                    /* y data read in from data file */
static double *datt;                    /* generated data */
static int Ndat;                        /* number of data points */
static char firsttime;
static double **Ptr;                    /* pointer to parameter pointer array */

struct info_line
   {
   char *line;
   };

/*****************************************************************************/
/* Do not modify any code not enclosed by this and the end banner            */
/*****************************************************************************/

/* the next two values are compared to values passed by the program as check */
/* NADD is the number of additional data columns required for this function  */
/* NPARAM is the number of parameters required by this function              */

#define NADD 4              
#define NPARAM 5

/* edit the following lines or add and delete new lines */ 
static struct info_line dll_info[] = { 
   "This function calculates the resistivity in GaAs",
   "to compare with measured rho values obtained fron Hall effect",
   "This function Uses Newton's method to solve a cubic.",
   "The first data column is the sample number",
   "The second column is the measured resistivity",
   "4 Additional data columns are read in from the data file", 
   "The first is the C concentration, the second the Zn conc",
   "the third the estimated donor concentration, last the mobility.",
   "It requires 5 parameters as follows:",
   "P1 is the EL2 density.",
   "P2 is the EL2 activation energy.",
   "P3 is the temperature.", 
   "P4 is the donor concentration multiplier",
   "P5 is the Zn acceptor concentration multiplier",
   NULL};
   
VOID EXPENTRY fitfunc(void)
   {
   int i;      
   double nets;
   double net;
   double Nc,Nv;
   double ni2;                  /* the intrinsic density squared */
   double n;                    /* the electron density */
   double n2;
   double Nel2;                 /* El2 concentration */
   double A,D,a,b,c,fn,dfn;     /* temp values */
   double ncexp;
   double T,Eg,Ea;              /* temperature, Gap and EL2 activation energy */
   double toler = 1e-8;
   double nseed = 1e19;
   double K = 8.617335e-5;      /* the Boltzmann constant in eV/K */
   double Dgx,Dgl;              /* the gamma,x and gamma,l separations in eV */
   double DS,ZS;                /* donor and zinc intensity divisors */

   /* set pointers to additional read in data columns (if any) */
   double *dat1,*dat2,*dat3,*dat4;
   dat1 = datx + Ndat * (ASLOT + 0); 
   dat2 = datx + Ndat * (ASLOT + 1); 
   dat3 = datx + Ndat * (ASLOT + 2); 
   dat4 = datx + Ndat * (ASLOT + 3); 

   /* The next section is executed the first time through only */
   /* so put any initialization code in here                   */
   if (firsttime) 
      {
      firsttime = FALSE;
      /* Mask off all floating point exceptions */
      _control87(MCW_EM,MCW_EM);
      }

   /* give parameters more meaningful names  */
   Nel2 = *Ptr[1];              /* The El2 concentration */
   Ea = *Ptr[2];                /* The EL2 activation energy */
   T = *Ptr[3];                 /* The temperature */
   DS = *Ptr[4];                /* divide measured donor intensity by this */
   ZS = *Ptr[5];                /* divide measured Zn intensity by this */

   /* fix up the activation energy to acount for temperature variation */
   Ea -= 2.37e-4 * T;
   /* calculate the band gap and gamma,x and gamma,l separations*/
   Eg = 1.519 - 5.405e-4 * T*T / (T + 204);
   Dgx = 0.462 + 8.05e-5 * T*T / (T + 204);
   Dgl = 0.296 - 6.45e-5 * T*T / (T + 204);
   /* calculate the effective conduction band density of states Blakemore R164 */
   Nc = 8.63e13 * sqrt(T*T*T) * ((1 - 1.93e-4 * T - 4.19e-8 * T*T) +
		21 * exp(-Dgx / K / T) + 44 * exp(-Dgl / K / T));
   Nv = 1.83e15 * sqrt(T*T*T);
   ncexp = Nc * exp(-Ea / K / T);
   /* calculate the intrinsic electron density squared */
   ni2 = Nc * Nv * exp(-Eg / K / T);
   c = ni2 * ncexp;

   /* This is the main loop that calculates the function datt[] */
   /* for all data points. Put your function in the loop body   */ 
   for (i = 0;i < Ndat;i++)
      { 
      n = nseed;
      A = dat1[i] + dat2[i] / ZS;             /* total acceptors */
      D = dat3[i] * A / DS;                   /* shallow donors */
      nets = A - D;
      Nel2 = Nel2 + nets;                     /* total EL2 concentration */
      net = nets - Nel2;
      a = nets - ncexp;
      b = net * ncexp - ni2;
      while (TRUE)
	 {
	 n2 = n * n;
	 fn = n * n2 + a * n2 + b * n - c;
	 dfn = 3 * n2 + 2 * a * n + b;
	 if ((fn / dfn / n) < toler) break;
	 n = n - fn / dfn;
	 }
      datt[i] = n;                 /*the electron density */
      }
   }

/*****************************************************************************/
/* Do not modify any code not enclosed by this and the start banner          */
/*****************************************************************************/



char EXPENTRY init_dll_globals(int N,int mb,int NP,double *dat,double *p[])
   {
   /* set globals */
   firsttime = TRUE;
   Ndat = N;
   if (NADD != mb - STDCOL) return -1;
   if (NPARAM != NP) return -2;
   /* set array pointers */
   data = dat;
   daty = dat + N;
   datx = dat;
   datt = dat + (mb + CALCSLOT) * N;     /* the calculated data goes into the third last slot */
   /* set parameter pointers */
   Ptr = p;
   return 0;
   }

struct info_line * EXPENTRY get_dll_info(void)
   {
   /* return pointer to array of info strings set up above */
   return dll_info; 
   }

