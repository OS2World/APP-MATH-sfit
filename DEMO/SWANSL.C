#define INCL_WIN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <os2.h>
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

#define NADD 2              
#define NPARAM 6

#define PI 3.14159

/* edit the following lines or add and delete new lines */ 
static struct info_line dll_info[] = { 
   "This function calculates magnetization ",
   "of the liquid pool in magnetization transfer.",
   "P1 is the solid spin relaxation time",
   "P2 is the solid transverse relaxation time",
   "P3 is the exchange rate",
   "P4 is RM = RM0b/Ra",
   "P5 is RT = 1/RaT2a",
   "P6 is solid dipolar relaxation time",
   NULL};
   
VOID EXPENTRY fitfunc(void)
   {
   int i, j;
   double w1, D, Rrfb, num, den1, den2, Rb, T2b, R, RM, RT, T1d;
   double Ds, T1b, M2, y, g, r, M2r;

   /* set pointers to additional read in data columns (if any) */
   double *dat1,*dat2;
   dat1 = datx + Ndat * (ASLOT + 0); 
   dat2 = datx + Ndat * (ASLOT + 1); 

   /* The next section is executed the first time through only */
   /* so put any initialization code in here                   */
   if (firsttime) 
      {
      firsttime = FALSE;
      }

   /* give parameters more meaningful names  */
   Rb   = *Ptr[1];          
   T2b  = *Ptr[2];
   R    = *Ptr[3];
   RM   = *Ptr[4];          
   RT   = *Ptr[5];          
   T1d  = *Ptr[6];

   Ds  = 1/T2b;
   T1b = 1/Rb;
   M2  = 1/(T2b*T2b);
   
   for (i = 0; i < Ndat;i++)
      {
      g = 0;
      w1 = 2.0*PI * dat1[i];
      D = 1000.0*2.0*PI*datx[i];
      for (j=1;j<=91;j++)
	 {
	 r = (j-1)* PI/180.0;
	 M2r = M2 * (3.0*cos(r)*cos(r)-1.0) * (3.0*cos(r)*cos(r)-1) / 4.0;
	 y = sin(r) * 1/sqrt(2*PI*M2r) * exp(-D*D/(2.0*M2r));
	 g = g + y * PI /180.0;
	 }
      Rrfb = w1*w1*PI*g;
      num = (Ds*Ds+D*D*Rrfb*T1d)*(1+RM+R*T1b) + Ds*Ds*Rrfb*T1b;
      den1 = (Ds*Ds+D*D*Rrfb*T1d)*(1+RM+R*T1b+(w1/D)*(w1/D)*RT*(1+T1b*R));
      den2 = Ds*Ds*Rrfb*T1b*(1+(w1/D)*(w1/D)*RT+RM);
      datt[i] = num/(den1+den2); 
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

