#define INCL_WIN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <os2.h>
#include <float.h>
#define STDCOL 5
#define XSLOT 0
#define YSLOT 1
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

#define NADD 0              
#define NPARAM  4

/* edit the following lines or add and delete new lines */ 
static struct info_line dll_info[] = { 
   "This function calculates a quadratic.",
   "It requires 4 parameters as follows:",
   "P1 is the shift along the x axis.",
   "P2 is the shift along the y value.",
   "P3 is the magnitude of the linear term.", 
   "P4 is the magnitude of the quadratic term",
   NULL};
   
VOID EXPENTRY fitfunc(void)
   {
   int i;
   double b, a0, a1, a2;
   double x;

   /* The next section is executed the first time through only */
   /* so put any initialization code in here                   */
   if (firsttime) 
      {
      firsttime = FALSE;
      /* Mask off all floating point exceptions */
      _control87(MCW_EM,MCW_EM);
      }

   /* give parameters more meaningful names  */
   b = *Ptr[1];
   a0 = *Ptr[2];
   a1 = *Ptr[3];
   a2 = *Ptr[4];

   /* This is the main loop that calculates the function datt[] */
   /* for all data points. Put your function in the loop body   */ 
   for (i = 0;i < Ndat;i++)
      {
      x = datx[i] - b;
      datt[i] = a0 + x * (a1 + x * a2);
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

