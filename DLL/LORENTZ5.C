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
#define NPARAM 16

/* edit the following lines or add and delete new lines */ 
static struct info_line dll_info[] = { 
   "This function calculates five lorentzians.",
   "It requires 16 parameters as follows:",
   "P1 is the constant background",
   "P2 is the center posiion of lorentzian 1 along the x axis.",
   "P3 is the amplitude of lorentzian 1.",
   "P4 is the width of lorentzian 1.", 
   "P5 is the center posiion of lorentzian 2 along the x axis.",
   "P6 is the amplitude of lorentzian 2.",
   "P7 is the width of lorentzian 2.", 
   "P8 is the center posiion of lorentzian 3 along the x axis.",
   "P9 is the amplitude of lorentzian 3.",
   "P10is the width of lorentzian 3.", 
   "P11  is the center posiion of lorentzian 4 along the x axis.",
   "P12 is the amplitude of lorentzian 4.",
   "P13 is the width of lorentzian 4.", 
   "P14 is the center posiion of lorentzian 5 along the x axis.",
   "P15 is the amplitude of lorentzian 5.",
   "P16 is the width of lorentzian 5.", 
   NULL};
   
VOID EXPENTRY fitfunc(void)
   {
   int i;
   double wid1, wid2;
   double wid3, wid4, wid5;
   double amp1, amp2;
   double amp3, amp4, amp5;
   double pos1, pos2;
   double pos3, pos4, pos5;
   double E1, E2;
   double E3, E4, E5;
   double lor1, lor2;
   double lor3, lor4, lor5;
   double offset;

   /* The next section is executed the first time through only */
   /* so put any initialization code in here                   */
   if (firsttime) 
      {
      firsttime = FALSE;
      /* Mask off all floating point exceptions */
      _control87(MCW_EM,MCW_EM);
      }

   /* give parameters more meaningful names  */
   offset=*Ptr[1];
   pos1 = *Ptr[2];
   amp1 = *Ptr[3];
   wid1 = *Ptr[4];
   pos2 = *Ptr[5];
   amp2 = *Ptr[6];
   wid2 = *Ptr[7];
   pos3 = *Ptr[8];
   amp3 = *Ptr[9];
   wid3 = *Ptr[10];
   pos4 = *Ptr[11];
   amp4 = *Ptr[12];
   wid4 = *Ptr[13];
   pos5 = *Ptr[14];
   amp5 = *Ptr[15];
   wid5 = *Ptr[16];

   /* This is the main loop that calculates the function datt[] */
   /* for all data points. Put your function in the loop body   */ 
   for (i = 0;i < Ndat;i++)
      {
      E1 = datx[i] - pos1;
      E2 = datx[i] - pos2;
      E3 = datx[i] - pos3;
      E4 = datx[i] - pos4;
      E5 = datx[i] - pos5;
      lor1 = amp1 / (1 + E1*E1/wid1/wid1) ;   /*Lorentzian */
      lor2 = amp2 / (1 + E2*E2/wid2/wid2) ;  
      lor3 = amp3 / (1 + E3*E3/wid3/wid3) ;  
      lor4 = amp4 / (1 + E4*E4/wid4/wid4) ;  
      lor5 = amp5 / (1 + E5*E5/wid5/wid5) ;  
      datt[i] = lor1 + lor2 + lor3 + lor4 + lor5 + offset;
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

