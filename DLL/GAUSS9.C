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
#define NPARAM  28

/* edit the following lines or add and delete new lines */ 
static struct info_line dll_info[] = { 
   "This function calculates nine gaussians.",
   "It requires 28 parameters as follows:",
   "P1 is the constant background.",
   "P2 is the center posiion of gaussian 1 along the x axis.",
   "P3 is the amplitude of gaussian 1.",
   "P4 is the width of gaussian 1.", 
   "P5 is the center posiion of gaussian 2 along the x axis.",
   "P6 is the amplitude of gaussian 2.",
   "P7 is the width of gaussian 2.", 
   "P8 is the center posiion of gaussian 3 along the x axis.",
   "P9 is the amplitude of gaussian 3.",
   "P10is the width of gaussian 3.", 
   "P11  is the center posiion of gaussian 4 along the x axis.",
   "P12 is the amplitude of gaussian 4.",
   "P13 is the width of gaussian 4.", 
   "P14 is the center posiion of gaussian 5 along the x axis.",
   "P15 is the amplitude of gaussian 5.",
   "P16 is the width of gaussian 5.", 
   "P17 is the center posiion of gaussian 6 along the x axis.",
   "P18 is the amplitude of gaussian 6.",
   "P19 is the width of gaussian 6.", 
   "P20 is the center posiion of gaussian 7 along the x axis.",
   "P21 is the amplitude of gaussian 7.",
   "P22 is the width of gaussian 7.", 
   "P23 is the center posiion of gaussian 8 along the x axis.",
   "P24 is the amplitude of gaussian 8.",
   "P25 is the width of gaussian 8.", 
   "P26 is the center posiion of gaussian 9 along the x axis.",
   "P27 is the amplitude of gaussian 9.",
   "P28 is the width of gaussian 9.", 
   NULL};
   
VOID EXPENTRY fitfunc(void)
   {
   int i;
   double wid1, wid2, wid6, wid7; 
   double wid3, wid4, wid5, wid8, wid9;
   double amp1, amp2 , amp6, amp7; 
   double amp3, amp4, amp5, amp8, amp9;
   double pos1, pos2, pos6, pos7;
   double pos3, pos4, pos5, pos8, pos9;
   double E1, E2, E6, E7;
   double E3, E4, E5, E8, E9;
   double gau1, gau2, gau6, gau7;
   double gau3, gau4, gau5, gau8, gau9;
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
   pos6 = *Ptr[17];
   amp6 = *Ptr[18];
   wid6 = *Ptr[19];
   pos7 = *Ptr[20];
   amp7 = *Ptr[21];
   wid7 = *Ptr[22];
   pos8 = *Ptr[23];
   amp8 = *Ptr[24];
   wid8 = *Ptr[25];
   pos9 = *Ptr[26];
   amp9 = *Ptr[27];
   wid9 = *Ptr[28];

   /* This is the main loop that calculates the function datt[] */
   /* for all data points. Put your function in the loop body   */ 
   for (i = 0;i < Ndat;i++)
      {
      E1 = datx[i] - pos1;
      E2 = datx[i] - pos2;
      E3 = datx[i] - pos3;
      E4 = datx[i] - pos4;
      E5 = datx[i] - pos5;
      E6 = datx[i] - pos6;
      E7 = datx[i] - pos7;
      E8 = datx[i] - pos8;
      E9 = datx[i] - pos9;
      gau1 = amp1 * exp(-E1*E1/wid1/wid1);        /*Gaussian */
      gau2 = amp2 * exp(-E2*E2/wid2/wid2); 
      gau3 = amp3 * exp(-E3*E3/wid3/wid3); 
      gau4 = amp4 * exp(-E4*E4/wid4/wid4); 
      gau5 = amp5 * exp(-E5*E5/wid5/wid5); 
      gau6 = amp6 * exp(-E6*E6/wid6/wid6); 
      gau7 = amp7 * exp(-E7*E7/wid7/wid7); 
      gau8 = amp8 * exp(-E8*E8/wid8/wid8); 
      gau9 = amp9 * exp(-E9*E9/wid9/wid9); 
      datt[i] = gau1 + gau2 + gau3 + gau4 + gau5 + gau6 + gau7 + gau8 + gau9 +offset;
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

