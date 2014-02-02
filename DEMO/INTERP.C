#define NPARAM 5
#define NADD 0

/* Example Interpreted function */
char firsttime = 1;
int Ndat;
/* consider the data arrays to be datx[Ndat] etc. */
double *datx,*daty,*datt;

/* The parameters  */
double *P1;     // x shift  
double *P2;     // y constant 
double *P3;     // y ^ 1 
double *P4;     // y ^ 2 
double *P5;     // y ^ 3 

main()
   {
   if (!verify(Ndat,NPARAM,NADD)) exit(1);
   /* do not put anything else in here */
   }

fitfunc()
   {
   int i;
   double x[Ndat];
   if (firsttime)
      {
      /* start up code goes in here */
      firsttime = 0;
      }
   for (i = 0;i < Ndat;i++)
      {
      /* edit the contents of this loop as necessary */
      /* to calculate datt as a function of the data */
      /* and parameters */
      x[i] = datx[i] - *P1;
      datt[i] = x[i] * (x[i] * (x[i] * *P5 + *P4) + *P3) + *P2;
      }
   }

 
