/*--------------------------------------------------------------------------
  Program name    : hetero.c

  Author          : Lars S Jermiin

  Institutions    : School of Biological Sciences
                    Heydon-Lawrence Building A08
                    University of Sydney
                    Sydney, NSW 2006, Australia
 
                  : CSIRO Entomology
                    GPO Box 1700
                    Canberra, ACT 2601, Augstralia

  Date begun      : 22 August, 1999

  Date modified   : 03 September, 2009

  Copyright       : 2003-2009. University of Sydney & CSIRO.
  
  Summary         : Hetero is designed to generate four nucleotide sequences
                    using a phylogeny and six nucleotide substitution models.
                    
                    The ancestral sequence is generated at random from a user-
                    specified frequency distribution.
                    
                    Nucleotide substitution model Ra applies to edge a.
                    Nucleotide substitution model Rb applies to edge b.
                    Nucleotide substitution model Rc applies to edge c.
                    Nucleotide substitution model Rd applies to edge d.
                    Nucleotide substitution model Re applies to edge e.
                    Nucleotide substitution model Rf applies to edge f.
                    
                    Edges e and f connects the root to the ancestors of the
                    lineages a, b, c and d.
                                        
                    The models of nucleotide substitution rates employ all
                    12 rate options, for each of the six edges.
                    
                    The equilibrium nucleotide contentis specified by the
                    user for each substitution model.
                    
                    The results are printed to two files; on contains details
                    about the simulation, and the other contains the aligned
                    nucleotide sequences in the sequential PHYLIP format.
                                        
                    The program is described in the following papers:
                    
                       Jermiin LS, Ho SYW, Ababneh F, Robinson J & Larkum AWD
                       (2003). Hetero: a program to simulate the evolution of
                       DNA on a four-taxon tree -- Applied Bioinformatics 2,
                       159-163.
                                           
                       Ababneh F, Jermiin LS, Robinson J (2006). Generation of
                       the exact distribution and simulation of matched
                       nucleotide sequences on a phylogenetic tree. Journal of
                       Mathematical Modelling and Algorithms 5, 291-308.
                       
Debugging          : Error in print out to file with information found and 
                     corrected [LSJ November 9 2004]
                   : Order of entry of conditional rates revised [LSJ August
                     3 2009]
----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define FALSE        0        /*Used for boolean decisions*/
#define TRUE         1        /*Used for boolean decisions*/
#define MAX_HITS    20        /*Limit of the number of multiple substitutions*/
#define RATE         0.40     /*Default conditional rate of change = 0.200*/
#define UNIFORM      0.25     /*Default nucleotide content = 0.250*/
#define LONG_TIME    0.95     /*Default length of edges a, b, c & d = 0.950*/
#define SHORT_TIME   0.05     /*Default length of edges e & f = 0.050*/
#define array_size  100000    /*Maximum length of sequences*/
#define IM1  2147483563       /*Constant from Numerical Recipes*/
#define IM2  2147483399       /*Constant from Numerical Recipes*/
#define AM   (1.0/IM1)        /*Constant from Numerical Recipes*/
#define IMM1 (IM1 - 1)        /*Constant from Numerical Recipes*/
#define IA1  40014            /*Constant from Numerical Recipes*/
#define IA2  40692            /*Constant from Numerical Recipes*/
#define IQ1  53668            /*Constant from Numerical Recipes*/
#define IQ2  52774            /*Constant from Numerical Recipes*/
#define IR1  12211            /*Constant from Numerical Recipes*/
#define IR2  3791             /*Constant from Numerical Recipes*/
#define NTAB 32               /*Constant from Numerical Recipes*/
#define NDIV (1 + IMM1/NTAB)  /*Constant from Numerical Recipes*/
#define EPS  1.2e-7           /*Constant from Numerical Recipes*/
#define RNMX (1.0 - EPS)      /*Constant from Numerical Recipes*/

/*
 *
 *---------------- Declaration of function Prototype -----------------------
 *
 */

void  Information(void);
void  Buffer_flush(void);
void  Screen_flush(void);
void  Tree_drawing(int feature);
float ran2(long *idum);

/*
 *
 *------------------- Declaration of External Variables --------------------
 *
 */

char  copyright[50] = "2003-2009, University of Sydney & CSIRO";
char  programName[7] = "Hetero";
char  version[10] = "1.1";

/*
 *
 *--------------------------- Main Program ---------------------------------
 *
 */

int main()
{
   char           c, d, choice = 'Y', use_time = 'Y';
                  /*Variables to pick up single letter input from the user*/

   char           outName1[75], outName2[75], line[120];
                  /*String variables to hold file names, etc.*/

   static char    seqA[array_size], seqB[array_size];
   static char    seqC[array_size], seqD[array_size];
                  /*Arrays to hold the four sequences during the evolution*/

   static int     hits[array_size];
                  /*Array to hold the number of times that a site has changed*/

   int            i, j, k, SUCCESS, length, cycles, pos, linelength;
                  /*Variables used predominantly in for-loops*/

   long           seed;
                  /*Variable needed by the random number generator*/

   unsigned long  sum_hits[MAX_HITS], tot_sum_hits[MAX_HITS];
                  /*Arrays with the number sites that have changed X times*/

   unsigned long  counter, counter_1, counter_2;
                  /*Variables used predominantly in for-loops*/
                  
   float          div_a, div_b, div_c, div_d, div_e, div_f, number, sum;
                  /*Variables holding mainly the edge lengths*/

   float          sub, splitA, splitB, splitC, splitD, splitE, splitF, splitG;
   float          constant_sites, hypervariable_sites;
                  /*Variables holding the number of different splits in the data*/

   float          Ra[4][4], Rb[4][4], Rc[4][4], Rd[4][4], Re[4][4], Rf[4][4];
                  /*Arrays holding the six rate matrices*/

   float          P0[4], P1[4], P2[4], P3[4], P4[4], P5[4], P6[4];
                  /*Arrays holding the nucleotide frequencies*/

   float          ac, ag, at, cg, ct, gt, ca, ga, gc, ta, tc, tg;
                  /*Variables holding the conditional rates of change*/

   float          pA, pC, pG, pT;
                  /*Variables holding the nucleotide frequencies*/

   float          average, percentage;
                  /*Your guess is as good as mine...*/

   float          GC_pair_AB, GC_pair_AC, GC_pair_AD;
   float          GC_pair_CD, GC_pair_BD, GC_pair_BC;
                  /*Variables holding the GC content of sequence pairs*/

   float          sumDif_GC_AB = 0.0, sumDif_GC_AC = 0.0, sumDif_GC_AD = 0.0;
   float          sumDif_GC_BC = 0.0, sumDif_GC_BD = 0.0, sumDif_GC_CD = 0.0;
                  /*Variables holding the sum of GD-content differences*/

   float          sum_splitA = 0.0, sum_splitB = 0.0, sum_splitC = 0.0;
   float          sum_splitD = 0.0, sum_splitE = 0.0, sum_splitF = 0.0;
   float          sum_splitG = 0.0, sum_constant_sites = 0.0;
   float          sum_hypervariable_sites = 0.0;
                  /*Variables holding the sum of different splits*/

   float          rate1, rate2, rate3, rate4, rate5, rate6;
                  /*Variables holding the average rate of change for the models*/

   float          T1[4][4], T2[4][4], T3[4][4], T4[4][4], T5[4][4], T6[4][4];
                  /*Arrays holding the threshold matrices*/

   float          obsMat[4][4];
                  /*Array holding the 4 x 4 divergence matrix*/

   time_t         systime;
                  /*Define systime as variable holding computer time*/

   FILE           *outFile1, *outFile2, *fopen();
                  /*Pointers to different files*/


/*
 *
 *Initilize R-matrices and PI-matrices
 *
 */

   for(i = 0; i < 4; i++) {
      for(j = 0; j < 4; j++) {
         if (i != j) {
            Ra[i][j] = Rb[i][j] = Rc[i][j] = RATE;
            Rd[i][j] = Re[i][j] = Rf[i][j] = RATE;
         }
      }
   }
   for(i = 0; i < 4; i++) {
      P0[i] = P1[i] = P2[i] = P3[i] = P4[i] = P5[i] = P6[i] = UNIFORM;
   }
   div_a = div_b = div_c = div_d = LONG_TIME;
   div_e = div_f = SHORT_TIME;

/*
 *
 *Welcome and information for users
 *
 */

   counter_1 = 0;
   SUCCESS = FALSE;
   while (SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nInformation [Y/N] ? ");
      c = getchar();
      d = getchar();
      if (d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         if (c != 'Y' && c != 'y' && c != 'N' && c != 'n')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if (c == 'Y' || c == 'y') {
               SUCCESS = TRUE;
               Information();
            }
            else
               SUCCESS = TRUE;
         }
      }
   }
   if (SUCCESS == FALSE) {
      fprintf(stderr,"\n\nSorry - Program aborted\n\n");
      exit(1);
   }


/*
 *
 *Data entry
 *
 */

   printf("\n_____________________ Data Entry _____________________\n");
   printf("\nThis program considers differences between an ancestor");
   printf("\nand its two decendants in terms of time or the average");
   printf("\nnumber of substitutions per site.\n");
   counter_1 = 0;
   SUCCESS = FALSE;
   while (SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use time [y/n] ? ");
      c = getchar();
      d = getchar();
      if (d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         if (c != 'Y' && c != 'y' && c != 'N' && c != 'n')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            SUCCESS = TRUE;
            use_time = toupper(c);
         }
      }
   }
   if (SUCCESS == FALSE) {
      fprintf(stderr,"\n\nSorry - Program aborted\n\n");
      exit(1);
   }
   Screen_flush();
   printf("\nPlease consider the following labelled phylogeny.  Six");
   printf("\nnucleotide substitution models (Ra, Rb, ..., Rf) apply");
   printf("\nto the edges (a, b, ..., f) of this rooted, four-taxon");
   printf("\ntree:\n");
   Tree_drawing(-1);
   printf("\n\nTo allow the Monte Carlo simulation to occur, you must");
   printf("\nenter the parameters that control the evolution of the");
   printf("\nsequence from the origin towards the four descendants.");

/*
 *
 *Entry of the tree's edge lengths
 *
 */

   printf("\n\nDefault edge lengths of a, b, c, d, e and f\n");
   printf("\n%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",div_a, div_b, div_c, div_d, div_e, div_f);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these values ........ [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(0);
               printf("\n\nEnter the length of the six edges a, b, c, d, e and f,");
               printf("\nwith blanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (6 == sscanf(line,"%f %f %f %f %f %f",&div_a, &div_b, &div_c, &div_d, &div_e, &div_f))
                        if(div_a <= 10 && div_b <= 10 && div_c <= 10 && div_d <= 10 && div_e <= 10 && div_f <= 10)
                           SUCCESS = TRUE;
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     else 
                        printf("\n\nToo few values -- please re-enter...\n\n");
                   }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of the nucleotide content in the anscestral sequence
 *
 */

   Screen_flush();
   Tree_drawing(7);
   printf("\n\nThe nucleotide content (A, C, G and T) of the last");
   printf("\ncommon ancestral sequence (the origin)\n");
   printf("\n%10f %10f %10f %10f\n",P0[0], P0[1], P0[2], P0[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(7);
               printf("\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P0[0] = pA;
                           P0[1] = pC;
                           P0[2] = pG;
                           P0[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of nucleotide frequencies for model Re
 *
 */

   Screen_flush();
   Tree_drawing(5);
   printf("\n\n\nThe frequencies of A, C, G and T in model Re\n");
   printf("\n%10f %10f %10f %10f\n",P5[0], P5[1], P5[2], P5[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(5);
               printf("\n\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P5[0] = pA;
                           P5[1] = pC;
                           P5[2] = pG;
                           P5[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of conditional rates of change for model Re
 *
 */

   Screen_flush();
   Tree_drawing(5);
   printf("\n\nConditional rates of change under model Re (listed in order\n");
   printf("A->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C T->G)\n");
   printf("\n%4.2f %4.2f %4.2f %4.2f",Re[0][1], Re[0][2], Re[0][3], Re[1][2]);
   printf(" %4.2f %4.2f %4.2f %4.2f",Re[1][3], Re[2][3], Re[1][0], Re[2][0]);
   printf(" %4.2f %4.2f %4.2f %4.2f\n",Re[3][0] ,Re[2][1], Re[3][1], Re[3][2]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these conditional rates ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(5);
               printf("\n\nEnter conditional rates of change for model Re (list in order");
               printf("\nA->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C & T->G");
               printf("\nwith blanks separating each value) .............. [0.0 - 1.0]\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (12 == sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&ac,&ag,&at,&cg,&ct,&gt,&ca,&ga,&ta,&gc,&tc,&tg)) {
                        if(ac<=1 && ag<=1 && at<=1 && cg<=1 && ct<=1 && gt<=1 && ca<=1 && ga<=1 && ta<=1 && gc<=1 && tc<=1 && tg<=1) {
                           Re[0][1] = ac;
                           Re[0][2] = ag;
                           Re[0][3] = at;
                           Re[1][2] = cg;
                           Re[1][3] = ct;
                           Re[2][3] = gt;
                           Re[1][0] = ca;
                           Re[2][0] = ga;
                           Re[3][0] = ta;
                           Re[2][1] = gc;
                           Re[3][1] = tc;
                           Re[3][2] = tg;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Generating and testing the full rate matrix for model Re
 *
 */

   Re[0][1] = Re[0][1] * P5[1];
   Re[0][2] = Re[0][2] * P5[2];
   Re[0][3] = Re[0][3] * P5[3];
   Re[1][0] = Re[1][0] * P5[0];
   Re[1][2] = Re[1][2] * P5[2];
   Re[1][3] = Re[1][3] * P5[3];
   Re[2][0] = Re[2][0] * P5[0];
   Re[2][1] = Re[2][1] * P5[1];
   Re[2][3] = Re[2][3] * P5[3];
   Re[3][0] = Re[3][0] * P5[0];
   Re[3][1] = Re[3][1] * P5[1];
   Re[3][2] = Re[3][2] * P5[2];
   Re[0][0] = 0.0 - Re[0][1] - Re[0][2] - Re[0][3];
   Re[1][1] = 0.0 - Re[1][0] - Re[1][2] - Re[1][3];
   Re[2][2] = 0.0 - Re[2][0] - Re[2][1] - Re[2][3];
   Re[3][3] = 0.0 - Re[3][0] - Re[3][1] - Re[3][2];
   SUCCESS = TRUE;
   for(i = 0; i < 4; i++)
      if(Re[i][i] < -1.0)
         SUCCESS = FALSE;
   if(SUCCESS == FALSE) {
     printf("\n\nProblem -- diagonal elements in model Re");
     printf("\nmay not be less than -1.0.\n\n");
     for(i = 0; i < 4; i++) {
       for(j = 0; j < 4; j++)
          printf("%10.6f",Re[i][j]);
       printf("\n");
     }
     printf("\nSorry - Program aborted\n\n");    
     exit(1);
   }

/*
 *
 *Estimating the average rate of change for model Re
 *
 */

   rate5 = -(Re[0][0]*P5[0] + Re[1][1]*P5[1] + Re[2][2]*P5[2] + Re[3][3]*P5[3]);

/*
 *
 *Entry of nucleotide frequencies for model Rf
 *
 */

   Screen_flush();
   Tree_drawing(6);
   printf("\n\n\nThe frequencies of A, C, G and T in model Rf\n");
   printf("\n%10f %10f %10f %10f\n",P6[0], P6[1], P6[2], P6[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(6);
               printf("\n\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P6[0] = pA;
                           P6[1] = pC;
                           P6[2] = pG;
                           P6[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of conditional rates of change for model Rf
 *
 */

   Screen_flush();
   Tree_drawing(6);
   printf("\n\nConditional rates of change under model Rf (listed in order\n");
   printf("A->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C T->G)\n");
   printf("\n%4.2f %4.2f %4.2f %4.2f",Rf[0][1], Rf[0][2], Rf[0][3], Rf[1][2]);
   printf(" %4.2f %4.2f %4.2f %4.2f",Rf[1][3], Rf[2][3], Rf[1][0], Rf[2][0]);
   printf(" %4.2f %4.2f %4.2f %4.2f\n", Rf[3][0],Rf[2][1], Rf[3][1], Rf[3][2]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these conditional rates ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(6);
               printf("\n\n\nEnter conditional rates of change for model Rf (list in order");
               printf("\nA->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C & T->G");
               printf("\nwith blanks separating each value) .............. [0.0 - 1.0]\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (12 == sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&ac,&ag,&at,&cg,&ct,&gt,&ca,&ga,&ta,&gc,&tc,&tg)) {
                        if(ac<=1 && ag<=1 && at<=1 && cg<=1 && ct<=1 && gt<=1 && ca<=1 && ga<=1 && ta<=1 && gc<=1 && tc<=1 && tg<=1) {
                           Rf[0][1] = ac;
                           Rf[0][2] = ag;
                           Rf[0][3] = at;
                           Rf[1][2] = cg;
                           Rf[1][3] = ct;
                           Rf[2][3] = gt;
                           Rf[1][0] = ca;
                           Rf[2][0] = ga;
                           Rf[3][0] = ta;
                           Rf[2][1] = gc;
                           Rf[3][1] = tc;
                           Rf[3][2] = tg;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Generating and testing the full rate matrix for model Rf
 *
 */

   Rf[0][1] = Rf[0][1] * P6[1];
   Rf[0][2] = Rf[0][2] * P6[2];
   Rf[0][3] = Rf[0][3] * P6[3];
   Rf[1][0] = Rf[1][0] * P6[0];
   Rf[1][2] = Rf[1][2] * P6[2];
   Rf[1][3] = Rf[1][3] * P6[3];
   Rf[2][0] = Rf[2][0] * P6[0];
   Rf[2][1] = Rf[2][1] * P6[1];
   Rf[2][3] = Rf[2][3] * P6[3];
   Rf[3][0] = Rf[3][0] * P6[0];
   Rf[3][1] = Rf[3][1] * P6[1];
   Rf[3][2] = Rf[3][2] * P6[2];
   Rf[0][0] = 0.0 - Rf[0][1] - Rf[0][2] - Rf[0][3];
   Rf[1][1] = 0.0 - Rf[1][0] - Rf[1][2] - Rf[1][3];
   Rf[2][2] = 0.0 - Rf[2][0] - Rf[2][1] - Rf[2][3];
   Rf[3][3] = 0.0 - Rf[3][0] - Rf[3][1] - Rf[3][2];
   SUCCESS = TRUE;
   for(i = 0; i < 4; i++)
      if(Re[i][i] < -1.0)
         SUCCESS = FALSE;
   if(SUCCESS == FALSE) {
     printf("\n\nProblem -- diagonal elements in model Rf");
     printf("\nmay not be less than -1.0.\n\n");
     for(i = 0; i < 4; i++) {
       for(j = 0; j < 4; j++)
          printf("%10.6f",Rf[i][j]);
       printf("\n");
     }
     printf("\nSorry - Program aborted\n\n");    
     exit(1);
   }

/*
 *
 *Estimating the average rate of change for model Rf
 *
 */

   rate6 = -(Rf[0][0]*P6[0] + Rf[1][1]*P6[1] + Rf[2][2]*P6[2] + Rf[3][3]*P6[3]);

/*
 *
 *Entry of nucleotide frequencies for model Ra
 *
 */

   Screen_flush();
   Tree_drawing(1);
   printf("\n\n\nThe frequencies of A, C, G and T in model Ra\n");
   printf("\n%10f %10f %10f %10f\n",P1[0], P1[1], P1[2], P1[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(1);
               printf("\n\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P1[0] = pA;
                           P1[1] = pC;
                           P1[2] = pG;
                           P1[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of conditional rates of change for model Ra
 *
 */

   Screen_flush();
   Tree_drawing(1);
   printf("\n\nConditional rates of change under model Ra (listed in order\n");
   printf("A->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C T->G)\n");
   printf("\n%4.2f %4.2f %4.2f %4.2f",Ra[0][1], Ra[0][2], Ra[0][3], Ra[1][2]);
   printf(" %4.2f %4.2f %4.2f %4.2f",Ra[1][3], Ra[2][3], Ra[1][0], Ra[2][0]);
   printf(" %4.2f %4.2f %4.2f %4.2f\n",Ra[3][0], Ra[2][1], Ra[3][1], Ra[3][2]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these conditional rates ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(1);
               printf("\n\nEnter conditional rates of change for model Ra (list in order");
               printf("\nA->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C & T->G");
               printf("\nwith blanks separating each value) .............. [0.0 - 1.0]\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (12 == sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&ac,&ag,&at,&cg,&ct,&gt,&ca,&ga,&ta,&gc,&tc,&tg)) {
                        if(ac<=1 && ag<=1 && at<=1 && cg<=1 && ct<=1 && gt<=1 && ca<=1 && ga<=1 && ta<=1 && gc<=1 && tc<=1 && tg<=1) {
                           Ra[0][1] = ac;
                           Ra[0][2] = ag;
                           Ra[0][3] = at;
                           Ra[1][2] = cg;
                           Ra[1][3] = ct;
                           Ra[2][3] = gt;
                           Ra[1][0] = ca;
                           Ra[2][0] = ga;
                           Ra[3][0] = ta;
                           Ra[2][1] = gc;
                           Ra[3][1] = tc;
                           Ra[3][2] = tg;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Generating and testing the full rate matrix for model Ra
 *
 */

   Ra[0][1] = Ra[0][1] * P1[1];
   Ra[0][2] = Ra[0][2] * P1[2];
   Ra[0][3] = Ra[0][3] * P1[3];
   Ra[1][0] = Ra[1][0] * P1[0];
   Ra[1][2] = Ra[1][2] * P1[2];
   Ra[1][3] = Ra[1][3] * P1[3];
   Ra[2][0] = Ra[2][0] * P1[0];
   Ra[2][1] = Ra[2][1] * P1[1];
   Ra[2][3] = Ra[2][3] * P1[3];
   Ra[3][0] = Ra[3][0] * P1[0];
   Ra[3][1] = Ra[3][1] * P1[1];
   Ra[3][2] = Ra[3][2] * P1[2];
   Ra[0][0] = 0.0 - Ra[0][1] - Ra[0][2] - Ra[0][3];
   Ra[1][1] = 0.0 - Ra[1][0] - Ra[1][2] - Ra[1][3];
   Ra[2][2] = 0.0 - Ra[2][0] - Ra[2][1] - Ra[2][3];
   Ra[3][3] = 0.0 - Ra[3][0] - Ra[3][1] - Ra[3][2];
   SUCCESS = TRUE;
   for(i = 0; i < 4; i++)
      if(Re[i][i] < -1.0)
         SUCCESS = FALSE;
   if(SUCCESS == FALSE) {
     printf("\n\nProblem -- diagonal elements in model Ra");
     printf("\nmay not be less than -1.0.\n\n");
     for(i = 0; i < 4; i++) {
       for(j = 0; j < 4; j++)
          printf("%10.6f",Ra[i][j]);
       printf("\n");
     }
     printf("\nSorry - Program aborted\n\n");    
     exit(1);
   }

/*
 *
 *Estimating the average rate of change for model Ra
 *
 */

   rate1 = -(Ra[0][0]*P1[0] + Ra[1][1]*P1[1] + Ra[2][2]*P1[2] + Ra[3][3]*P1[3]);

/*
 *
 *Entry of nucleotide frequencies for model Rb
 *
 */

   Screen_flush();
   Tree_drawing(2);
   printf("\n\n\nThe frequencies of A, C, G and T in model Rb\n");
   printf("\n%10f %10f %10f %10f\n",P2[0], P2[1], P2[2], P2[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(2);
               printf("\n\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P2[0] = pA;
                           P2[1] = pC;
                           P2[2] = pG;
                           P2[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of conditional rates of change for model Rb
 *
 */

   Screen_flush();
   Tree_drawing(2);
   printf("\n\nConditional rates of change under model Rb (listed in order\n");
   printf("A->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C T->G)\n");
   printf("\n%4.2f %4.2f %4.2f %4.2f",Rb[0][1], Rb[0][2], Rb[0][3], Rb[1][2]);
   printf(" %4.2f %4.2f %4.2f %4.2f",Rb[1][3], Rb[2][3], Rb[1][0], Rb[2][0]);
   printf(" %4.2f %4.2f %4.2f %4.2f\n",Rb[3][0], Rb[2][1], Rb[3][1], Rb[3][2]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these conditional rates ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(2);
               printf("\n\nEnter conditional rates of change for model Rb (list in order");
               printf("\nA->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C & T->G");
               printf("\nwith blanks separating each value) .............. [0.0 - 1.0]\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (12 == sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&ac,&ag,&at,&cg,&ct,&gt,&ca,&ga,&ta,&gc,&tc,&tg)) {
                        if(ac<=1 && ag<=1 && at<=1 && cg<=1 && ct<=1 && gt<=1 && ca<=1 && ga<=1 && ta<=1 && gc<=1 && tc<=1 && tg<=1) {
                           Rb[0][1] = ac;
                           Rb[0][2] = ag;
                           Rb[0][3] = at;
                           Rb[1][2] = cg;
                           Rb[1][3] = ct;
                           Rb[2][3] = gt;
                           Rb[1][0] = ca;
                           Rb[2][0] = ga;
                           Rb[3][0] = ta;
                           Rb[2][1] = gc;
                           Rb[3][1] = tc;
                           Rb[3][2] = tg;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Generating and testing the full rate matrix for model Rb
 *
 */

   Rb[0][1] = Rb[0][1] * P2[1];
   Rb[0][2] = Rb[0][2] * P2[2];
   Rb[0][3] = Rb[0][3] * P2[3];
   Rb[1][0] = Rb[1][0] * P2[0];
   Rb[1][2] = Rb[1][2] * P2[2];
   Rb[1][3] = Rb[1][3] * P2[3];
   Rb[2][0] = Rb[2][0] * P2[0];
   Rb[2][1] = Rb[2][1] * P2[1];
   Rb[2][3] = Rb[2][3] * P2[3];
   Rb[3][0] = Rb[3][0] * P2[0];
   Rb[3][1] = Rb[3][1] * P2[1];
   Rb[3][2] = Rb[3][2] * P2[2];
   Rb[0][0] = 0.0 - Rb[0][1] - Rb[0][2] - Rb[0][3];
   Rb[1][1] = 0.0 - Rb[1][0] - Rb[1][2] - Rb[1][3];
   Rb[2][2] = 0.0 - Rb[2][0] - Rb[2][1] - Rb[2][3];
   Rb[3][3] = 0.0 - Rb[3][0] - Rb[3][1] - Rb[3][2];
   SUCCESS = TRUE;
   for(i = 0; i < 4; i++)
      if(Re[i][i] < -1.0)
         SUCCESS = FALSE;
   if(SUCCESS == FALSE) {
     printf("\n\nProblem -- diagonal elements in model Rb");
     printf("\nmay not be less than -1.0.\n\n");
     for(i = 0; i < 4; i++) {
       for(j = 0; j < 4; j++)
          printf("%10.6f",Rb[i][j]);
       printf("\n");
     }
     printf("\nSorry - Program aborted\n\n");    
     exit(1);
   }

/*
 *
 *Estimating the average rate of change for model Rb
 *
 */

   rate2 = -(Rb[0][0]*P2[0] + Rb[1][1]*P2[1] + Rb[2][2]*P2[2] + Rb[3][3]*P2[3]);

/*
 *
 *Entry of nucleotide frequencies for model Rc
 *
 */

   Screen_flush();
   Tree_drawing(3);
   printf("\n\n\nThe frequencies of A, C, G and T in model Rc\n");
   printf("\n%10f %10f %10f %10f\n",P3[0], P3[1], P3[2], P3[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(3);
               printf("\n\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P3[0] = pA;
                           P3[1] = pC;
                           P3[2] = pG;
                           P3[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of conditional rates of change for model Rc
 *
 */

   Screen_flush();
   Tree_drawing(3);
   printf("\n\nConditional rates of change under model Rc (listed in order\n");
   printf("A->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C T->G)\n");
   printf("\n%4.2f %4.2f %4.2f %4.2f",Rc[0][1], Rc[0][2], Rc[0][3], Rc[1][2]);
   printf(" %4.2f %4.2f %4.2f %4.2f",Rc[1][3], Rc[2][3], Rc[1][0], Rc[2][0]);
   printf(" %4.2f %4.2f %4.2f %4.2f\n",Rc[3][0], Rc[2][1], Rc[3][1], Rc[3][2]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these conditional rates ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(3);
               printf("\n\nEnter conditional rates of change for model Rc (list in order");
               printf("\nA->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C & T->G");
               printf("\nwith blanks separating each value) .............. [0.0 - 1.0]\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (12 == sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&ac,&ag,&at,&cg,&ct,&gt,&ca,&ga,&ta,&gc,&tc,&tg)) {
                        if(ac<=1 && ag<=1 && at<=1 && cg<=1 && ct<=1 && gt<=1 && ca<=1 && ga<=1 && ta<=1 && gc<=1 && tc<=1 && tg<=1) {
                           Rc[0][1] = ac;
                           Rc[0][2] = ag;
                           Rc[0][3] = at;
                           Rc[1][2] = cg;
                           Rc[1][3] = ct;
                           Rc[2][3] = gt;
                           Rc[1][0] = ca;
                           Rc[2][0] = ga;
                           Rc[3][0] = ta;
                           Rc[2][1] = gc;
                           Rc[3][1] = tc;
                           Rc[3][2] = tg;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Generating and testing the full rate matrix for model Rc
 *
 */

   Rc[0][1] = Rc[0][1] * P3[1];
   Rc[0][2] = Rc[0][2] * P3[2];
   Rc[0][3] = Rc[0][3] * P3[3];
   Rc[1][0] = Rc[1][0] * P3[0];
   Rc[1][2] = Rc[1][2] * P3[2];
   Rc[1][3] = Rc[1][3] * P3[3];
   Rc[2][0] = Rc[2][0] * P3[0];
   Rc[2][1] = Rc[2][1] * P3[1];
   Rc[2][3] = Rc[2][3] * P3[3];
   Rc[3][0] = Rc[3][0] * P3[0];
   Rc[3][1] = Rc[3][1] * P3[1];
   Rc[3][2] = Rc[3][2] * P3[2];
   Rc[0][0] = 0.0 - Rc[0][1] - Rc[0][2] - Rc[0][3];
   Rc[1][1] = 0.0 - Rc[1][0] - Rc[1][2] - Rc[1][3];
   Rc[2][2] = 0.0 - Rc[2][0] - Rc[2][1] - Rc[2][3];
   Rc[3][3] = 0.0 - Rc[3][0] - Rc[3][1] - Rc[3][2];
   SUCCESS = TRUE;
   for(i = 0; i < 4; i++)
      if(Re[i][i] < -1.0)
         SUCCESS = FALSE;
   if(SUCCESS == FALSE) {
     printf("\n\nProblem -- diagonal elements in model Rc");
     printf("\nmay not be less than -1.0.\n\n");
     for(i = 0; i < 4; i++) {
       for(j = 0; j < 4; j++)
          printf("%10.6f",Rc[i][j]);
       printf("\n");
     }
     printf("\nSorry - Program aborted\n\n");    
     exit(1);
   }

/*
 *
 *Estimating the average rate of change for model Rc
 *
 */

   rate3 = -(Rc[0][0]*P3[0] + Rc[1][1]*P3[1] + Rc[2][2]*P3[2] + Rc[3][3]*P3[3]);

/*
 *
 *Entry of nucleotide frequencies for model Rd
 *
 */

   Screen_flush();
   Tree_drawing(4);
   printf("\n\n\nThe frequencies of A, C, G and T in model Rd\n");
   printf("\n%10f %10f %10f %10f\n",P4[0], P4[1], P4[2], P4[3]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these frequencies ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(4);
               printf("\n\n\nEnter nucleotide frequency (A, C, G and T), with");
               printf("\nblanks separating each value, then press RETURN.\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (4 == sscanf(line,"%f %f %f %f", &pA, &pC, &pG, &pT)) {
                        sum = pA + pC + pG + pT;
                        if(0.9999 < sum && sum < 1.0001) {
                           P4[0] = pA;
                           P4[1] = pC;
                           P4[2] = pG;
                           P4[3] = pT;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues do not add up to 1.0 -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Entry of conditional rates of change for model Rd
 *
 */

   Screen_flush();
   Tree_drawing(4);
   printf("\n\nConditional rates of change under model Rd (listed in order\n");
   printf("A->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C T->G)\n");
   printf("\n%4.2f %4.2f %4.2f %4.2f",Rd[0][1], Rd[0][2], Rd[0][3], Rd[1][2]);
   printf(" %4.2f %4.2f %4.2f %4.2f",Rd[1][3], Rd[2][3], Rd[1][0], Rd[2][0]);
   printf(" %4.2f %4.2f %4.2f %4.2f\n",Rd[3][0], Rd[2][1], Rd[3][1], Rd[3][2]);
   SUCCESS = FALSE;
   counter_1 = 0;
   while(SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      printf("\nDo you want to use these conditional rates ... [y/n] : ");
      c = getchar();
      d = getchar();
      if(d != '\n') {
         printf("\n\nERROR - This character is not an option!\n");
         Buffer_flush();
      }
      else {
         c = toupper(c);
         if(c != 'Y' && c != 'N')
            printf("\n\nERROR - This character is not an option!\n");
         else {
            if(c == 'Y') {
               SUCCESS = TRUE;
            }
            else {
               Screen_flush();
               Tree_drawing(4);
               printf("\n\nEnter conditional rates of change for model Rd (list in order");
               printf("\nA->C A->G A->T C->G C->T G->T C->A G->A T->A G->C T->C & T->G");
               printf("\nwith blanks separating each value) .............. [0.0 - 1.0]\n\n");
               counter_2 = 0;
               while (SUCCESS == FALSE && counter_2 < 3) {
                  ++counter_2;
                  i = 0;
                  while((c = getchar()) != '\n')
                     line[i++] = c;
                  line[i] = '\0';
                  linelength = i;
                  for(i = 0; i < linelength; i++) {
                     c = line[i];
                     if(c == '0' || c == '1' || c == '2' || c == '3' ||
                        c == '4' || c == '5' || c == '6' || c == '7' ||
                        c == '8' || c == '9' || c == '.' || c == ' ')
                        ;
                     else
                        SUCCESS = TRUE;
                  }
                  if (SUCCESS == TRUE) {
                     SUCCESS = FALSE;
                     printf("\n\nInvalid characters used -- please re-enter...\n\n");
                  }
                  else {
                     if (12 == sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f",&ac,&ag,&at,&cg,&ct,&gt,&ca,&ga,&ta,&gc,&tc,&tg)) {
                        if(ac<=1 && ag<=1 && at<=1 && cg<=1 && ct<=1 && gt<=1 && ca<=1 && ga<=1 && ta<=1 && gc<=1 && tc<=1 && tg<=1) {
                           Rd[0][1] = ac;
                           Rd[0][2] = ag;
                           Rd[0][3] = at;
                           Rd[1][2] = cg;
                           Rd[1][3] = ct;
                           Rd[2][3] = gt;
                           Rd[1][0] = ca;
                           Rd[2][0] = ga;
                           Rd[3][0] = ta;
                           Rd[2][1] = gc;
                           Rd[3][1] = tc;
                           Rd[3][2] = tg;
                           SUCCESS = TRUE;
                        }
                        else
                           printf("\n\nValues outside the range -- please re-enter...\n\n");
                     }
                     else
                        printf("\n\nToo few values -- please re-enter...\n\n");
                  }
               }
               if(SUCCESS == FALSE) {
                  printf("\n\nSorry - Program aborted\n\n");
                  exit(1);
               }
            }
         }
      }
   }
   if(SUCCESS == FALSE) {
     printf("\n\nSorry - Program aborted\n\n");
     exit(1);
   }

/*
 *
 *Generating and testing the full rate matrix for model Rd
 *
 */

   Rd[0][1] = Rd[0][1] * P4[1];
   Rd[0][2] = Rd[0][2] * P4[2];
   Rd[0][3] = Rd[0][3] * P4[3];
   Rd[1][0] = Rd[1][0] * P4[0];
   Rd[1][2] = Rd[1][2] * P4[2];
   Rd[1][3] = Rd[1][3] * P4[3];
   Rd[2][0] = Rd[2][0] * P4[0];
   Rd[2][1] = Rd[2][1] * P4[1];
   Rd[2][3] = Rd[2][3] * P4[3];
   Rd[3][0] = Rd[3][0] * P4[0];
   Rd[3][1] = Rd[3][1] * P4[1];
   Rd[3][2] = Rd[3][2] * P4[2];
   Rd[0][0] = 0.0 - Rd[0][1] - Rd[0][2] - Rd[0][3];
   Rd[1][1] = 0.0 - Rd[1][0] - Rd[1][2] - Rd[1][3];
   Rd[2][2] = 0.0 - Rd[2][0] - Rd[2][1] - Rd[2][3];
   Rd[3][3] = 0.0 - Rd[3][0] - Rd[3][1] - Rd[3][2];
   SUCCESS = TRUE;
   for(i = 0; i < 4; i++)
      if(Re[i][i] < -1.0)
         SUCCESS = FALSE;
   if(SUCCESS == FALSE) {
     printf("\n\nProblem -- diagonal elements in model Rd");
     printf("\nmay not be less than -1.0.\n\n");
     for(i = 0; i < 4; i++) {
       for(j = 0; j < 4; j++)
          printf("%10.6f",Rd[i][j]);
       printf("\n");
     }
     printf("\nSorry - Program aborted\n\n");    
     exit(1);
   }

/*
 *
 *Estimating the average rate of change for model Rd
 *
 */

   rate4 = -(Rd[0][0]*P4[0] + Rd[1][1]*P4[1] + Rd[2][2]*P4[2] + Rd[3][3]*P4[3]);

/*
 *
 *Entry of number of simulations
 *
 */

   Screen_flush();
   printf("\n\nEnter other details for simulation:\n");
   printf("\nSequence length .................. [10 - 100000] : ");
   counter_1 = 0;
   SUCCESS = FALSE;
   while (SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      i = 0;
      while((c = getchar()) != '\n')
         line[i++] = c;
      line[i] = '\0';
      linelength = i;
      for(i = 0; i < linelength; i++) {
         c = line[i];
         if(c == '0' || c == '1' || c == '2' || c == '3' ||
            c == '4' || c == '5' || c == '6' || c == '7' ||
            c == '8' || c == '9')
            ;
         else
            SUCCESS = TRUE;
      }
      if (SUCCESS == TRUE) {
         SUCCESS = FALSE;
         printf("\n\nInvalid characters used -- please re-enter...\n\n");
      }
      else {
         if (1 == sscanf(line,"%d",&length))
            if(10 <= length && length <= 100000)
               SUCCESS = TRUE;
            else
               printf("\n\nValues outside the range -- please re-enter...\n\n");
         else 
            printf("\n\nToo few values -- please re-enter...\n\n");
      }
   }
   if(SUCCESS == FALSE) {
      printf("\n\nSorry - Program aborted\n\n");
      exit(1);
   }

/*
 *
 *Entry of number of simulations
 *
 */

   printf("\nNumber of cycles in simulation ..... [1 - 10000] : ");
   counter_1 = 0;
   SUCCESS = FALSE;
   while (SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      i = 0;
      while((c = getchar()) != '\n')
         line[i++] = c;
      line[i] = '\0';
      linelength = i;
      for(i = 0; i < linelength; i++) {
         c = line[i];
         if(c == '0' || c == '1' || c == '2' || c == '3' ||
            c == '4' || c == '5' || c == '6' || c == '7' ||
            c == '8' || c == '9')
            ;
         else
            SUCCESS = TRUE;
      }
      if (SUCCESS == TRUE) {
         SUCCESS = FALSE;
         printf("\n\nInvalid characters used -- please re-enter...\n\n");
      }
      else {
         if (1 == sscanf(line,"%d",&cycles))
            if(cycles <= 10000)
               SUCCESS = TRUE;
            else
               printf("\n\nValues outside the range -- please re-enter...\n\n");
         else 
            printf("\n\nToo few values -- please re-enter...\n\n");
      }
   }
   if(SUCCESS == FALSE) {
      printf("\n\nSorry - Program aborted\n\n");
      exit(1);
   }

/*
 *
 *Entry of seed for simulations
 *
 */

   printf("\nSeed (needed for simulation) ....... [1 - 32767] : ");
   counter_1 = 0;
   SUCCESS = FALSE;
   while (SUCCESS == FALSE && counter_1 < 3) {
      ++counter_1;
      i = 0;
      while((c = getchar()) != '\n')
         line[i++] = c;
      line[i] = '\0';
      linelength = i;
      for(i = 0; i < linelength; i++) {
         c = line[i];
         if(c == '0' || c == '1' || c == '2' || c == '3' ||
            c == '4' || c == '5' || c == '6' || c == '7' ||
            c == '8' || c == '9')
            ;
         else
            SUCCESS = TRUE;
      }
      if (SUCCESS == TRUE) {
         SUCCESS = FALSE;
         printf("\n\nInvalid characters used -- please re-enter...\n\n");
      }
      else {
         if (1 == sscanf(line,"%ld",&seed))
            if(seed <= 32767)
               SUCCESS = TRUE;
            else
               printf("\n\nValues outside the range -- please re-enter...\n\n");
         else 
            printf("\n\nToo few values -- please re-enter...\n\n");
      }
   }
   if(SUCCESS == FALSE) {
      printf("\n\nSorry - Program aborted\n\n");
      exit(1);
   }
   if (seed >= 0)
     seed = -seed;

/*
 *
 *Determining whether multiple hits should be allowed
 *
 */

   if(use_time == 'Y') {
      sum = div_a * rate1 + div_b * rate2 + div_e * rate5 +
            div_c * rate3 + div_d * rate4 + div_f * rate6;
   }
   else {
      sum = div_a + div_b + div_c + div_d + div_e + div_f;
   }
   if(sum < 1.0) {
      counter_1 = 0;
      SUCCESS = FALSE;
      while (SUCCESS == FALSE && counter_1 < 3) {
         ++counter_1;
         printf("\nAllow multiple hits ...................... [y/n] : ");
         c = getchar();
         d = getchar();
         if (d != '\n') {
            printf("\n\nERROR - This character is not an option!\n");
            Buffer_flush();
         }
         else {
            if (c != 'Y' && c != 'y' && c != 'N' && c != 'n')
               printf("\n\nERROR - This character is not an option!\n");
            else {
               SUCCESS = TRUE;
               choice = toupper(c);
            }
         }
      }
      if (SUCCESS == FALSE) {
         fprintf(stderr,"\n\nSorry - Program aborted\n\n");
         exit(1);
      }
   }

   printf("\n\nEnter name of output file:\n\n");

/*
 *
 *Entering name of first output file
 *
 */

   printf("Store information on simulation in file ........ : ");
   i = 0;
   while((c = getchar()) != '\n')
      outName1[i++] = c;
   outName1[i] = '\0';
   linelength = i;
   SUCCESS = FALSE;
   counter_1 = 0;
   while (counter_1 < 3 && SUCCESS == FALSE) {
      ++counter_1;
      if (linelength == 0) {
         printf("\nNo name entered - please re-enter file name .... : ");
         i = 0;
         while((c = getchar()) != '\n')
            outName1[i++] = c;
         outName1[i] = '\0';
         linelength = i;
      }
      else {
         for(i = 0; i < linelength; i++) {
            c = line[i];
            if(c != ' ')
               SUCCESS = TRUE;
         }
      }
      if (SUCCESS == TRUE) {
         if ((outFile1 = fopen(outName1, "w")) == NULL)
            SUCCESS = FALSE;
      }
   }
   if (SUCCESS == FALSE) {
      fprintf(stderr,"Program aborted - could not open file named");
      fprintf(stderr," %s.\n",outName1);
      exit(1);
   }

/*
 *
 *Entering name of second output file
 *
 */

   printf("\nStore sequences in file ........................ : ");
   j = 0;
   while((c = getchar()) != '\n')
      outName2[j++] = c;
   outName2[j] = '\0';
   linelength = j;
   SUCCESS = FALSE;
   counter_1 = 0;
   while (counter_1 < 3 && SUCCESS == FALSE) {
      ++counter_1;
      if (linelength == 0) {
         printf("\nNo name entered - please re-enter file name .... : ");
         j = 0;
         while((c = getchar()) != '\n')
            outName2[j++] = c;
         outName2[j] = '\0';
         linelength = j;
      }
      else {
         for(j = 0; j < linelength; j++) {
            c = line[j];
            if(c != ' ')
               SUCCESS = TRUE;
         }
      }
      if (SUCCESS == TRUE) {
         if (strcmp(outName1, outName2))
            SUCCESS = TRUE;
         else {
            SUCCESS = FALSE;
            printf("\nFile names identical - please re-enter file name : ");
            j = 0;
            while((c = getchar()) != '\n')
               outName2[j++] = c;
            outName2[j] = '\0';
            linelength = j;
         }
      }
      if (SUCCESS == TRUE) {
         if ((outFile2 = fopen(outName2, "w")) == NULL)
            SUCCESS = FALSE;
      }
   }
   if (SUCCESS == FALSE) {
      fprintf(stderr,"Program aborted - could not open file named");
      fprintf(stderr," %s.\n",outName2);
      exit(1);
   }

/*
 *
 *-------------------- Printing simulation details ----------------------
 *
 */

   systime = time(NULL);
   fprintf(outFile1,"DETAILS OF PROGRAM\n\n");
   fprintf(outFile1,"Program        %s %s\n",programName, version);
   fprintf(outFile1,"Copyright      %s\n",copyright);
   fprintf(outFile1,"Output file    %s\n",outName2);
   fprintf(outFile1,"Time           %s\n",ctime(&systime));
   fprintf(outFile1,"\nPROPERTIES OF THE TREE\n");
   fprintf(outFile1,"\nEdges entered in terms of ");
   if(use_time == 'Y')
      fprintf(outFile1,"time (units)\n");
   else
      fprintf(outFile1,"rate (substitutions per site)\n");
   fprintf(outFile1,"\nLength of a   %8.5f", div_a);
   fprintf(outFile1,"\nLength of b   %8.5f", div_b);
   fprintf(outFile1,"\nLength of c   %8.5f", div_c);
   fprintf(outFile1,"\nLength of d   %8.5f", div_d);
   fprintf(outFile1,"\nLength of e   %8.5f", div_e);
   fprintf(outFile1,"\nLength of f   %8.5f", div_f);
   fprintf(outFile1,"\n\nThe tree used with the Monte Carlo simulation,");
   fprintf(outFile1,"\nwritten in the Newick format with edge lengths");
   fprintf(outFile1,"\ngiven in terms of (1) time or (2) average rate");
   fprintf(outFile1,"\nof nucleotide substitution per site:\n");
   if(use_time == 'Y') {
      fprintf(outFile1,"\n(1) - ");
      fprintf(outFile1,"((SeqA:%f,SeqB:%f):%f,(SeqC:%f,SeqD:%f):%f);\n",div_a, div_b, div_e, div_c, div_d, div_f);
      fprintf(outFile1,"\n(2) - ");
      fprintf(outFile1,"((SeqA:%f,SeqB:%f):%f,",div_a * rate1, div_b * rate2, div_e * rate5);
      fprintf(outFile1,"(SeqC:%f,SeqD:%f):%f);",div_c * rate3, div_d * rate4, div_f * rate6);
   }
   else {
      fprintf(outFile1,"\n(1) - ");
      fprintf(outFile1,"((SeqA:%f,SeqB:%f):%f,",div_a / rate1, div_b / rate2, div_e / rate5);
      fprintf(outFile1,"(SeqC:%f,SeqD:%f):%f);\n",div_c / rate3, div_d / rate4, div_f / rate6);

      fprintf(outFile1,"\n(2) - ");
      fprintf(outFile1,"((SeqA:%f,SeqB:%f):%f,(SeqC:%f,SeqD:%f):%f);",div_a, div_b, div_e, div_c, div_d, div_f);
   }
   fprintf(outFile1,"\n\n\nAVERAGE RATES OF CHANGE ALONG THE EDGES\n");
   fprintf(outFile1,"\nRate along edge a     : %8.5f",rate1);
   fprintf(outFile1,"\nRate along edge b     : %8.5f",rate2);
   fprintf(outFile1,"\nRate along edge c     : %8.5f",rate3);
   fprintf(outFile1,"\nRate along edge d     : %8.5f",rate4);
   fprintf(outFile1,"\nRate along edge e     : %8.5f",rate5);
   fprintf(outFile1,"\nRate along edge f     : %8.5f",rate6);

   fprintf(outFile1,"\n\n\nOTHER RELEVANT INFORMATION\n");
   fprintf(outFile1,"\nSequence length       : %5d", length);
   fprintf(outFile1,"\nNumber of cycles      : %5d", cycles);
   fprintf(outFile1,"\nSeed                  : %5ld", seed);
   fprintf(outFile1,"\nMultiple hits         : ");
   if (choice == 'Y')
      fprintf(outFile1,"yes");
   else
      fprintf(outFile1,"no");
   fprintf(outFile1,"\nOrder of nucleotides  : A, C, G & T");

   fprintf(outFile1,"\n\n\nPROPERTIES OF THE ANCESTRAL SEQUENCE\n");
   fprintf(outFile1,"\nFrequency of A   %8.5f", P0[0]);
   fprintf(outFile1,"\nFrequency of T   %8.5f", P0[1]);
   fprintf(outFile1,"\nFrequency of C   %8.5f", P0[2]);
   fprintf(outFile1,"\nFrequency of G   %8.5f", P0[3]);

   fprintf(outFile1,"\n\n\nPROPERTIES OF THE SUBSTITUTION MODELS");

   fprintf(outFile1,"\n\nModel Ra -- Nucleotide frequencies");
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f%10.5f",P1[0],P1[1],P1[2],P1[3]);

   fprintf(outFile1,"\n\nModel Ra -- Conditional rates");
   fprintf(outFile1,"\n   -------%10.5f%10.5f%10.5f",Ra[0][1]/P1[1],Ra[0][2]/P1[2],Ra[0][3]/P1[3]);
   fprintf(outFile1,"\n%10.5f   -------%10.5f%10.5f",Ra[1][0]/P1[0],Ra[1][2]/P1[2],Ra[1][3]/P1[3]);
   fprintf(outFile1,"\n%10.5f%10.5f   -------%10.5f",Ra[2][0]/P1[0],Ra[2][1]/P1[1],Ra[2][3]/P1[3]);
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f   -------",Ra[3][0]/P1[0],Ra[3][1]/P1[1],Ra[3][2]/P1[2]);
   
   fprintf(outFile1,"\n\nModel Ra -- Rates of change along edge a\n\n");
   for(i = 0; i < 4; ++i) {
     for(j = 0; j < 4; ++j) {
       fprintf(outFile1,"%10.5f",Ra[i][j]);
     }
     fprintf(outFile1,"\n");
   }

   fprintf(outFile1,"\n\nModel Rb -- Nucleotide frequencies");
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f%10.5f",P2[0],P2[1],P2[2],P2[3]);

   fprintf(outFile1,"\n\nModel Rb -- Conditional rates");
   fprintf(outFile1,"\n   -------%10.5f%10.5f%10.5f",Rb[0][1]/P2[1],Rb[0][2]/P2[2],Rb[0][3]/P2[3]);
   fprintf(outFile1,"\n%10.5f   -------%10.5f%10.5f",Rb[1][0]/P2[0],Rb[1][2]/P2[2],Rb[1][3]/P2[3]);
   fprintf(outFile1,"\n%10.5f%10.5f   -------%10.5f",Rb[2][0]/P2[0],Rb[2][1]/P2[1],Rb[2][3]/P2[3]);
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f   -------",Rb[3][0]/P2[0],Rb[3][1]/P2[1],Rb[3][2]/P2[2]);

   fprintf(outFile1,"\n\nModel Rb -- Rates of change along edge b\n\n");
   for(i = 0; i < 4; ++i) {
     for(j = 0; j < 4; ++j) {
       fprintf(outFile1,"%10.5f",Rb[i][j]);
     }
     fprintf(outFile1,"\n");
   }

   fprintf(outFile1,"\n\nModel Rc -- Nucleotide frequencies");
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f%10.5f",P3[0],P3[1],P3[2],P3[3]);

   fprintf(outFile1,"\n\nModel Rc -- Conditional rates");
   fprintf(outFile1,"\n   -------%10.5f%10.5f%10.5f",Rc[0][1]/P3[1],Rc[0][2]/P3[2],Rc[0][3]/P3[3]);
   fprintf(outFile1,"\n%10.5f   -------%10.5f%10.5f",Rc[1][0]/P3[0],Rc[1][2]/P3[2],Rc[1][3]/P3[3]);
   fprintf(outFile1,"\n%10.5f%10.5f   -------%10.5f",Rc[2][0]/P3[0],Rc[2][1]/P3[1],Rc[2][3]/P3[3]);
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f   -------",Rc[3][0]/P3[0],Rc[3][1]/P3[1],Rc[3][2]/P3[2]);

   fprintf(outFile1,"\n\nModel Rc -- Rates of change along edge c\n\n");
   for(i = 0; i < 4; ++i) {
     for(j = 0; j < 4; ++j) {
       fprintf(outFile1,"%10.5f",Rc[i][j]);
     }
     fprintf(outFile1,"\n");
   }

   fprintf(outFile1,"\n\nModel Rd -- Nucleotide frequencies");
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f%10.5f",P4[0],P4[1],P4[2],P4[3]);

   fprintf(outFile1,"\n\nModel Rd -- Conditional rates");
   fprintf(outFile1,"\n   -------%10.5f%10.5f%10.5f",Rd[0][1]/P4[1],Rd[0][2]/P4[2],Rd[0][3]/P4[3]);
   fprintf(outFile1,"\n%10.5f   -------%10.5f%10.5f",Rd[1][0]/P4[0],Rd[1][2]/P4[2],Rd[1][3]/P4[3]);
   fprintf(outFile1,"\n%10.5f%10.5f   -------%10.5f",Rd[2][0]/P4[0],Rd[2][1]/P4[1],Rd[2][3]/P4[3]);
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f   -------",Rd[3][0]/P4[0],Rd[3][1]/P4[1],Rd[3][2]/P4[2]);

   fprintf(outFile1,"\n\nModel Rd -- Rates of change along edge d\n\n");
   for(i = 0; i < 4; ++i) {
     for(j = 0; j < 4; ++j) {
       fprintf(outFile1,"%10.5f",Rd[i][j]);
     }
     fprintf(outFile1,"\n");
   }

   fprintf(outFile1,"\n\nModel Re -- Nucleotide frequencies");
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f%10.5f",P5[0],P5[1],P5[2],P5[3]);

   fprintf(outFile1,"\n\nModel Re -- Conditional rates");
   fprintf(outFile1,"\n   -------%10.5f%10.5f%10.5f",Re[0][1]/P5[1],Re[0][2]/P5[2],Re[0][3]/P5[3]);
   fprintf(outFile1,"\n%10.5f   -------%10.5f%10.5f",Re[1][0]/P5[0],Re[1][2]/P5[2],Re[1][3]/P5[3]);
   fprintf(outFile1,"\n%10.5f%10.5f   -------%10.5f",Re[2][0]/P5[0],Re[2][1]/P5[1],Re[2][3]/P5[3]);
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f   -------",Re[3][0]/P5[0],Re[3][1]/P5[1],Re[3][2]/P5[2]);

   fprintf(outFile1,"\n\nModel Re -- Rates of change along edge e\n\n");
   for(i = 0; i < 4; ++i) {
     for(j = 0; j < 4; ++j) {
       fprintf(outFile1,"%10.5f",Re[i][j]);
     }
     fprintf(outFile1,"\n");
   }
   fprintf(outFile1,"\n\nModel Rf -- Nucleotide frequencies");
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f%10.5f",P6[0],P6[1],P6[2],P6[3]);

   fprintf(outFile1,"\n\nModel Rf -- Conditional rates");
   fprintf(outFile1,"\n   -------%10.5f%10.5f%10.5f",Rf[0][1]/P6[1],Rf[0][2]/P6[2],Rf[0][3]/P6[3]);
   fprintf(outFile1,"\n%10.5f   -------%10.5f%10.5f",Rf[1][0]/P6[0],Rf[1][2]/P6[2],Rf[1][3]/P6[3]);
   fprintf(outFile1,"\n%10.5f%10.5f   -------%10.5f",Rf[2][0]/P6[0],Rf[2][1]/P6[1],Rf[2][3]/P6[3]);
   fprintf(outFile1,"\n%10.5f%10.5f%10.5f   -------",Rf[3][0]/P6[0],Rf[3][1]/P6[1],Rf[3][2]/P6[2]);

   fprintf(outFile1,"\n\nModel Rf -- Rates of change along edge f\n\n");
   for(i = 0; i < 4; ++i) {
     for(j = 0; j < 4; ++j) {
       fprintf(outFile1,"%10.5f",Rf[i][j]);
     }
     fprintf(outFile1,"\n");
   }

/*
 *
 *------- Converts rate matrices to threshold matrices --------
 *
 */

   for(i = 0; i < 4; ++i) {
     Ra[i][i] = 1.0 + Ra[i][i];
     Rb[i][i] = 1.0 + Rb[i][i];
     Rc[i][i] = 1.0 + Rc[i][i];
     Rd[i][i] = 1.0 + Rd[i][i];
     Re[i][i] = 1.0 + Re[i][i];
     Rf[i][i] = 1.0 + Rf[i][i];
   }
   for(i = 0; i < 4; ++i) {
     T1[i][0] = Ra[i][0];
     T1[i][1] = Ra[i][0] + Ra[i][1];
     T1[i][2] = Ra[i][0] + Ra[i][1] + Ra[i][2];
     T1[i][3] = Ra[i][0] + Ra[i][1] + Ra[i][2] + Ra[i][3];
     T2[i][0] = Rb[i][0];
     T2[i][1] = Rb[i][0] + Rb[i][1];
     T2[i][2] = Rb[i][0] + Rb[i][1] + Rb[i][2];
     T2[i][3] = Rb[i][0] + Rb[i][1] + Rb[i][2] + Rb[i][3];
     T3[i][0] = Rc[i][0];
     T3[i][1] = Rc[i][0] + Rc[i][1];
     T3[i][2] = Rc[i][0] + Rc[i][1] + Rc[i][2];
     T3[i][3] = Rc[i][0] + Rc[i][1] + Rc[i][2] + Rc[i][3];
     T4[i][0] = Rd[i][0];
     T4[i][1] = Rd[i][0] + Rd[i][1];
     T4[i][2] = Rd[i][0] + Rd[i][1] + Rd[i][2];
     T4[i][3] = Rd[i][0] + Rd[i][1] + Rd[i][2] + Rd[i][3];
     T5[i][0] = Re[i][0];
     T5[i][1] = Re[i][0] + Re[i][1];
     T5[i][2] = Re[i][0] + Re[i][1] + Re[i][2];
     T5[i][3] = Re[i][0] + Re[i][1] + Re[i][2] + Re[i][3];
     T6[i][0] = Rf[i][0];
     T6[i][1] = Rf[i][0] + Rf[i][1];
     T6[i][2] = Rf[i][0] + Rf[i][1] + Rf[i][2];
     T6[i][3] = Rf[i][0] + Rf[i][1] + Rf[i][2] + Rf[i][3];

   }

/*
 *
 *The threshold matrices are needed because random numbers are used to
 *determine whether a given nucleotide remains the same or changes.
 *
 *The value of the threshold matrix becomes clear when a random number
 *is provided.  If the random number is below T[x][0] the nucleotide
 *at the site becomes an A; if it is larger than T[x][0] and smaller
 *than T[x][1] the nucleotide becomes a C; and so forth.
 *
 *The combined use of threshold matrices and random numbers ensure the
 *process of change is stochastic.
 *
 */

/*
 *
 *-------------------- Printing simulation details ----------------------
 *
 */

   fprintf(outFile1,"\n\nTHRESSHOLD MATRICES\n");
   fprintf(outFile1,"\nThresholds used under Model Ra (edge a)\n\n");
   for(i = 0; i < 4; ++i) {
      for(j = 0; j < 4; ++j) {
         fprintf(outFile1,"%10.5f",T1[i][j]);
      }
      fprintf(outFile1,"\n");
   }
   fprintf(outFile1,"\nThresholds used under Model Rb (edge b)\n\n");
   for(i = 0; i < 4; ++i) {
      for(j = 0; j < 4; ++j) {
         fprintf(outFile1,"%10.5f",T2[i][j]);
      }
      fprintf(outFile1,"\n");
   }
   fprintf(outFile1,"\nThresholds used under Model Rc (edge c)\n\n");
   for(i = 0; i < 4; ++i) {
      for(j = 0; j < 4; ++j) {
         fprintf(outFile1,"%10.5f",T3[i][j]);
      }
      fprintf(outFile1,"\n");
   }
   fprintf(outFile1,"\nThresholds used under Model Rd (edge d)\n\n");
   for(i = 0; i < 4; ++i) {
      for(j = 0; j < 4; ++j) {
         fprintf(outFile1,"%10.5f",T4[i][j]);
      }
      fprintf(outFile1,"\n");
   }
   fprintf(outFile1,"\nThresholds used under Model Re (edge e)\n\n");
   for(i = 0; i < 4; ++i) {
      for(j = 0; j < 4; ++j) {
         fprintf(outFile1,"%10.5f",T5[i][j]);
      }
      fprintf(outFile1,"\n");
   }
   fprintf(outFile1,"\nThresholds used under Model Rf (edge f)\n\n");
   for(i = 0; i < 4; ++i) {
      for(j = 0; j < 4; ++j) {
         fprintf(outFile1,"%10.5f",T6[i][j]);
      }
      fprintf(outFile1,"\n");
   }
   fflush(outFile1);

/*
 *
 *-------------------- Printing simulation details ----------------------
 *
 */

   fprintf(outFile1,"\n\n\nRESULTS OF SIMULATION");
   fprintf(outFile1,"\n\nColumn ..... [1]   Dif. in GC content (SeqA vs SeqB)");
   fprintf(outFile1,"\n             [2]   Dif. in GC content (SeqA vs SeqC)");
   fprintf(outFile1,"\n             [3]   Dif. in GC content (SeqA vs SeqD)");
   fprintf(outFile1,"\n             [4]   Dif. in GC content (SeqB vs SeqC)");
   fprintf(outFile1,"\n             [5]   Dif. in GC content (SeqB vs SeqD)");
   fprintf(outFile1,"\n             [6]   Dif. in GC content (SeqC vs SeqD)");
   fprintf(outFile1,"\n             [7]   Constant sites");
   fprintf(outFile1,"\n             [8]   Split A (A|BCD)");
   fprintf(outFile1,"\n             [9]   Split B (B|ACD)");
   fprintf(outFile1,"\n            [10]   Split C (C|ABD)");
   fprintf(outFile1,"\n            [11]   Split D (D|ABC)");
   fprintf(outFile1,"\n            [12]   Split E (AB|CD)");
   fprintf(outFile1,"\n            [13]   Split F (AC|BD)");
   fprintf(outFile1,"\n            [14]   Split G (AD|BC)");
   fprintf(outFile1,"\n            [15]   Hypervariable sites");
   
   fprintf(outFile1,"\n\n     [1]       [2]       [3]       [4]       [5]       [6]       [7]       [8]");
   fprintf(outFile1,"       [9]      [10]      [11]      [12]      [13]      [14]      [15]\n\n");


   printf("\n\n_______________ Simulation in Progress _______________\n\n");

   srand(seed);
   sumDif_GC_AB = 0.0;
   sumDif_GC_AC = 0.0;
   sumDif_GC_AD = 0.0;
   for(i = 0; i < MAX_HITS; i++) {
      tot_sum_hits[i] = 0;
   }

   
/*
 *
 *The loop starts the generation of sequences.
 *
 *First, a sequence of ancestral nucleotides is drawn at random from
 *values given by the user.  The ancestral sequence is placed in the
 *arrays called seqA and seqC.
 *Initially, SeqA and SeqC hold the sequences that evolve along the
 *edges e and f, respectively; later, after the second divergences,
 *they hold the sequences evolving along the edges a and c.  At the
 *point of the second divergences, the nucleotide sequence held in
 *SeqA is copied into the array called SeqB and the sequence held in
 *SeqC is copied into the array SeqD, and the final evolution along
 *edges b and d is completed.
 *
 */

   for (k = 0; k < cycles; ++k) {
      for(i = 0; i < MAX_HITS; i++) {
         sum_hits[i] = 0;
      }
      i = 0;
      while (i < length) {
         SUCCESS = FALSE;
         hits[i] = 0;
         sub = ran2(&seed);
         if(SUCCESS == FALSE && sub <= P0[0]) {
            seqA[i] = seqC[i] = 'A';
            SUCCESS = TRUE;
            i++;
         }
         if(SUCCESS == FALSE && P0[0] < sub && sub <= (P0[0] + P0[1])) {
            seqA[i] = seqC[i] = 'C';
            SUCCESS = TRUE;
            i++;
         }
         if(SUCCESS == FALSE && (P0[0] + P0[1]) < sub && sub <= (P0[0] + P0[1] + P0[2])) {
            seqA[i] = seqC[i] = 'G';
            SUCCESS = TRUE;
            i++;
         }
         if(SUCCESS == FALSE && (P0[0] + P0[1] + P0[2]) < sub && sub <= 1.0) {
            seqA[i] = seqC[i] = 'T';
            SUCCESS = TRUE;
            i++;
         }
      }
      
/*
 *
 *Evolution from the Origin along edge e towards A and B
 *
 */      

      counter = 0;
      if(use_time == 'Y')
         number = length * div_e;
      else
         number = length * div_e / rate5;
      while (counter < number) {
         if(choice == 'Y') {
            pos = length * ran2(&seed);
         }
         else {
            do {
               pos = length * ran2(&seed);
            } while (hits[pos] >= 1);
         }
         sub = ran2(&seed);
         ++counter;
         if('A' == seqA[pos]) {
            if(T5[0][0] < sub && sub <= T5[0][1]) {
               seqA[pos] = 'C';
               ++hits[pos];
            }
            if(T5[0][1] < sub && sub <= T5[0][2]) {
               seqA[pos] = 'G';
               ++hits[pos];
            }
            if(T5[0][2] < sub && sub <= T5[0][3]) {
               seqA[pos] = 'T';
               ++hits[pos];
            }
         }
         if('C' == seqA[pos]) {
            if(sub <= T5[1][0]) {
               seqA[pos] = 'A';
               ++hits[pos];
            }
            if(T5[1][1] < sub && sub <= T5[1][2]) {
               seqA[pos] = 'G';
               ++hits[pos];
            }
            if(T5[1][2] < sub && sub <= T5[1][3]) {
               seqA[pos] = 'T';
               ++hits[pos];
            }
         }
         if('G' == seqA[pos]) {
            if(sub <= T5[2][0]) {
               seqA[pos] = 'A';
               ++hits[pos];
            }
            if(T5[2][0] < sub && sub <= T5[2][1]) {
               seqA[pos] = 'C';
               ++hits[pos];
            }
            if(T5[2][2] < sub && sub <= T5[2][3]) {
               seqA[pos] = 'T';
               ++hits[pos];
            }
         }
         if('T' == seqA[pos]) {
            if(sub <= T5[3][0]) {
               seqA[pos] = 'A';
               ++hits[pos];
            }
            if(T5[3][0] < sub && sub <= T5[3][1]) {
               seqA[pos] = 'C';
               ++hits[pos];
            }
            if(T5[3][1] < sub && sub <= T5[3][2]) {
               seqA[pos] = 'G';
               ++hits[pos];
            }
         }
      }

/*
 *
 *Copy seqA into seqB before further simulation
 *
 */

      for(i = 0; i < length; ++i)
         seqB[i] = seqA[i];
      
      
/*
 *
 *Evolution along edge a towards A
 *
 */      

      counter = 0;
      if(use_time == 'Y')
         number = length * div_a;
      else
         number = length * div_a / rate1;
      while (counter < number) {
         if(choice == 'Y') {
            pos = length * ran2(&seed);
         }
         else {
            do {
               pos = length * ran2(&seed);
            } while (hits[pos] >= 1);
         }
         sub = ran2(&seed);
         ++counter;
         if('A' == seqA[pos]) {
            if(T1[0][0] < sub && sub <= T1[0][1]) {
               seqA[pos] = 'C';
               ++hits[pos];
            }
            if(T1[0][1] < sub && sub <= T1[0][2]) {
               seqA[pos] = 'G';
               ++hits[pos];
            }
            if(T1[0][2] < sub && sub <= T1[0][3]) {
               seqA[pos] = 'T';
               ++hits[pos];
            }
         }
         if('C' == seqA[pos]) {
            if(sub <= T1[1][0]) {
               seqA[pos] = 'A';
               ++hits[pos];
            }
            if(T1[1][1] < sub && sub <= T1[1][2]) {
               seqA[pos] = 'G';
               ++hits[pos];
            }
            if(T1[1][2] < sub && sub <= T1[1][3]) {
               seqA[pos] = 'T';
               ++hits[pos];
            }
         }
         if('G' == seqA[pos]) {
            if(sub <= T1[2][0]) {
               seqA[pos] = 'A';
               ++hits[pos];
            }
            if(T1[2][0] < sub && sub <= T1[2][1]) {
               seqA[pos] = 'C';
               ++hits[pos];
            }
            if(T1[2][2] < sub && sub <= T1[2][3]) {
               seqA[pos] = 'T';
               ++hits[pos];
            }
         }
         if('T' == seqA[pos]) {
            if(sub <= T1[3][0]) {
               seqA[pos] = 'A';
               ++hits[pos];
            }
            if(T1[3][0] < sub && sub <= T1[3][1]) {
               seqA[pos] = 'C';
               ++hits[pos];
            }
            if(T1[3][1] < sub && sub <= T1[3][2]) {
               seqA[pos] = 'G';
               ++hits[pos];
            }
         }
      }
      
      
/*
 *
 *Evolution along edge b towards B
 *
 */      

      counter = 0;
      if(use_time == 'Y')
         number = length * div_b;
      else
         number = length * div_b / rate2;
      while (counter < number) {
         if(choice == 'Y') {
            pos = length * ran2(&seed);
         }
         else {
            do {
               pos = length * ran2(&seed);
            } while (hits[pos] >= 1);
         }
         sub = ran2(&seed);
         ++counter;
         if('A' == seqB[pos]) {
            if(T2[0][0] < sub && sub <= T2[0][1]) {
               seqB[pos] = 'C';
               ++hits[pos];
            }
            if(T2[0][1] < sub && sub <= T2[0][2]) {
               seqB[pos] = 'G';
               ++hits[pos];
            }
            if(T2[0][2] < sub && sub <= T2[0][3]) {
               seqB[pos] = 'T';
               ++hits[pos];
            }
         }
         if('C' == seqB[pos]) {
            if(sub <= T2[1][0]) {
               seqB[pos] = 'A';
               ++hits[pos];
            }
            if(T2[1][1] < sub && sub <= T2[1][2]) {
               seqB[pos] = 'G';
               ++hits[pos];
            }
            if(T2[1][2] < sub && sub <= T2[1][3]) {
               seqB[pos] = 'T';
               ++hits[pos];
            }
         }
         if('G' == seqB[pos]) {
            if(sub <= T2[2][0]) {
               seqB[pos] = 'A';
               ++hits[pos];
            }
            if(T2[2][0] < sub && sub <= T2[2][1]) {
               seqB[pos] = 'C';
               ++hits[pos];
            }
            if(T2[2][2] < sub && sub <= T2[2][3]) {
               seqB[pos] = 'T';
               ++hits[pos];
            }
         }
         if('T' == seqB[pos]) {
            if(sub <= T2[3][0]) {
               seqB[pos] = 'A';
               ++hits[pos];
            }
            if(T2[3][0] < sub && sub <= T2[3][1]) {
               seqB[pos] = 'C';
               ++hits[pos];
            }
            if(T2[3][1] < sub && sub <= T2[3][2]) {
               seqB[pos] = 'G';
               ++hits[pos];
            }
         }
      }


/*
 *
 *Evolution from the Origin along edge f towards C and D
 *
 */      

      counter = 0;
      if(use_time == 'Y')
         number = length * div_f;
      else
         number = length * div_f / rate6;
      while (counter < number) {
         if(choice == 'Y') {
            pos = length * ran2(&seed);
         }
         else {
            do {
               pos = length * ran2(&seed);
            } while (hits[pos] >= 1);
         }
         sub = ran2(&seed);
         ++counter;
         if('A' == seqC[pos]) {
            if(T6[0][0] < sub && sub <= T6[0][1]) {
               seqC[pos] = 'C';
               ++hits[pos];
            }
            if(T6[0][1] < sub && sub <= T6[0][2]) {
               seqC[pos] = 'G';
               ++hits[pos];
            }
            if(T6[0][2] < sub && sub <= T6[0][3]) {
               seqC[pos] = 'T';
               ++hits[pos];
            }
         }
         if('C' == seqC[pos]) {
            if(sub <= T6[1][0]) {
               seqC[pos] = 'A';
               ++hits[pos];
            }
            if(T6[1][1] < sub && sub <= T6[1][2]) {
               seqC[pos] = 'G';
               ++hits[pos];
            }
            if(T6[1][2] < sub && sub <= T6[1][3]) {
               seqC[pos] = 'T';
               ++hits[pos];
            }
         }
         if('G' == seqC[pos]) {
            if(sub <= T6[2][0]) {
               seqC[pos] = 'A';
               ++hits[pos];
            }
            if(T6[2][0] < sub && sub <= T6[2][1]) {
               seqC[pos] = 'C';
               ++hits[pos];
            }
            if(T6[2][2] < sub && sub <= T6[2][3]) {
               seqC[pos] = 'T';
               ++hits[pos];
            }
         }
         if('T' == seqC[pos]) {
            if(sub <= T6[3][0]) {
               seqC[pos] = 'A';
               ++hits[pos];
            }
            if(T6[3][0] < sub && sub <= T6[3][1]) {
               seqC[pos] = 'C';
               ++hits[pos];
            }
            if(T6[3][1] < sub && sub <= T6[3][2]) {
               seqC[pos] = 'G';
               ++hits[pos];
            }
         }
      }
      
/*
 *
 *Copy seqD into seqC before further simulation
 *
 */

      for(i = 0; i < length; ++i)
         seqD[i] = seqC[i];

/*
 *
 *Evolution from edge c towards C
 *
 */      

      counter = 0;
      if(use_time == 'Y')
         number = length * div_c;
      else
         number = length * div_c / rate3;
      while (counter < number) {
         if(choice == 'Y') {
            pos = length * ran2(&seed);
         }
         else {
            do {
               pos = length * ran2(&seed);
            } while (hits[pos] >= 1);
         }
         sub = ran2(&seed);
         ++counter;
         if('A' == seqC[pos]) {
            if(T3[0][0] < sub && sub <= T3[0][1]) {
               seqC[pos] = 'C';
               ++hits[pos];
            }
            if(T3[0][1] < sub && sub <= T3[0][2]) {
               seqC[pos] = 'G';
               ++hits[pos];
            }
            if(T3[0][2] < sub && sub <= T3[0][3]) {
               seqC[pos] = 'T';
               ++hits[pos];
            }
         }
         if('C' == seqC[pos]) {
            if(sub <= T3[1][0]) {
               seqC[pos] = 'A';
               ++hits[pos];
            }
            if(T3[1][1] < sub && sub <= T3[1][2]) {
               seqC[pos] = 'G';
               ++hits[pos];
            }
            if(T3[1][2] < sub && sub <= T3[1][3]) {
               seqC[pos] = 'T';
               ++hits[pos];
            }
         }
         if('G' == seqC[pos]) {
            if(sub <= T3[2][0]) {
               seqC[pos] = 'A';
               ++hits[pos];
            }
            if(T3[2][0] < sub && sub <= T3[2][1]) {
               seqC[pos] = 'C';
               ++hits[pos];
               SUCCESS = TRUE;
            }
            if(T3[2][2] < sub && sub <= T3[2][3]) {
               seqC[pos] = 'T';
               ++hits[pos];
            }
         }
         if('T' == seqC[pos]) {
            if(sub <= T3[3][0]) {
               seqC[pos] = 'A';
               ++hits[pos];
            }
            if(T3[3][0] < sub && sub <= T3[3][1]) {
               seqC[pos] = 'C';
               ++hits[pos];
            }
            if(T3[3][1] < sub && sub <= T3[3][2]) {
               seqC[pos] = 'G';
               ++hits[pos];
            }
         }
      }

/*
 *
 *Evolution along edge d towards D
 *
 */      

      counter = 0;
      if(use_time == 'Y')
         number = length * div_d;
      else
         number = length * div_d / rate4;
      while (counter < number) {
         if(choice == 'Y') {
            pos = length * ran2(&seed);
         }
         else {
            do {
               pos = length * ran2(&seed);
            } while (hits[pos] >= 1);
         }
         sub = ran2(&seed);
         ++counter;
         if('A' == seqD[pos]) {
            if(T4[0][0] < sub && sub <= T4[0][1]) {
               seqD[pos] = 'C';
               ++hits[pos];
            }
            if(T4[0][1] < sub && sub <= T4[0][2]) {
               seqD[pos] = 'G';
               ++hits[pos];
            }
            if(T4[0][2] < sub && sub <= T4[0][3]) {
               seqD[pos] = 'T';
               ++hits[pos];
            }
         }
         if('C' == seqD[pos]) {
            if(sub <= T4[1][0]) {
               seqD[pos] = 'A';
               ++hits[pos];
            }
            if(T4[1][1] < sub && sub <= T4[1][2]) {
               seqD[pos] = 'G';
               ++hits[pos];
            }
            if(T4[1][2] < sub && sub <= T4[1][3]) {
               seqD[pos] = 'T';
               ++hits[pos];
            }
         }
         if('G' == seqD[pos]) {
            if(sub <= T4[2][0]) {
               seqD[pos] = 'A';
               ++hits[pos];
            }
            if(T4[2][0] < sub && sub <= T4[2][1]) {
               seqD[pos] = 'C';
               ++hits[pos];
            }
            if(T4[2][2] < sub && sub <= T4[2][3]) {
               seqD[pos] = 'T';
               ++hits[pos];
            }
         }
         if('T' == seqD[pos]) {
            if(sub <= T4[3][0]) {
               seqD[pos] = 'A';
               ++hits[pos];
            }
            if(T4[3][0] < sub && sub <= T4[3][1]) {
               seqD[pos] = 'C';
               ++hits[pos];
            }
            if(T4[3][1] < sub && sub <= T4[3][2]) {
               seqD[pos] = 'G';
               ++hits[pos];
            }
         }
      }

/*
 *
 *Simulated sequences will now be printed to outFile2
 *
 */

      fprintf(outFile2,"   4  %d\n",length);
      fprintf(outFile2,"Seq1      ");
      for (i = 0; i < length; i++)
        fprintf(outFile2,"%c",seqA[i]);
      fprintf(outFile2,"\n");
      fflush(outFile2);
      fprintf(outFile2,"Seq2      ");
      for (i = 0; i < length; i++)
        fprintf(outFile2,"%c",seqB[i]);
      fprintf(outFile2,"\n");
      fflush(outFile2);
      fprintf(outFile2,"Seq3      ");
      for (i = 0; i < length; i++)
        fprintf(outFile2,"%c",seqC[i]);
      fprintf(outFile2,"\n");
      fflush(outFile2);
      fprintf(outFile2,"Seq4      ");
      for (i = 0; i < length; i++)
        fprintf(outFile2,"%c",seqD[i]);
      fprintf(outFile2,"\n");
      fflush(outFile2);

/*
 *
 *The output format is designed to allow the data to be analysed using
 *the PHYLIP program packages
 *
 */
 

/*
 *
 *From here the program proceeds to compare the nucleotide content of
 *the four sequences -- the comparison is done in a pairwise fashion.
 *
 */

/*
 *
 *Set elements in obsMat to 0.0
 *
 */ 

      for(j = 0; j < 4; j++) {
         for(i = 0; i < 4; i++)
            obsMat[i][j] = 0.0;
      }

/*
 *
 *Summing up matrix with nucleotides
 *
 */

      for(i = 0; i < length; i++) {
         if(seqA[i] == 'A')
            obsMat[0][0] = obsMat[0][0] + 1.0;
         if(seqA[i] == 'C')
            obsMat[0][1] = obsMat[0][1] + 1.0;
         if(seqA[i] == 'G')
            obsMat[0][2] = obsMat[0][2] + 1.0;
         if(seqA[i] == 'T')
            obsMat[0][3] = obsMat[0][3] + 1.0;
         if(seqB[i] == 'A')
            obsMat[1][0] = obsMat[1][0] + 1.0;
         if(seqB[i] == 'C')
            obsMat[1][1] = obsMat[1][1] + 1.0;
         if(seqB[i] == 'G')
            obsMat[1][2] = obsMat[1][2] + 1.0;
         if(seqB[i] == 'T')
            obsMat[1][3] = obsMat[1][3] + 1.0;
         if(seqC[i] == 'A')
            obsMat[2][0] = obsMat[2][0] + 1.0;
         if(seqC[i] == 'C')
            obsMat[2][1] = obsMat[2][1] + 1.0;
         if(seqC[i] == 'G')
            obsMat[2][2] = obsMat[2][2] + 1.0;
         if(seqC[i] == 'T')
            obsMat[2][3] = obsMat[2][3] + 1.0;
         if(seqD[i] == 'A')
            obsMat[3][0] = obsMat[3][0] + 1.0;
         if(seqD[i] == 'C')
            obsMat[3][1] = obsMat[3][1] + 1.0;
         if(seqD[i] == 'G')
            obsMat[3][2] = obsMat[3][2] + 1.0;
         if(seqD[i] == 'T')
            obsMat[3][3] = obsMat[3][3] + 1.0;
      }
/*
 *
 *Count the frequency of the three splits
 *
 */

      splitA = splitB = splitC = splitD = splitE = splitF = splitG = 0.0;
      constant_sites = hypervariable_sites= 0.0;
      for(i = 0; i < length; i++) {
         SUCCESS = FALSE;
         if(seqA[i] == seqB[i] && seqB[i] == seqC[i] && seqC[i] == seqD[i]) {
            constant_sites = constant_sites + 1.0;
            SUCCESS = TRUE;
         }
         if(seqA[i] != seqB[i] && seqB[i] == seqC[i] && seqC[i] == seqD[i]) {
            splitA = splitA + 1.0;
            SUCCESS = TRUE;
         }
         if(seqB[i] != seqA[i] && seqA[i] == seqC[i] && seqC[i] == seqD[i]) {
            splitB = splitB + 1.0;
            SUCCESS = TRUE;
         }
         if(seqC[i] != seqD[i] && seqD[i] == seqA[i] && seqA[i] == seqB[i]) {
            splitC = splitC + 1.0;
            SUCCESS = TRUE;
         }
         if(seqD[i] != seqC[i] && seqC[i] == seqA[i] && seqA[i] == seqB[i]) {
            splitD = splitD + 1.0;
            SUCCESS = TRUE;
         }
         if(seqA[i] == seqB[i] && seqB[i] != seqC[i] && seqC[i] == seqD[i]) {
            splitE = splitE + 1.0;
            SUCCESS = TRUE;
          }
         if(seqA[i] == seqC[i] && seqC[i] != seqD[i] && seqD[i] == seqB[i]) {
            splitF = splitF + 1.0;
            SUCCESS = TRUE;
         }
         if(seqA[i] == seqD[i] && seqD[i] != seqC[i] && seqC[i] == seqB[i]) {
            splitG = splitG + 1.0;
            SUCCESS = TRUE;
         }
         if(SUCCESS == FALSE)
            hypervariable_sites = hypervariable_sites + 1.0;
      }

/*
 *
 *Counter of number of simulations
 *
 */

      if(cycles > 1) {
         printf("\rNumber of cycles = %5d",k);
         fflush(NULL);
      }
      

/*
 *
 *Printing out differences in GC contents for all sequence pairs
 *
 */

      GC_pair_AB = (obsMat[0][1] + obsMat[0][2] - obsMat[1][1] - obsMat[1][2])/length;
      GC_pair_AC = (obsMat[0][1] + obsMat[0][2] - obsMat[2][1] - obsMat[2][2])/length;
      GC_pair_AD = (obsMat[0][1] + obsMat[0][2] - obsMat[3][1] - obsMat[3][2])/length;
      GC_pair_BC = (obsMat[1][1] + obsMat[1][2] - obsMat[2][1] - obsMat[2][2])/length;
      GC_pair_BD = (obsMat[1][1] + obsMat[1][2] - obsMat[3][1] - obsMat[3][2])/length;
      GC_pair_CD = (obsMat[2][1] + obsMat[2][2] - obsMat[3][1] - obsMat[3][2])/length;
      
      fprintf(outFile1,"%10.5f",GC_pair_AB);   
      fprintf(outFile1,"%10.5f",GC_pair_AC);   
      fprintf(outFile1,"%10.5f",GC_pair_AD);
      fprintf(outFile1,"%10.5f",GC_pair_BC);   
      fprintf(outFile1,"%10.5f",GC_pair_BD);   
      fprintf(outFile1,"%10.5f",GC_pair_CD);
      fprintf(outFile1,"%10.0f",constant_sites);
      fprintf(outFile1,"%10.0f",splitA);
      fprintf(outFile1,"%10.0f",splitB);
      fprintf(outFile1,"%10.0f",splitC);
      fprintf(outFile1,"%10.0f",splitD);
      fprintf(outFile1,"%10.0f",splitE);
      fprintf(outFile1,"%10.0f",splitF);
      fprintf(outFile1,"%10.0f",splitG);
      fprintf(outFile1,"%10.0f\n",hypervariable_sites);

      sumDif_GC_AB = sumDif_GC_AB + GC_pair_AB;
      sumDif_GC_AC = sumDif_GC_AC + GC_pair_AC;
      sumDif_GC_AD = sumDif_GC_AD + GC_pair_AD;
      sumDif_GC_BC = sumDif_GC_BC + GC_pair_BC;
      sumDif_GC_BD = sumDif_GC_BD + GC_pair_BD;
      sumDif_GC_CD = sumDif_GC_CD + GC_pair_CD;
      sum_constant_sites = sum_constant_sites + constant_sites;
      sum_splitA = sum_splitA + splitA;
      sum_splitB = sum_splitB + splitB;
      sum_splitC = sum_splitC + splitC;
      sum_splitD = sum_splitD + splitD;
      sum_splitE = sum_splitE + splitE;
      sum_splitF = sum_splitF + splitF;
      sum_splitG = sum_splitG + splitG;
      sum_hypervariable_sites = sum_hypervariable_sites + hypervariable_sites;
      for (i = 0; i < length; i++) {
         ++sum_hits[hits[i]];
      }
      for (j = 0; j < MAX_HITS; j++) {
         tot_sum_hits[j] = tot_sum_hits[j] + sum_hits[j];
      }
      
   } /* end of cycles loop */


/*
 *
 *Finishing off by printing summary statistics
 *
 */
   
   if(cycles > 1) {
      printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
      printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
      printf("\n\n_________________ Summary Statistics _________________\n\n");
      printf("Repeats                    = %5d\n",cycles);
      printf("\nMean dif. in GC content (Seq_A vs Seq_B) [1]    = %12.6f\n",sumDif_GC_AB / cycles);
      printf("\nMean dif. in GC content (Seq_A vs Seq_C) [2]    = %12.6f\n",sumDif_GC_AC / cycles);
      printf("\nMean dif. in GC content (Seq_A vs Seq_D) [3]    = %12.6f\n",sumDif_GC_AD / cycles);
      printf("\nMean dif. in GC content (Seq_B vs Seq_C) [4]    = %12.6f\n",sumDif_GC_BC / cycles);
      printf("\nMean dif. in GC content (Seq_B vs Seq_D) [5]    = %12.6f\n",sumDif_GC_BD / cycles);
      printf("\nMean dif. in GC content (Seq_C vs Seq_D) [6]    = %12.6f\n",sumDif_GC_CD / cycles);
      printf("\nMean number of constant sites                   = %12.6f\n",sum_constant_sites / cycles);
      printf("\nMean number of sites supporting split D (A|BCD) = %12.6f\n",sum_splitA / cycles);
      printf("\nMean number of sites supporting split E (B|ABD) = %12.6f\n",sum_splitB / cycles);
      printf("\nMean number of sites supporting split F (C|ABD) = %12.6f\n",sum_splitC / cycles);
      printf("\nMean number of sites supporting split G (D|ABC) = %12.6f\n",sum_splitD / cycles);
      printf("\nMean number of sites supporting split A (AB|CD) = %12.6f\n",sum_splitE / cycles);
      printf("\nMean number of sites supporting split B (AC|BD) = %12.6f\n",sum_splitF / cycles);
      printf("\nMean number of sites supporting split C (AD|BC) = %12.6f\n",sum_splitG / cycles);
      printf("\nMean number of hypervariable sites              = %12.6f\n",sum_hypervariable_sites / cycles);
      fprintf(outFile1,"\n\nAVERAGE VALUES\n");
      fprintf(outFile1,"\n\n     [1]       [2]       [3]       [4]       [5]       [6]       [7]       [8]");
      fprintf(outFile1,"       [9]      [10]      [11]      [12]      [13]      [14]      [15]\n\n");
      fprintf(outFile1,"%10.5f",sumDif_GC_AB / cycles);   
      fprintf(outFile1,"%10.5f",sumDif_GC_AC / cycles);   
      fprintf(outFile1,"%10.5f",sumDif_GC_AD / cycles);
      fprintf(outFile1,"%10.5f",sumDif_GC_BC / cycles);   
      fprintf(outFile1,"%10.5f",sumDif_GC_BD / cycles);   
      fprintf(outFile1,"%10.5f",sumDif_GC_CD / cycles);
      fprintf(outFile1,"%10.2f",sum_constant_sites / cycles);
      fprintf(outFile1,"%10.2f",sum_splitA / cycles);
      fprintf(outFile1,"%10.2f",sum_splitB / cycles);
      fprintf(outFile1,"%10.2f",sum_splitC / cycles);
      fprintf(outFile1,"%10.2f",sum_splitD / cycles);
      fprintf(outFile1,"%10.2f",sum_splitE / cycles);
      fprintf(outFile1,"%10.2f",sum_splitF / cycles);
      fprintf(outFile1,"%10.2f",sum_splitG / cycles);
      fprintf(outFile1,"%10.2f\n",sum_hypervariable_sites / cycles);

      printf("\n\nAverage number of sites with X hits\n\n");
      fprintf(outFile1,"\n\nAverage number of sites with X hits\n\n");
      printf("    X      # Sites      Percentage\n");
      fprintf(outFile1,"    X      # Sites      Percentage\n");
      for(j = 0; j < MAX_HITS; j++) {
         average = (float) tot_sum_hits[j] / cycles;
         percentage = 100 * (average / length);
         printf("   %2d\t%10.3f\t%10.3f\n", j, average, percentage);
         fprintf(outFile1,"   %2d\t%10.3f\t%10.3f\n", j, average, percentage);
      }


   }
   printf("\n\n________________ Simulation Completed ________________\n\n");
   printf("\a\a");
   fclose(outFile1);
   return (0);
}


/*
 *
 *---------------------- Declaration of Functions ---------------------------
 *
 */

/*
 *
 *Information and instructions provided
 *
 */

void Information(void)
{
   int   c;

   printf("\nProgram   - %s Version %s",programName, version);
   printf("\n\nAuthor    - Lars S. Jermiin");
 	printf("\n            CSIRO Entomology");
	printf("\n            GPO Box 1700");
	printf("\n            Canberra, ACT 2601, Australia.");
	printf("\n\n            School of Biological Sciences, A08");
   printf("\n            University of Sydney");
   printf("\n            Sydney, NSW 2006, Australia.");
   printf("\n\n            E-mail lars.jermiin@csiro.au");
   printf("\n\nCopyright - %s",copyright);
   printf("\n\nDetails   - www.bio.usyd.edu.au/~jermiin/hetero.htm."); 
   printf("\n\n\nType RETURN to continue...");
   c = getchar();
   Screen_flush();
}

/*
 *
 *Flushing input buffer
 *
 */

void Buffer_flush(void)
{
   int c;

   while ((c = getchar()) != '\n')
      ;
}

/*
 *
 *Flushing screen
 *
 */

void Screen_flush(void)
{
  int  c;
  
  for(c = 0; c < 60; c++) 
    printf("\n");
}

/*
 *
 *Screen drawing of tree with models and edge lengths
 *
 */
 
void Tree_drawing(int feature)
{
   printf("\n              A        B  C        D");
   printf("\n               \\      /    \\      /");
   printf("\n              ");
   if(feature == -1 || feature == 0) printf("a");
   else printf(" ");
   printf(" \\    / ");
   if(feature == -1 || feature == 0) printf("b");
   else printf(" ");
   printf("  ");
   if(feature == -1 || feature == 0) printf("c");
   else printf(" ");
   printf(" \\    / ");
   if(feature == -1 || feature == 0) printf("d");
   printf("\n              ");
   if(feature == -1 || feature == 1) printf("Ra");
   else printf("  ");
   printf(" \\  / ");
   if(feature == -1 || feature == 2) printf("Rb");
   else printf("  ");
   printf("  ");
   if(feature == -1 || feature == 3) printf("Rc");
   else printf("  ");
   printf(" \\  / ");
   if(feature == -1 || feature == 4) printf("Rd");
   printf("\n                  \\/          \\/");
   printf("\n                   \\          /");
   printf("\n                    \\        /");
   printf("\n                   ");
   if(feature == -1 || feature == 0) printf("e");
   else printf(" ");
   printf(" \\      / ");
   if(feature == -1 || feature == 0) printf("f");
   printf("\n                   ");
   if(feature == -1 || feature == 5) printf("Re");
   else printf("  ");
   printf(" \\    / ");
   if(feature == -1 || feature == 6) printf("Rf");
   printf("\n                       \\  /");
   printf("\n                        \\/");
   printf("\n                      Origin");
}

/*
 *
 *Random number generator from Numerical recipes (ran2)
 *
 */

float ran2(long *idum)
{
	int	j;
	long	k;
	static long	idum2 = 123456789;
	static long	iy = 0;
	static long	iv[NTAB];
	float	temp;
	
	if(*idum <= 0) {
		if(-(*idum) < 1)
			*idum = 1;
		else
			*idum = -(*idum);
		idum2 = (*idum);
		for(j = NTAB + 7; j >= 0; j--) {
			k = (*idum)/IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;
			if(*idum < 0)
				*idum += IM1;
			if(j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum = IA1 * (*idum - k * IQ1) - k * IR1;
	if (*idum < 0)
		*idum += IM1;
	k = idum2/IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
	if (idum2 < 0)
		idum2 += IM2;
	j = iy/NDIV;
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if(iy < 1)
		iy += IMM1;
	if((temp = AM * iy) > RNMX)
		return RNMX;
	else
		return temp;
}
