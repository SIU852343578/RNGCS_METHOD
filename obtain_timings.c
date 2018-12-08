//program to measure the actual running times of original TSP, original AUCTSP and reduced-size TSP, reduced-size AUCTSP, where the reduction in gene size is found by RNGCS

#include<stdio.h>
#include<math.h>
#include <errno.h>
#include<strings.h>
#include<stdlib.h>
#include <time.h>
#include <sys/resource.h>


#define GENECNT 13000
#define SUBJECTCNT 103
#define OVERFIT_MARGIN 0.15

//the following are macros to avoid problems with float inaccuracies

#define GREATER(A,B)    (A - B >= 0.000001) 
#define LESS(A,B)    ((A - B) < 0.000001) 
#define ABS(A)    ((A >= 0.000001) ? A : -A) 
#define EQUAL(A, B) ((A - B < 0.000001) && (B - A < 0.000001)) 


#define BETA 500 //Number of bins

#define FACTOR_INT_PART 1000 // Multiplicative factor for turning the significant decimal digits of AUC value to an integer, to avoid rounding errors.
                            //for example,  0.005/0.005 when truncated to an integer returns 0 whereas it should be 1
                            // Bin width is W = 0.5/BETA = 0.5/500 =  0.001 so FACTOR_INT_PART = 1000



typedef struct binstruc { int geneindex; float geneauc; struct binstruc *next;} binelement;

struct binarraystruc { binelement *list; int count;} BINARRAY[BETA];

binelement *p;


struct rusage timing0, timing1;


float geneauc;
int seed;

int PAIR1_TSP, PAIR2_TSP;


char GENE_NAME[GENECNT+1][50];
int SELECTED[GENECNT+1];//holds the indices of the genes selected

float DATA_HEALTHY[GENECNT+1][SUBJECTCNT], DATA_SICK[GENECNT+1][SUBJECTCNT];
int TOPPAIR1_TSP, TOPPAIR2_TSP;
int TOPPAIR1_AUC, TOPPAIR2_AUC;


int SIMULTRIES;






float TSP_PR_HEALTHYval, TSP_PR_SICKval;
float AUC_PR_HEALTHYval, AUC_PR_SICKval;

float TOPSCORE_TSP;
float TOPSCORE_AUC;

float TSP_SCOREval;
float AUC_SCOREval;
float ftemp;




int NGENESALLOWED;
int NHEALTHY, NSICK;



char geneid[100];

main(argc, argv)
int argc;
char *argv[];
{


int i, j, genecount, binindex, bintosplit, selected_i, try, gene, gi, gj, k, ik, originalNGENES;
float theta;

setbuf(stdout,NULL);


sscanf(argv[1],"%d", &originalNGENES);
sscanf(argv[2],"%d", &NGENESALLOWED);
sscanf(argv[3],"%d", &NHEALTHY);
sscanf(argv[4],"%d", &NSICK);
sscanf(argv[5],"%d", &TRIES);
sscanf(argv[6],"%d", &seed);



for(gi=1;gi<=originalNGENES;gi++){

scanf("%s", geneid);

        for(i=1;i <= NHEALTHY;i++){scanf("%f", &DATA_HEALTHY[gi][i]); }
        for(i=1;i <= NSICK;i++){scanf("%f", &DATA_SICK[gi][i]); }
	    


} //genes have been read in  


if(NGENESALLOWED >  originalNGENES) NGENESALLOWED = originalNGENES;


getrusage(RUSAGE_SELF, &timing0);//snapshot of the timer before the TRIES start

for(try=1; try<= TRIES; try++){//measure the time over a number of tries for better estimation, and report the average


	
//take into account  the cost of RNCGS i.e., the cost of the AUC filtering that would need to be done in order to obtain the allowed number of genes NGENESALLOWED
//This filtering is needed only when NGENESALLOWED < originalNGENES

if(NGENESALLOWED >=  originalNGENES){
	for(gi=1;gi<=NGENESALLOWED;gi++) SELECTED[gi] = gi; // this stores the gene indices selected.
}else{ // (NGENESALLOWED < originalNGENES)

//initialize each entry of BINARRAY to empty
for(i=0; i< BETA; i++) {BINARRAY[i].list = NULL; BINARRAY[i].count = 0;}

for(gi=1; gi < originalNGENES; gi++){

      geneauc= 0;
	for(k=1; k<= NHEALTHY; k++)
	for(j=1; j<= NSICK; j++){
	if(GREATER(DATA_HEALTHY[gi][k], DATA_SICK[gi][j])) geneauc = geneauc + 1;
	if(EQUAL(DATA_HEALTHY[gi][k], DATA_SICK[gi][j])) geneauc = geneauc + 0.5;
	}
	geneauc = geneauc /(NHEALTHY*NSICK);
	if(geneauc < 0.5) geneauc = 1 -geneauc;


//find the bin index to which this geneauc value belongs
binindex  = geneauc  * FACTOR_INT_PART - FACTOR_INT_PART/2 ; 
//note: the straightforward expression binindex = (geneauc -0.5)*FACTOR_INT_PART may give rounding errors!
//
//printf("binindex of gene %d with AUC= %f is %d\n", gi, geneauc, binindex);

p = malloc(sizeof(binelement)) ;
if( (p  == NULL)){ printf("MALLOC memory problem\n"); exit(1);}

//insert the new member to the top (convenient for speed: order does not matter!) of the list pointed to by BINARRAY[binindex] 
p->geneindex = gi;
p->geneauc = geneauc;
p->next = BINARRAY[binindex].list;
BINARRAY[binindex].list = p;
BINARRAY[binindex].count++;

}

//find the binindex to stop at, based on how many genes can be afforded, by starting from the larger bin indices

genecount = 0;

for(i=BETA-1; i>= 0; i--){
 genecount += BINARRAY[i].count;
 if(genecount > NGENESALLOWED) break;

} 

  bintosplit = i;
//the excess amount is genecount - BINARRAY[bintosplit].count
//therefore, go down the list of BINARRAY[bintosplit].list and pick only the allowed number of genes

//select all genes before the bintosplit
//mark also the minumun AUC value (theta) of all genes selected
//Note: all genes selected have AUC values >= theta 
//but there may be genes beloging to bin with index bintosplit that have AUC value >= theta 
//but they are not selected, since the genes in each bin are unsorted.

selected_i = 0;
theta = 1;
for(i=BETA-1; i>bintosplit; i--) {//print the members of list BINARRAY[i].list
   p = BINARRAY[i].list; 
   while(p != NULL){ 
      SELECTED[++selected_i] = p->geneindex; 
      if(p->geneauc < theta) theta = p->geneauc; 
      p = p->next;}
}

//select  now only the allowed number of genes from the list of BINARRAY[bintosplit].list

genecount -= BINARRAY[bintosplit].count; //the excess amount is genecount - BINARRAY[bintosplit].count

//go down the list of BINARRAY[bintosplit].list and pick only the allowed number of genes


p = BINARRAY[bintosplit].list; 
while(p != NULL){ 
   if(genecount == NGENESALLOWED) break;
   genecount++;
   SELECTED[++selected_i] = p->geneindex; 
   if(p->geneauc < theta) theta = p->geneauc;
   p = p->next;
}

}

//COMPUTE THE SCORES OF the pairs (keeping only the top one (ignoring any ties))


    //COMPUTE THE TSP SCORE 


	TOPSCORE_TSP = 0;

for(gi=1; gi < NGENESALLOWED; gi++){
   for(gj= gi+1; gj<= NGENESALLOWED; gj++){
	   



	 
     TSP_PR_HEALTHYval= 0;
     for(k=1; k<= NHEALTHY; k++){
         if(GREATER(DATA_HEALTHY[SELECTED[gi]][k], DATA_HEALTHY[SELECTED[gj]][k])) TSP_PR_HEALTHYval = TSP_PR_HEALTHYval + 1;
	 }
     TSP_PR_HEALTHYval= TSP_PR_HEALTHYval/NHEALTHY;
    
     TSP_PR_SICKval = 0;
     for(k=1; k<= NSICK; k++){
         if(GREATER(DATA_SICK[SELECTED[gi]][k] , DATA_SICK[SELECTED[gj]][k])) TSP_PR_SICKval = TSP_PR_SICKval + 1;
	 }
     TSP_PR_SICKval = TSP_PR_SICKval/NSICK;
         

     TSP_SCOREval = TSP_PR_HEALTHYval - TSP_PR_SICKval; if (TSP_SCOREval < 0) TSP_SCOREval = - TSP_SCOREval;




    if(GREATER(TSP_SCOREval , TOPSCORE_TSP)){
        TOPSCORE_TSP = TSP_SCOREval; 
        TOPPAIR1_TSP = SELECTED[gi];
        TOPPAIR2_TSP = SELECTED[gj];
	}
}
}
   printf("Top TSP is %d and %d with score %f\n", TOPPAIR1_TSP, TOPPAIR2_TSP, TOPSCORE_TSP); 

 

//Compute the AUCTSP score
//This is commented in/out to obtain the timings of AUCTSP instead of TSP and vice-versa

/*

	TOPSCORE_AUC = 0;

for(gi=1; gi < NGENESALLOWED; gi++){
   for(gj= gi+1; gj <= NGENESALLOWED; gj++){
	   

      AUC_PR_HEALTHYval= 0;
	for(k=1; k<= NHEALTHY; k++)
	for(ik=1; ik<= NHEALTHY; ik++){
	if(k == ik) continue; 
	if(GREATER(DATA_HEALTHY[gi][k], DATA_HEALTHY[gj][ik])) AUC_PR_HEALTHYval = AUC_PR_HEALTHYval + 1;
	if(EQUAL(DATA_HEALTHY[gi][k], DATA_HEALTHY[gj][ik])) AUC_PR_HEALTHYval = AUC_PR_HEALTHYval + 0.5;
	}
	AUC_PR_HEALTHYval= AUC_PR_HEALTHYval/(NHEALTHY*NHEALTHY);

      AUC_PR_SICKval= 0;
	for(k=1; k<= NSICK; k++)
	for(ik=1; ik<= NSICK; ik++){
	if(k == ik) continue; 
	if(GREATER(DATA_SICK[gi][k], DATA_SICK[gj][ik])) AUC_PR_SICKval = AUC_PR_SICKval + 1;
	if(EQUAL(DATA_SICK[gi][k], DATA_SICK[gj][ik])) AUC_PR_SICKval = AUC_PR_SICKval + 0.5;
	}
	AUC_PR_SICKval= AUC_PR_SICKval/(NSICK*NSICK);


AUC_SCOREval = AUC_PR_HEALTHYval - AUC_PR_SICKval; if (LESS(AUC_SCOREval , 0)) AUC_SCOREval = - AUC_SCOREval;



    if(GREATER(AUC_SCOREval , TOPSCORE_AUC)){
        TOPSCORE_AUC = AUC_SCOREval; 
        TOPPAIR1_AUC = gi;
        TOPPAIR2_AUC = gj;

   }
}
}
   printf("Top AUCTSP is %d and %d with score %f\n", TOPPAIR1_AUC, TOPPAIR2_AUC, TOPSCORE_AUC); 

*/

} //end tries

getrusage(RUSAGE_SELF, &timing1);//snapsot of the timer after the TRIES, so as to better compute the average duration
printf(" EXEC TIME  %ld usec\n\n", ((timing1.ru_utime.tv_sec - timing0.ru_utime.tv_sec)*1000000L +timing1.ru_utime.tv_usec - timing0.ru_utime.tv_usec)/TRIES);

}
