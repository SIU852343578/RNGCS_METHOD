#include <stdio.h>
#include <stdlib.h>
#include<string.h>

#define BETA 500 //Number of bins

#define FACTOR_INT_PART 1000 // Multiplicative factor for turning the significant decimal digits of AUC value to an integer, to avoid rounding errors.
                            //for example,  0.005/0.005 when truncated to an integer returns 0 whereas it should be 1
                            // Bin width is W = 0.5/BETA = 0.5/500 =  0.001 so FACTOR_INT_PART = 1000

#define MAXGENEIDSIZE 100

typedef struct binstruc { char geneid[MAXGENEIDSIZE]; float geneauc; struct binstruc *next;} binelement;

struct binarraystruc { binelement *list; int count;} BINARRAY[BETA];

binelement *p;

//the program expects the original number of genes and the target number of genes as  program arguments 
//and the file of the AUC values of the genes by redirection to the standard input

main(argc, argv)
int argc;
char *argv[];
{

float geneauc;
float theta;
int i, bintosplit, binindex, genecount, numberofgenes, allowedgenecount;
char geneid[MAXGENEIDSIZE]; 

sscanf(argv[1],"%d", &numberofgenes);
sscanf(argv[2],"%d", &allowedgenecount);

//initialize each entry of BINARRAY to empty
for(i=0; i< BETA; i++) {BINARRAY[i].list = NULL; BINARRAY[i].count = 0;}

//read in the AUC values of the genes (assumed to be read in by redirection to the standard input)
for(i=1; i<= numberofgenes; i++){

scanf("%s %f", geneid, &geneauc);

//find the bin index to which this geneauc value belongs

binindex  = geneauc  * FACTOR_INT_PART - FACTOR_INT_PART/2 ; 
//note: the straightforward expression binindex = (geneauc -0.5)*FACTOR_INT_PART may give rounding errors!
//
//
printf("binindex of gene %s with AUC= %f is %d\n", geneid, geneauc, binindex);

p = malloc(sizeof(binelement)) ;
if( (p  == NULL)){ printf("MALLOC memory problem\n"); exit(1);}

//insert the new member to the top (convenient for speed: order does not matter!) of the list pointed to by BINARRAY[binindex] 
strcpy(p->geneid, geneid);
p->geneauc = geneauc;
p->next = BINARRAY[binindex].list;
BINARRAY[binindex].list = p;
BINARRAY[binindex].count++;
}


//find the binindex to stop, based on how many genes can be afforded, by starting from the larger bin indices

genecount = 0;

for(i=BETA-1; i>= 0; i--){
 genecount += BINARRAY[i].count;
 if(genecount > allowedgenecount) break;

} 
 
printf("\n\n");

 if(genecount <= allowedgenecount){ // this means that all genes in the dataset can be afforded!
    
	 printf("ALL %d GENES IN THE DATASET CAN BE AFFORDED!\n", allowedgenecount);
	 
}else{
  bintosplit = i;
//the excess amount is genecount - BINARRAY[bintosplit].count
//therefore, go down the list of BINARRAY[bintosplit].list and pick only the allowed number of genes

//select (i.e., print out) all genes before the bintosplit
//mark also the minumun AUC value (theta) of all genes selected
//Note: all genes selected have AUC values >= theta 
//but there may be genes beloging to bin with index bintosplit that have AUC value >= theta 
//but they are not selected, since the genes in each bin are unsorted.
theta = 1;
for(i=BETA-1; i>bintosplit; i--) {//print the members of list BINARRAY[i].list
   p = BINARRAY[i].list; 
   while(p != NULL){ 
      printf("%s %f\n", p->geneid, p->geneauc); 
      if(p->geneauc < theta) theta = p->geneauc; 
      p = p->next;}
}

//select (print out) now only the allowed number of genes from the list of BINARRAY[bintosplit].list

genecount -= BINARRAY[bintosplit].count; //the excess amount is genecount - BINARRAY[bintosplit].count

//go down the list of BINARRAY[bintosplit].list and pick only the allowed number of genes


p = BINARRAY[bintosplit].list; 
while(p != NULL){ 
   if(genecount == allowedgenecount) break;
   genecount++;
   printf("%s %f\n", p->geneid, p->geneauc); 
   if(p->geneauc < theta) theta = p->geneauc;
   p = p->next;
}


printf("\n\n TOTAL NUMBER OF GENES SELECTED : %d\n", genecount);
printf(" MINIMUM AUC OF SELECTED GENES (THETA)  : %f\n", theta);

}

}
