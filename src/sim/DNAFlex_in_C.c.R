/* --------------- CHANGE IN THIS PROGRAM TO ORIGINAL ---------------- */

/* 

1) The arrays seq1,stif and geom are adjusted so that to each single bp-step a stiffness matrix can be assigned and not just an average stiffness matrix as before 

2) Another output file "table_all_(helpar).dat" is provided which gives all (helpar) parameters for each MC optimization (each line of the table equals the (helpar) parameters of the examined sequence)

*/




/* Headers from DNA Flex */
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Headers from SCHNArP */

#include "schna_ar_p.h"

#include <time.h>

/* End Headers */


/* -------------------Functions DNA Flex------------------------ */

void order(char **names,int *ord)
{
int i;
for(i=0;i<6;i++){
	 if (names[i] =="shif"){
		ord[i]=0;}
	  else if (names[i]=="slid") {
		ord[i]=1;}
	  else if (names[i]=="rise") {
		ord[i]=2;}
	  else if (names[i]=="tilt") {
		ord[i]=3;}
	  else if (names[i]=="roll") {
		ord[i]=4;}
	  else if (names[i]=="twis") {
		ord[i]=5;}
	}

	
}


/*----------------Read sequence parameters from seq.dat in 'ist' --------------------------*/

void readsq_di(int *pointer_itot,char ist[1000][10])
{
FILE *myfile;
  char line[800];
  int i;
  int len;
  char buffer[3];


  myfile=fopen("my_seq_56mer.dat", "r");


while(fgets(line, 800, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
 {
len = strlen(line);
*pointer_itot = len-2;

for(i=0;i<*pointer_itot;i++)
   {
        sprintf(ist[i], "%c%c", line[i],line[i+1]);  
  
   }
 }

 fclose(myfile);


}

void readsq_tetra(int *pointer_itot,char ist[1000][10])
{
FILE *myfile;
  char line[800];
  int i;
  int len;



  myfile=fopen("my_seq_56mer.dat", "r");


while(fgets(line, 800, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
 {
len = strlen(line);
*pointer_itot = len-4;

for(i=0;i<*pointer_itot;i++)
   {
        sprintf(ist[i], "%c%c%c%c", line[i],line[i+1],line[i+2],line[i+3]);  
  
   }
 }

 fclose(myfile);


}


/*---------------Read stiffness parameters from tetranucleotide step from stif.dat---------*/





void readss_ind(double stif[1000][6][6],char seq1[1000][10],double geom[1000][6],int ord[6])
{
FILE *myfile;
FILE *myseq_to_file;
   int i,j,k;
   double buffer;
   char line[800];
   int len;
   

   myseq_to_file = fopen ("my_seq_56mer.dat", "rt");  	
   myfile = fopen ("fcte_56merSPCE.MCbsc1seqdep.dat", "rt");  

while(fgets(line, 800, myseq_to_file) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
 {
len = strlen(line);
j=0;
//len = (len-1)/2;
len = len - 2;
//printf("len %d \n",len);


for(i=0;i<len;i++)
   {
	
        sprintf(seq1[i], "%c%c", line[j],line[j+1]);  
	j=j+1;
  
   }
 }




 for(i=0;i<len;i++){

 for(j=0;j<6;j++){
	for(k=0;k<6;k++){
		fscanf(myfile, "%lf", &buffer);		// Get stif parameters 6x6 (line 2-7)
		stif[i][ord[j]][ord[k]] = buffer;
		}
		 }

for(k=0;k<6;k++){
	fscanf(myfile, "%lf", &buffer);			// Get geom parameters 1x6 (line 8)
	geom[i][ord[k]] = buffer;
		}

	}

   fclose(myseq_to_file);
   fclose(myfile);  /* close the file prior to exiting the routine */

}


void readss_di(double stif[1000][6][6],char seq1[1000][10],double geom[1000][6],int ord[6])
{
FILE *myfile;
   int i,j,k;
   double buffer;
   char line[800];
   char tetra[10];
   int len;
   

 	
   myfile = fopen ("stif_bsc0_dimer.dat", "rt");  



 for(i=0;i<16;i++){	

fscanf(myfile, "%s", tetra);	
strcpy(seq1[i],tetra); 

 for(j=0;j<6;j++){
	for(k=0;k<6;k++){
		fscanf(myfile, "%lf", &buffer);		// Get stif parameters 6x6 (line 2-7)
		stif[i][j][k] = buffer;
		}
		 }

for(k=0;k<6;k++){
	fscanf(myfile, "%lf", &buffer);			// Get geom parameters 1x6 (line 8)
	geom[i][k] = buffer;
		}

	}


   fclose(myfile);  /* close the file prior to exiting the routine */

}




void readss_tetra(double stif[1000][6][6],char seq1[1000][10],double geom[1000][6],int ord[6])
{
FILE *myfile;
   int i,j,k;
   double buffer;
   char line[800];
   char tetra[10];
   int len;
   

   myfile = fopen ("output_reordered.tab", "rt");  



 for(i=0;i<256;i++){	

fscanf(myfile, "%s", tetra);	
strcpy(seq1[i],tetra); 

 for(j=0;j<6;j++){
	for(k=0;k<6;k++){
		fscanf(myfile, "%lf", &buffer);		// Get stif parameters 6x6 (line 2-7)
		stif[i][j][k] = buffer;
		}
		 }

for(k=0;k<6;k++){
	fscanf(myfile, "%lf", &buffer);			// Get geom parameters 1x6 (line 8)
	geom[i][k] = buffer;
		}

	}

   fclose(myfile);  /* close the file prior to exiting the routine */

}


void assign_ind(int itot,char ist[1000][10],char seq1[1000][10],double stif[1000][6][6],double stiff[1000][6][6],double geom[1000][6],double xconf0[6][1000])
{
int i,j,k,k1,k2;
int f=0;
for (i=0; i<itot; i++) {
				

				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff[i][k1][k2] = stif[i][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0[k][i] = geom[i][k]; //xconf0[helical values][# of tetranucleotide step]
							}
						
			 }


}


void assign_di(int itot,char ist[1000][10],char seq1[1000][10],double stif[1000][6][6],double stiff[1000][6][6],double geom[1000][6],double xconf0[6][1000])
{
int i,j,k,k1,k2;
int f=0;
for (i=0; i<itot; i++) {
	for (j=0; j<16; j++) {
				if (ist[i][0] == seq1[j][0] && ist[i][1] == seq1[j][1]) 
				{ 
				f++;
				
				
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff[i][k1][k2] = stif[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0[k][i] = geom[j][k]; 
							}
				   }

			      }							
			 }

}


void assign_tetra(int itot,char ist[1000][10],char seq1[1000][10],double stif[1000][6][6],double stiff[1000][6][6],double geom[1000][6],double xconf0[6][1000])
{
int i,j,k,k1,k2;
int f=0;
for (i=0; i<itot; i++) {
	for (j=0; j<256; j++) {
				if (strcmp(ist[i],seq1[j]) == 0) 
				{ 
				f++;
				
				
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff[i][k1][k2] = stif[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0[k][i] = geom[j][k]; 
							}
				   }

			      }							
			 }

}




void readini(int itot,double xconfi[6][1000])
{
FILE *myfile;
   int i,j,k;
   double buffer;

myfile = fopen ("conf.dat", "rt");  

for (i=0; i<itot; i++){
	for(j=0;j<6;j++){
	fscanf(myfile, "%lf", &buffer);			
	xconfi[j][i] = buffer;		
			}
		}

  fclose(myfile);
}


void diff(int itot,double xconf[6][1000],double xconf0[6][1000],double dx[6][1000])
{
int i,j;

for (i=0; i<itot; i++) {
	for (j=0; j<6; j++) {
				dx[j][i] = xconf[j][i]-xconf0[j][i];
			    }
			}

}


void energy(int itot,double xconf[6][1000],double xconf0[6][1000],double dx[6][1000],double stiff[1000][6][6],double ener[1000])
{
int i,j,k;
double xt=0.0;

diff(itot,xconf,xconf0,dx);

for (i=0; i<itot; i++) {
	for (j=0; j<6; j++) {
				xt=xt+0.5*(stiff[i][j][j]*dx[j][i]*dx[j][i]);
				for (k=j+1; k<6; k++) {xt=xt+(stiff[i][j][k]*dx[j][i]*dx[k][i]); //i atom number, j and k summation over different helical parameters
						     }
				
    			     }
				
			 ener[i] = xt;

			}

}




void putini(int itot,double xconf[6][1000],double enertot,char ist[1000][10])
{
char *seq[60];
int i,j;
j=0;
printf("\n\n************************************************\n \n");
printf("MC-DNA RUN \n");
printf("************************************************\n \n");
printf("************************************************\n");
printf("Sequence \n");
 //WRITE SEQUENCE
for(i=0; i<itot; i++){printf("%s",ist[j]);
			j=j+2;}

printf("\n\n************************************************ \n");
printf("Starting coordinates \n");
for(i=0; i<itot; i++){
			for(j=0; j<6; j++){
printf("%f   ",xconf[j][i]);
					    }
printf("\n ");
		      }
printf("\n\n************************************************ \n\n");
printf("Starting energy kcal/mol %f \n",enertot);
printf("\n************************************************ \n\n");

}


void energ(int it,double xconf[6][1000],double xconf0[6][1000],double stiff[1000][6][6],double *enex, double dd[6])
{
int i,j,k;
double xt=0.0;

for(j=0;j<6;j++){dd[j] = xconf[j][it]-xconf0[j][it];}

for (j=0; j<6; j++) {
				xt=xt+0.5*(stiff[it][j][j]*dd[j]*dd[j]);
				for (k=j+1; k<6; k++) {
						xt=xt+(stiff[it][j][k]*dd[j]*dd[k]); //i atom number, j and k summation over different helical parameters
						     }
				
    			     }

	*enex = xt;
}

void metropolis(double boltz, double temp,double enertot, double enerp, int *iflag)
{
double r,delta,beta,emet;


delta=enerp-enertot;

//printf("delta %lf \n",delta);

if(delta <= 0.0) {*iflag=1;}	// Configuration accepted

//printf("iflag %d \n",*iflag);
beta = 1.0/(boltz*temp);
r = ((double) rand() / (RAND_MAX)); 

//printf("r %lf \n",r);
//printf("beta %f \n",beta);

emet = exp(-beta*delta);
if(emet<r) {*iflag=2;}
	   else{*iflag=1;}

}



// void putinter(itot,xconf,ist,icon,idiv)

void output_str_di(char output[1000][100],int itot, char ist[1000][10], double xconf[6][1000], int phi)
{
FILE *myfile;
//char str_out[1000];
char seq1;
int i,j;
j=2;
int len_seq;
len_seq = itot + 1;

//char c[10];
//sprintf(c, "%d", phi); //Make phi as a string
char filename[64];

sprintf(filename, "output_dnaflex/output_coordinates_%.6d.dat", phi);

myfile=fopen(filename, "w");

//myfile=fopen("output_coordinates.dat", "w");

sprintf(output[0], "%d base-pairs \n", len_seq);
sprintf(output[1], "Shift  Slide    Rise    Tilt    Roll   Twist \n");
for(i=0;i<itot;i++){
		
        if (ist[i][0]=='A'){ 
		seq1='T';}
	else if (ist[i][0]=='T') {
		seq1='A';}
	else if (ist[i][0]=='C') {
		seq1='G';}
	else if (ist[i][0]=='G') {
		seq1='C';}
       
sprintf(output[j], "%c-%c	%lf	%lf	%lf	%lf	%lf	%lf \n", ist[i][0],seq1,xconf[0][i-1],xconf[1][i-1],xconf[2][i-1],xconf[3][i-1],xconf[4][i-1],xconf[5][i-1]);

			j=j+1;

	}

        if (ist[itot-1][1]=='A'){ 		// For the last basepair in sequence
		seq1='T';}
	else if (ist[itot-1][1]=='T') {
		seq1='A';}
	else if (ist[itot-1][1]=='C') {
		seq1='G';}
	else if (ist[itot-1][1]=='G') {
		seq1='C';}

sprintf(output[j+1], "%c-%c	%lf	%lf	%lf	%lf	%lf	%lf \n", ist[itot-1][1],seq1,xconf[0][itot-1],xconf[1][itot-1],xconf[2][itot-1],xconf[3][itot-1],xconf[4][itot-1],xconf[5][itot-1]);

for(i=0;i<itot+4;i++){		//Put strings in file
fputs(output[i],myfile);}

	

fclose(myfile);

}



void output_str_tetra(char output[1000][100],int itot, char ist[1000][10], double xconf[6][1000], int phi)
{

// Tetramer output is 'length of sequence - 2' because of comparing tetramers for stiffness matrices

FILE *myfile;
//char str_out[1000];
char seq1;
int i,j;
j=2;
int len_seq;
len_seq = itot + 1;

//char c[10];
//sprintf(c, "%d", phi); //Make phi as a string
char filename[64];

sprintf(filename, "output_dnaflex/output_coordinates_%.6d.dat", phi);

myfile=fopen(filename, "w");

//myfile=fopen("output_coordinates.dat", "w");

sprintf(output[0], "%d base-pairs \n", len_seq);
sprintf(output[1], "Shift  Slide    Rise    Tilt    Roll   Twist \n");
for(i=0;i<itot;i++){
		
        if (ist[i][1]=='A'){ 
		seq1='T';}
	else if (ist[i][1]=='T') {
		seq1='A';}
	else if (ist[i][1]=='C') {
		seq1='G';}
	else if (ist[i][1]=='G') {
		seq1='C';}
       
sprintf(output[j], "%c-%c	%lf	%lf	%lf	%lf	%lf	%lf \n", ist[i][1],seq1,xconf[0][i-1],xconf[1][i-1],xconf[2][i-1],xconf[3][i-1],xconf[4][i-1],xconf[5][i-1]);

			j=j+1;

	}


        if (ist[itot-1][2]=='A'){ 		// For the last basepair in sequence
		seq1='T';}
	else if (ist[itot-1][2]=='T') {
		seq1='A';}
	else if (ist[itot-1][2]=='C') {
		seq1='G';}
	else if (ist[itot-1][2]=='G') {
		seq1='C';}

sprintf(output[j+1], "%c-%c	%lf	%lf	%lf	%lf	%lf	%lf \n", ist[itot-1][2],seq1,xconf[0][itot-1],xconf[1][itot-1],xconf[2][itot-1],xconf[3][itot-1],xconf[4][itot-1],xconf[5][itot-1]);

for(i=0;i<itot+4;i++){		//Put strings in file
fputs(output[i],myfile);}

	

fclose(myfile);

}



void write_hel_coord_in_table(char table_heli[1000][100],int itot, int ord[6], double xconf[6][1000], int phi, int ind)
{

//,char table_heli_rise[1000][100],char table_heli_tilt[1000][100],char table_heli_roll[1000][100],char table_heli_twis[1000][100],

int i;


for(i=0;i<itot-1;i++){
sprintf(table_heli[i], "%lf ", xconf[ind][i]);
		    }
sprintf(table_heli[itot-1], "%lf", xconf[ind][itot-1]);
sprintf(table_heli[itot], "\n\n");

}


/*--------------------End of Functions------------------------ */





/*--------------------------Beginning of DNAFlex-------------------------*/


int main(void)
{

/* Initialize all variables at beginning of code!!*/ 

int nmc,i0,ibuff2,juega,last;
int iflag; 	//it should be a bool, but I dont know how to assign true and false to a bool in c
char seq1[1000][10], ist[1000][10];
double stif[1000][6][6],stiff[1000][6][6],geom[1000][6];
int ord[6];
double xconfi[6][1000],xconf0[6][1000],xconf[6][1000];
double xconft[6][1000],ener[1000],sc[6],dd[6];
double dx[6][1000];
char *names[6];   				// No maximum size of char in C
char line;
double boltz=0.00198717;
double enex;
int itot;
int i,j,k;	// Iteration parameters
int test;	// Test parameter for # MC moves
char output[1000][100];			// Writing helical coordinates in this string
int test1,phi;
double avg_energy;
double *test_array;
float sum_energy;
int selection;

/* Initialize arrays and files for huge tables */

char table_heli_shif[1000][100];
char table_heli_slid[1000][100];
char table_heli_rise[1000][100];
char table_heli_tilt[1000][100];
char table_heli_roll[1000][100];
char table_heli_twis[1000][100];

FILE *file_shif;
file_shif=fopen("output_tables_helpar/table_all_shif.dat", "w");

FILE *file_slid;
file_slid=fopen("output_tables_helpar/table_all_slid.dat", "w");

FILE *file_rise;
file_rise=fopen("output_tables_helpar/table_all_rise.dat", "w");

FILE *file_tilt;
file_tilt=fopen("output_tables_helpar/table_all_tilt.dat", "w");

FILE *file_roll;
file_roll=fopen("output_tables_helpar/table_all_roll.dat", "w");

FILE *file_twis;
file_twis=fopen("output_tables_helpar/table_all_twis.dat", "w");


/* ------- End Initialze for tables and files --------- */







srand (time(NULL));		// initialize random seed (the random seed is initialized to a value representing the current time (calling time) to generate a different value every time the program is run.)


clock_t start_all = clock(), diff_all; //Timer

/* icuan = 1;
icnt = 1; */

//printf("This is the variable %i \n", nmc);

/* This is the part which should be read from a txt file (see fortran code) PROBLEM!!! */

nmc=100000;
double agr=1.0;
int ibuff=10000;
ibuff2=100;
double temp=298.0,fact=0.3;
last=1000000;
sc[0]=5.148;
sc[1]=3.558;
sc[2]=2.052;
sc[3]=29.472;
sc[4]=34.266;
sc[5]=53.862;
names[0]="shif";
names[1]="slid";
names[2]="rise";
names[3]="tilt";
names[4]="roll";
names[5]="twis";

/* ----------------------------------End of part--------------------------------------- */





/* ------------------------------- Looping and selection variables ------------------------------- */


selection = 2;			// 0 = individual, 1 = dimer, 2 = tetramer


test1 = 200;			// number of executing all the MC and reconstruction part
test = 100000;			// number of MC steps


/* ------------------------------- END Looping and selection variables ------------------------------- */






order(names,ord);



if(selection == 0){		// Select individual bsc1 stiffness matrices
readsq_di(&itot,ist);

readss_ind(stif,seq1,geom,ord);

assign_ind(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 1){		// Select average bsc0 dimer stiffness matrices
readsq_di(&itot,ist);

readss_di(stif,seq1,geom,ord);

assign_di(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 2){		// Select average bsc1 tetramer stiffness matrices
readsq_tetra(&itot,ist);

readss_tetra(stif,seq1,geom,ord);

assign_tetra(itot,ist,seq1,stif,stiff,geom,xconf0);
}





readini(itot,xconfi);










/* ------------------------------- Looping starts ------------------------------- */



test_array = malloc(test1*sizeof(double));
for(phi=0;phi<test1;phi++){			// do optimization test1 times to check if MC converges


clock_t start = clock(), diff; //Timer


// Compute energy
energy(itot,xconfi,xconf0,dx,stiff,ener);



double ener0=0.0;
double enertot;

for(i=0; i<itot; i++){ener0 = ener0 + ener[i];}

enertot = ener0;



// Set working coordinates and energies equal to the starting ones



for(i=0; i<6; i++){
for(j=0; j<itot; j++){
			xconf[i][j] = xconfi[i][j];  
		      }	
		  }




// Write initial conformation in output

//putini(itot,xconf,enertot,ist);

















/*---------------- Start of MC -----------------*/



//Loop to iterate MC algorithm 'test' times
for(i0=0;i0<test;i0++){


double r;
r = ((double) rand() / (RAND_MAX)); 	//Generates a random number between 0 and 1






// Select a random step (itt) to change

double xtot = (double) itot+1e-9;
int itt = 0 + (int) (xtot*r);		// the step to change
//printf("itt %d \n", itt);

double enerbase = enertot - ener[itt];	// reference values for calculating energies


for(i=0;i<6;i++){xconft[i][itt] = xconf[i][itt];} 

double r1 = ((double) rand() / (RAND_MAX)); 
int it1 = 1+ (int) 3.0*r1;			// it1 is the number of variables to change
int jf1=0;
double r2,r3;
for(i=0;i<it1;i++){

/* --- Change coordinates randomly --- */

	r2 = ((double) rand() / (RAND_MAX)); 	


		jf1 = 0 + (int) 6.0*r2;			// jf1 is the coordinate to change
			
		sc[jf1] = 12.0*sqrt(0.6/(6.0 * stiff[itt][jf1][jf1]));	//new entry of scaling matrix

		r3 = ((double) rand() / (RAND_MAX)); 
		r3 = (r3*2.0)-1.0;
		r3 = r3*agr*sc[jf1]*fact;

		xconft[jf1][itt] = xconft[jf1][itt] + r3;

		}


enex=0.0;
energ(itt,xconft,xconf0,stiff,&enex,dd);



double enerp = enerbase + enex;





// Call Metropolis algorithm

metropolis(boltz,temp,enertot,enerp,&iflag);



if(iflag==1) {	for(i=0;i<6;i++){xconf[i][itt] = xconft[i][itt];}	// What to do when move accepted
		enertot = enerp;					// New total energy is energy from accepted move
		ener[itt] = enex;
	      }





} 

/*---------------------------------End of MC algorithm-----------------------------*/




test_array[phi] = enertot;		// Save Energy after MC in array to compute average later


//printf("Energy after optimizing: %f \n \n", test_array[phi]);	// printing output energy








	if(selection == 2){
	output_str_tetra(output,itot,ist,xconf,phi);
			  }
	else{
	output_str_di(output,itot,ist,xconf,phi);	// Output hel coord in file, for input to SCHNArP
	}


/* ------------ End of DNA Flex Program -------------- */










diff = clock() - start;

int msec = diff * 1000 / CLOCKS_PER_SEC;
printf("Time taken for MC: %d seconds %d milliseconds\n", msec/1000, msec%1000);











clock_t start1 = clock(), diff1; // Start next timer





/* -------------------- Beginning of SCHNArP main schnarp.c -------------- */

  long ich;
    char str[BUF512];

    system("clear");


/* ----- Commented to probe for efficiency ----------- */

ich = 1;

/*
    printf("Choose one the following rebuilding methods:\n\n");
    printf("1. Use LOCAL CEHS base-pair and base-step parameters\n");
    printf("2. Use GLOBAL helical parameters\n");
    printf("3. Quit\n");
    printf("Your choice(1-3, Dft 1): ");
    read_stdin_str(str);
    ich = atoi(str);

*/

    if ((ich < 1) || (ich > 3))
        ich = 1;

    if (ich == 1)
        CEHS_build(phi);  /* Local CEHS */
    else if (ich == 2)
        GLH_build();  /* GLOBAL */
    else
        exit(0);  /* return to OS */


/*---------------End of SCHNArP main schnarp.c --------------------*/






diff1 = clock() - start1;

int msec1 = diff1 * 1000 / CLOCKS_PER_SEC;
printf("Time taken for reconstruction: %d seconds %d milliseconds\n", msec1/1000, msec1%1000);

printf("Percentage of computation time for SCHNArP %f \n",(double) msec1/(msec+msec1));








/* ---- Write values for output tables ---- */

int ind;

ind=ord[0];						// assigns which hel par is printed in table
write_hel_coord_in_table(table_heli_shif, itot, ord, xconf, phi, ind);

ind=ord[1];
write_hel_coord_in_table(table_heli_slid, itot, ord, xconf, phi, ind);

ind=ord[2];
write_hel_coord_in_table(table_heli_rise, itot, ord, xconf, phi, ind);

ind=ord[3];
write_hel_coord_in_table(table_heli_tilt, itot, ord, xconf, phi, ind);

ind=ord[4];
write_hel_coord_in_table(table_heli_roll, itot, ord, xconf, phi, ind);

ind=ord[5];
write_hel_coord_in_table(table_heli_twis, itot, ord, xconf, phi, ind);




for(i=0;i<itot+1;i++){		//Put strings in file = one line in the file
fputs(table_heli_shif[i],file_shif);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_slid[i],file_slid);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_rise[i],file_rise);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_tilt[i],file_tilt);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_roll[i],file_roll);}

for(i=0;i<itot+1;i++){		
fputs(table_heli_twis[i],file_twis);}







} 
/* ------------------- End of for loop to execute the whole program test1 times ------------------- */


fclose(file_shif);
fclose(file_slid);
fclose(file_rise);
fclose(file_tilt);
fclose(file_roll);
fclose(file_twis);



sum_energy=0.0;
for(i=0;i<test1;i++){
sum_energy = sum_energy + test_array[i];
//printf("%f \n", test_array[i]);
}


free (test_array);				//free memory which was allocated before
//free (table_heli_shif);			//free memory which was allocated before


avg_energy = sum_energy/(double) test1;


printf("Mean energy after %d steps and optimizing %d times: %f \n \n", test, test1, avg_energy);



diff_all = clock() - start_all;

int msec_all = diff_all * 1000 / CLOCKS_PER_SEC;
printf("Time taken for MC: %d seconds %d milliseconds\n", msec_all/1000, msec_all%1000);

}
