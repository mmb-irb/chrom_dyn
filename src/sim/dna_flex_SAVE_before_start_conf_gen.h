/* -------------------- (De)Allocate memory ----------------------- */

int **int_aloc2d(int x, int y)
{
  int i;
  int **array = malloc(sizeof(int *) * x + sizeof(int) * x * y);
  for (i = 0; i < x; i++)
  {
    array[i] = ((int *)(array + x)) + y * i;
  }
  return array;
}

double **dbl_aloc2d(int x, int y)
{
  int i;
  double **array = malloc(sizeof(double *) * x + sizeof(double) * x * y);
  for (i = 0; i < x; i++)
  {
    array[i] = ((double *)(array + x)) + y * i;
  }
  return array;
}


double ***dbl_aloc3d(int x,int y,int z)
{
    int i,j;
    double *allElements = malloc(x * y * z * sizeof(double));
    double ***array3D = malloc(x * sizeof(double **));

    for(i = 0; i < x; i++)
    {
        array3D[i] = malloc(y * sizeof(double *));

        for(j = 0; j < y; j++)
        {
            array3D[i][j] = allElements + (i * y * z) + (j * z);
        }
    }
	return array3D;
}


char **char_aloc2d(int x, int y)
{
  int i;
  char **array = malloc(sizeof(char *) * x + sizeof(char) * x * y);
  for (i = 0; i < x; i++)
  {
    array[i] = ((char *)(array + x)) + y * i;
  }
  return array;
}


/*
double **dbl_aloc2d(int rows,int cols)
{
    double **arr2D;
    int i;

    arr2D = (double**)malloc(rows*sizeof(double*));
    for(i=0;i<rows;i++)
    {
    	arr2D[i] = (double*)malloc(cols*sizeof(double));
    }
   return arr2D;
}


double ***dbl_aloc3d(int l,int m,int n)
{
double ***arr3D;
int i,j,k;

arr3D = (double***)malloc(l * sizeof(double **));

for(i=0;i<l;i++)
{
	arr3D[i] = (double**)malloc(m * sizeof(double*));
	for(j=0;j<m;j++)
	{
		arr3D[i][j] = (double*)malloc(n*sizeof(double));
	}
}

return arr3D;
}


char **char_aloc2d(int rows,int cols)
{
    char **arr2D;
    int i;

    arr2D = (char**)malloc(rows*sizeof(char*));
    for(i=0;i<rows;i++)
    {
    	arr2D[i] = (char*)malloc(cols*sizeof(char));
    }
   return arr2D;
}
*/

/*
void free_dbl2d(double **array_2d, int x)
{
  int i;
  for (i = 0; i < x; i++)
  {
    free(array_2d[i]);
  }
  free(array_2d);
}
*/


void free_dbl3d(double ***array3D,int x)
{
  int i;
    for(i = 0; i < x; i++)
    {
        free(array3D[i]);
    }
    free (array3D);
}

/*
void free_dbl3d(double ***array_3d, int x, int y)
{
   int i,j,k;
    for(i = 0; i < x; i++)
    {
        for(j = 0; j < y; j++)
        {
            free(array_3d[i][j]);
        }
	free(array_3d[i]);
    }
  free(array_3d);
}


void free_char2d(char **array_2d, int x)
{
  int i;
  for (i = 0; i < x; i++)
  {
    free(array_2d[i]);
  }
  free(array_2d);
}
*/

/* ------------------ My added functions ----------------------- */

void get_seq_len(int *pointer_len_seq, int *pointer_itot, char seq_file[])
{
FILE *myfile;
  char line[10000000];
  int i;
  int len;
  char buffer[3];


//  myfile=fopen("nucl3_linker4_56mer_seq.dat", "r");
  // myfile=fopen("seq_ddd_160bp.dat", "r");
   myfile=fopen(seq_file, "r");


while(fgets(line, 10000000, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
 {
len = strlen(line);
*pointer_len_seq = len+5;		// +5 just to overcome possible edge effects in arrays
*pointer_itot = len-2;
 }
 fclose(myfile);

}


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


void readsq_ind(int *pointer_itot,char **ist,char seq_file[])
{
FILE *myfile;
  char line[10000000];
  int i;
  int len;
  char buffer[3];


//myfile=fopen("nucl3_linker4_56mer_seq.dat", "r");
//myfile=fopen("seq_56merTIP3P.dat", "r");
//myfile=fopen("seq_ddd_160bp.dat", "r");
  myfile=fopen(seq_file, "r");


while(fgets(line, 10000000, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
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

void readsq_di(int *pointer_itot,char **ist,char seq_file[])
{
FILE *myfile;
  char line[10000000];
  int i;
  int len;
  char buffer[3];


 //myfile=fopen("my_seq_56mer.dat", "r");
   myfile=fopen(seq_file, "r");


while(fgets(line, 10000000, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
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

void readsq_tetra(int *pointer_itot,char **ist,char seq_file[])
{
FILE *myfile;
  char line[10000000];
  int i;
  int len;



//myfile=fopen("my_seq_56mer.dat", "r");
  myfile=fopen(seq_file, "r");


while(fgets(line, 10000000, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
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


void readsq_tetra_di(int *pointer_itot,char **ist,char seq_file[])
{
FILE *myfile;
  char line[10000000];
  int i;
  int len;
  int len_tetra;


//myfile=fopen("my_seq_56mer.dat", "r");
  myfile=fopen(seq_file, "r");


while(fgets(line, 10000000, myfile) != NULL) // get a line, up to 800 chars from fr.  done if NULL
 {
len = strlen(line);
*pointer_itot = len-2;	// Let itot be the number of all tetramer steps with the fictive steps at both ends (this is why I put -2 and not -4)
len_tetra = len-4;

sprintf(ist[0], "C%c%c%c", line[0],line[1],line[2]);		// Assign a tetramer sequence for first bp ('C' is just for completeness)

for(i=1;i<len_tetra+1;i++)
   {
        sprintf(ist[i], "%c%c%c%c", line[i-1],line[i],line[i+1],line[i+2]);  
  
   }

sprintf(ist[len_tetra+1], "%c%c%cC", line[len_tetra],line[len_tetra+1],line[len_tetra+2]);	// Assign a tetramer sequence for last bp ('C' is just for completeness)	

 }

 fclose(myfile);
/*
for(i=0;i<*pointer_itot;i++){
printf("%s\n",ist[i]); }

printf("%d\n", i);
printf("%s\n",ist[*pointer_itot-1]);
*/

}


/*---------------Read stiffness parameters from tetranucleotide step from stif.dat---------*/





void readss_ind(double ***stif,char **seq1,double **geom,int ord[6],char seq_file[],char stif_file[])
{
FILE *myfile;
FILE *myseq_to_file;
   int i,j,k;
   double buffer;
   char line[10000000];
   int len;
   

   //myseq_to_file = fopen ("nucl3_linker4_56mer_seq.dat", "rt");  	
   //myfile = fopen ("nucl3_linker4_56mer_stiff.dat", "rt");  

   //myseq_to_file = fopen ("seq_56merTIP3P.dat", "rt");  	
   //myfile = fopen ("stiff_56merTIP3P.dat", "rt"); 

   //myseq_to_file = fopen ("seq_ddd_160bp.dat", "rt"); 
   //myfile = fopen ("DSpcCheaII_FC_160bp.dat", "rt");

   myseq_to_file = fopen(seq_file, "rt");
   myfile = fopen(stif_file, "rt");

while(fgets(line, 10000000, myseq_to_file) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
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


void readss_di(double ***stif,char **seq1,double **geom,int ord[6])
{
FILE *myfile;
   int i,j,k;
   double buffer;
   char tetra[10];
   int len;
   

 	
   myfile = fopen ("stif_bsc1_di_k_avg_miniabc.dat", "rt");  



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




void readss_tetra(double ***stif,char **seq1,double **geom,int ord[6])
{
FILE *myfile;
   int i,j,k;
   double buffer;
   char tetra[10];
   int len;
   

   myfile = fopen ("stif_bsc1_k_avg_miniabc.dat", "rt");  



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



void readss_tetra_di(double ***stif_tetra,char **seq1_tetra,double **geom_tetra, double ***stif_di,char **seq1_di,double **geom_di,int ord[6])
{
readss_tetra(stif_tetra,seq1_tetra,geom_tetra,ord);
readss_di(stif_di,seq1_di,geom_di,ord);

}




void assign_ind(int itot,char **ist,char **seq1,double ***stif,double ***stiff,double **geom,double **xconf0)
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


void assign_di(int itot,char **ist,char **seq1,double ***stif,double ***stiff,double **geom,double **xconf0)
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


void assign_tetra(int itot,char **ist,char **seq1,double ***stif,double ***stiff,double **geom,double **xconf0)
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


void assign_tetra_di(int itot,char **ist,double ***stif_tetra,char **seq1_tetra,double **geom_tetra,double ***stif_di,char **seq1_di,double **geom_di,double ***stiff,double **xconf0)
{
int i,j,k,k1,k2;
int f=0;

/*
printf("%d \n",itot);
printf("%c%c \n",ist[0][1],ist[0][2]);
printf("%c%c \n",seq1_di[0][0],seq1_di[0][1]);
*/

for (j=0; j<16; j++) {											// i=0
if (ist[0][1] == seq1_di[j][0] && ist[0][2] == seq1_di[j][1]) { 
								for (k1=0; k1<6; k1++) {
								  for (k2=0; k2<6; k2++) {
								     stiff[0][k1][k2] = stif_di[j][k1][k2];
												}
												}
								  for (k=0; k<6; k++) {
									xconf0[k][0] = geom_di[j][k]; 
											}
							}
			}

/*printf("%lf\n", xconf0[5][0]);
printf("%lf\n", stiff[0][4][4]);
printf("%lf\n", geom_di[5][3]);*/

for (i=1; i<itot-1; i++) {
	for (j=0; j<256; j++) {
				if (strcmp(ist[i],seq1_tetra[j]) == 0) 
				{ 
				f++;
								
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff[i][k1][k2] = stif_tetra[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0[k][i] = geom_tetra[j][k]; 
							}
				   }

			      }							
			 }


for (j=0; j<16; j++) {											// i=itot-1
if (ist[itot-1][1] == seq1_di[j][0] && ist[itot-1][2] == seq1_di[j][1]) { 
								for (k1=0; k1<6; k1++) {
								  for (k2=0; k2<6; k2++) {
								     stiff[itot-1][k1][k2] = stif_di[j][k1][k2];
												}
												}
								  for (k=0; k<6; k++) {
									xconf0[k][itot-1] = geom_di[j][k]; 
											}
							}
			}


/*printf("%lf\n", xconf0[5][itot-1]);
printf("%lf\n", stiff[itot-1][4][4]);
printf("%lf\n", geom_di[5][3]);
*/

}



void readini(int itot,double **xconfi)
{
FILE *myfile;
   int i,j,k;
   double buffer;

//myfile = fopen ("conf.dat", "rt");			// EDIT for every new sequence!!
//myfile = fopen ("conf_straight.dat", "rt");			// For straight starting sequence
//myfile = fopen ("conf_straight_mod_tw36.dat", "rt");  	// for nucleosome
myfile = fopen ("conf_fiber.dat", "rt");		// Reconstructing X-ray structure

for (i=0; i<itot; i++){
	for(j=0;j<6;j++){
	fscanf(myfile, "%lf", &buffer);			
	xconfi[j][i] = buffer;		
			}
		}

  fclose(myfile);
}


void diff(int itot,double **xconf,double **xconf0,double **dx)
{
int i,j;

for (i=0; i<itot; i++) {
	for (j=0; j<6; j++) {
				dx[j][i] = xconf[j][i]-xconf0[j][i];
			    }
			}

}


void energy(int itot,double **xconf,double **xconf0,double **dx,double ***stiff,double *ener)
{
int i,j,k;
double xt=0.0;

diff(itot,xconf,xconf0,dx);

for (i=0; i<itot; i++) {
	xt = 0.0;
	for (j=0; j<6; j++) {
				xt=xt+0.5*(stiff[i][j][j]*dx[j][i]*dx[j][i]);
				for (k=j+1; k<6; k++) {xt=xt+(stiff[i][j][k]*dx[j][i]*dx[k][i]); //i atom number, j and k summation over different helical parameters
						     }
				
    			     }
				
			 ener[i] = xt;

			}

}




void putini(int itot,double **xconf,double enertot,char **ist)
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


void energ(int it,double **xconf,double **xconf0,double ***stiff,double *enex, double dd[6])
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


void metropolis(double boltz, double temp, double enertot, double enerp, int *iflag)
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


void monte_carlo(int it_mc, int itot, double enertot, double *ener, double **xconft, double **xconf, double sc[6], double **xconf0, double ***stiff, double enex, double dd[6], double boltz, double temp, double *endenergy, int metr_on, int *no_mc_move)
{

int i0,i,buf_no_mc_move,iflag;
buf_no_mc_move = *no_mc_move;
iflag = 2;


//Loop to iterate MC algorithm 'it_mc' times
for(i0=0;i0<it_mc;i0++){


	double r;
	r = ((double) rand() / (RAND_MAX)); 	// Generates a random number between 0 and 1
	if(r>0.99999){r=0.99999;}
	// Select a random step (itt) to change

	double xtot = (double) itot+1e-9;
	int itt = 0 + (int) (xtot*r);		// the step to change
	//printf("itt %d \n", itt);
	
	double enerbase = enertot - ener[itt];	// reference values for calculating energies


	for(i=0;i<6;i++){xconft[i][itt] = xconf[i][itt];} 

	double r1 = ((double) rand() / (RAND_MAX)); 
	int it1 = 1+ (int) 3.0*r1;			// it1 is the number of variables to change
	int jf1=0;
	double agr=1.0;
	double fact=0.3;
	double r2,r3;
	for(i=0;i<it1;i++){

	/* --- Change coordinates randomly --- */

		r2 = ((double) rand() / (RAND_MAX)); 	

			if(r2>0.98){r2=0.98;}	// You have to put this so that (int) never rounds it to 6
			
			jf1 = 0 + (int) 6.0*r2;			// jf1 is the coordinate to change
				
			sc[jf1] = 12.0*sqrt(0.6/(6.0 * stiff[itt][jf1][jf1]));	//new entry of scaling matrix
	
			r3 = ((double) rand() / (RAND_MAX)); 
			r3 = (r3*2.0)-1.0;
			r3 = r3*agr*sc[jf1]*fact;
	
			xconft[jf1][itt] = xconft[jf1][itt] + r3;
	
			}

	
	enex=0.0;
	energ(itt,xconft,xconf0,stiff,&enex,dd);
	if(enex < 0){printf("enex %lf\n",enex);}


	double enerp = enerbase + enex;

	// Call Metropolis algorithm

if(metr_on==1){metropolis(boltz,temp,enertot,enerp,&iflag);}
	else{iflag = 1;}



	if(iflag==1) {	for(i=0;i<6;i++){xconf[i][itt] = xconft[i][itt];}	// What to do when move accepted
			enertot = enerp;					// New total energy is energy from 	accepted move
			ener[itt] = enex;
		      }
		//else {	*no_mc_move = buf_no_mc_move + 1; }

	} 
//printf("enerp %lf \n", enerp);
if(enex < 0){printf("enertot %lf \n", enertot);}
*endenergy = enertot*0.043424;	// kcal/mol recalculated in eV

}


void sim_ann(double tmax, double it, double *t, double factor, double boltz, double temp, double enertot, double enerp, int *iflag)
{
double r,delta,beta,emet,buf_t;

buf_t = *t;

delta=enerp-enertot;

if(it == 0){*t = tmax;}
	else{*t = buf_t*factor;}

//printf("t %lf \n",*t);
//printf("psi %lf \n",it);

if(delta <= 0.0) {*iflag=1;}	// Configuration accepted

//printf("iflag %d \n",*iflag);
beta = 1.0/((*t+boltz*temp));
//printf("beta %lf \n",beta);
r = ((double) rand() / (RAND_MAX)); 

//printf("r %lf \n",r);
//printf("beta %f \n",beta);

emet = exp(-beta*delta);
if(emet<r) {*iflag=2;}
	   else{*iflag=1;}

}



void monte_carlo_ind(int it_mc, int itot, double enertot, double *ener, double **xconft, double **xconf, double sc[6], double **xconf0, double ***stiff, double enex, double dd[6], double boltz, double temp, double *endenergy, int metr_on, int *get_itt)
{

int i0,i;
int itt,iflag;
double xtot, r;
double r1 = ((double) rand() / (RAND_MAX)); 
int it1;			// it1 is the number of variables to change
int jf1=0;
double agr=1.0;
double fact=0.3;
double r2,r3;
double enerbase;
double enerp;

iflag = 2;

do{
//Loop to iterate MC algorithm 'it_mc' times
for(i0=0;i0<it_mc;i0++){


	
	r = ((double) rand() / (RAND_MAX)); 	// Generates a random number between 0 and 1
	if(r>0.99999){r=0.99999;}
	// Select a random step (itt) to change

	xtot = (double) itot+1e-9;
	itt = 0 + (int) (xtot*r);		// the step to change
	//printf("itt %d \n", itt);
	
	enerbase = enertot - ener[itt];	// reference values for calculating energies


	for(i=0;i<6;i++){xconft[i][itt] = xconf[i][itt];} 

	r1 = ((double) rand() / (RAND_MAX)); 
	it1 = 1+ (int) 3.0*r1;
	jf1=0;
	
	for(i=0;i<it1;i++){

	/* --- Change coordinates randomly --- */

		r2 = ((double) rand() / (RAND_MAX)); 	

			if(r2>0.98){r2=0.98;}	// You have to put this so that (int) never rounds it to 6
			
			jf1 = 0 + (int) 6.0*r2;			// jf1 is the coordinate to change
				
			sc[jf1] = 12.0*sqrt(0.6/(6.0 * stiff[itt][jf1][jf1]));	//new entry of scaling matrix
	
			r3 = ((double) rand() / (RAND_MAX)); 
			r3 = (r3*2.0)-1.0;
			r3 = r3*agr*sc[jf1]*fact;
	
			xconft[jf1][itt] = xconft[jf1][itt] + r3;
	
			}

	
	enex=0.0;
	energ(itt,xconft,xconf0,stiff,&enex,dd);
	if(enex < 0){printf("enex %lf\n",enex);}


	enerp = enerbase + enex;

	// Call Metropolis algorithm

if(metr_on==1){metropolis(boltz,temp,enertot,enerp,&iflag);}
	else{iflag = 1;}



	if(iflag==1) {	for(i=0;i<6;i++){xconf[i][itt] = xconft[i][itt];}	// What to do when move accepted
			enertot = enerp;					// New total energy is energy from 	accepted move
			ener[itt] = enex;
		      }

	} 
} while(iflag == 2);					// Do this loop to always have a change in helical parameters even when metr_on = 1

//printf("enerp %lf \n", enerp);
//printf("itt %d\n",itt);
if(enex < 0){printf("enertot %lf \n", enertot);}
*endenergy = enertot/23;	// kcal/mol recalculated in eV
*get_itt = itt;
}


// void putinter(itot,xconf,ist,icon,idiv)

void output_str_di(char **output,int itot, char **ist, double **xconf, int phi, char out_folder[])
{
FILE *myfile;
//char str_out[1000];
char seq1;
double buf;
buf = 0.0;
int i,j,k;
j=2;
int len_seq;
len_seq = itot + 1;

//char c[10];
//sprintf(c, "%d", phi); //Make phi as a string
char filename[1000];

//sprintf(filename, "%s/output_dnaflex/output_coordinates_%.6d.dat", out_folder, phi);

//myfile=fopen(filename, "w");

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


if(i==0){sprintf(output[j], "%c-%c %lf %lf %lf %lf %lf %lf \n", ist[i][0],seq1,buf,buf,buf,buf,buf,buf);	
	}
else{
sprintf(output[j], "%c-%c %f %f %f %f %f %f \n", ist[i][0],seq1,xconf[0][i-1],xconf[1][i-1],xconf[2][i-1],xconf[3][i-1],xconf[4][i-1],xconf[5][i-1]);
     }
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

sprintf(output[j+1], "%c-%c %lf %lf %lf %lf %lf %lf \n", ist[itot-1][1],seq1,xconf[0][itot-1],xconf[1][itot-1],xconf[2][itot-1],xconf[3][itot-1],xconf[4][itot-1],xconf[5][itot-1]);

/*for(i=0;i<itot+4;i++){		//Put strings in file
fputs(output[i],myfile);}

//printf("%s \n", output[j+1]);
//printf("%d \n", j+1);		

fclose(myfile);*/

}



void output_str_tetra(char **output,int itot, char **ist, double **xconf, int phi, char out_folder[])
{

// Tetramer output is 'length of sequence - 2' because of comparing tetramers for stiffness matrices

FILE *myfile;
//char str_out[1000];
char seq1;
int i,j;
j=2;
int len_seq;
double zero = 0.0;
len_seq = itot + 1;

//char c[10];
//sprintf(c, "%d", phi); //Make phi as a string
char filename[1000];

//sprintf(filename, "%s/output_dnaflex/output_coordinates_%.6d.dat", out_folder, phi);

//myfile=fopen(filename, "w");

//myfile=fopen("output_coordinates.dat", "w");

sprintf(output[0], "%d base-pairs \n", len_seq);
sprintf(output[1], "Shift  Slide    Rise    Tilt    Roll   Twist \n");

        if (ist[0][1]=='A'){ 
		seq1='T';}
	else if (ist[0][1]=='T') {
		seq1='A';}
	else if (ist[0][1]=='C') {
		seq1='G';}
	else if (ist[0][1]=='G') {
		seq1='C';}
       
sprintf(output[2], "%c-%c %lf %lf %lf %lf %lf %lf \n", ist[0][1],seq1,zero,zero,zero,zero,zero,zero);

			j=j+1;

for(i=1;i<itot;i++){
		
        if (ist[i][1]=='A'){ 
		seq1='T';}
	else if (ist[i][1]=='T') {
		seq1='A';}
	else if (ist[i][1]=='C') {
		seq1='G';}
	else if (ist[i][1]=='G') {
		seq1='C';}
       
sprintf(output[j], "%c-%c %lf %lf %lf %lf %lf %lf \n", ist[i][1],seq1,xconf[0][i-1],xconf[1][i-1],xconf[2][i-1],xconf[3][i-1],xconf[4][i-1],xconf[5][i-1]);

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

sprintf(output[j+1], "%c-%c %lf %lf %lf %lf %lf %lf \n", ist[itot-1][2],seq1,xconf[0][itot-1],xconf[1][itot-1],xconf[2][itot-1],xconf[3][itot-1],xconf[4][itot-1],xconf[5][itot-1]);

/*for(i=0;i<itot+4;i++){		//Put strings in file
fputs(output[i],myfile);}


fclose(myfile);*/

}


void write_hel_coord_in_table(char **table_heli,int itot, int ord[6], double **xconf, int phi, int ind)
{

//,char table_heli_rise[1000][100],char table_heli_tilt[1000][100],char table_heli_roll[1000][100],char table_heli_twis[1000][100],

int i;


for(i=0;i<itot-1;i++){
sprintf(table_heli[i], "%lf ", xconf[ind][i]);
		    }
sprintf(table_heli[itot-1], "%lf", xconf[ind][itot-1]);
sprintf(table_heli[itot], "\n\n");

}



/* ----- Functions for Chromatin model ----- */

void read_nucl(int itot,double **xconfi,char **ist, int *pointer_itot)
{
FILE *myfile;
   int i,j,k;
   double buf;

myfile = fopen ("conf_nucl_merge_3dna_curves.dat", "rt");
//myfile = fopen ("avg_helpar_f1kx5.dat", "rt");

for (i=0; i<itot; i++){
	for(j=0;j<6;j++){
	fscanf(myfile, "%lf", &buf);			
	xconfi[j][i] = buf;		
			}
		}

  fclose(myfile);

/* Read sequence */
  char line[10000000];
  int len;
  char buffer[3];


 myfile=fopen("seq_nucl_1kx5.dat", "rt");


while(fgets(line, 10000000, myfile) != NULL) /* get a line, up to 800 chars from fr.  done if NULL */
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


void transfer_hel_coord(int start, int itot, double **xconf, double **xconfi)
{
int i,j;

for(i=0; i<6; i++){
for(j=start; j<itot; j++){
			xconf[i][j] = xconfi[i][j];  
		      }	
		  }

}


void transfer_hel_coord_gen(int start1, int stop1, double **xconf, int start2, double **xconfi)
{
int i,j,k;

for(i=0; i<6; i++){
		k=start2;
for(j=start1; j<stop1; j++){
			xconf[i][j] = xconfi[i][k];
			k++;  
		      }	
		  }

}


void transfer_stif_gen(int start1, int stop1, double ***dest, int start2, double ***src)
{
int i,j,k,k1,k2;

j=start2;

for(i=start1;i<stop1;i++){
	  for (k1=0; k1<6; k1++) {
	  for (k2=0; k2<6; k2++) {
	     dest[i][k1][k2] = src[j][k1][k2];
				}
				}
		j++;
			}
}


void transfer_seq_gen(int start1, int stop1, char **dest, int start2, char **src)
{
int i,j,k;


		k=start2;
for(j=start1; j<stop1; j++){
for(i=0; i<2; i++){
			dest[j][i] = src[k][i];
			
		   }	
			dest[j][2] = '\0';
			k++;
		  	    }

}

void cmb_link_nucl(double **xconf_mc, double **xconf_nucl, double **xconf_all, char **ist_mc, char **ist_nucl, char **ist_all, int ins_nucl[], int no_nucl, int len_bp_nucl, int len_seq_nucl, int itot_mc, int len_all, int get_itt, int *real_itt)
{
int i,j,k,l,k1,k2;
int buf,buf_mc,buf_nucl,buf_itt,bufbuf_mc,bufbuf;
//char ist_lin_nucl[4], ist_nucl_lin[4];
transfer_hel_coord(0,ins_nucl[0],xconf_all,xconf_mc);	// first linker
transfer_seq_gen(0,ins_nucl[0],ist_all,0,ist_mc);

buf_mc = ins_nucl[0];
buf = ins_nucl[0]+len_bp_nucl;
buf_nucl = ins_nucl[0];

if(buf_mc >= get_itt){buf_itt = get_itt; }

transfer_hel_coord_gen(buf_mc,buf,xconf_all,0,xconf_nucl);	// first nucleosome
transfer_seq_gen(buf_mc,buf_nucl+len_bp_nucl,ist_all,0,ist_nucl);	

for(i=1;i<no_nucl;i++){

bufbuf_mc = buf_mc;
bufbuf = buf;

transfer_hel_coord_gen(buf,buf+ins_nucl[i]-ins_nucl[i-1]+1,xconf_all,buf_mc,xconf_mc);	//+2 because of bp step at end of first and start of second nucl
transfer_seq_gen(buf,buf+ins_nucl[i]-ins_nucl[i-1]+1,ist_all,buf_mc,ist_mc);

buf_mc = buf_mc + ins_nucl[i]-ins_nucl[i-1]+1;
buf = buf+ins_nucl[i]-ins_nucl[i-1]+1;


if(buf_mc > get_itt && bufbuf_mc <= get_itt){buf_itt = (get_itt - bufbuf_mc) + bufbuf; }

transfer_hel_coord_gen(buf,buf+len_bp_nucl,xconf_all,0,xconf_nucl);	// nucleosome
transfer_seq_gen(buf,buf+len_bp_nucl,ist_all,0,ist_nucl);

buf = buf+len_bp_nucl;
}


if(buf_mc <= get_itt){buf_itt = buf + get_itt - buf_mc; }
transfer_hel_coord_gen(buf,len_all,xconf_all,buf_mc,xconf_mc);	// for last linker part
transfer_seq_gen(buf,len_all,ist_all,buf_mc,ist_mc);

//printf("get_itt %d\n",get_itt);
//printf("buf_itt %d\n",buf_itt);

*real_itt = buf_itt;
}



void prepare_mc(double **xconf0, double **xconf0_mc, char **ist_mc, char **ist, char **ist_nucl, int ins_nucl[], int no_nucl, int len_seq_nucl, double **geom_tetra_gen, char **seq_tetra_gen, int itot, int itot_mc, double ***stiff, double ***stiff_mc, double ***stif_tetra_gen)
{
int i,j,k,l,buf,k1,k2,it_conf;
transfer_hel_coord_gen(0,ins_nucl[0]-1,xconf0_mc,0,xconf0);
transfer_stif_gen(0,ins_nucl[0]-1,stiff_mc,0,stiff);
//transfer_hel_coord_gen(ins_nucl[0]-1,ins_nucl[0]+len_bp_nucl,xconf_all,0,xconf_nucl);
transfer_seq_gen(0,ins_nucl[0]-1,ist_mc,0,ist);

buf = ins_nucl[0]-1;

	char ist_lin_nucl[5] = {ist[buf-1][0],ist[buf-1][1],ist_nucl[0][0],ist_nucl[0][1],'\0'};
	char ist_nucl_lin[5] = {ist_nucl[len_seq_nucl-1][0],ist_nucl[len_seq_nucl-1][1],ist[buf+1][0],ist[buf+1][1],'\0'};

	//printf("ist_lin_nucl %s\n",ist_lin_nucl);
	//printf("ist_nucl_lin %s\n",ist_nucl_lin);
	ist_mc[buf][0] = ist[buf-1][1]; ist_mc[buf][1] = ist_nucl[0][0]; ist_mc[buf][2] = '\0';
	ist_mc[buf+1][0] = ist_nucl[len_seq_nucl-1][1]; ist_mc[buf+1][1] = ist[buf+1][0]; ist_mc[buf+1][2] = '\0';
		
		for (j=0; j<256; j++) {
				if (strcmp(ist_lin_nucl,seq_tetra_gen[j]) == 0) 
				{ 	
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff_mc[buf][k1][k2] = stif_tetra_gen[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0_mc[k][buf] = geom_tetra_gen[j][k]; 
							}
				   }


				if (strcmp(ist_nucl_lin,seq_tetra_gen[j]) == 0) 
				{ 	
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff_mc[buf+1][k1][k2] = stif_tetra_gen[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0_mc[k][buf+1] = geom_tetra_gen[j][k]; 
							}
				   }

			      }		

if(no_nucl == 1){
transfer_hel_coord_gen(buf+2,itot_mc,xconf0_mc,buf+1,xconf0);
transfer_stif_gen(buf+2,itot_mc,stiff_mc,buf+1,stiff);
transfer_seq_gen(buf+2,itot_mc,ist_mc,buf+1,ist);
}
else{	it_conf = buf;
	for(i=1;i<no_nucl;i++){
				//printf("buf %d \n",buf);
				//printf("it_conf %d \n",it_conf);
				//buf = ins_nucl[i-1]-1;
				transfer_hel_coord_gen(buf+2,buf+2+ins_nucl[i]-ins_nucl[i-1]-1,xconf0_mc,it_conf+1,xconf0);
				transfer_stif_gen(buf+2,buf+2+ins_nucl[i]-ins_nucl[i-1]-1,stiff_mc,it_conf+1,stiff);
				transfer_seq_gen(buf+2,buf+2+ins_nucl[i]-ins_nucl[i-1]-1,ist_mc,it_conf+1,ist);
				buf = buf+2+ins_nucl[i]-ins_nucl[i-1]-1;
				it_conf = it_conf+ins_nucl[i]-ins_nucl[i-1]+0;

				//printf("buf %d\n",buf);

				char ist_lin_nucl[5] = {ist[ins_nucl[i]-2][0],ist[ins_nucl[i]-2][1],ist_nucl[0][0],ist_nucl[0][1],'\0'};
				char ist_nucl_lin[5] = {ist_nucl[len_seq_nucl-1][0],ist_nucl[len_seq_nucl-1][1],ist[ins_nucl[i]][0],ist[ins_nucl[i]][1],'\0'};

	//printf("ist_lin_nucl %s\n",ist_lin_nucl);
	//printf("ist_nucl_lin %s\n",ist_nucl_lin);

				ist_mc[buf][0] = ist[ins_nucl[i]-2][1]; ist_mc[buf][1] = ist_nucl[0][0]; ist_mc[buf][2] = '\0';
				ist_mc[buf+1][0] = ist_nucl[len_seq_nucl-1][1]; ist_mc[buf+1][1] = ist[ins_nucl[i]][0]; ist_mc[buf+1][2] = '\0';
		
		for (j=0; j<256; j++) {
				if (strcmp(ist_lin_nucl,seq_tetra_gen[j]) == 0) 
				{ 	
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff_mc[buf][k1][k2] = stif_tetra_gen[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0_mc[k][buf] = geom_tetra_gen[j][k]; 
							}
				   }


				if (strcmp(ist_nucl_lin,seq_tetra_gen[j]) == 0) 
				{ 	
				  for (k1=0; k1<6; k1++) {
				  for (k2=0; k2<6; k2++) {
				     stiff_mc[buf+1][k1][k2] = stif_tetra_gen[j][k1][k2];
								}
								}
				  for (k=0; k<6; k++) {
					xconf0_mc[k][buf+1] = geom_tetra_gen[j][k]; 
							}
				   }

			      }	
		

			}
	transfer_hel_coord_gen(buf+2,itot_mc,xconf0_mc,it_conf+1,xconf0);
	transfer_stif_gen(buf+2,itot_mc,stiff_mc,it_conf+1,stiff);
	transfer_seq_gen(buf+2,itot_mc,ist_mc,it_conf+1,ist);

      }


}

void write_bp_pos_in_table(char **table_bp,int itot, int phi, double **pos_store, char out_folder[], int dyn)
{

int i,j,l;
long k;
double zero = 0.0;
char filename[300];

if(dyn==-1){
	    sprintf(filename,"%s/output_cart/cart_%.6d.xyz",out_folder,phi);
	   }
else{
	sprintf(filename,"%s/output_cart/cart_%.6d_%.6d.xyz",out_folder,phi,dyn);
	}


FILE *file_cart;
file_cart=fopen(filename, "w");

//if(phi==0){k=0;}else{k=13*phi;}

	k=0;
	for(i=0;i<3;i++){sprintf(table_bp[k], "%d\n",itot+1);}
	k++;
	for(i=0;i<3;i++){sprintf(table_bp[k], "generated by Juergen Walther\n");}
	k++;
	for(i=0;i<3;i++){sprintf(table_bp[k], "C %.3lf %.3lf %.3lf\n",zero,zero,zero);}
	k++;
	for(i=0;i<itot;i++){
			sprintf(table_bp[k], "C %.3lf %.3lf %.3lf\n",pos_store[i][0],pos_store[i][1],pos_store[i][2]);		
			k++;
			    }
			
for(i=0;i<itot+3;i++){		
fputs(table_bp[i],file_cart);}

fclose(file_cart);

}


void write_bp_triad_in_table(char **table_bp_triad,int itot, int phi, double ***store_orient_t_i, char out_folder[], int dyn)
{

int i,j,l;
long k;
//if(phi==0){k=0;}else{k=13*phi;}
char filename[300];

if(dyn==-1){
	    sprintf(filename,"%s/output_triad/triad_%.6d.dat",out_folder,phi);
	   }
else{
	sprintf(filename,"%s/output_triad/triad_%.6d_%.6d.dat",out_folder,phi,dyn);
	}

FILE *file_triad;
file_triad=fopen(filename, "w");


k=0;

			sprintf(table_bp_triad[k], "%.3lf %.3lf %.3lf\n",1.0,0.0,0.0);
			k++;
			sprintf(table_bp_triad[k], "%.3lf %.3lf %.3lf\n",0.0,1.0,0.0);
			k++;
			sprintf(table_bp_triad[k], "%.3lf %.3lf %.3lf\n",0.0,0.0,1.0);
			k++;


	l=0;
	for(l=0;l<itot;l++){
	for(i=0;i<3;i++){
			sprintf(table_bp_triad[k], "%.3lf %.3lf %.3lf\n",store_orient_t_i[l][i][0], store_orient_t_i[l][i][1], store_orient_t_i[l][i][2]);
			k++;
			}
			//k++;
			}
			//sprintf(table_bp_triad[k], "\n");



for(i=0;i<k;i++){		
fputs(table_bp_triad[i],file_triad);}

fclose(file_triad);

}


void write_pars_in_file(char seq_file[], int ins_nucl[], int no_nucl, int selection, int it_mcsteps, int it_prog, int it_mc, int metr_on, double tot_en, int num_red_mc, int acc_vol, int acc_mc, int do_sim_ann, double tmax, double factor, double lj_cut_up, double dh_cut_up, int rnd_start, int acc, int acc_ph, double e_lj, double e_dh, char out_folder[])
{

int i,j,l;
long k;
char filename[300];
char **table_bp;
table_bp = char_aloc2d(10000,200);
sprintf(filename,"%s/output_schnarp/pars_mcsteps_%d_mc_%d_metr_%d_simann_%d_tot_en_%.4f.dat", out_folder,it_mcsteps,it_mc,metr_on,do_sim_ann,tot_en);
FILE *file_cart;
file_cart=fopen(filename, "w");

	k=0;
	sprintf(table_bp[k], "Sequence %s \n",seq_file);
	k++;
	sprintf(table_bp[k], "Nucleosome positions %d ",ins_nucl[0]);
	k++;
	for(i=1;i<no_nucl;i++){sprintf(table_bp[k], "%d ",ins_nucl[i]);
				k++;
				}
	sprintf(table_bp[k], "\n");
	k++;
	sprintf(table_bp[k], "Selection %d \n",selection);
	k++;
	sprintf(table_bp[k], "it_mcsteps %d \n",it_mcsteps);
	k++;
	sprintf(table_bp[k], "it_prog %d \n",it_prog);
	k++;
	sprintf(table_bp[k], "it_mc %d \n",it_mc);
	k++;
	sprintf(table_bp[k], "metr_on %d \n",metr_on);
	k++;
	sprintf(table_bp[k], "do_sim_ann %d \n",do_sim_ann);
	k++;
	sprintf(table_bp[k], "tmax %lf \n",tmax);
	k++;
	sprintf(table_bp[k], "factor %lf \n",factor);
	k++;
	sprintf(table_bp[k], "lj_cut_up %lf \n",lj_cut_up);
	k++;
	sprintf(table_bp[k], "dh_cut_up %lf \n",dh_cut_up);
	k++;
	sprintf(table_bp[k], "rnd_start %d \n",rnd_start);
	k++;
	sprintf(table_bp[k], "acc %d \n",acc);
	k++;
	sprintf(table_bp[k], "acc_ph %d \n",acc_ph);
	k++;
	sprintf(table_bp[k], "e_lj %lf \n",e_lj);
	k++;
	sprintf(table_bp[k], "e_dh %lf \n",e_dh);
	k++;
	sprintf(table_bp[k], "Final total energy tot_en %lf \n",tot_en);
	k++;
	sprintf(table_bp[k], "Acceptance rate for Volume overlap: %lf\n", (double) acc_vol/ (double) num_red_mc);
	k++;
	sprintf(table_bp[k], "Acceptance rate for overall MC: %lf\n", (double) acc_mc/ (double) num_red_mc);


	
for(i=0;i<k+1;i++){		
fputs(table_bp[i],file_cart);}

free(table_bp);
fclose(file_cart);

}


void write_acc_mc_moves_in_table(int sel, int *acc_arr, int accepted,int it_mcsteps, int phi, int it_mc, int metr_on, int do_sim_ann, char stif_file[], char seq_file[], char out_folder[])
{

char **table_acc;
int i,j,l,n,m;
long k;
double zero = 0.0;
char filename[300], seq[100], stif[100];
table_acc = char_aloc2d(accepted+1,10);
n = (int)strlen(seq_file);
m = (int)strlen(stif_file);
//printf("n %d m %d\n",n,m); 
for(i=0;i<(n-4);i++){seq[i] = seq_file[i];} seq[n-4] = '\0';
for(i=0;i<(m-4);i++){stif[i] = stif_file[i];} stif[m-4] = '\0';

if(sel == 0){sprintf(filename,"%s/output_acc/acc_dna_mc_mcsteps_%d_mc_%d_metr_%d_simann_%d_%d.dat", out_folder,it_mcsteps,it_mc,metr_on,do_sim_ann,phi);}
if(sel == 1){sprintf(filename,"%s/output_acc/acc_vol_mcsteps_%d_mc_%d_metr_%d_simann_%d_%d.dat", out_folder,it_mcsteps,it_mc,metr_on,do_sim_ann,phi);}
if(sel == 2){sprintf(filename,"%s/output_acc/acc_wholemc_mcsteps_%d_mc_%d_metr_%d_simann_%d_%d.dat", out_folder,it_mcsteps,it_mc,metr_on,do_sim_ann,phi);}
FILE *file_acc;
file_acc=fopen(filename, "w");

	for(i=0;i<accepted;i++){
		sprintf(table_acc[i], "%d ",acc_arr[i]);		
			    }
	sprintf(table_acc[accepted], "\n");
			
for(i=0;i<(accepted+1);i++){		
fputs(table_acc[i],file_acc);}

fclose(file_acc);
free(table_acc);

}


void write_nucl_pos_in_table(char **table_nucl,int no_nucl, int phi, double **nucl_pos, char out_folder[], int dyn)
{

int i,j,l;
long k;
double zero = 0.0;
char filename[300];

if(dyn==-1){
	    sprintf(filename,"%s/output_nucl/nucl_%.6d.xyz",out_folder,phi);
	   }
else{
	sprintf(filename,"%s/output_nucl/nucl_%.6d_%.6d.xyz",out_folder,phi,dyn);
	}
FILE *file_cart;
file_cart=fopen(filename, "w");

	k=0;
	for(i=0;i<3;i++){sprintf(table_nucl[k], "%d\n",no_nucl);}
	k++;
	for(i=0;i<3;i++){sprintf(table_nucl[k], "generated by Juergen Walther\n");}
	k++;
	for(i=0;i<no_nucl;i++){
			sprintf(table_nucl[k], "N %.3lf %.3lf %.3lf\n",nucl_pos[i][0],nucl_pos[i][1],nucl_pos[i][2]);		
			k++;
			    }
		
for(i=0;i<no_nucl+2;i++){		
fputs(table_nucl[i],file_cart);}

fclose(file_cart);

}


void write_helpars_in_table(int len_all, int phi, double **xconf_all, char out_folder[], int dyn)
{

int i,j,l;
long k;
double zero = 0.0;
char filename[300];
char **table_ener;

table_ener  = char_aloc2d(len_all+1,200);

if(dyn==-1){
	    sprintf(filename,"%s/output_dnaflex/helpars_%.6d.dat",out_folder,phi);
	   }
else{
	sprintf(filename,"%s/output_dnaflex/helpars_%.6d_%.6d.dat",out_folder,phi,dyn);
	}

FILE *file_cart;
file_cart=fopen(filename, "w");

	k=0;
	sprintf(table_ener[k],"shift slide rise tilt roll twist \n");
	k++;
	for(i=0;i<len_all;i++){
			sprintf(table_ener[k], "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",xconf_all[0][i], xconf_all[1][i], xconf_all[2][i], xconf_all[3][i], xconf_all[4][i], xconf_all[5][i]);
			k++;
			    }
		
for(i=0;i<(len_all+1);i++){		
fputs(table_ener[i],file_cart);}

fclose(file_cart);

free(table_ener);
}


void write_energy_in_table(int it_en, int phi, double **save_en, char out_folder[])
{

int i,j,l;
long k;
double zero = 0.0;
char filename[300];
char **table_ener;

table_ener  = char_aloc2d(it_en+1,200);

sprintf(filename,"%s/output_schnarp/energies_%.6d.dat",out_folder,phi);
FILE *file_cart;
file_cart=fopen(filename, "w");

	k=0;
	sprintf(table_ener[k],"psi acc_mc e_lj_bp_bpfl e_dh_ph_phfl e_nucl_lj e_nucl_dh endenergy/acc tot_en \n");
	k++;
	for(i=0;i<it_en;i++){
			sprintf(table_ener[k], "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",save_en[i][0],save_en[i][1],save_en[i][2],save_en[i][3],save_en[i][4],save_en[i][5], save_en[i][6], save_en[i][7]);
			k++;
			    }
		
for(i=0;i<(it_en+1);i++){		
fputs(table_ener[i],file_cart);}

fclose(file_cart);

free(table_ener);
}


/*--------------------End of Functions------------------------ */

