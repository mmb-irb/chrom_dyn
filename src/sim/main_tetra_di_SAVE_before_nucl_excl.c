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

#include "dna_flex.h"	// Functions used in DNA Flex
#include "umbrella.h"	// Functions needed for Cartesian reconstruction
#include "energy.h"	// Functions for calculating int energies

/* Headers from SCHNArP */

#include "schna_ar_p.h"

#include <time.h>

/* End Headers */







/*--------------------------Beginning of DNAFlex-------------------------*/


int main(void)
{

char out_folder[] = "/orozco/homes/pluto/jwalther/Programs/Chromatin/src_test";

char seq_file[] = "seq_150mer.dat";							
char stif_file[] = "stif_150mer.dat";
int ins_nucl[] = {20,100};

/* Initialize all variables at beginning of code!!*/ 

int nmc,i0,ibuff2,juega,last,psi;
int iflag; 	//it should be a bool, but I dont know how to assign true and false to a bool in c
char **seq1, **ist; //seq1[1000][10], ist[1000][10];
double ***stif, ***stiff, **geom; // stif[1000][6][6],stiff[1000][6][6],geom[1000][6];
int ord[6];
double **xconfi, **xconf0, **xconf; // xconfi[6][1000],xconf0[6][1000],xconf[6][1000];
double **xconft, *ener; // xconft[6][1000],ener[1000]
double sc[6],dd[6];
double **dx; //dx[6][1000];
char *names[6];   				// No maximum size of char in C
char line;
double boltz=0.00198717;
double enex;
int itot;
int i,j,k;	// Iteration parameters
int it_mc;	// Test parameter for # MC moves
char **output;  //output[1000][100];			// Writing helical coordinates in this string
int it_prog,phi;
double ind_energy;
double *ind_ener;
float sum_energy;
int selection;
double ***stif_tetra, ***stif_di, **geom_tetra, **geom_di; //stif_tetra[1000][6][6],stif_di[1000][6][6],geom_tetra[1000][6],geom_di[1000][6];
char **seq1_tetra, **seq1_di;  //seq1_tetra[1000][10], seq1_di[1000][10];
int it_snapshots, it_mcsteps, alpha, it_mod;
char **out_part;
double endenergy;

double **xyz_float;
int msec,msec1,msec2,msec3;
double **center_bp;
double r_bp;
char **pdb_data;
char **tmp_at;
char *baset;
char *strand;
int *base_num;
double **tmp_xyz;
double *q_fl;
int num_bp;
double **d_bp;
double **d_bp_fl;
double **lj_bp, **lj_bp_fl;
double **mid_bp;
double **dir_bp;
double **dir_ph;
double **pos_ph;
double **dh_ph;
double **dh_ph_fl;
double q_ph;
double **d_ph1, **d_ph2;
double **d_ph_fl;
double tot_dh_lj, tot_lj, tot_dh;
int k1,k2,k3,k4;
char **pdb_tot;
int tot_num_at;
int num_fl;	// Number of floating patricles
double shift;
double shift_in;
int acc;
int acc_ph;
double bp_excl_vol;
double bp_fl_excl_vol;
double lj_k,lj_cut_up,lj_cut_low,lj_a_bp,lj_a_fl;
double eps,kappa,f_pi_eps0,dh_cut_low,dh_cut_up;
double **xconf_buf;
double tot_en,prev_tot_en;
double d_ph_midbp;
double **xyz_float_buf;
int count_no_excl;


double ***orient_old, ***orient_new, **pos_new, ***orient_store, **pos_store, ***rot_store, ***rot_new;
double **pos_old;

int no_rot_pars = 3;
int vol_olap_start;

char **ist_nucl, **ist_all, **ist_mc;
char **seq_tetra, **seq_tetra_gen;
double **geom_tetra_gen, ***stif_tetra_gen, **xconf_all;
double **xconf_mc, **xconf0_mc, ***stiff_mc, **xconfi_mc, **xconft_mc, **xconf_nucl;

int len_seq;	// Variable to determine dynamic size of all arrays
int seq_pcs;	// Length of sequence pieces to save (i.e. 4 for tetramer seq)
int no_helpars;
int len_out_line;
int len_bp_nucl,len_seq_nucl;
int len_all,bp_all;
int itot_mc;
int no_nucl = sizeof(ins_nucl)/sizeof(ins_nucl[0]);
double buf_tot_en = 1000000;
double buf_enertot = 1000000;

len_bp_nucl = 146;
seq_pcs = 10;
no_helpars = 6;
len_out_line = 100;

num_fl = 0;			// Number of floating particles


get_seq_len(&len_seq,&itot,seq_file);

printf("len_seq %d \n",len_seq);

itot_mc = itot + no_nucl;
printf("itot_mc %d\n",itot_mc);

len_all = no_nucl*len_bp_nucl+itot_mc;
printf("len_all %d\n",len_all);

bp_all = len_all + 1;

int alc_space;
if(len_seq < 256){alc_space = 256;}else{alc_space = len_seq;}

/* Allocate memory for arrays */
seq1 = char_aloc2d(alc_space,seq_pcs);
ist  = char_aloc2d(len_seq,seq_pcs);
stif = dbl_aloc3d(alc_space,no_helpars,no_helpars);
stiff = dbl_aloc3d(len_seq,no_helpars,no_helpars);
geom = dbl_aloc2d(alc_space,no_helpars);
xconfi = dbl_aloc2d(no_helpars,len_seq);
xconf0 = dbl_aloc2d(no_helpars,len_seq);
xconf = dbl_aloc2d(no_helpars,len_seq);
xconft = dbl_aloc2d(no_helpars,len_seq);
ener = malloc(itot_mc*sizeof(double));
dx = dbl_aloc2d(no_helpars,itot_mc);
output  = char_aloc2d(len_all+5,len_out_line);

seq1_tetra = char_aloc2d(256,seq_pcs);
seq1_di = char_aloc2d(16,seq_pcs);
stif_tetra = dbl_aloc3d(256,no_helpars,no_helpars);
stif_di = dbl_aloc3d(16,no_helpars,no_helpars);
geom_tetra = dbl_aloc2d(256,no_helpars);
geom_di = dbl_aloc2d(16,no_helpars);

out_part  = char_aloc2d(len_seq,len_out_line);
xyz_float = dbl_aloc2d(num_fl,no_helpars);
q_fl = malloc(len_seq*sizeof(double));
center_bp = dbl_aloc2d(bp_all,no_helpars);
/*pdb_data = char_aloc2d(2*len_seq,len_out_line);
tmp_at = char_aloc2d(2*len_seq,no_helpars);
baset = malloc(2*len_seq*sizeof(char));
strand = malloc(2*len_seq*sizeof(char));
base_num = malloc(2*len_seq*sizeof(int));
tmp_xyz = dbl_aloc2d(2*len_seq,no_helpars);*/
d_bp = dbl_aloc2d(bp_all,bp_all);
d_bp_fl = dbl_aloc2d(bp_all,bp_all);
lj_bp = dbl_aloc2d(2*bp_all,no_helpars);
lj_bp_fl = dbl_aloc2d(2*bp_all,no_helpars);
mid_bp = dbl_aloc2d(bp_all,no_helpars);
dir_bp = dbl_aloc2d(bp_all,no_helpars);
dir_ph = dbl_aloc2d(bp_all,no_helpars);
pos_ph = dbl_aloc2d(2*bp_all,no_helpars);
d_ph1 = dbl_aloc2d(2*bp_all,2*bp_all);
d_ph2 = dbl_aloc2d(2*bp_all,2*bp_all);
d_ph_fl = dbl_aloc2d(2*bp_all,bp_all);
dh_ph = dbl_aloc2d(4*bp_all,no_helpars);
dh_ph_fl = dbl_aloc2d(2*bp_all,no_helpars);

xconf_buf = dbl_aloc2d(no_helpars,itot_mc);
xyz_float_buf = dbl_aloc2d(num_fl,no_helpars);

pdb_tot = char_aloc2d(20*len_seq,len_out_line);

orient_old = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
orient_new = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
orient_store = dbl_aloc3d(len_all,no_rot_pars,no_rot_pars);
pos_old = dbl_aloc2d(len_all,no_rot_pars);
pos_new = dbl_aloc2d(len_all,no_rot_pars);
pos_store = dbl_aloc2d(len_all,no_rot_pars);
rot_new = dbl_aloc3d(2*len_all,no_rot_pars,no_rot_pars);
rot_store = dbl_aloc3d(2*len_all,no_rot_pars,no_rot_pars);

xconf_nucl = dbl_aloc2d(no_helpars,len_bp_nucl);
ist_nucl  = char_aloc2d(len_bp_nucl+1,seq_pcs);
ist_all  = char_aloc2d(len_all,seq_pcs);
xconf_all  = dbl_aloc2d(no_helpars,len_all);
seq_tetra_gen = char_aloc2d(256,seq_pcs);
stif_tetra_gen = dbl_aloc3d(256,no_helpars,no_helpars);
geom_tetra_gen = dbl_aloc2d(256,no_helpars);


xconf_mc = dbl_aloc2d(no_helpars,itot_mc);
xconf0_mc = dbl_aloc2d(no_helpars,itot_mc);
stiff_mc = dbl_aloc3d(itot_mc,no_helpars,no_helpars);
xconfi_mc = dbl_aloc2d(no_helpars,itot_mc);
xconft_mc = dbl_aloc2d(no_helpars,itot_mc);
ist_mc  = char_aloc2d(itot_mc,seq_pcs);


/* Initialize arrays and files for huge tables */

char **table_heli_shif;
char **table_heli_slid;
char **table_heli_rise;
char **table_heli_tilt;
char **table_heli_roll;
char **table_heli_twis;

table_heli_shif  = char_aloc2d(len_seq,len_out_line);
table_heli_slid  = char_aloc2d(len_seq,len_out_line);
table_heli_rise  = char_aloc2d(len_seq,len_out_line);
table_heli_tilt  = char_aloc2d(len_seq,len_out_line);
table_heli_roll  = char_aloc2d(len_seq,len_out_line);
table_heli_twis  = char_aloc2d(len_seq,len_out_line);

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


selection = 3;			// 0 = individual, 1 = dimer, 2 = tetramer, 3 = tetramer with dimer end


it_mcsteps = 10000;		// number of MC steps in floating process before looking at position of particle
it_snapshots = 1;		// number of pdb's generated during the floating process
it_prog = 1;			// number of executing all the MC and reconstruction part (THE WHOLE PROGRAM)
it_mc = 5*it_mcsteps;		// number of MC steps to bring structure to minimum energy before particle inserted
				/* 100000 -> 56mer
				   10000000 -> 1000mer
				   2000000 -> 450mer
				   5000000 -> 750mer
				*/


/* ------------------------------- END Looping and selection variables ------------------------------- */


/* ------------------------------- Parameters for dynamical part ------------------------------------- */


shift_in = 600.0;		// Volume occupied by floating particles
shift = 5.0;			// Range of Shift of floating particle per MC move

/* ---- Parameters for excluded volume ---- */
acc = 10;			// bp-steps for calculating distance
acc_ph = 20;			// phosphate = bp steps for calculating distance
bp_excl_vol = 20;		// minimum distance between two bp centers
bp_fl_excl_vol = 10;		// minimum distance between bp center and floating particle

d_ph_midbp = 10.25;		// Distance of Phosphate of helical axis



/* ---- Parameters for LJ ---- */	
lj_k = 4*0.0257; 		// in eV (0.0257 eV = 1kT, adapted from Arya 2014)
lj_cut_up = 1.5*bp_excl_vol;	
lj_cut_low = 0.1;

lj_a_bp = bp_excl_vol;		// = 10A as excluded volume (radius of excluded volume around center_bp pos)
lj_a_fl = bp_fl_excl_vol;	// 



/* ---- Parameters for DH ---- */
eps = 80.0;
kappa = 0.033; 			// unit: 1/A	taken from Arya paper
f_pi_eps0 = 1.1126*pow(10,-10);
q_ph = -1.0;

dh_cut_low = 0.1;
dh_cut_up = bp_excl_vol;	// 20.0


/* ------------------------------- END Parameters for dynamical part ----------------------------------- */


order(names,ord);



if(selection == 0){		// Select individual bsc1 stiffness matrices
readsq_ind(&itot,ist,seq_file);

readss_ind(stif,seq1,geom,ord,seq_file,stif_file);

assign_ind(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 1){		// Select average bsc0 dimer stiffness matrices
readsq_di(&itot,ist,seq_file);

readss_di(stif,seq1,geom,ord);

assign_di(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 2){		// Select average bsc0 tetramer stiffness matrices
readsq_tetra(&itot,ist,seq_file);

readss_tetra(stif,seq1,geom,ord);

assign_tetra(itot,ist,seq1,stif,stiff,geom,xconf0);
}


if(selection == 3){		// Select average bsc0 tetramer stiffness matrices with bsc0 dimer at end
readsq_tetra_di(&itot,ist,seq_file);

readss_tetra_di(stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,ord);

assign_tetra_di(itot,ist,stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,stiff,xconf0);
}



// Get stiffness tetramer values
readss_tetra(stif_tetra_gen,seq_tetra_gen,geom_tetra_gen,ord);
// Get coordinates for nucleosome
read_nucl(len_bp_nucl,xconf_nucl,ist_nucl,&len_seq_nucl);

readini(itot,xconfi);
// Set working coordinates and energies equal to the starting ones
transfer_hel_coord(0,itot,xconf,xconfi);
// Prepare MC sequence xconf0_mc, stiff_mc
prepare_mc(xconf0,xconf0_mc,ist_mc,ist,ist_nucl,ins_nucl,no_nucl,len_seq_nucl,geom_tetra_gen,seq_tetra_gen,itot,itot_mc,stiff,stiff_mc,stif_tetra_gen);
// Initialize starting linker configuration
readini(itot_mc,xconfi_mc);
// Set working coordinates and energies equal to the starting ones
transfer_hel_coord(0,itot_mc,xconf_mc,xconfi_mc);

//for(i=0;i<itot_mc;i++){for(j=0;j<6;j++){printf("%f ",xconf0_mc[j][i]);} printf(" \n");}
//for(i=0;i<len_all;i++){for(j=0;j<6;j++){printf("%f ",xconf_all[j][i]);} printf(" \n");}
//for(i=0;i<len_all;i++){printf("%s \n",ist_all[i]);}
//for(i=0;i<itot_mc;i++){printf("%s \n",ist_mc[i]);}


num_bp = itot+1;

// Compute energy
	energy(itot_mc,xconfi_mc,xconf0_mc,dx,stiff_mc,ener);


	double ener0=0.0;
	double enertot;

	for(i=0; i<itot_mc; i++){ener0 = ener0 + ener[i];}

	enertot = ener0;

	buf_enertot = enertot;	

	int num_red_mc = it_mc/it_mcsteps;

	
	printf("num_red_mc %d\n",num_red_mc);

	/* --- Introduce floating particle --- */
	double br_mot[3];	// Initial values floating particle
	br_mot[0] = 0.0;
	br_mot[1] = 0.0;
	br_mot[2] = 0.0;
	
	double rx,ry,rz;	

/* ------------------------------- Looping starts ------------------------------- */

ind_ener = malloc(it_snapshots*sizeof(double));

for(phi=0;phi<it_prog;phi++){			// do optimization it_prog times to check if MC converges


	clock_t start = clock(), diff; //Timer


	/*---------------------------------MC algorithm-----------------------------*/

	/* --- Initialize floating particle --- */

	initial_fl_pos(num_fl,shift_in,br_mot,xyz_float);

	// Charge configuration of floating particles

	for(i=0;i<num_fl;i++){
			q_fl[i] = 1.0;
			      }

    for(psi=0;psi<num_red_mc;psi++){
	printf("psi %d\n",psi);

	/*for(j=0; j<(itot_mc); j++){
	for(i=0; i<6; i++){
			printf("%lf ",xconf_buf[i][j]);  
		      }	
			printf("\n");
		  }
printf("\n");*/


	transfer_hel_coord(0,itot_mc,xconf_buf,xconf_mc);
	
	


		/*for(j=0; j<itot_mc; j++){
	for(i=0; i<6; i++){
			printf("%lf ",xconf_mc[i][j]);  
		      }	
			printf("\n");
		  }*/

	// Do Monte Carlo
	monte_carlo(it_mcsteps,itot_mc,enertot,ener,xconft_mc,xconf_mc,sc,xconf0_mc,stiff_mc,enex,dd,boltz,temp,iflag, &endenergy);

	enertot = endenergy/0.043424;

	// Combine linker and nucleosomal DNA
	cmb_link_nucl(xconf_mc,xconf_nucl,xconf_all,ist_mc,ist_nucl,ist_all,ins_nucl,no_nucl,len_bp_nucl,len_seq_nucl,itot_mc,len_all);
	
	// Cartesian reconstruction of whole structure used for MC
	cart_rec_old(len_all,no_rot_pars,xconf_all,pos_old,orient_old,rot_store);

	/*---------------------------------End of MC algorithm-----------------------------*/

	//printf("enertot %lf \n", enertot);
	//printf("endenergy %lf \n", endenergy);	
	//printf("buf_enertot %lf \n", buf_enertot);

	diff = clock() - start;

	msec = diff * 1000 / CLOCKS_PER_SEC;

	/* --------------------------- Check if created structure overlaps --------------------------- */

	/* ---- Position of bp and floating particles and calculate if overlap ---- */

	check_olap(bp_all,num_fl,len_all,acc,center_bp,pos_old,d_bp,d_bp_fl,xyz_float,bp_excl_vol,bp_fl_excl_vol,&vol_olap_start);

	printf("Overlap starting structure1 %d\n",vol_olap_start);

	if(vol_olap_start == 1){if(psi>0){enertot = buf_enertot;}
				psi--;
				transfer_hel_coord(0,itot_mc,xconf_mc,xconf_buf);
				printf("Start new\n");
				}
     	else{
	
		buf_enertot = enertot;

	/* ---- Get DH + LJ for starting configuration ---- */

	
	/* -------------- Calculate Phosphate position --------------*/

	ph_pos(bp_all,len_all,d_ph_midbp,mid_bp,center_bp,dir_bp,orient_old,dir_ph,pos_ph);
	
	// Get distances between ph and between ph and floating particles
	
	dist_ph(bp_all,num_fl,acc_ph,pos_ph,d_ph1,d_ph2,d_ph_fl,xyz_float);

	//for(i=0;i<2*bp_all;i++){for(j=0;j<2*bp_all;j++){if(d_ph1[i][j] > 0){printf("%lf ",d_ph1[i][j]);}} printf("\n");}

	/* ---- Calculate Lenard-Jones ---- */

	calc_lj(&k1,&k2,bp_all,d_bp,lj_cut_low,lj_cut_up,lj_bp,lj_k,lj_a_bp,acc,num_fl,lj_bp_fl,d_bp_fl,lj_a_fl);

	/* ---------- Calculate DH potential ---------- */

	calc_dh(&k3,&k4,bp_all,dh_ph,d_ph1,d_ph2,dh_cut_low,dh_cut_up,q_ph,eps,kappa,f_pi_eps0,acc_ph,num_fl,dh_ph_fl,d_ph_fl,q_fl);

	printf("LJ %d %d\n",k1,k2);
	printf("DH %d %d\n",k3,k4);
	/*----Total energy (DH ph-ph/ph-fl, LJ bp-bp/bp_fl)----*/

	tot_dh = 0;
	tot_lj = 0;
	tot_dh_lj = 0;
	for(i=0;i<k1;i++){tot_dh_lj = tot_dh_lj + lj_bp[i][0];
			  tot_lj = tot_lj + lj_bp[i][0];}
	for(i=0;i<k2;i++){tot_dh_lj = tot_dh_lj + lj_bp_fl[i][0];
			  tot_lj = tot_lj + lj_bp_fl[i][0];}
	for(i=0;i<k3;i++){tot_dh_lj = tot_dh_lj + dh_ph[i][0];
			  tot_dh = tot_dh + dh_ph[i][0];}
	for(i=0;i<k4;i++){tot_dh_lj = tot_dh_lj + dh_ph_fl[i][0];
			  tot_dh = tot_dh + dh_ph_fl[i][0];}

	tot_en = tot_dh_lj + endenergy;

	printf("Total energy LJ + DH + elastic = %lf + %lf + %lf = %lf eV\n",tot_lj,tot_dh,endenergy,tot_en);	


	metropolis(boltz,temp,buf_tot_en,tot_en,&iflag);	// Metropolis for total energy
	printf("iflag %d\n",iflag);
	printf("tot_en %lf\n",tot_en);
	printf("buf_tot_en %lf\n",buf_tot_en);
		if(iflag==2){	// rejected
				psi--;
				transfer_hel_coord(0,itot_mc,xconf_mc,xconf_buf);
				enertot = buf_enertot;
	      		     }
		else{buf_tot_en = tot_en;}	// accepted
	}
	printf("End if statement vol_olap_start\n");
    }
    printf("End loop with psi\n");
	/* ----- Output pdb ----- */

	if(selection == 0 || selection == 1){
	output_str_di(output,len_all,ist_all,xconf_all,phi,out_folder);
	}
	if(selection == 2 || selection == 3){
	output_str_tetra(output,itot,ist,xconf,phi,out_folder);
			  }

	CEHS_build(len_all,ist_all,xconf_all,phi,output,out_folder);
	

/* ---------- Loop the manipulation part DELETE THIS PART BECAUSE WE ONLY WANT TO GET END STRUCTURES (KEEP CALCULATION OF POTENTIALS!!)---------- */

	count_no_excl = 0;
	it_mod = 0;
	alpha = 0;




} 
/* ------------------- End of for loop to execute the whole program it_prog times ------------------- */


fclose(file_shif);
fclose(file_slid);
fclose(file_rise);
fclose(file_tilt);
fclose(file_roll);
fclose(file_twis);



sum_energy=0.0;
for(i=0;i<it_snapshots;i++){
sum_energy = sum_energy + ind_ener[i];
//printf("%f \n", ind_ener[i]);
}
free (ind_ener);


ind_energy = sum_energy/(double) it_snapshots;



printf("Time taken for stabilizing MC: %d seconds %d milliseconds\n", msec/1000, msec%1000);
//printf("%d steps needed to get %d structures that are %d times energy-accepted by Metropolis = %lf efficiency\n", it_mod, count_no_excl, alpha ,(double) alpha/ (double) it_mod);
printf("Time taken for one of the %d variations with %d MC steps: %d seconds %d milliseconds\n", it_snapshots, it_mcsteps, msec2/1000, msec2%1000);
printf("Time taken for reconstruction: %d seconds %d milliseconds\n", msec1/1000, msec1%1000);
printf("Mean energy after %d stabilization steps and modifying %d times by %d MC steps: %f eV\n", it_mc, it_snapshots, it_mcsteps, ind_energy);



diff_all = clock() - start_all;

int msec_all = diff_all * 1000 / CLOCKS_PER_SEC;
printf("Time taken for whole algorithm: %d seconds %d milliseconds\n", msec_all/1000, msec_all%1000);


/* Free memory of arrays */

free(seq1);
free(ist);
free_dbl3d(stif,alc_space);
free_dbl3d(stiff,len_seq);
free(geom);
free(xconfi);
free(xconf0);
free(xconf);
free(xconft);
free(ener);
free(dx);
free(output);

free(seq1_tetra);
free(seq1_di);
free_dbl3d(stif_tetra,256);
free_dbl3d(stif_di,16);
free(geom_tetra);
free(geom_di);
//printf("Test11\n");
free(table_heli_shif);
free(table_heli_slid);
free(table_heli_rise);
free(table_heli_tilt);
free(table_heli_roll);
free(table_heli_twis);
//printf("Test12\n");
free(out_part);
free(xyz_float);
free(q_fl);
free(center_bp);
/*free(pdb_data);
free(tmp_at);
free(baset);
free(strand);
free(base_num);
free(tmp_xyz);*/
free(d_bp);
free(d_bp_fl);
free(lj_bp);
free(lj_bp_fl);
free(mid_bp);
free(dir_bp);
free(dir_ph);
free(pos_ph);
free(d_ph1);
free(d_ph2);
free(d_ph_fl);
free(dh_ph);
free(dh_ph_fl);
//printf("Test13\n");
free(xconf_buf);
free(xyz_float_buf);

free(pdb_tot);
//printf("Test14\n");
free_dbl3d(orient_old,len_all);
free_dbl3d(orient_new,len_all);
free_dbl3d(orient_store,len_all);
free(pos_old);
free(pos_new);
free(pos_store);
free_dbl3d(rot_new,2*len_all);
free_dbl3d(rot_store,2*len_all);
//printf("Test15\n");
free(xconf_nucl);
free(ist_nucl);
free(ist_all);
free(xconf_all);
free(seq_tetra_gen);
free_dbl3d(stif_tetra_gen,256);
free(geom_tetra_gen);
//printf("Test16\n");
free(xconf_mc);
free(xconf0_mc);
free_dbl3d(stiff_mc,itot_mc);
free(xconfi_mc);
free(xconft_mc);
free(ist_mc);


/* End free memory */



}
