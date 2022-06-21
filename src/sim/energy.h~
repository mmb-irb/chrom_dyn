/* ---- Functions for calculating interaction energies ---- */

double pow2(double x)
{
  return x*x;
}


double dist (double *at1, double *at2)
{
int i;
double dist,buf;
dist = 0; buf = 0;
for(i=0;i<3;i++){
		buf = pow2(at1[i]-at2[i]);
		dist = dist + buf;
		}

		dist = sqrt(dist);
return dist;
}

/*double vec_len (double *vec)
{
int i;
double len=0.0;

for(i=0;i<3;i++){
		len = len + pow2(vec[i]);
		}
	len = sqrt(len);

}*/

double lj (double r, double k, double a, double b)
{
double lj;

lj = 4*k*(pow((a/r),12) - pow((b/r),6));

return lj;
}

double dh (double r, double q1, double q2, double eps, double kappa, double f_pi_eps0)
{
double dh;

dh = q1*1.6*pow(10,-19) * q2*1.6*pow(10,-19) / (f_pi_eps0*eps) * 1/r * pow(10,10) * exp(-kappa*r) /(1.6*pow(10,-19));	// 10^10 wegen r in Angstrom

return dh;
}

void positions(int num_bp, int num_fl, int itot, int acc, double **center_bp, double **pos, double **d_bp, double **d_bp_fl, double **xyz_float, int red_bp_all)
{
int i,j,k,m;
	
for(i=0;i<3;i++){center_bp[0][i] == 0.0;}
for(i=1;i<num_bp;i++){		
for(j=0;j<3;j++){
		center_bp[i][j] = pos[i-1][j];
		}
		     }
				

	// Get distances between bp and between bp and floating particles
	
	//m = (int) (double) (num_bp-1)/(double)acc;
	k=1;
	for(i=0;i<red_bp_all;i++){
	for(j=i+1;j<red_bp_all;j++){ //num_bp-acc-i
				if(j*acc < num_bp){
					d_bp[i][j] = dist(center_bp[i*acc],center_bp[j*acc]);
					//if(d_bp[i][k-1] < 20){printf("%lf %d %d %d\n",d_bp[i][k-1],i,k-1,i+k*acc);}
					k++;
						     }
			     }
				k=1;
				}



	// Get distances between bp and between bp and floating particles (OLD WAY)

	/*k=1;
	for(i=0;i<num_bp;i++){
	for(j=0;j<num_bp;j++){ //num_bp-acc-i
				if(i+k*acc < num_bp){
					d_bp[i][k-1] = dist(center_bp[i],center_bp[i+k*acc]);
					//if(d_bp[i][k-1] < 20){printf("%lf %d %d %d\n",d_bp[i][k-1],i,k-1,i+k*acc);}
					k++;
						     }
			     }
				k=1;
				}*/
	

	for(i=0;i<num_bp;i++){
	for(j=0;j<num_fl;j++){
				d_bp_fl[i][j] = dist(center_bp[i],xyz_float[j]);
			     }
				}

}

void ph_pos(int num_bp, int itot, double d_ph_midbp, double **mid_bp, double **center_bp, double **dir_bp, double ***orient_old, double **dir_ph, double **pos_ph)
{
int i,j,k;
const float pi = 3.14159265358979323846;

double **x_axis_bp, **y_axis_bp, **mid_bp_x, **mid_bp_y;
x_axis_bp = dbl_aloc2d(num_bp,3);
y_axis_bp = dbl_aloc2d(num_bp,3);
mid_bp_x = dbl_aloc2d(num_bp,3);
mid_bp_y = dbl_aloc2d(num_bp,3);


	x_axis_bp[0][0] = 1; x_axis_bp[0][1] = 0; x_axis_bp[0][2] = 0;
	y_axis_bp[0][0] = 0; y_axis_bp[0][1] = -1; y_axis_bp[0][2] = 0;	
	
	for(i=1;i<num_bp;i++){
	for(j=0;j<3;j++){
	x_axis_bp[i][j] = orient_old[i-1][0][j];		// normalized
	y_axis_bp[i][j] = -1 * orient_old[i-1][1][j];		// -1 to get right phosphate position
			}
				}
	

	for(i=0;i<num_bp-1;i++){
			mid_bp[i][0] = 0.5*(center_bp[i+1][0] + center_bp[i][0]);
			mid_bp[i][1] = 0.5*(center_bp[i+1][1] + center_bp[i][1]);
			mid_bp[i][2] = 0.5*(center_bp[i+1][2] + center_bp[i][2]);		
				}


	double cos_x,cos_y;	
	double buf_x, buf_y;	

	for(i=0;i<num_bp-1;i++){
	cos_x = dot_product(x_axis_bp[i],x_axis_bp[i+1]);
	cos_y = dot_product(y_axis_bp[i],y_axis_bp[i+1]);

	//printf("%lf %lf %lf %lf \n", cos_x, cos(acos(cos_x)), cos_y, cos(acos(cos_y)));

	

	for(j=0;j<3;j++){
	mid_bp_x[i][j] = 0.5*cos_x*x_axis_bp[i][j] + 0.5*x_axis_bp[i+1][j];
	mid_bp_y[i][j] = 0.5*cos_y*y_axis_bp[i][j] + 0.5*y_axis_bp[i+1][j];
			}
	
	buf_x = vec_len(mid_bp_x[i],3);
	buf_y = vec_len(mid_bp_y[i],3);

	for(j=0;j<3;j++){
	mid_bp_x[i][j] = mid_bp_x[i][j] / buf_x;
	mid_bp_y[i][j] = mid_bp_y[i][j] / buf_y;
			}

				}


	double a1,a2,a_adj;
	a_adj = 12;		// Adjust angle values
	a1 = 33+90-a_adj;	
	a2 = 147+90+a_adj;

	//printf("%lf %lf \n",cos(a1*PI/180),sin(a1*PI/180));

	k=0;
	for(i=0;i<num_bp-1;i++){
	for(j=0;j<3;j++){
	pos_ph[k][j] = mid_bp[i][j] + cos(a2*pi/180)*d_ph_midbp*mid_bp_x[i][j] + sin(a2*pi/180)*d_ph_midbp*mid_bp_y[i][j];
			}
	k++;

	for(j=0;j<3;j++){
	pos_ph[k][j] = mid_bp[i][j] + cos(a1*pi/180)*d_ph_midbp*mid_bp_x[i][j] + sin(a1*pi/180)*d_ph_midbp*mid_bp_y[i][j];
			}
	k++;
				}


free(x_axis_bp);
free(y_axis_bp);
free(mid_bp_x);
free(mid_bp_y);

}


void partial_ph_pos(int num_bp, int itot, double d_ph_midbp, double **mid_bp, double **center_bp, double **dir_bp, double ***orient_old, double **dir_ph, double **pos_ph, int real_itt)
{

if(real_itt == 0){ph_pos(num_bp,itot,d_ph_midbp,mid_bp,center_bp,dir_bp,orient_old,dir_ph,pos_ph);}
else{

int i,j,k;
const float pi = 3.14159265358979323846;

double **x_axis_bp, **y_axis_bp, **mid_bp_x, **mid_bp_y;
x_axis_bp = dbl_aloc2d(num_bp,3);
y_axis_bp = dbl_aloc2d(num_bp,3);
mid_bp_x = dbl_aloc2d(num_bp,3);
mid_bp_y = dbl_aloc2d(num_bp,3);


	//x_axis_bp[0][0] = 1; x_axis_bp[0][1] = 0; x_axis_bp[0][2] = 0;
	//y_axis_bp[0][0] = 0; y_axis_bp[0][1] = -1; y_axis_bp[0][2] = 0;	
	
	for(i=real_itt;i<num_bp;i++){
	for(j=0;j<3;j++){
	x_axis_bp[i][j] = orient_old[i-1][0][j];		// normalized
	y_axis_bp[i][j] = -1 * orient_old[i-1][1][j];		// -1 to get right phosphate position
			}
				}
	

	for(i=real_itt;i<num_bp-1;i++){
			mid_bp[i][0] = 0.5*(center_bp[i+1][0] + center_bp[i][0]);
			mid_bp[i][1] = 0.5*(center_bp[i+1][1] + center_bp[i][1]);
			mid_bp[i][2] = 0.5*(center_bp[i+1][2] + center_bp[i][2]);		
				}


	double cos_x,cos_y;	
	double buf_x, buf_y;	

	for(i=real_itt;i<num_bp-1;i++){
	cos_x = dot_product(x_axis_bp[i],x_axis_bp[i+1]);
	cos_y = dot_product(y_axis_bp[i],y_axis_bp[i+1]);

	//printf("%lf %lf %lf %lf \n", cos_x, cos(acos(cos_x)), cos_y, cos(acos(cos_y)));

	

	for(j=0;j<3;j++){
	mid_bp_x[i][j] = 0.5*cos_x*x_axis_bp[i][j] + 0.5*x_axis_bp[i+1][j];
	mid_bp_y[i][j] = 0.5*cos_y*y_axis_bp[i][j] + 0.5*y_axis_bp[i+1][j];
			}
	
	buf_x = vec_len(mid_bp_x[i],3);
	buf_y = vec_len(mid_bp_y[i],3);

	for(j=0;j<3;j++){
	mid_bp_x[i][j] = mid_bp_x[i][j] / buf_x;
	mid_bp_y[i][j] = mid_bp_y[i][j] / buf_y;
			}

				}


	double a1,a2,a_adj;
	a_adj = 12;		// Adjust angle values
	a1 = 33+90-a_adj;	
	a2 = 147+90+a_adj;

	//printf("%lf %lf \n",cos(a1*PI/180),sin(a1*PI/180));

	k=2*(real_itt);
	for(i=real_itt;i<num_bp-1;i++){
	for(j=0;j<3;j++){
	pos_ph[k][j] = mid_bp[i][j] + cos(a2*pi/180)*d_ph_midbp*mid_bp_x[i][j] + sin(a2*pi/180)*d_ph_midbp*mid_bp_y[i][j];
			}
	k++;

	for(j=0;j<3;j++){
	pos_ph[k][j] = mid_bp[i][j] + cos(a1*pi/180)*d_ph_midbp*mid_bp_x[i][j] + sin(a1*pi/180)*d_ph_midbp*mid_bp_y[i][j];
			}
	k++;
				}


free(x_axis_bp);
free(y_axis_bp);
free(mid_bp_x);
free(mid_bp_y);

} //End of else that real_itt > 0
}



void nucl_dist(int num_bp, int num_nucl, int num_fl, int itot, int acc, double **center_bp, double **nucl_pos, int ins_nucl[], double **pos, double **d_nucl, double **d_nucl_bp, double **d_nucl_fl, double **xyz_float, double d, int *get_m)
{
int i,j,k,l,m;
double r;
	
for(i=0;i<3;i++){center_bp[0][i] == 0.0;}
for(i=1;i<num_bp;i++){		
for(j=0;j<3;j++){
		center_bp[i][j] = pos[i-1][j];
		}
		     }


k=0;
l=1;
m = (int) (double)(num_bp-1) / (double)acc;
				//printf("m %d\n",m);
				
	// Nucleosomes with nucleosomes

for(i=0;i<num_nucl;i++){
	for(j=i+1;j<num_nucl;j++){
				d_nucl[i][j] = dist(nucl_pos[i],nucl_pos[j]);
				//if(d_nucl[i][j] < 310){printf("d_nucl %lf %d %d\n", d_nucl[i][j], i,j);}
			        }



	// Nucleosomes with bp's




				r = ((double) rand() / (RAND_MAX));
				k = (int) (r*(double)acc);
				k=0;
				//printf("k %d\n",k);
	for(j=0;j<m;j++){  //num_bp-acc-i 
				
				if((k+j*acc) < num_bp && (k+j*acc < (int)((ins_nucl[i]-1+(i*147))-d) || k+j*acc > (int)((ins_nucl[i]-1+(i*147))+147+d)) ) 					// Exclude DNA of nucleosome and d bp before and after nucleosome
						     {
					d_nucl_bp[i][l-1] = dist(nucl_pos[i],center_bp[k+j*acc]);
					//printf("d_nucl_bp %lf %d %d %d\n",d_nucl_bp[i][l-1],i,k-1,k+j*acc);
					//if(d_nucl_bp[i][l-1] < 65){printf("d_nucl_bp %lf %d %d %d\n",d_nucl_bp[i][l-1],i,k-1,k+j*acc);}
					l++;
						     }
			     }
				l=1;
				
	
	// 

	// Nucleosomes with floating particles


	for(j=0;j<num_fl;j++){
				d_nucl_fl[i][j] = dist(nucl_pos[i],xyz_float[j]);
			     }
			}


	*get_m = m;

}


void check_olap(int num_bp, int num_fl, int itot, int acc, double **center_bp, double **pos, double **d_bp, double **d_bp_fl, double **xyz_float, double bp_excl_vol, double bp_fl_excl_vol, int no_nucl, double **d_nucl, double **d_nucl_bp, double **d_nucl_fl, double nucl_nucl_excl_vol, double nucl_bp_excl_vol, double nucl_fl_excl_vol, int red_bp_all, int *vol_olap)
{
int i,j,k;

	/* ---- Calculate position and distances of bp's and floating particles ----*/
	
	positions(num_bp,num_fl,itot,acc,center_bp,pos,d_bp,d_bp_fl,xyz_float,red_bp_all);	

	/* ---- Condition to start calculating LJ and DH potentials ---- */

	int overlap;
	overlap = 0;

	//double bp_excl_vol = 10;
	//double bp_fl_excl_vol = 5;	

	// Overlap bp/bp and bp/float
	for(i=0;i<red_bp_all;i++){
	for(j=i+2;j<red_bp_all;j++){	// i+2 to allow DNA to explore full elasticity (i+1 would restrict: with acc=7 distance between two DNA would be 20nm while average 7*3nm = 21nm)
				if(d_bp[i][j] > 0.1 && d_bp[i][j] < bp_excl_vol){overlap = 1; 
										//printf("i %d j %d d_bp %lf \n",i,j,d_bp[i][j]);
										}
				if(d_bp_fl[i][j] > 0.1 && d_bp_fl[i][j] < bp_fl_excl_vol){overlap = 1; 
											//printf("i %d j %d d_bp_fl %lf \n",i,j,d_bp_fl[i][j]);
											}	// be careful when I introduce floating particles
			      }
			      }

	// Overlap nucl/bp and nucl/float
	for(i=0;i<no_nucl;i++){
	for(j=0;j<red_bp_all;j++){
				if(d_nucl_bp[i][j] > 0.1 && d_nucl_bp[i][j] < nucl_bp_excl_vol){overlap = 1; 
												//printf("i %d j %d d_nucl_bp %lf \n",i,j,d_nucl_bp[i][j]);
												}
				if(d_nucl_fl[i][j] > 0.1 && d_nucl_fl[i][j] < nucl_fl_excl_vol){overlap = 1; 
												//printf("i %d j %d d_nucl_fl %lf \n",i,j,d_nucl_fl[i][j]);
												}
			      }
			      }


	// Overlap nucl/nucl
	for(i=0;i<no_nucl;i++){
	for(j=0;j<no_nucl;j++){
				if(d_nucl[i][j] > 0.1 && d_nucl[i][j] < nucl_nucl_excl_vol){overlap = 1; 
												//printf("i %d j %d d_nucl %lf \n",i,j,d_nucl[i][j]);
												}
			      }
			      }



	*vol_olap = overlap;

}

/*void transfer_hel_coord(int itot, double **xconf, double **xconfi)
{
int i,j;

for(i=0; i<6; i++){
for(j=0; j<itot; j++){
			xconf[i][j] = xconfi[i][j];  
		      }	
		  }

}*/

void initial_fl_pos(int num_fl, double shift, double br_mot[3], double **xyz_float)
{
int i;
double rx,ry,rz;

for(i=0;i<num_fl;i++){
			rx = ((double) rand() / (RAND_MAX)); 
			ry = ((double) rand() / (RAND_MAX)); 
			rz = ((double) rand() / (RAND_MAX)); 			
			xyz_float[i][0] = br_mot[0] + shift*(rx - 0.5);	// Position particles between -200 and 200
			xyz_float[i][1] = br_mot[1] + shift*(ry - 0.5);
			xyz_float[i][2] = br_mot[2] + shift*(rz - 0.5);
		}

}


void dist_ph(int num_bp, int num_fl, int acc_ph, double **pos_ph, double **d_ph1, double **d_ph2, double **d_ph_fl, double **xyz_float, int num_nucl, int ins_nucl[], double **d_nucl_ph, double **nucl_pos, int ***in_q_ph_ph, int **in_q_nucl_ph, double d, int *get_m, int red_ph)
{
int k,ki,kj,i,j,l,m,num_ph;
double r;
ki=0;
kj=0;
num_ph = 2*num_bp-2;
//m = (int) (double)(num_ph) / (double)acc_ph;
//printf("num ph %d red_ph %d\n",num_ph,red_ph);

	// Phosphates with Phosphates
	for(i=0;i<red_ph;i++){
	
	for(j=i+1;j<red_ph;j++){ //num_bp-acc-i 
				if((j*acc_ph+1) < num_ph){

					//printf("kj %d\n",kj);
					d_ph1[ki][kj] = dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph]);	
					in_q_ph_ph[ki][kj][0] = i*acc_ph;
					in_q_ph_ph[ki][kj][1] = j*acc_ph;
					d_ph1[ki+1][kj] = dist(pos_ph[i*acc_ph+1],pos_ph[j*acc_ph]);
					in_q_ph_ph[ki+1][kj][0] = i*acc_ph+1;
					in_q_ph_ph[ki+1][kj][1] = j*acc_ph;
					d_ph1[ki][kj+1] = dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph+1]);
					in_q_ph_ph[ki][kj+1][0] = i*acc_ph;
					in_q_ph_ph[ki][kj+1][1] = j*acc_ph+1;
					d_ph1[ki+1][kj+1] = dist(pos_ph[i*acc_ph+1],pos_ph[j*acc_ph+1]);
					in_q_ph_ph[ki+1][kj+1][0] = i*acc_ph+1;
					in_q_ph_ph[ki+1][kj+1][1] = j*acc_ph+1;
					//if(d_ph1[i][j] < 20){printf("%d %d %d %lf\n",num_bp,i*acc_ph,j*acc_ph,d_ph1[i][j]);}
					kj=kj+2;
				
						     }
			     }
		
	//printf("ki %d\n",ki);
	ki=ki+2;
	kj=ki;		
	


	// Phosphates with Phosphates (take minimum of ph_distance)
	/*for(i=0;i<red_ph;i++){
	for(j=i+1;j<red_ph;j++){ //num_bp-acc-i 
				if(j*acc_ph < 2*num_bp-2){
					d_ph1[i][j] = min(dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph]),dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph+1]));
					if(dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph]) < dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph+1])){in_q_ph_ph[i][k-1] = j*acc_ph;}
						else{in_q_ph_ph[i][j] = j*acc_ph+1;}	// INDEX FUNKTIONIERT NOCH NICHT!!!!
k++;
				
						     }
			     }
				k=1;*/


	// Phosphates with Phosphates (OLD WAY)
	/*for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){ //num_bp-acc-i 
				if(i+1+k*acc_ph < 2*num_bp-2){
					d_ph1[i][k-1] = min(dist(pos_ph[i],pos_ph[i+k*acc_ph]),dist(pos_ph[i],pos_ph[i+1+k*acc_ph]));
					if(dist(pos_ph[i],pos_ph[i+k*acc_ph]) < dist(pos_ph[i],pos_ph[i+1+k*acc_ph])){in_q_ph_ph[i][k-1] = i+k*acc_ph;}
						else{in_q_ph_ph[i][k-1] = i+1+k*acc_ph;}
					//d_ph1[i][k-1] = dist(pos_ph[i],pos_ph[i+k*acc_ph]);
					d_ph2[i][k-1] = dist(pos_ph[i],pos_ph[i+1+k*acc_ph]);
					//printf("%d %d %d %lf\n",num_bp,i,i+1+k*acc_ph,pos_ph[i+1+k*acc_ph][0]);
					k++;
				
						     }
			     }
				k=1;*/
				


	// Phosphates with floating particle

	for(j=0;j<num_fl;j++){
				d_ph_fl[i][j] = dist(pos_ph[i],xyz_float[j]);
			     }

				} //End loop



	// Nucleosome with Phosphates
	r = ((double) rand() / (RAND_MAX));
	k = (int) (r*(double)acc_ph);
	k=0;
	//printf("k %d\n",k);
	//printf("m %d\n",m);
	l=1;


	for(i=0;i<num_nucl;i++){
	for(j=0;j<red_ph;j++){  //num_bp-acc-i 
				
				if((j*acc_ph+1) < num_ph && (j*acc_ph < (int)2*((ins_nucl[i]-1+(i*147))-d) || j*acc_ph > (int)2*((ins_nucl[i]-1+(i*147))+147+d)) ) 					// Exclude phosphates of nucleosome and d bp before and after nucleosome
						     {
					d_nucl_ph[i][l] = dist(nucl_pos[i],pos_ph[j*acc_ph]);
					in_q_nucl_ph[i][l] = j*acc_ph;
					d_nucl_ph[i][l+1] = dist(nucl_pos[i],pos_ph[1+j*acc_ph]);
					in_q_nucl_ph[i][l+1] = 1+j*acc_ph;
					//printf("d_nucl_ph %lf %d %d %d\n",d_nucl_ph[i][l-1],i,l-1,k+j*acc_ph);
					//if(d_nucl_ph[i][l-1] < 65){printf("d_nucl_ph %lf %d %d %d\n",d_nucl_ph[i][l-1],i,l-1,k+j*acc_ph);}
					l=l+2;
						     }
			     }
				l=0;
	
				}


// OLD WAY phosphate-nucleosome
/*	for(i=0;i<num_nucl;i++){
	for(j=0;j<red_ph;j++){  //num_bp-acc-i 
				
				if((k+j*acc_ph) < num_ph && (k+j*acc_ph < (int)2*((ins_nucl[i]-1+(i*147))-d) || k+j*acc_ph > (int)2*((ins_nucl[i]-1+(i*147))+147+d)) ) 					// Exclude phosphates of nucleosome and d bp before and after nucleosome
						     {
					d_nucl_ph[i][l-1] = min(dist(nucl_pos[i],pos_ph[k+j*acc_ph]),dist(nucl_pos[i],pos_ph[k+1+j*acc_ph]));
					if(dist(nucl_pos[i],pos_ph[k+j*acc_ph]) < dist(nucl_pos[i],pos_ph[k+1+j*acc_ph])){in_q_nucl_ph[i][l-1] = k+j*acc_ph;}
						else{in_q_nucl_ph[i][l-1] = k+1+j*acc_ph;}
					//printf("d_nucl_ph %lf %d %d %d\n",d_nucl_ph[i][l-1],i,l-1,k+j*acc_ph);
					//if(d_nucl_ph[i][l-1] < 65){printf("d_nucl_ph %lf %d %d %d\n",d_nucl_ph[i][l-1],i,l-1,k+j*acc_ph);}
					l++;
						     }
			     }
				l=1;
	
				}*/
				
				

	*get_m = m;

	/*for(i=500;i<510;i++){	
	printf("%lf %lf %lf %lf\n",d_ph_fl[i][0],d_ph_fl[i][1],d_ph_fl[i][2],d_ph_fl[i][3]);
	}*/

}


void calc_lj(int *k1, int *k2, int num_bp, double **d_bp, double lj_cut_low, double lj_cut_up, double **lj_bp, double lj_k, double lj_a_bp, int acc, int num_fl, double **lj_bp_fl, double **d_bp_fl, double lj_a_fl, int red_bp_all, double *e_lj_bp_bpfl)
{
int i,j,k1_buf,k2_buf;
double buf_energy;
buf_energy = 0;
	// Calculate inter-bp LJ pot (for lj_cut_low < r < lj_cut_up)


	k1_buf=0;
	k2_buf=0;
	for(i=0;i<red_bp_all;i++){
	for(j=i+1;j<red_bp_all;j++){		
				if(d_bp[i][j] > lj_cut_low && d_bp[i][j] < lj_cut_up){ //printf("yes1");
						lj_bp[i][j] = lj(d_bp[i][j],lj_k,lj_a_bp,lj_a_bp);
						buf_energy = buf_energy + lj_bp[i][j];
						k1_buf++;
					/*	lj_bp[k1_buf][0] = lj(d_bp[i][j],lj_k,lj_a_bp,lj_a_bp);
						lj_bp[k1_buf][1] = i;
						lj_bp[k1_buf][2] = i+j*acc;
						lj_bp[k1_buf][3] = d_bp[i][j];
						k1_buf++;   */
												
						
						    }
				  }
			   

	

	
	// Calculate bp - floating patricle LJ pot (for lj_cut_low < r < lj_cut_up)
	
	
	for(j=0;j<num_fl;j++){		
				if(d_bp_fl[i][j] > lj_cut_low && d_bp_fl[i][j] < lj_cut_up){ //printf("yes2");
						lj_bp_fl[i][j] = lj(d_bp_fl[i][j],lj_k,lj_a_fl,lj_a_fl);
						buf_energy = buf_energy + lj_bp_fl[i][j];
						k2_buf++;
						/*lj_bp_fl[k2_buf][0] = lj(d_bp_fl[i][j],lj_k,lj_a_fl,lj_a_fl);
						lj_bp_fl[k2_buf][1] = i;
						lj_bp_fl[k2_buf][2] = j;
						lj_bp_fl[k2_buf][3] = d_bp_fl[i][j];
						k2_buf++;	*/
						    }
				  }
			     }
	
	*k1 = k1_buf;
	*k2 = k2_buf;

	*e_lj_bp_bpfl = buf_energy;

}

void calc_dh(int *k3, int *k4, int num_bp, double **dh_ph, double **d_ph1, double **d_ph2, double dh_cut_low, double dh_cut_up, double *q_ph, double eps, double kappa, double f_pi_eps0, int acc_ph, int num_fl, double **dh_ph_fl, double **d_ph_fl, double *q_fl, int ***in_q_ph_ph, int red_ph, double *e_dh_ph_phfl)
{
int i,j,m,num_ph,k3_buf,k4_buf,k;
double buf_energy;
buf_energy = 0;
k=0;
//num_ph = 2*num_bp-2;
//m = (int) (double)(num_ph) / (double)acc_ph;


	// Calculate inter-ph DH pot for dh_cut_low < r < dh_cut_up


	k3_buf=0;
	k4_buf=0;

	for(i=0;i<2*red_ph;i++){
	for(j=0;j<2*red_ph;j++){		
				if(d_ph1[i][j] > dh_cut_low && d_ph1[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[i][j] = dh(d_ph1[i][j],q_ph[in_q_ph_ph[i][j][0]],q_ph[in_q_ph_ph[i][j][1]],eps,kappa,f_pi_eps0);  //in ev
						buf_energy = buf_energy + dh_ph[i][j];
						k3_buf++;
						    }
				  }
			

	// With minimum of ph_distance
	/*for(i=0;i<red_ph;i++){
	for(j=i+1;j<red_ph;j++){		
				if(d_ph1[i][j] > dh_cut_low && d_ph1[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[i][j] = dh(d_ph1[i][j],q_ph[i*acc_ph],q_ph[in_q_ph_ph[i][j]],eps,kappa,f_pi_eps0);  //in ev
						buf_energy = buf_energy + dh_ph[i][j];
						k3_buf++;

						    }
				  }*/
			     

	/*for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){		
				if(d_ph2[i][j] > dh_cut_low && d_ph2[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[k3_buf][0] = dh(d_ph2[i][j],q_ph,q_ph,eps,kappa,f_pi_eps0);
						dh_ph[k3_buf][1] = i;
						dh_ph[k3_buf][2] = i+1+j*acc_ph;
						dh_ph[k3_buf][3] = d_ph2[i][j];
						k3_buf++;
						    }
				  }
			     }*/
		
	

	// Calculate bp - floating patricle LJ pot for lj_cut_low < r < lj_cut_up 
	
	//for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<num_fl;j++){		
				if(d_ph_fl[i][j] > dh_cut_low && d_ph_fl[i][j] < dh_cut_up){ //printf("yes4");
						dh_ph_fl[i][j] = dh(d_ph_fl[i][j],q_ph[i],q_fl[j],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_ph_fl[i][j];
						k4_buf++;
						/*dh_ph_fl[k4_buf][0] = dh(d_ph_fl[i][j],q_ph[i],q_fl[j],eps,kappa,f_pi_eps0);
						dh_ph_fl[k4_buf][1] = i;
						dh_ph_fl[k4_buf][2] = j;
						dh_ph_fl[k4_buf][3] = d_ph_fl[i][j];
						k4_buf++;*/
						    }
				  }
			     }

	*k3 = k3_buf;
	*k4 = k4_buf;

	*e_dh_ph_phfl = buf_energy;
}



void nucl_positions(double **nucl_pos, double ***orient_old, double **pos_old, double *center_nucl, int ins_nucl[], int no_nucl)
{
int i,j,k;
double *buf_vec;
buf_vec = malloc(3*sizeof(double));

for(i=0;i<no_nucl;i++){
		mat_vec_mult_n(buf_vec,orient_old[ins_nucl[i]-1+(i*147)],center_nucl,3);
		//for(i=0;i<3;i++){printf("%lf ",buf_vec[i]);} printf("\n");
		for(j=0;j<3;j++){nucl_pos[i][j] = pos_old[ins_nucl[i]-1+(i*147)][j] + buf_vec[j];}
		      }

free(buf_vec);
}


void sep_nucl_dna(double **linker_pos, double ***linker_orient, double **nucl_dna_pos, double ***nucl_dna_orient, double ***orient_old, double **pos_old, int ins_nucl[], int no_nucl, int len_all)
{
int h,i,j,k,start;
int count_link,count_nucl;
i=0;
count_link = 0;
count_nucl = 0;



for(j=0;j<ins_nucl[0];j++){	// linker DNA until first bp of nucleosome
		copy_2darray(linker_orient[j],orient_old[j],3);
		copy_1darray(linker_pos[j],pos_old[j],3);
		i++; 
		count_link++;
			  }

//printf("link-nucl i %d\n",i);

for(j=0;j<146;j++){	// first nucleosome
		copy_2darray(nucl_dna_orient[j],orient_old[i],3);
		copy_1darray(nucl_dna_pos[j],pos_old[i],3);
		i++; 
		count_nucl++;
			  }

//printf("nucl-link i %d\n",i);

for(h=1;h<no_nucl;h++){
		//printf("h %d\n",h);
		start = count_link;
		//printf("start %d\n",start);
		for(j=i;j<(ins_nucl[h]+(h*147));j++){	// linker DNA
		copy_2darray(linker_orient[count_link],orient_old[i],3);
		copy_1darray(linker_pos[count_link],pos_old[i],3);
					i++;
					count_link++;
		
		//printf("link i %d\n",i);
					   }
		//printf("link-nucl i %d\n",i);
		//printf("first_bp_nucl %d\n",ins_nucl[h]+(h*147));

		start = count_nucl;
		for(j=0;j<146;j++){	// nucleosomal DNA
		copy_2darray(nucl_dna_orient[count_nucl],orient_old[i],3);
		copy_1darray(nucl_dna_pos[count_nucl],pos_old[i],3);
				i++;	
				count_nucl++;
		//printf("nucl i %d\n",i);
				   }
		
		//printf("nucl-link i %d\n",i);
		//for(i=0;i<3;i++){printf("%lf ",buf_vec[i]);} printf("\n");
		//for(j=0;j<3;j++){nucl_pos[i][j] = pos_old[ins_nucl[i]-1+(i*147)][j] + buf_vec[j];}
		      }

				start = count_link;
				//printf("start %d\n",start);
	for(j=i;j<len_all;j++){copy_1darray(linker_pos[count_link],pos_old[j],3);
			       copy_2darray(linker_orient[count_link],orient_old[j],3);
			       count_link++;
			       i++;
				//printf("link i %d\n",j);
			     }

				//printf("link-end i %d\n",i);				
		
//printf("count_link %d\n", count_link);
//printf("count_nucl %d\n", count_nucl);


}


void index_nucl_dna(int *linker_nucl_index, int ins_nucl[], int no_nucl, int len_all)
{
int h,i,j,k,start;
int count_link,count_nucl;
i=0;
count_link = 0;
count_nucl = 0;
int len_nucl = 146;

linker_nucl_index[0] = 0;
linker_nucl_index[1] = ins_nucl[0];
linker_nucl_index[2] = ins_nucl[0]+len_nucl;

i=3;
for(h=1;h<no_nucl;h++){
		linker_nucl_index[i] = ins_nucl[h]+(h*147);		
		linker_nucl_index[i+1] = linker_nucl_index[i]+len_nucl;	
		i=i+2;
		      }

linker_nucl_index[i] = len_all;

}



void sep_nucl_linker_ph(double **linker_ph, double **nucl_ph, double **pos_ph, int ins_nucl[], int no_nucl, int len_all)
{
int h,i,j,k,start,m;
int count_link,count_nucl;
i=0;
count_link = 0;
count_nucl = 0;
m=2;


for(j=0;j<m*ins_nucl[0];j++){	// linker DNA until first bp of nucleosome	
		copy_1darray(linker_ph[j],pos_ph[j],3);
		i++; 
		count_link++;
			  }


for(j=0;j<m*146;j++){	// first nucleosome
		copy_1darray(nucl_ph[j],pos_ph[i],3);
		i++; 
		count_nucl++;
			  }


for(h=1;h<no_nucl;h++){
		//printf("h %d\n",h);
		start = count_link;
		//printf("start %d\n",start);
		for(j=i;j<m*(ins_nucl[h]+(h*147));j++){	// linker DNA
		copy_1darray(linker_ph[count_link],pos_ph[i],3);
					i++;
					count_link++;
		
		//printf("link i %d\n",i);
					   }
		
		//printf("first_bp_nucl %d\n",ins_nucl[h]+(h*147));

		start = count_nucl;
		for(j=0;j<m*146;j++){	// nucleosomal DNA
		copy_1darray(nucl_ph[count_nucl],pos_ph[i],3);
				i++;	
				count_nucl++;
		//printf("nucl i %d\n",i);
				   }
		
		//for(i=0;i<3;i++){printf("%lf ",buf_vec[i]);} printf("\n");
		//for(j=0;j<3;j++){nucl_pos[i][j] = pos_old[ins_nucl[i]-1+(i*147)][j] + buf_vec[j];}
		      }

				start = count_link;
				//printf("start %d\n",start);
	for(j=i;j<m*len_all;j++){copy_1darray(linker_ph[count_link],pos_ph[j],3);	   
			       count_link++;
			       
				//printf("link i %d\n",j);
			     }

//printf("count_link %d\n", count_link);
//printf("count_nucl %d\n", count_nucl);


}


void charges_ph(double *arr_q_ph, double q_link, double q_nucl_dna, int *ph_link_nucl_index, int ins_nucl[], int no_nucl, int len_all)
{
int i,j,k;
int m = 2;
int a = 2*no_nucl+2;
for(j=0;j<(a-1);j++){
		for(i=ph_link_nucl_index[j];i<ph_link_nucl_index[j+1];i++){
							if(j % 2 == 0){arr_q_ph[i] = q_link;}
								else{arr_q_ph[i] = q_nucl_dna;}
							
									}
		//printf("i %d\n",i);
		     }



}



void calc_lj_nucl(int *k1, int *k2, int *k3, int num_bp, int num_nucl, double **d_nucl, double **d_nucl_bp, double **d_nucl_fl, double lj_nucl_cut_low, double lj_nucl_cut_up, double **lj_nucl, double **lj_nucl_bp, double **lj_nucl_fl, double lj_k_nucl, double lj_k_nucl_bp, double lj_k_nucl_fl, double lj_a_nucl, double lj_a_bp_nucl, int acc, int num_fl, double lj_a_nucl_fl, double *energy)
{
int i,j,k1_buf,k2_buf,k3_buf;
double buf_energy, nucl_nucl_e, nucl_ph_e;
buf_energy = 0;
nucl_nucl_e = 0;
nucl_ph_e = 0;

	// Calculate inter-nucl LJ pot (for lj_cut_low < r < lj_cut_up)


	k1_buf=0;
	k2_buf=0;
	k3_buf=0;
	for(i=0;i<num_nucl;i++){

	for(j=i;j<num_nucl;j++){	//printf("%lf\n",d_nucl[i][j]);
				if(d_nucl[i][j] > lj_nucl_cut_low && d_nucl[i][j] < lj_nucl_cut_up){ //printf("yes\n");
						lj_nucl[i][j] = lj(d_nucl[i][j],lj_k_nucl,lj_a_nucl,lj_a_nucl);
						buf_energy = buf_energy + lj_nucl[i][j];
						nucl_nucl_e = nucl_nucl_e + lj_nucl[i][j];
						k1_buf++;
						/*lj_nucl[k1_buf][0] = lj(d_nucl[i][j],lj_k_nucl,lj_a_nucl,lj_a_nucl);
						lj_nucl[k1_buf][1] = i;
						lj_nucl[k1_buf][2] = j;
						lj_nucl[k1_buf][3] = d_nucl[i][j];
						buf_energy = buf_energy + lj_nucl[k1_buf][0];
						nucl_nucl_e = nucl_nucl_e + lj_nucl[k1_buf][0];
						k1_buf++;*/
						    }
				  }



	// Calculate nucl-bp LJ pot (for lj_cut_low < r < lj_cut_up)


	
	for(j=0;j<num_bp;j++){		//printf("%lf\n",d_nucl_bp[i][j]);
				if(d_nucl_bp[i][j] > lj_nucl_cut_low && d_nucl_bp[i][j] < lj_nucl_cut_up){ //printf("yes1 j %d\n",j);
						lj_nucl_bp[i][j] = lj(d_nucl_bp[i][j],lj_k_nucl_bp,lj_a_bp_nucl,lj_a_bp_nucl);
						buf_energy = buf_energy + lj_nucl_bp[i][j];
						nucl_ph_e = nucl_ph_e + lj_nucl_bp[i][j];
						k2_buf++;
						/*lj_nucl_bp[k2_buf][0] = lj(d_nucl_bp[i][j],lj_k_nucl_bp,lj_a_bp_nucl,lj_a_bp_nucl);
						lj_nucl_bp[k2_buf][1] = i;
						lj_nucl_bp[k2_buf][2] = i+j*acc;
						lj_nucl_bp[k2_buf][3] = d_nucl_bp[i][j];
						buf_energy = buf_energy + lj_nucl_bp[k2_buf][0];
						nucl_ph_e = nucl_ph_e + lj_nucl_bp[k2_buf][0];
						k2_buf++;*/
						    }
				  }


	
	// Calculate nucl - floating patricle LJ pot (for lj_cut_low < r < lj_cut_up)
	

	for(j=0;j<num_fl;j++){		
				if(d_nucl_fl[i][j] > lj_nucl_cut_low && d_nucl_fl[i][j] < lj_nucl_cut_up){ //printf("yes2");
						lj_nucl_fl[i][j] = lj(d_nucl_fl[i][j],lj_k_nucl_fl,lj_a_nucl_fl,lj_a_nucl_fl);
						buf_energy = buf_energy + lj_nucl_fl[i][j];
						k3_buf++;
						/*lj_nucl_fl[k1_buf][0] = lj(d_nucl_fl[i][j],lj_k_nucl_fl,lj_a_nucl_fl,lj_a_nucl_fl);
						lj_nucl_fl[k1_buf][1] = i;
						lj_nucl_fl[k1_buf][2] = j;
						lj_nucl_fl[k1_buf][3] = d_nucl_fl[i][j];
						buf_energy = buf_energy + lj_nucl_fl[k3_buf][0];
						k3_buf++;*/
						    }
				  }

			     }
	
	*k1 = k1_buf;
	*k2 = k2_buf;
	*k3 = k3_buf;

	*energy = buf_energy;

//printf("LJ nucl_nucl_e %lf\n",nucl_nucl_e);
//printf("LJ nucl_bp_e %lf\n",nucl_ph_e);
}


void calc_dh_nucl(int *k1, int *k2, int *k3, int num_bp, int num_nucl, double **d_nucl, double **d_nucl_ph, double **d_nucl_fl, double dh_nucl_cut_low, double dh_nucl_cut_up, double **dh_nucl, double **dh_nucl_ph, double **dh_nucl_fl, double *q_nucl, double *q_ph, double *q_fl, int **in_q_nucl_ph, int acc_ph, int num_fl, double eps, double kappa, double f_pi_eps0, int red_ph, double *energy)
{
int i,j,k1_buf,k2_buf,k3_buf;
double buf_energy, nucl_nucl_e, nucl_ph_e;
buf_energy = 0;
nucl_nucl_e = 0;
nucl_ph_e = 0;

	// Calculate inter-nucl dh pot (for dh_cut_low < r < dh_cut_up)


	k1_buf=0;
	k2_buf=0;
	k3_buf=0;
	for(i=0;i<num_nucl;i++){

	for(j=i;j<num_nucl;j++){	//printf("%lf\n",d_nucl[i][j]);
				if(d_nucl[i][j] > dh_nucl_cut_low && d_nucl[i][j] < dh_nucl_cut_up){ //printf("yes\n");
						dh_nucl[i][j] = dh(d_nucl[i][j],q_nucl[i],q_nucl[j],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_nucl[i][j];
						nucl_nucl_e = nucl_nucl_e + dh_nucl[i][j];
						k1_buf++;
						/*dh_nucl[k1_buf][0] = dh(d_nucl[i][j],q_nucl[i],q_nucl[j],eps,kappa,f_pi_eps0);
						dh_nucl[k1_buf][1] = i;
						dh_nucl[k1_buf][2] = j;
						dh_nucl[k1_buf][3] = d_nucl[i][j];
						buf_energy = buf_energy + dh_nucl[k1_buf][0];
						nucl_nucl_e = nucl_nucl_e + dh_nucl[k1_buf][0];
						k1_buf++;*/
						    }
				  }



	// Calculate nucl-bp dh pot (for dh_cut_low < r < dh_cut_up)


	
	for(j=0;j<2*red_ph;j++){	//printf("%lf\n",d_nucl_bp[i][j]);
				if(d_nucl_ph[i][j] > dh_nucl_cut_low && d_nucl_ph[i][j] < dh_nucl_cut_up){ //printf("yes1 j %d\n",j);
						dh_nucl_ph[i][j] = dh(d_nucl_ph[i][j],q_nucl[i],q_ph[in_q_nucl_ph[i][j]],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_nucl_ph[i][j];
						nucl_ph_e = nucl_ph_e + dh_nucl_ph[i][j];
						k2_buf++;
						/*dh_nucl_ph[k2_buf][0] = dh(d_nucl_ph[i][j],q_nucl[i],q_ph[in_q_nucl_ph[i][j]],eps,kappa,f_pi_eps0);
						dh_nucl_ph[k2_buf][1] = i;
						dh_nucl_ph[k2_buf][2] = i+j*acc_ph;
						dh_nucl_ph[k2_buf][3] = d_nucl_ph[i][j];
						buf_energy = buf_energy + dh_nucl_ph[k2_buf][0];
						nucl_ph_e = nucl_ph_e + dh_nucl_ph[k2_buf][0];
						k2_buf++;*/
						    }
				  }


	
	// Calculate nucl - floating patricle dh pot (for dh_cut_low < r < dh_cut_up)
	

	for(j=0;j<num_fl;j++){		
				if(d_nucl_fl[i][j] > dh_nucl_cut_low && d_nucl_fl[i][j] < dh_nucl_cut_up){ //printf("yes2");
						dh_nucl_fl[i][j] = dh(d_nucl_fl[i][j],q_nucl[i],q_fl[j],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_nucl_fl[i][j];
						k3_buf++;
						/*dh_nucl_fl[k1_buf][0] = dh(d_nucl_fl[i][j],q_nucl[i],q_fl[j],eps,kappa,f_pi_eps0);
						dh_nucl_fl[k1_buf][1] = i;
						dh_nucl_fl[k1_buf][2] = j;
						dh_nucl_fl[k1_buf][3] = d_nucl_fl[i][j];
						buf_energy = buf_energy + dh_nucl_fl[k3_buf][0];
						k3_buf++;*/
						    }
				  }

			     }
	
	*k1 = k1_buf;
	*k2 = k2_buf;
	*k3 = k3_buf;

	*energy = buf_energy;

//printf("DH nucl_nucl_e %lf\n",nucl_nucl_e);
//printf("DH nucl_ph_e %lf\n",nucl_ph_e);
}




/* Functions to calculate with minimum of distance of phosphates */
void dist_ph_MIN(int num_bp, int num_fl, int acc_ph, double **pos_ph, double **d_ph1, double **d_ph2, double **d_ph_fl, double **xyz_float, int num_nucl, int ins_nucl[], double **d_nucl_ph, double **nucl_pos, int **in_q_ph_ph, int **in_q_nucl_ph, double d, int *get_m, int red_ph)
{
int k,i,j,l,m,num_ph;
double r;
k=1;
num_ph = 2*num_bp-2;
//m = (int) (double)(num_ph) / (double)acc_ph;
//printf("m %d\n",m);

	// Phosphates with Phosphates
	for(i=0;i<red_ph;i++){
	for(j=i+1;j<red_ph;j++){ //num_bp-acc-i 
				if(j*acc_ph < 2*num_bp-2){
					d_ph1[i][j] = min(dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph]),dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph+1]));
					if(dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph]) < dist(pos_ph[i*acc_ph],pos_ph[j*acc_ph+1])){in_q_ph_ph[i][k-1] = j*acc_ph;}
						else{in_q_ph_ph[i][j] = j*acc_ph+1;}	// INDEX FUNKTIONIERT NOCH NICHT!!!!
					//d_ph1[i][k-1] = dist(pos_ph[i],pos_ph[i+k*acc_ph]);
					//d_ph2[i][k-1] = dist(pos_ph[i],pos_ph[i+1+k*acc_ph]);
					//if(d_ph1[i][j] < 20){printf("%d %d %d %lf\n",num_bp,i*acc_ph,j*acc_ph,d_ph1[i][j]);}
					k++;
				
						     }
			     }
				k=1;

	// Phosphates with Phosphates (OLD WAY)
	/*for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){ //num_bp-acc-i 
				if(i+1+k*acc_ph < 2*num_bp-2){
					d_ph1[i][k-1] = min(dist(pos_ph[i],pos_ph[i+k*acc_ph]),dist(pos_ph[i],pos_ph[i+1+k*acc_ph]));
					if(dist(pos_ph[i],pos_ph[i+k*acc_ph]) < dist(pos_ph[i],pos_ph[i+1+k*acc_ph])){in_q_ph_ph[i][k-1] = i+k*acc_ph;}
						else{in_q_ph_ph[i][k-1] = i+1+k*acc_ph;}
					//d_ph1[i][k-1] = dist(pos_ph[i],pos_ph[i+k*acc_ph]);
					d_ph2[i][k-1] = dist(pos_ph[i],pos_ph[i+1+k*acc_ph]);
					//printf("%d %d %d %lf\n",num_bp,i,i+1+k*acc_ph,pos_ph[i+1+k*acc_ph][0]);
					k++;
				
						     }
			     }
				k=1;*/
				


	// Phosphates with floating particle

	for(j=0;j<num_fl;j++){
				d_ph_fl[i][j] = dist(pos_ph[i],xyz_float[j]);
			     }

				} //End loop



	// Nucleosome with Phosphates
	r = ((double) rand() / (RAND_MAX));
	k = (int) (r*(double)acc_ph);
	k=0;
	//printf("k %d\n",k);
	//printf("m %d\n",m);
	l=1;

	for(i=0;i<num_nucl;i++){
	for(j=0;j<red_ph;j++){  //num_bp-acc-i 
				
				if((k+j*acc_ph) < num_ph && (k+j*acc_ph < (int)2*((ins_nucl[i]-1+(i*147))-d) || k+j*acc_ph > (int)2*((ins_nucl[i]-1+(i*147))+147+d)) ) 					// Exclude phosphates of nucleosome and d bp before and after nucleosome
						     {
					d_nucl_ph[i][l-1] = min(dist(nucl_pos[i],pos_ph[k+j*acc_ph]),dist(nucl_pos[i],pos_ph[k+1+j*acc_ph]));
					if(dist(nucl_pos[i],pos_ph[k+j*acc_ph]) < dist(nucl_pos[i],pos_ph[k+1+j*acc_ph])){in_q_nucl_ph[i][l-1] = k+j*acc_ph;}
						else{in_q_nucl_ph[i][l-1] = k+1+j*acc_ph;}
					//printf("d_nucl_ph %lf %d %d %d\n",d_nucl_ph[i][l-1],i,l-1,k+j*acc_ph);
					//if(d_nucl_ph[i][l-1] < 65){printf("d_nucl_ph %lf %d %d %d\n",d_nucl_ph[i][l-1],i,l-1,k+j*acc_ph);}
					l++;
						     }
			     }
				l=1;
	
				}
				
				

	*get_m = m;

	/*for(i=500;i<510;i++){	
	printf("%lf %lf %lf %lf\n",d_ph_fl[i][0],d_ph_fl[i][1],d_ph_fl[i][2],d_ph_fl[i][3]);
	}*/

}




void calc_dh_MIN(int *k3, int *k4, int num_bp, double **dh_ph, double **d_ph1, double **d_ph2, double dh_cut_low, double dh_cut_up, double *q_ph, double eps, double kappa, double f_pi_eps0, int acc_ph, int num_fl, double **dh_ph_fl, double **d_ph_fl, double *q_fl, int **in_q_ph_ph, int red_ph, double *e_dh_ph_phfl)
{
int i,j,m,num_ph,k3_buf,k4_buf;
double buf_energy;
buf_energy = 0;
//num_ph = 2*num_bp-2;
//m = (int) (double)(num_ph) / (double)acc_ph;


	// Calculate inter-ph DH pot for dh_cut_low < r < dh_cut_up


	k3_buf=0;
	k4_buf=0;

	for(i=0;i<red_ph;i++){
	for(j=i+1;j<red_ph;j++){		
				if(d_ph1[i][j] > dh_cut_low && d_ph1[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[i][j] = dh(d_ph1[i][j],q_ph[i*acc_ph],q_ph[in_q_ph_ph[i][j]],eps,kappa,f_pi_eps0);  //in ev
						buf_energy = buf_energy + dh_ph[i][j];
						k3_buf++;
						/*dh_ph[k3_buf][0] = dh(d_ph1[i][j],q_ph[i*acc_ph],q_ph[in_q_ph_ph[i][j]],eps,kappa,f_pi_eps0);  //in ev
						dh_ph[k3_buf][1] = i*acc_ph;
						dh_ph[k3_buf][2] = j*acc_ph;
						dh_ph[k3_buf][3] = d_ph1[i][j];
						k3_buf++;*/
						    }
				  }
			
			     

	/*for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){		
				if(d_ph2[i][j] > dh_cut_low && d_ph2[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[k3_buf][0] = dh(d_ph2[i][j],q_ph,q_ph,eps,kappa,f_pi_eps0);
						dh_ph[k3_buf][1] = i;
						dh_ph[k3_buf][2] = i+1+j*acc_ph;
						dh_ph[k3_buf][3] = d_ph2[i][j];
						k3_buf++;
						    }
				  }
			     }*/
		
	

	// Calculate bp - floating patricle LJ pot for lj_cut_low < r < lj_cut_up 
	
	//for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<num_fl;j++){		
				if(d_ph_fl[i][j] > dh_cut_low && d_ph_fl[i][j] < dh_cut_up){ //printf("yes4");
						dh_ph_fl[i][j] = dh(d_ph_fl[i][j],q_ph[i],q_fl[j],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_ph_fl[i][j];
						k4_buf++;
						/*dh_ph_fl[k4_buf][0] = dh(d_ph_fl[i][j],q_ph[i],q_fl[j],eps,kappa,f_pi_eps0);
						dh_ph_fl[k4_buf][1] = i;
						dh_ph_fl[k4_buf][2] = j;
						dh_ph_fl[k4_buf][3] = d_ph_fl[i][j];
						k4_buf++;*/
						    }
				  }
			     }

	*k3 = k3_buf;
	*k4 = k4_buf;

	*e_dh_ph_phfl = buf_energy;
}



void calc_dh_nucl_MIN(int *k1, int *k2, int *k3, int num_bp, int num_nucl, double **d_nucl, double **d_nucl_ph, double **d_nucl_fl, double dh_nucl_cut_low, double dh_nucl_cut_up, double **dh_nucl, double **dh_nucl_ph, double **dh_nucl_fl, double *q_nucl, double *q_ph, double *q_fl, int **in_q_nucl_ph, int acc_ph, int num_fl, double eps, double kappa, double f_pi_eps0, int red_ph, double *energy)
{
int i,j,k1_buf,k2_buf,k3_buf;
double buf_energy, nucl_nucl_e, nucl_ph_e;
buf_energy = 0;
nucl_nucl_e = 0;
nucl_ph_e = 0;

	// Calculate inter-nucl dh pot (for dh_cut_low < r < dh_cut_up)


	k1_buf=0;
	k2_buf=0;
	k3_buf=0;
	for(i=0;i<num_nucl;i++){

	for(j=i;j<num_nucl;j++){	//printf("%lf\n",d_nucl[i][j]);
				if(d_nucl[i][j] > dh_nucl_cut_low && d_nucl[i][j] < dh_nucl_cut_up){ //printf("yes\n");
						dh_nucl[i][j] = dh(d_nucl[i][j],q_nucl[i],q_nucl[j],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_nucl[i][j];
						nucl_nucl_e = nucl_nucl_e + dh_nucl[i][j];
						k1_buf++;
						/*dh_nucl[k1_buf][0] = dh(d_nucl[i][j],q_nucl[i],q_nucl[j],eps,kappa,f_pi_eps0);
						dh_nucl[k1_buf][1] = i;
						dh_nucl[k1_buf][2] = j;
						dh_nucl[k1_buf][3] = d_nucl[i][j];
						buf_energy = buf_energy + dh_nucl[k1_buf][0];
						nucl_nucl_e = nucl_nucl_e + dh_nucl[k1_buf][0];
						k1_buf++;*/
						    }
				  }



	// Calculate nucl-bp dh pot (for dh_cut_low < r < dh_cut_up)


	
	for(j=0;j<red_ph;j++){	//printf("%lf\n",d_nucl_bp[i][j]);
				if(d_nucl_ph[i][j] > dh_nucl_cut_low && d_nucl_ph[i][j] < dh_nucl_cut_up){ //printf("yes1 j %d\n",j);
						dh_nucl_ph[i][j] = dh(d_nucl_ph[i][j],q_nucl[i],q_ph[in_q_nucl_ph[i][j]],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_nucl_ph[i][j];
						nucl_ph_e = nucl_ph_e + dh_nucl_ph[i][j];
						k2_buf++;
						/*dh_nucl_ph[k2_buf][0] = dh(d_nucl_ph[i][j],q_nucl[i],q_ph[in_q_nucl_ph[i][j]],eps,kappa,f_pi_eps0);
						dh_nucl_ph[k2_buf][1] = i;
						dh_nucl_ph[k2_buf][2] = i+j*acc_ph;
						dh_nucl_ph[k2_buf][3] = d_nucl_ph[i][j];
						buf_energy = buf_energy + dh_nucl_ph[k2_buf][0];
						nucl_ph_e = nucl_ph_e + dh_nucl_ph[k2_buf][0];
						k2_buf++;*/
						    }
				  }


	
	// Calculate nucl - floating patricle dh pot (for dh_cut_low < r < dh_cut_up)
	

	for(j=0;j<num_fl;j++){		
				if(d_nucl_fl[i][j] > dh_nucl_cut_low && d_nucl_fl[i][j] < dh_nucl_cut_up){ //printf("yes2");
						dh_nucl_fl[i][j] = dh(d_nucl_fl[i][j],q_nucl[i],q_fl[j],eps,kappa,f_pi_eps0);
						buf_energy = buf_energy + dh_nucl_fl[i][j];
						k3_buf++;
						/*dh_nucl_fl[k1_buf][0] = dh(d_nucl_fl[i][j],q_nucl[i],q_fl[j],eps,kappa,f_pi_eps0);
						dh_nucl_fl[k1_buf][1] = i;
						dh_nucl_fl[k1_buf][2] = j;
						dh_nucl_fl[k1_buf][3] = d_nucl_fl[i][j];
						buf_energy = buf_energy + dh_nucl_fl[k3_buf][0];
						k3_buf++;*/
						    }
				  }

			     }
	
	*k1 = k1_buf;
	*k2 = k2_buf;
	*k3 = k3_buf;

	*energy = buf_energy;

//printf("DH nucl_nucl_e %lf\n",nucl_nucl_e);
//printf("DH nucl_ph_e %lf\n",nucl_ph_e);
}


double calc_sed_coeff(int no_nucl, double **d_nucl)
{
int i,j;
double inv_d,sed_co;
inv_d = 0;
for(i=0;i<no_nucl;i++){
for(j=i+1;j<no_nucl;j++){
			inv_d = inv_d + 1/d_nucl[i][j];
			}
		       }
sed_co = 11.1*(1 + 2*54.6/no_nucl * inv_d);

return(sed_co);
}


double calc_rg(int no_nucl, double **nucl_pos)
{
int i,j;
double mean_dist[3],rg_sq[3],rg,rg_buf[no_nucl],rg_sum;
rg_sum=0;
for(i=0;i<3;i++) mean_dist[i] = 0;

for(j=0;j<3;j++){
for(i=0;i<no_nucl;i++){
			mean_dist[j] = mean_dist[j] + nucl_pos[i][j];
		       }
		}

for(i=0;i<3;i++) mean_dist[i] = mean_dist[i]/no_nucl;
//printf("mean %lf %lf %lf\n", mean_dist[0],mean_dist[1],mean_dist[2]);

for(i=0;i<no_nucl;i++){
for(j=0;j<3;j++){
			rg_sq[j] = (nucl_pos[i][j] - mean_dist[j])*(nucl_pos[i][j] - mean_dist[j]);
		}
		rg_buf[i] = rg_sq[0]+rg_sq[1]+rg_sq[2];
		rg_sum = rg_sum + rg_buf[i];
		       }


rg = sqrt(1/(double) no_nucl * rg_sum);
return(rg);
}



