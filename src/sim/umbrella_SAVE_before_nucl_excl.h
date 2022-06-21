// Functions needed for Umbrella sampling

void unit_matrix3(double **matrix)
{
int i,j;
for(i=0;i<3;i++){
for(j=0;j<3;j++){
		if(i==j){matrix[i][j] = 1.0;}
		else {matrix[i][j] = 0.0;}
		}
		}

}

void zero_matrix_n(double **matrix, int n)
{
int i,j;
for(i=0;i<n;i++){
for(j=0;j<n;j++){
		matrix[i][j] = 0.0;
		}
		}

}

void rot_x(double **rot, double phi)
{
int i,j;
rot[0][1] = rot[1][0] = rot[0][2] = rot[2][0] = 0;
rot[0][0] = 1.0;
rot[1][1] = cos(phi);
rot[2][2] = cos(phi);
rot[2][1] = sin(phi);
rot[1][2] = -sin(phi);
}

void rot_y(double **rot, double phi)
{
int i,j;
rot[0][1] = rot[1][0] = rot[1][2] = rot[2][1] = 0;
rot[1][1] = 1.0;
rot[0][0] = cos(phi);
rot[2][2] = cos(phi);
rot[0][2] = sin(phi);
rot[2][0] = -sin(phi);
}


void rot_z(double **rot, double phi)
{
int i,j;
rot[0][2] = rot[2][0] = rot[1][2] = rot[2][1] = 0;
rot[2][2] = 1.0;
rot[0][0] = cos(phi);
rot[1][1] = cos(phi);
rot[1][0] = sin(phi);
rot[0][1] = -sin(phi);
}

void mat_mult_n(double **result, double **a, double **b, int n)
{
int i,j,k;

zero_matrix_n(result,n);

for(i=0;i<n;i++){
for(j=0;j<n;j++){
for(k=0;k<n;k++){
result[i][j] = result[i][j] + a[i][k]*b[k][j];
		}
		}
		}

}

void mat_vec_mult_n(double *result, double **matrix, double *vec, int n)
{
int i,j,k;

for(i=0;i<n;i++){result[i]=0;}

for(i=0;i<n;i++){
for(k=0;k<n;k++){
result[i] = result[i] + matrix[i][k]*vec[k];
		}
		}

}

void copy_2darray(double **a, double **b, int n)
{
int i,j;
for(i=0;i<n;i++){
for(j=0;j<n;j++){
		a[i][j] = b[i][j];
		}
		}

}

void transp_matr_n(double **tr_matrix, double **matrix, int n)
{
int i,j;
for(i=0;i<n;i++){
for(j=0;j<n;j++){tr_matrix[i][j] = matrix[j][i];}}

}

double vec_len(double *vec, int n)
{
int i;
double len,res;
len = 0;
for(i=0;i<n;i++){len = len + vec[i]*vec[i];}
res = sqrt(len);

return(res);
}

double dot_product (double *vec1, double *vec2)
{
double dot_pr;
dot_pr = 0;
int i;

for(i=0;i<3;i++){dot_pr = dot_pr + vec1[i]*vec2[i];}

return dot_pr;
}

void end_to_end_dist(int itot, int no_rot_pars, double **xconf, double **store_r_mst, double ***store_orient_t_i, double *end_r_mst, double *ete_dist)
{
int i,j,k;
const float pi = 3.14159265358979323846;
double **rotz1, **roty, **rotz2;
double g,p,tw;	//gamma,phi,twist (for end-to-end distance)
double **res_buf1, **res_buf2, **result, **res_prev;
double *transl;
double **t_mst, *r_mst, *r_bppar;
double ete_dist_buf;


rotz1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
rotz2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
roty = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
result = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_prev = dbl_aloc2d(no_rot_pars,no_rot_pars);
transl = malloc(no_rot_pars*sizeof(double));
t_mst = dbl_aloc2d(no_rot_pars,no_rot_pars);
r_mst = malloc(no_rot_pars*sizeof(double));
r_bppar = malloc(no_rot_pars*sizeof(double));


	unit_matrix3(res_prev);		// initialize matrix as unit matrix (= T_i)
	for(i=0;i<3;i++){end_r_mst[i]=0;}


	for(i=0;i<itot;i++){	// result is multiplication of all

			
			p = atan(xconf[3][i]/xconf[4][i]); //in rad (mixture of tilt and roll from Hassan&Calladine 95)
			g = xconf[4][i]/cos(p)*pi/180;	//in rad (mixture of tilt and roll from Hassan&Calladine 95)
			tw = xconf[5][i]*pi/180;	//in rad
		

			//Get Mid-step triad
			rot_y(roty,g/2);
			rot_z(rotz1,p);
			rot_z(rotz2,tw/2-p);

			

			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);


			mat_mult_n(t_mst,res_prev,res_buf2,3);	//get mid-step matrix out of T_i (=T_i * R_z*R_y*R_z)


			//Get r-vector i+1 from mid-step triad
		
			r_bppar[0] = xconf[0][i];	//shift
			r_bppar[1] = xconf[1][i];	//slide	
			r_bppar[2] = xconf[2][i];	//rise

			mat_vec_mult_n(r_mst,t_mst,r_bppar,3);

			double buf;
			for(j=0;j<3;j++){buf = end_r_mst[j];
					end_r_mst[j] = buf + r_mst[j];
					store_r_mst[i][j] = end_r_mst[j];
					}	// This gives the end r-vector


				
			//Get T_i+1 out of T_i
			rot_y(roty,g);
			rot_z(rotz1,tw/2+p);
			rot_z(rotz2,tw/2-p);

			
			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);

			
			
			mat_mult_n(result,res_prev,res_buf2,3);	//get T_i+1 out of T_i

			copy_2darray(res_prev,result,3);	//Copy T_i+1 in res_prev

			for(j=0;j<3;j++){
			for(k=0;k<3;k++){
					store_orient_t_i[i][j][k] = res_prev[j][k]; // store orientation matrix T_i+1
					}
					}

				}


	ete_dist_buf = vec_len(end_r_mst,3);	// End-to-end distance

	*ete_dist = ete_dist_buf;


free(rotz1);
free(rotz2);
free(roty);
free(res_buf1);
free(res_buf2);
free(res_prev);
free(result);
free(transl);
free(t_mst);
free(r_mst);
free(r_bppar);

}


void partial_ete(int itt, int itot, int no_rot_pars, double **xconf, double **store_r_mst, double ***store_orient_t_i, double *end_r_mst, double *ete_dist)
{
int i,j,k;
const double pi = 3.14159265358979323846;
double **rotz1, **roty, **rotz2;
double g,p,tw;	//gamma,phi,twist (for end-to-end distance)
double **res_buf1, **res_buf2, **result, **res_prev;
double *transl;
double **t_mst, *r_mst, *r_bppar;
double ete_dist_buf;


rotz1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
rotz2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
roty = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
result = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_prev = dbl_aloc2d(no_rot_pars,no_rot_pars);
transl = malloc(no_rot_pars*sizeof(double));
t_mst = dbl_aloc2d(no_rot_pars,no_rot_pars);
r_mst = malloc(no_rot_pars*sizeof(double));
r_bppar = malloc(no_rot_pars*sizeof(double));

	if(itt == 0){
	unit_matrix3(res_prev);		// initialize matrix as unit matrix (= T_i)
	for(i=0;i<3;i++){end_r_mst[i]=0;}
	}
	else{
	for(i=0;i<3;i++){end_r_mst[i]=store_r_mst[itt-1][i];}
	copy_2darray(res_prev,store_orient_t_i[itt-1],3);
	}
	


	for(i=itt;i<itot;i++){	// result is multiplication of all

			
			p = atan(xconf[3][i]/xconf[4][i]); //in rad (mixture of tilt and roll from Hassan&Calladine 95)
			g = xconf[4][i]/cos(p)*pi/180;	//in rad (mixture of tilt and roll from Hassan&Calladine 95)
			tw = xconf[5][i]*pi/180;	//in rad
		

			//Get Mid-step triad
			rot_y(roty,g/2);
			rot_z(rotz1,p);
			rot_z(rotz2,tw/2-p);

			

			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);


			mat_mult_n(t_mst,res_prev,res_buf2,3);	//get mid-step matrix out of T_i (=T_i * R_z*R_y*R_z)


			//Get r-vector i+1 from mid-step triad
		
			r_bppar[0] = xconf[0][i];	//shift
			r_bppar[1] = xconf[1][i];	//slide	
			r_bppar[2] = xconf[2][i];	//rise

			mat_vec_mult_n(r_mst,t_mst,r_bppar,3);

			double buf;
			for(j=0;j<3;j++){buf = end_r_mst[j];
					end_r_mst[j] = buf + r_mst[j];
					store_r_mst[i][j] = end_r_mst[j];
					}	// This gives the end r-vector


				
			//Get T_i+1 out of T_i
			rot_y(roty,g);
			rot_z(rotz1,tw/2+p);
			rot_z(rotz2,tw/2-p);

			
			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);

			
			
			mat_mult_n(result,res_prev,res_buf2,3);	//get T_i+1 out of T_i

			copy_2darray(res_prev,result,3);	//Copy T_i+1 in res_prev

			for(j=0;j<3;j++){
			for(k=0;k<3;k++){
					store_orient_t_i[i][j][k] = res_prev[j][k]; // store orientation matrix T_i+1
					}
					}

				}


	ete_dist_buf = vec_len(end_r_mst,3);	// End-to-end distance

	*ete_dist = ete_dist_buf;


free(rotz1);
free(rotz2);
free(roty);
free(res_buf1);
free(res_buf2);
free(res_prev);
free(result);
free(transl);
free(t_mst);
free(r_mst);
free(r_bppar);

}



void rot_matr(int itt, double **xconf, double **rot_mat, int half_step)
{
int i,j,k;
const float pi = 3.14159265358979323846;
int no_rot_pars = 3;
double **rotz1, **roty, **rotz2;
double g,p,tw;	//gamma,phi,twist (for defining rotation matrices)
double **res_buf1, **res_buf2;



rotz1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
rotz2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
roty = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf2 = dbl_aloc2d(no_rot_pars,no_rot_pars);


p = atan(xconf[3][itt]/xconf[4][itt]); //in rad (mixture of tilt and roll from Hassan&Calladine 95)
g = xconf[4][itt]/cos(p)*pi/180;	//in rad (mixture of tilt and roll from Hassan&Calladine 95)
//g = sqrt(xconf[3][itt]*xconf[3][itt]+xconf[4][itt]*xconf[4][itt])*pi/180;
//p = acos(xconf[4][itt]*(pi/180)/g);
tw = xconf[5][itt]*pi/180;	//in rad


			if(half_step == 1){ //Get Mid-step triad
			
			rot_y(roty,g/2);
			rot_z(rotz1,p);
			rot_z(rotz2,tw/2-p);
			
					}

			else{	// Rotation matrix for T_i+1 out of T_i


			rot_y(roty,g);
			rot_z(rotz1,tw/2+p);
			rot_z(rotz2,tw/2-p);
				}

			
			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);
	

			copy_2darray(rot_mat,res_buf2,3);	//Copy rotation matrix for mid-step triad
			

free(rotz1);
free(rotz2);
free(roty);
free(res_buf1);
free(res_buf2);

}



void smart_rec(double **xconf, double **xconft, int no_rot_pars, int itt, int itot, double ***smart_t_i, double **smart_r_i, double ***store_orient_t_i, double **store_r_mst)
{
int i,j,flag;
double **rot_mat;
double **t_mst;
double internal[6],r_mst[3];
double buf_r[3], buf;
double **tr_store_orient_t_i, **Q;
for(i=0;i<6;i++){internal[i] = xconft[i][itt];}

rot_mat = dbl_aloc2d(no_rot_pars,no_rot_pars);
t_mst = dbl_aloc2d(no_rot_pars,no_rot_pars);
tr_store_orient_t_i = dbl_aloc2d(no_rot_pars,no_rot_pars);
Q = dbl_aloc2d(no_rot_pars,no_rot_pars);

flag=0;
for(i=3;i<6;i++){if(xconf[i][itt] != xconft[i][itt]){flag = 1;}}	// flag=1 when at least one rotational coordinate changes

//use the Hassan & Calladine method for basepair i when change bp-steps at position i
if(itt == 0){

if(flag == 1){
rot_matr(itt,xconft,rot_mat,0);
for(i=0;i<3;i++){
for(j=0;j<3;j++){
smart_t_i[itt][i][j] = rot_mat[i][j];
		}
		}
	      }

rot_matr(itt,xconft,rot_mat,1);
for(i=0;i<3;i++){
for(j=0;j<3;j++){
		t_mst[i][j] = rot_mat[i][j];	
		}
		}

mat_vec_mult_n(smart_r_i[itt],t_mst,internal,3);


	    }
else{
if(flag == 1){
rot_matr(itt,xconft,rot_mat,0);
mat_mult_n(smart_t_i[itt],store_orient_t_i[itt-1],rot_mat,3);	// Get T_i
	      }

rot_matr(itt,xconft,rot_mat,1);
mat_mult_n(t_mst,store_orient_t_i[itt-1],rot_mat,3);	// Get T_mst
mat_vec_mult_n(r_mst,t_mst,internal,3);

for(j=0;j<3;j++){
	smart_r_i[itt][j] = smart_r_i[itt-1][j] + r_mst[j];			// Get r_itt
		}


     }


if(flag == 1){
// Obtain Transformation matrix Q
transp_matr_n(tr_store_orient_t_i,store_orient_t_i[itt],3);
mat_mult_n(Q,smart_t_i[itt],tr_store_orient_t_i,3);	// Get T_mst


// Apply Q to all bp triads from itt+1 to end
for(i=itt+1;i<itot;i++){
	mat_mult_n(smart_t_i[i],Q,store_orient_t_i[i],3);
			}
	    

// Apply Q to all position vectors from itt+1 to end
double r_dif[3];
for(i=itt+1;i<itot;i++){
	for(j=0;j<3;j++){r_dif[j] = store_r_mst[i][j]-store_r_mst[i-1][j];}
		mat_vec_mult_n(buf_r,Q,r_dif,3);
	for(j=0;j<3;j++){smart_r_i[i][j] = smart_r_i[i-1][j] + buf_r[j];}
			}


	      }
	      else{
			for(i=itt+1;i<itot;i++){
				for(j=0;j<3;j++){smart_r_i[i][j] = smart_r_i[i-1][j] + store_r_mst[i][j]-store_r_mst[i-1][j];}
						}
		   }

free(rot_mat);
free(t_mst);
free(tr_store_orient_t_i);
free(Q);

}

void transfer_smart_to_store(int itt, int itot, double ***smart_t_i, double **smart_r_i, double ***store_orient_t_i, double **store_r_mst)
{
int i,j,k;

for(i=itt;i<itot;i++){
			for(j=0;j<3;j++){
			for(k=0;k<3;k++){
					store_orient_t_i[i][j][k] = smart_t_i[i][j][k]; // store orientation matrix T_i+1
					}
					}


			for(j=0;j<3;j++){
					store_r_mst[i][j] = smart_r_i[i][j];
					}	

	  	    }

}

void transfer_store_to_smart(int itt, int itot, double ***smart_t_i, double **smart_r_i, double ***store_orient_t_i, double **store_r_mst)
{
int i,j,k;

for(i=itt;i<itot;i++){
			for(j=0;j<3;j++){
			for(k=0;k<3;k++){
					smart_t_i[i][j][k] = store_orient_t_i[i][j][k]; // store orientation matrix T_i+1
					}
					}


			for(j=0;j<3;j++){
					smart_r_i[i][j] = store_r_mst[i][j];
					}	

	  	    }

}


double Determinant(double **a,int n)
{
    int i,j,j1,j2 ;                    // general loop and matrix subscripts
    double det = 0 ;                   // init determinant
    double **m = NULL ;                // pointer to pointers to implement 2d
                                       // square array

    if (n < 1)    {   }                // error condition, should never get here

    else if (n == 1) {                 // should not get here
        det = a[0][0] ;
        }

    else if (n == 2)  {                // basic 2X2 sub-matrix determinate
                                       // definition. When n==2, this ends the
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1] ;// the recursion series
        }


                                       // recursion continues, solve next sub-matrix
    else {                             // solve the next minor by building a
                                       // sub matrix
        det = 0 ;                      // initialize determinant of sub-matrix

                                           // for each column in sub-matrix
        for (j1 = 0 ; j1 < n ; j1++) {
                                           // get space for the pointer list
            m = (double **) malloc((n-1)* sizeof(double *)) ;

            for (i = 0 ; i < n-1 ; i++)
                m[i] = (double *) malloc((n-1)* sizeof(double)) ;
                       //     i[0][1][2][3]  first malloc
                       //  m -> +  +  +  +   space for 4 pointers
                       //       |  |  |  |          j  second malloc
                       //       |  |  |  +-> _ _ _ [0] pointers to
                       //       |  |  +----> _ _ _ [1] and memory for
                       //       |  +-------> _ a _ [2] 4 doubles
                       //       +----------> _ _ _ [3]
                       //
                       //                   a[1][2]
                      // build sub-matrix with minor elements excluded
            for (i = 1 ; i < n ; i++) {
                j2 = 0 ;               // start at first sum-matrix column position
                                       // loop to copy source matrix less one column
                for (j = 0 ; j < n ; j++) {
                    if (j == j1) continue ; // don't copy the minor column element

                    m[i-1][j2] = a[i][j] ;  // copy source element into new sub-matrix
                                            // i-1 because new sub-matrix is one row
                                            // (and column) smaller with excluded minors
                    j2++ ;                  // move to next sub-matrix column position
                    }
                }

            det += pow(-1.0,1.0 + j1 + 1.0) * a[0][j1] * Determinant(m,n-1) ;
                                            // sum x raised to y power
                                            // recursively get determinant of next
                                            // sub-matrix which is now one
                                            // row & column smaller

            for (i = 0 ; i < n-1 ; i++) free(m[i]) ;// free the storage allocated to
                                            // to this minor's set of pointers
            free(m) ;                       // free the storage for the original
                                            // pointer to pointer
        }
    }
    return(det) ;
}



void cart_rec_old(int itot, int no_rot_pars, double **xconf, double **store_r_mst, double ***store_orient_t_i, double ***rot_store)
{
int i,j,k,m;
const float pi = 3.14159265358979323846;
double **rotz1, **roty, **rotz2;
double g,p,tw;	//gamma,phi,twist (for end-to-end distance)
double **res_buf1, **res_buf2, **result, **res_prev;
double *transl;
double **t_mst, *r_mst, *r_bppar, *end_r_mst;
double ete_dist_buf;


rotz1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
rotz2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
roty = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
result = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_prev = dbl_aloc2d(no_rot_pars,no_rot_pars);
transl = malloc(no_rot_pars*sizeof(double));
t_mst = dbl_aloc2d(no_rot_pars,no_rot_pars);
r_mst = malloc(no_rot_pars*sizeof(double));
r_bppar = malloc(no_rot_pars*sizeof(double));
end_r_mst = malloc(no_rot_pars*sizeof(double));


	unit_matrix3(res_prev);		// initialize matrix as unit matrix (= T_i)
	for(i=0;i<3;i++){end_r_mst[i]=0;}

	m=0;
	for(i=0;i<itot;i++){	// result is multiplication of all

			
			p = atan(xconf[3][i]/xconf[4][i]); //in rad (mixture of tilt and roll from Hassan&Calladine 95)
			g = xconf[4][i]/cos(p)*pi/180;	//in rad (mixture of tilt and roll from Hassan&Calladine 95)
			tw = xconf[5][i]*pi/180;	//in rad
		

			//Get Mid-step triad
			rot_y(roty,g/2);
			rot_z(rotz1,p);
			rot_z(rotz2,tw/2-p);

			

			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);
			
			copy_2darray(rot_store[m],res_buf2,3);	// Store rot matrix for mid-step triad


			mat_mult_n(t_mst,res_prev,res_buf2,3);	//get mid-step matrix out of T_i (=T_i * R_z*R_y*R_z)


			//Get r-vector i+1 from mid-step triad
		
			r_bppar[0] = xconf[0][i];	//shift
			r_bppar[1] = xconf[1][i];	//slide	
			r_bppar[2] = xconf[2][i];	//rise

			mat_vec_mult_n(r_mst,t_mst,r_bppar,3);

			double buf;
			for(j=0;j<3;j++){buf = end_r_mst[j];
					end_r_mst[j] = buf + r_mst[j];	//printf("i %d j %d end_r_mst %f\n",i,j, end_r_mst[j]);
					store_r_mst[i][j] = end_r_mst[j]; //printf("2\n"); 
					}	// This gives the end r-vector


				
			//Get T_i+1 out of T_i
			rot_y(roty,g);
			rot_z(rotz1,tw/2+p);
			rot_z(rotz2,tw/2-p);

			
			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);

			
			copy_2darray(rot_store[m+1],res_buf2,3);	// Store rot matrix for T_i+1
			m=m+2;
			
			
			mat_mult_n(result,res_prev,res_buf2,3);	//get T_i+1 out of T_i

			copy_2darray(res_prev,result,3);	//Copy T_i+1 in res_prev

			for(j=0;j<3;j++){
			for(k=0;k<3;k++){
					store_orient_t_i[i][j][k] = res_prev[j][k]; // store orientation matrix T_i+1
					}
					}

				}


	ete_dist_buf = vec_len(end_r_mst,3);	// End-to-end distance



free(rotz1);
free(rotz2);
free(roty);
free(res_buf1);
free(res_buf2);
free(res_prev);
free(result);
free(transl);
free(t_mst);
free(r_mst);
free(r_bppar);
free(end_r_mst);

}


/* ---------- End of functions ---------- */
			
