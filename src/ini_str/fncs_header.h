
/* analysis.c */
void get_data(long *RY_idx, char **atom_id, char **bp_seq, double **xyz,
              char *base_type, char *strand, long *base_num,
              long *num_atom, long *num_bp, char *inpfil, FILE * fo);
void get_list(long *C1, long *N, long *P, long *O4, long *C68, long *chi,
              long *mc, long *nmc, long *sugar, long *base, long *n_base,
              long *RY_idx, char **atom_id, long num_bp, long num_atom, double **xyz);
void set_str(double **xyz, long *C1, long *N, long *C68, long *base,
             long *base_num, long n_base, long num_atom, long num_bp, FILE * fo);
void wrt_hel_dat(long num_atom, char **atom_id, char *base_type,
                 char *strand, long *base_num, double **xyz,
                 long *base, long n_base, char *filstr);
void sugar_ana(long *sugar, long num_bp, double **xyz, char **bp_seq, FILE * fo);
void polar_coord(long *alist, long midnum, char *asym, double **xyz,
                 char **bp_seq, FILE * fo);
void get_lambda(long *C1, long *N, long *C68, long num_bp, double **xyz,
                char **bp_seq, FILE * fo);
void groove_width(long *P, long midnum, long *O4, double **xyz, long num_bp,
                  char **bp_seq, FILE * fo);
void get_frame(double **bx, double **by, double **bz, double **borg, long *xdir,
               long *RY_idx, long *C68, long num_bp, long *base, long n_base,
               long *base_num, char **atom_id, double **xyz, char **bp_seq, FILE * fo);
void get_bp_par(double *bp_par, double **xyz_rb, double *org_rb,
                double **xyz_lb, double *org_lb);
void get_step_par(double *step_par, double **mstxyz, double **xyz_lbp,
                  double *org_lbp, double **xyz_ubp, double *org_ubp);
void CEHS_par(double **mxyz, double **morg, double **bp_par, double **P_mst,
              double *mtwist, double **bx, double **by, double **bz, double **borg,
              long *P, long numP1, long num_bp, char **bp_seq,
              double **xyz, char *filstr, FILE * fo);
void wrt_mst(char **atom_id, char *base_type, char *strand, long *base_num,
             double **xyz, char **bp_seq, long num_atom, long num_bp,
             double **mxyz, double **morg, char *filstr);
void GLH_par(double **bx, double **by, double **bz, double **borg,
             long num_bp, char **bp_seq, double **bp_par, char *filstr, FILE * fo);
void classify(double **P_mst, long num_bp, char **bp_seq, long xdir,
              double mtwist, FILE * fo);

/* cehs.c */
void get_data_lite(long *RY_idx, char **atom_id, char **bp_seq, double **xyz,
                   char *base_type, long *base_num, long *num_atom, long *num_bp,
                   char *inpfil, FILE * fo);
void get_list_lite(long *C68, long *base, long *n_base, long *RY_idx, char **atom_id,
                   long num_bp, long num_atom);
void get_frame_lite(double **bx, double **by, double **bz, double **borg,
                    long *RY_idx, long *C68, long num_bp, long *base, long n_base,
                    long *base_num, char **atom_id, double **xyz, char **bp_seq,
                    FILE * fo);
void CEHS_par_lite(double **bx, double **by, double **bz, double **borg,
                   long num_bp, char **bp_seq, char *filstr, FILE * fo);

/* cmn_fncs.c */
void jacobi(double **a, long n, double *d, double **v, long *nrot);
void eigsrt(double *d, double **v, long n);
void dludcmp(double **a, long n, long *indx, double *d);
void dlubksb(double **a, long n, long *indx, double *b);
double dist_ab(double a, double b);
void zero_dmatrix(double **a, long nr, long nc);
void one_dmatrix(double **a, long nr, long nc);
void zero_dvector(double *a, long n);
void zero_lvector(long *a, long n);
void identity_dmatrix(double **a, long n);
double deg2rad(double theta);
double rad2deg(double theta);
void rotx(double **rotmat, double theta);
void roty(double **rotmat, double theta);
void rotz(double **rotmat, double theta);
void mul_dmatrix(double **o, double **a, double **b,
                 long nr_a, long nc_a, long nr_b, long nc_b);
void tr_dmatrix(double **o, double **a, long nr, long nc);
void cov_dmatrix(double **cov_mtx, double **xyz, long nr, long nc);
void mean_dmatrix(double *o, double **a, long nr, long nc);
void std_dmatrix(double *o, double **a, long nr, long nc);
double mean_dvector(double *a, long n);
double std_dvector(double *a, long n);
double dot_dvector(double *a, double *b, long n);
void cross_dvector(double *o, double *a, double *b);
double mag_dvector(double *a, long n);
void norm_dvector(double *o, double *a, long n);
void lsplane(double *vec_normal, double *org_dist, double *pnts_dist,
             double **xyz, long nr);
void lsline(double *vec_line, double *xyzpoint, double *pnts_dist, double **xyz, long nr);
void alignz(double **rotmat, double *rot_axis);
void arbrot(double **rotmat, double *rot_axis, double theta);
double magang(double *veca, double *vecb);
double vec_ang(double *veca, double *vecb, double *vec_ref);
void bpfunc(double **orien, double *pos, double *param);
void stepfunc(double **orien, double **mst, double *pos, double *param);
void cirshift(long *o, long *a, long n, long m);
void vec_orth(double *o, double *veca, double *vecb);
double val_ang(double **abcxyz);
double torsion(double **abcdxyz);
void i_diff(long *o, long *a, long n);
void logical_not(long *o, long *a, long n);
void isort(long n, long *ra);
void dsort(long n, double *ra);
void dinverse(double **a, long n, double **y);
void dvec_x_dmtx(double *o, double *a, double **b, long na, long nr, long nc);
void duplicate_dmtx(double **o, double **a, long nr, long nc);
double max_dvector(double *a, long n);
double min_dvector(double *a, long n);
void upper(char *a, long n);
void read_stdin_str(char *strval);
void comp_base(char *o, char *a, long n);
FILE *open_file(char *filnam, char *filmod);
void get_sequence(char **bpseq, long *nbp);
void nrptx(char *o, char x, long n);
void prt_sep(FILE * f1, char x, long n);
void rdalc(long *num, long *nbond, char **asym, long **lkg, double **xyz, char *filnam);
void wrtalc(long num, long nbond, char **asym, long **lkg, double **xyz, char *filnam);
void rdpdb(long *num, char **asym, char *btype, char *strand,
           long *bnum, double **xyz, char *filnam);
void wrtpdb(long num, char **asym, char *btype, char *strand,
            long *bnum, double **xyz, char *filnam);
void str_idx(long *n_match, long *idx, char **matx, char *y, long nr);
void set_bp(long *num1, long *num2, char **asym1, char **asym2,
            double **xyz1, double **xyz2, char *bpname,
            double *param, long xdir, char *dir_name);

/* nrutil.c */
void nrerror(char error_text[]);
long *lvector(long nl, long nh);
char *cvector(long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
long **lmatrix(long nrl, long nrh, long ncl, long nch);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_lvector(long *v, long nl, long nh);
void free_cvector(char *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);

/* pdb2inp.c */

/* rebuild.c */
void align_helix(double **oxyz, double **xyz, long natom, double **org_xyz, long nbp);
long get_xdir(double twist);
void str_display(long imdl, char *filnam);
void set_bppar(double **bppar, long nbp, char **bpseq);
void alc_connect(char *filinp, char *filout);
void wrt_seq_CEHS(long nbp, char **bpseq, double **steppar);
void std_CEHS(long *nbp, char **bpseq, double **steppar);
void dimer_CEHS(long *nbp, char **bpseq, double **steppar);
void trimer_CEHS(long *nbp, char **bpseq, double **steppar);
void step_CEHS(long *nbp, char **bpseq, double **steppar,int itot, char **ist, double **xconf, int phi, char **output);
void exact_CEHS(long *nbp, char **bpseq, double **bppar, double **steppar);
void atomic_CEHS(long nbp, char **bpseq, double **bppar, double **steppar,
                 long xdir, char *filnam);
void block_CEHS(long nblk, double **steppar, long xdir, char *filnam);
void CEHS_build(int itot, char **ist, double **xconf, int phi, char **output, char out_folder[]);
void wrt_seq_GLH(long nbp, char **bpseq, double **glhpar);
void std_GLH(long *nbp, char **bpseq, double **glhpar);
void step_GLH(long *nbp, char **bpseq, double **glhpar);
void exact_GLH(long *nbp, char **bpseq, double **bppar, double **glhpar);
void atomic_GLH(long nbp, char **bpseq, double **bppar, double **glhpar,
                long xdir, char *filnam);
void block_GLH(long nblk, double **glhpar, long xdir, char *filnam);
void GLH_build(void);

/* schnaap.c */

/* schnarp.c */

/* set_base.c */

/* sr2inp.c */
