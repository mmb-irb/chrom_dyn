/*<<<
 * ------------------------------------------------------------------------
 * This file is part of the SCHNAarP software package for the analysis and
 * rebuilding of double helical nucleic acid structures. Detailed description
 * of the underlying algorithms can be found in the following three papers:
 * 
 *    (1) Xiang-Jun Lu, M. A. El Hassan and C. A. Hunter (1997), "Structure
 *        and Conformation of Helical Nucleic Acids: Analysis Program
 *        (SCHNAaP)", J. Mol. Biol. 273, 668--680.
 * 
 *    (2) Xiang-Jun Lu, M. A. El Hassan and C. A. Hunter (1997), "Structure
 *        and Conformation of Helical Nucleic Acids: Rebuilding Program
 *        (SCHNArP)", J. Mol. Biol. 273, 681--691.
 * 
 *    (3) M. A. El Hassan and C. R. Calladine (1995), "The Assessment of the
 *        Geometry of Dinucleotide Steps in Double-Helical DNA; a New Local
 *        Calculation Scheme", J. Mol. Biol. 251, 648--664.
 * 
 * Copyright (C) 1996-2002 Xiang-Jun Lu. Permission to use, copy, and modify
 * this software is hereby granted to all academic and not-for-profit
 * institutions without fee, provided that the above papers (or part of) are
 * properly cited in related publications.  However, the right to use this
 * software in conjunction with for profit activities, and the right to
 * distribute the software or modified or extended versions thereof for profit
 * are *NOT* granted except by prior arrangement and written consent of the
 * copyright holders.
 * 
 * The software is provided "AS-IS" and without warranty of any kind, express,
 * implied or otherwise, including without limitation, any warranty of
 * merchantability or fitness for a particular purpose.
 * 
 * The software was written by Xiang-Jun Lu (3dna.lu@gmail.com). It started
 * off his PhD thesis work at Dr. Chris Hunter Laboratory at University of
 * Sheffield, U.K. The source code should be valuable for understanding the
 * elegant CEHS algorithm, and many other aspects of nucleic acid structures.
 * For real world applications, however, you are strongly recommended to
 * switch to 3DNA (http://3dna.rutgers.edu).
 * 
 * Revision history -- please see file 'Change.log'
 * 
 * $LastChangedDate: 2009-06-01 18:34:25 -0400 (Mon, 01 Jun 2009) $
 * $LastChangedRevision: 1556 $
 * $LastChangedBy: xiangjun $
 * ------------------------------------------------------------------------
>>>*/

#include "schna_ar_p.h"

#define BASE_ATOM_NUM 100

/* Set a base in its reference triad.
   The input data file should be in standard PDB format */
int main(long argc, char *argv[])
{
    long i, j, num, *bnum, n_match, *idx;
    long posYnum, negYnum, *base_atm, n_base = 0;
    char **asym, *btype, *strand, inpfil[80], outfil[80];
    double **xyz, **xyz2, *y_axis, *tmp_vec;
    double *b_org, **base_xyz;
    double *b_normal, *z_axis, *x_axis;
    double dval, *pnts_dist, **xyz_axis;
    long N1, C6;

    if (argc == 1) {
        printf("Name of the original PDB data file: ");
        read_stdin_str(inpfil);
        printf("Name of the new PDB data file: ");
        read_stdin_str(outfil);
    } else if (argc == 2) {
        strcpy(inpfil, argv[1]);
        printf("Name of the new PDB data file: ");
        read_stdin_str(outfil);
    } else {
        strcpy(inpfil, argv[1]);
        strcpy(outfil, argv[2]);
    }

    asym = cmatrix(1, BASE_ATOM_NUM, 1, 4);
    btype = cvector(1, BASE_ATOM_NUM);
    strand = cvector(1, BASE_ATOM_NUM);
    bnum = lvector(1, BASE_ATOM_NUM);
    xyz = dmatrix(1, BASE_ATOM_NUM, 1, 3);
    idx = lvector(1, BASE_ATOM_NUM);

    rdpdb(&num, asym, btype, strand, bnum, xyz, inpfil);

    if (bnum[num] != bnum[1])
        nrerror("Too many (more than one) residues!");

    for (i = 1; i <= num; i++) {
        if (asym[i][2] == '*')
            asym[i][2] = '\'';
        if (strcmp(asym[i], "O1'") == 0)
            strcpy(asym[i], "O4'");
        if (strcmp(asym[i], "OL ") == 0)
            strcpy(asym[i], "O1P");
        if (strcmp(asym[i], "OR ") == 0)
            strcpy(asym[i], "O2P");

        strand[i] = 'A';
        bnum[i] = 1;
    }

    /* Find the proper atoms to define y-axis */

    str_idx(&n_match, idx, asym, "N9 ", num);  /* Check for N9 atom */
    if (n_match > 1)
        nrerror("More than on N9 atom for this base!!!");

    if (n_match == 1) {  /* R base */
        str_idx(&n_match, idx, asym, "C4 ", num);
        if (n_match != 1)
            nrerror("There should be only one C4 atom for purine!");
        posYnum = idx[1];

        str_idx(&n_match, idx, asym, "N1 ", num);
        if (n_match != 1)
            nrerror("There should be only one N1 atom for purine!");
        negYnum = idx[1];
    } else {  /* Y base */
        str_idx(&n_match, idx, asym, "C6 ", num);
        if (n_match != 1)
            nrerror("There should be only one C6 atom for pyrimidine!");
        posYnum = idx[1];

        str_idx(&n_match, idx, asym, "N3 ", num);
        if (n_match != 1)
            nrerror("There should be only one N3 atom for pyrimidine!");
        negYnum = idx[1];
    }

    /* Get y-axis & origin of the base */
    y_axis = dvector(1, 3);
    tmp_vec = dvector(1, 3);
    b_org = dvector(1, 3);

    for (i = 1; i <= 3; i++) {
        tmp_vec[i] = xyz[posYnum][i] - xyz[negYnum][i];
        b_org[i] = 0.5 * (xyz[posYnum][i] + xyz[negYnum][i]);
    }
    norm_dvector(y_axis, tmp_vec, 3);

    /* % Find the atom list to define the base normal
       %     Excluding Hs, sugar-phosphate backbone */
    base_atm = lvector(1, BASE_ATOM_NUM);

    for (i = 1; i <= num; i++) {
        if ((asym[i][0] != 'P') &&
            (asym[i][2] != 'P') && (asym[i][0] != 'H') && (asym[i][2] != '\'')) {
            n_base++;
            base_atm[n_base] = i;
        }
    }

    /* Get the base-pair normal & z-, x-axes */
    base_xyz = dmatrix(1, n_base, 1, 3);
    for (i = 1; i <= n_base; i++)
        for (j = 1; j <= 3; j++)
            base_xyz[i][j] = xyz[base_atm[i]][j];

    b_normal = dvector(1, 3);
    z_axis = dvector(1, 3);
    x_axis = dvector(1, 3);
    pnts_dist = dvector(1, n_base);

    lsplane(b_normal, &dval, pnts_dist, base_xyz, n_base);
    vec_orth(z_axis, b_normal, y_axis);
    cross_dvector(x_axis, y_axis, z_axis);

    /* % Check if +x-axis pointing from minor-groove to major-groove side
       %   by using its dot product with N1-->C6 vector.
       %   If NOT, reverse the sign of x- and z-axes */

    str_idx(&n_match, idx, asym, "N1 ", num);
    if (n_match != 1)
        nrerror("There are should only one N1 atom for a base!!!");
    N1 = idx[1];

    str_idx(&n_match, idx, asym, "C6 ", num);
    if (n_match != 1)
        nrerror("There are should only one C6 atom for a base!!!");
    C6 = idx[1];

    for (i = 1; i <= 3; i++)
        tmp_vec[i] = xyz[C6][i] - xyz[N1][i];
    dval = dot_dvector(x_axis, tmp_vec, 3);

    if (dval < 0.0)
        for (i = 1; i <= 3; i++) {
            x_axis[i] = -x_axis[i];
            z_axis[i] = -z_axis[i];
        }

    /* Set the base orientation */
    xyz_axis = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++) {
        xyz_axis[i][1] = x_axis[i];
        xyz_axis[i][2] = y_axis[i];
        xyz_axis[i][3] = z_axis[i];
    }

    for (i = 1; i <= num; i++)
        for (j = 1; j <= 3; j++)
            xyz[i][j] -= b_org[j];

    xyz2 = dmatrix(1, num, 1, 3);

    mul_dmatrix(xyz2, xyz, xyz_axis, num, 3, 3, 3);

    /* Write out the final structure in PDB format */
    wrtpdb(num, asym, btype, strand, bnum, xyz2, outfil);

    free_cvector(btype, 1, BASE_ATOM_NUM);
    free_cvector(strand, 1, BASE_ATOM_NUM);
    free_lvector(bnum, 1, BASE_ATOM_NUM);
    free_lvector(idx, 1, BASE_ATOM_NUM);
    free_dvector(y_axis, 1, 3);
    free_dvector(tmp_vec, 1, 3);
    free_dvector(b_org, 1, 3);
    free_lvector(base_atm, 1, BASE_ATOM_NUM);
    free_dvector(b_normal, 1, 3);
    free_dvector(z_axis, 1, 3);
    free_dvector(x_axis, 1, 3);
    free_dvector(pnts_dist, 1, n_base);
    free_cmatrix(asym, 1, BASE_ATOM_NUM, 1, 4);
    free_dmatrix(xyz, 1, BASE_ATOM_NUM, 1, 3);
    free_dmatrix(base_xyz, 1, n_base, 1, 3);
    free_dmatrix(xyz_axis, 1, 3, 1, 3);
    free_dmatrix(xyz2, 1, num, 1, 3);

    return 0;
}
