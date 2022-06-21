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

/*
   % PDB2INP processes an input PDb file to generate a simplified data
   % file as input to SCHNAaP. In Unix, you may simply use:
   %
   % grep '^ATOM' PDB_FILE | cut -c14-17,20,23-54 > INP_FILE
   %
   %       where PDB_FILE is the original PDB file and the INP_FILE
   %       is the input file to SCHNAaP.
   %
   % Note:
   %   (1) The original PDB file should be in Brookhaven PDB format, where
   %       residues on either strand are put along its 5'-->3' direction as
   %       showed schematically below:
   %
   %         Strand I residues          Strand II residues
   %       5'----------------->3' then 5'---------------->3'
   %         Required order of residues in the PDB file
   %
   %   (2) The two strands should be anti-parallel and of equal length
   %
   % For other complicated cases, please refer to the utility program
   % <sr2inp> for how to extract the paired residues and put them in a
   % suitable format for SCHNAaP to handle. Occasionally, you may need
   % to write a specific program to transform uncommon data format.
   %
   % Program written by Xiang-Jun Lu <xiangjun@rutchem.rutgers.edu>
   */

#include "schna_ar_p.h"

int main(long argc, char *argv[])
{
    long bnum, i;
    char asym[4], btype, inpfil[80], outfil[80], str[80], tmp[80];
    double x, y, z;
    FILE *finp, *fout;

    if (argc == 1) {
        printf("Name of the original PDB data file: ");
        read_stdin_str(inpfil);
        printf("Name of the input file to SCHNAaP: ");
        read_stdin_str(outfil);
    } else if (argc == 2) {
        strcpy(inpfil, argv[1]);
        printf("Name of the input file to SCHNAaP: ");
        read_stdin_str(outfil);
    } else {
        strcpy(inpfil, argv[1]);
        strcpy(outfil, argv[2]);
    }

    if ((finp = fopen(inpfil, "r")) == NULL) {
        printf("Failed to open file: %s\n", inpfil);
        exit(1);
    }

    if ((fout = fopen(outfil, "w")) == NULL) {
        printf("Failed to open file: %s\n", inpfil);
        exit(1);
    }

    while (fgets(str, sizeof(str), finp) != NULL) {
        for (i = 0; i < (long) strlen(str); i++)
            str[i] = toupper(str[i]);
        if (strncmp(str, "END", 3) == 0)
            break;
        if (strncmp(str, "ATOM", 4) == 0) {
            strncpy(asym, str + 13, 3);
            asym[3] = '\0';
            btype = str[19];
            tmp[0] = '\0';
            bnum = atoi(strncat(tmp, str + 23, 3));
            tmp[0] = '\0';
            x = atof(strncat(tmp, str + 30, 8));
            tmp[0] = '\0';
            y = atof(strncat(tmp, str + 38, 8));
            tmp[0] = '\0';
            z = atof(strncat(tmp, str + 46, 8));
            fprintf(fout, "%3s %c %3ld %8.3f %8.3f %8.3f\n", asym, btype, bnum, x, y, z);
        }
    }

    fclose(finp);
    fclose(fout);

    return 0;
}
