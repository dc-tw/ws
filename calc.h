/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _calc_h_
#define _calc_h_

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "util.h"
//#include "bisulfite_heap.h"
#include "vector.h"
/*-------------------------------------*/
#include <time.h>
static time_t now;
static struct tm *time_info;
#define print_status(M, ...)     \
    time(&now);                  \
    time_info = localtime(&now); \
    fprintf(stderr, M, ##__VA_ARGS__);
/*-------------------------------------*/

#define NT_CODES 21   // Size of nucleotide code table

/* Mapping table */
int seqnt_map[58];

/* Fastq quality to probability table */
double p_match[50], p_mismatch[50];

void init_seqnt_map(int *seqnt_map);
void init_q2p_table(double *p_match, double *p_mismatch, int size);

void set_prob_matrix(double *matrix, const read_t *read, const double *is_match, const double *no_match, const int *seqnt_map, const int bisulfite);
double calc_read_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *seqnt_map);
double calc_prob_region(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map);
double calc_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_mp);
double smith_waterman_gotoh(const double *matrix, int read_length, const char *seq, int seq_length, int start, int end, int gap_op, int gap_ex, int *seqnt_map);
double calc_prob_region_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int gap_op, int gap_ex, int *seqnt_map);
double calc_prob_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int gap_op, int gap_ex, int *seqnt_map);
void calc_prob_snps_region(double *prgu, double *prgv, vector_int_t *combo, variant_t **var_data, double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map);
void calc_prob_snps(double *prgu, double *prgv, vector_int_t *combo, variant_t **var_data, double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_map);

void set_prob_matrix_bisulfite(double *matrix, const read_t *read, const double *is_match, const double *no_match, const int *seqnt_map, const int bisulfite);
void set_prob_matrix_bisulfite2(double *matrix, const read_t *read, const double *is_match, const double *no_match, const int *seqnt_map, const int bisulfite);
double calc_read_prob_bisulfite(vector_t *var_set, double *matrix, int read_length, const char *seq, int seq_length, int pos, int *seqnt_map, int alt);
double calc_prob_region_bisulfite(vector_t *var_set, double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map, int alt);
double calc_prob_bisulfite(vector_t *var_set, double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_map, int alt);

#endif
