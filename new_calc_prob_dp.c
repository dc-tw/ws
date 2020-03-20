double new_smith_waterman_gotoh(const double *matrix, int read_length, const char *seq, int seq_length, 
                            int start, int end, int gap_op, int gap_ex, int *seqnt_map) { /* short in long version */
    int i, j;

    double prev[read_length + 1], curr[read_length + 1];
    double a_gap_curr[read_length + 1];
    double b_gap_prev[read_length + 1], b_gap_curr[read_length + 1];

    for (j = 0; j < read_length + 1; j++) prev[j] = 0;
    for (j = 0; j < read_length + 1; j++) b_gap_prev[j] = 0;

    double max_score = 0;
    //for (i = start; i < end; i++) {//這個是用來給不同起始點, 拉到calc_likelihood
        double row_max = 0;
        double upleft, open, extend;

        curr[0] = 0;
        a_gap_curr[0] = 0;
        b_gap_curr[0] = 0;
        for (j = 1; j <= read_length; j++) {
            int c = seq[i] - 'A';
            if (c < 0 || c > 57 || (c > 25 && c < 32)) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], i, seq_length); }

            upleft = prev[j - 1] + matrix[read_length * seqnt_map[c] + (j - 1)];

            open = curr[j - 1] - gap_op;
            extend = a_gap_curr[j - 1] - gap_ex;
            a_gap_curr[j] = (open >= extend) ? open : extend;

            open = prev[j] - gap_op;
            extend = b_gap_prev[j] - gap_ex;
            b_gap_curr[j] = (open >= extend) ? open : extend;

            curr[j] = upleft;
            if (a_gap_curr[j] >= curr[j]) curr[j] = a_gap_curr[j];
            if (b_gap_curr[j] >= curr[j]) curr[j] = b_gap_curr[j];
            if (curr[j] > row_max) row_max = curr[j];
        }
        if (row_max > max_score) max_score = row_max;

        memcpy(prev, curr, sizeof (prev));
        memcpy(b_gap_prev, b_gap_curr, sizeof (b_gap_prev));
    //}
    return max_score;
}

double new_calc_prob_region_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int gap_op, int gap_ex, int *seqnt_map) {
    if (start < 0) start = 0;
    else if (start >= seq_length) start = seq_length - 1;
    end += read_length;
    if (end < 0) end = 0;
    else if (end >= seq_length) end = seq_length - 1;
    return smith_waterman_gotoh(matrix, read_length, seq, seq_length, start, end, gap_op, gap_ex, seqnt_map);
}

double new_calc_prob_dp(const double *matrix, int read_length, const char *seq, 
                    int seq_length, int pos, int *splice_pos, int *splice_offset, 
                    int n_splice, int gap_op, int gap_ex, int *seqnt_map) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int start = pos - (read_length / 2);
    int end = pos + (read_length / 2);

    int i, j;
    double probability = 0;
    if (n_splice == 0) {
        probability = new_calc_prob_region_dp(matrix, read_length, seq, seq_length, pos, start, end, gap_op, gap_ex, seqnt_map);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        for (i = 0; i <= n_splice; i++) {
            int r_len = (i < n_splice) ? splice_pos[i] - r_pos + 1 : read_length - r_pos;
            int n = r_len / 2;
            start = g_pos - n;
            end = g_pos + n;

            double *submatrix = malloc(NT_CODES * r_len * sizeof (double));
            for (j = 0; j < NT_CODES; j++) memcpy(&submatrix[r_len * j], &matrix[read_length * j + r_pos], r_len * sizeof (double));
            probability += new_calc_prob_region_dp(submatrix, r_len, seq, seq_length, g_pos, start, end, gap_op, gap_ex, seqnt_map);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
    return probability;
}