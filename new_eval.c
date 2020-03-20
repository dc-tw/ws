typedef struct {
    vector_int_t *combo;
    vector_double_t *read_prgv;
    double ref, alt, het, mut;
    int ref_count, alt_count, seen;
} stats_t;

static char *new_evaluate(vector_t *var_set) {
    size_t i, readi, seti;

    variant_t **var_data = (variant_t **)var_set->data;

    /* Reference sequence */
    fasta_t *f = refseq_fetch(var_data[0]->chr, fa_file);
    if (f == NULL)
        return NULL;
    char *refseq = f->seq;
    int refseq_length = f->seq_length;

    /* Reads in variant region coordinates */
    vector_t *read_list = bam_fetch(bam_file, var_data[0]->chr, var_data[0]->pos, var_data[var_set->len - 1]->pos);
    if (read_list->len == 0)
    {
        vector_destroy(read_list);
        free(read_list);
        read_list = NULL;
        return NULL;
    }
    read_t **read_data = (read_t **)read_list->data;

    /* Variant combinations as a vector of vectors */
    //vector_t *combo = powerset(var_set->len, maxh);

    vector_t *stats = vector_create(var_set->len + 1, STATS_T);//result storage
    for(readi = 0; readi<read_list->len; ++readi){//for each read
        new_calc_likelihood(s, var_set, refseq, refseq_length, read_data, read_list->len, readi, seqnt_map);//Bam file=read_list
    }
    /*
    for (seti = 0; seti < combo->len; seti++) { // Print combinations
        fprintf(stderr, "%d\t", (int)seti); 
        for (i = 0; i < ((vector_int_t *)combo->data[seti])->len; i++) { fprintf(stderr, "%d;", ((vector_int_t *)combo->data[seti])->data[i]); } fprintf(stderr, "\t"); 
        for (i = 0; i < ((vector_int_t *)combo->data[seti])->len; i++) { variant_t *v = var_data[((vector_int_t *)combo->data[seti])->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", v->chr, v->pos, v->ref, v->alt); } fprintf(stderr, "\n"); 
    }
    */

    //vector_t *stats = vector_create(var_set->len + 1, STATS_T);

    /*for (seti = 0; seti < combo->len; seti++)
    { // all, singles
        stats_t *s = stats_create((vector_int_t *)combo->data[seti], read_list->len);
        vector_add(stats, s);
        calc_likelihood(s, var_set, refseq, refseq_length, read_data, read_list->len, seti, seqnt_map);
    }
    if (var_set->len > 1)
    { // doubles and beyond
        heap_t *h = heap_create(STATS_T);
        for (seti = 1; seti < combo->len; seti++)
            heap_push(h, ((stats_t *)stats->data[seti])->mut, stats->data[seti]);

        stats_t *s;
        while (s = heap_pop(h), s != NULL)
        {
            if ((int)stats->len - var_set->len - 1 >= maxh)
                break;

            vector_t *c = vector_create(8, VOID_T);
            derive_combo(c, s->combo, var_set->len);
            for (i = 0; i < c->len; i++)
            {
                stats_t *s = stats_create((vector_int_t *)c->data[i], read_list->len);
                vector_add(stats, s);
                calc_likelihood(s, var_set, refseq, refseq_length, read_data, read_list->len, stats->len - 1, seqnt_map);
                heap_push(h, s->mut, s);
            }
            vector_free(c); //combos in stat so don't destroy
        }
        heap_free(h);
    }
    vector_free(combo); //combos in stat so don't destroy*/

    stats_t **stat = (stats_t **)stats->data;

    /* Heterozygous non-reference haplotypes as mixture model hypotheses */
    int c[stats->len];
    memset(c, 0, sizeof(c));
    for (readi = 0; readi < read_list->len; readi++)
        c[read_data[readi]->index]++; // combinations, based on best combination in each read

    vector_int_t *haplotypes = vector_int_create(stats->len);
    for (i = 0; i < stats->len; i++)
    {
        if ((double)c[i] / (double)read_list->len >= 0.1)
            vector_int_add(haplotypes, i); // relevant combination if read count >= 10% of reads seen
    }
    combo = vector_create(haplotypes->len, VOID_T);
    if (haplotypes->len > 1)
        combinations(combo, 2, haplotypes->len); // combination pairs

    vector_double_t *prhap = vector_double_create(combo->len);
    for (seti = 0; seti < combo->len; seti++)
    { // mixture model probabilities of combination pairs
        int x = haplotypes->data[((vector_int_t *)combo->data[seti])->data[0]];
        int y = haplotypes->data[((vector_int_t *)combo->data[seti])->data[1]];
        vector_double_add(prhap, 0);
        for (readi = 0; readi < read_list->len; readi++)
        {
            if (stat[x]->read_prgv->data[readi] == -DBL_MAX && stat[y]->read_prgv->data[readi] == -DBL_MAX)
                continue;
            double phet = log_add_exp(LOG50 + stat[x]->read_prgv->data[readi], LOG50 + stat[y]->read_prgv->data[readi]);
            double phet10 = log_add_exp(LOG10 + stat[x]->read_prgv->data[readi], LOG90 + stat[y]->read_prgv->data[readi]);
            double phet90 = log_add_exp(LOG90 + stat[x]->read_prgv->data[readi], LOG10 + stat[y]->read_prgv->data[readi]);
            if (phet10 > phet)
                phet = phet10;
            if (phet90 > phet)
                phet = phet90;
            prhap->data[seti] += phet; // equal prior probability to ref since this assumes heterozygous non-reference variant
        }
    }
    if (debug >= 1)
    {
        for (seti = 0; seti < combo->len; seti++)
        {
            int x = haplotypes->data[((vector_int_t *)combo->data[seti])->data[0]];
            int y = haplotypes->data[((vector_int_t *)combo->data[seti])->data[1]];
            fprintf(stderr, "==\t%d, %d, %f\n", x, y, prhap->data[seti]);
        }
    }

    double total = log_add_exp(stat[0]->mut, stat[0]->ref);
    for (seti = 1; seti < stats->len; seti++)
    {
        total = log_add_exp(total, stat[seti]->mut);
        total = log_add_exp(total, stat[seti]->ref);
    }
    for (seti = 0; seti < combo->len; seti++)
        total = log_add_exp(total, prhap->data[seti]);

    char *output = malloc(sizeof(*output));
    output[0] = '\0';
    if (mvh)
    { /* Max likelihood variant hypothesis */
        size_t max_seti = 0;
        double r = stat[0]->mut - stat[0]->ref;
        double has_alt = stat[0]->mut;
        for (seti = 1; seti < stats->len; seti++)
        {
            if (stat[seti]->mut - stat[seti]->ref > r)
            {
                r = stat[seti]->mut - stat[seti]->ref;
                has_alt = stat[seti]->mut;
                max_seti = seti;
            }
        }
        vector_t *v = vector_create(var_set->len, VARIANT_T);
        for (i = 0; i < stat[max_seti]->combo->len; i++)
            vector_add(v, var_data[stat[max_seti]->combo->data[i]]);
        variant_print(&output, v, 0, stat[max_seti]->seen, stat[max_seti]->ref_count, stat[max_seti]->alt_count, log_add_exp(total, stat[max_seti]->ref), has_alt, stat[max_seti]->ref);
        vector_free(v); //variants in var_list so don't destroy
    }
    else
    { /* Marginal probabilities & likelihood ratios*/
        for (i = 0; i < var_set->len; i++)
        {
            double has_alt = 0;
            double not_alt = 0;
            int acount = -1;
            int rcount = -1;
            int seen = -1;
            for (seti = 0; seti < stats->len; seti++)
            {
                if (variant_find(stat[seti]->combo, i) != -1)
                { // if variant is in this combination
                    has_alt = (has_alt == 0) ? stat[seti]->mut : log_add_exp(has_alt, stat[seti]->mut);
                    not_alt = (not_alt == 0) ? stat[seti]->ref : log_add_exp(not_alt, stat[seti]->ref);
                    if (stat[seti]->seen > seen)
                        seen = stat[seti]->seen;
                    if (stat[seti]->alt_count > acount)
                    {
                        acount = stat[seti]->alt_count;
                        rcount = stat[seti]->ref_count;
                    }
                }
                else
                {
                    not_alt = (not_alt == 0) ? stat[seti]->mut : log_add_exp(not_alt, stat[seti]->mut);
                }
            }
            for (seti = 0; seti < combo->len; seti++)
            {
                int x = haplotypes->data[((vector_int_t *)combo->data[seti])->data[0]];
                int y = haplotypes->data[((vector_int_t *)combo->data[seti])->data[1]];
                if (variant_find(stat[x]->combo, i) != -1 || variant_find(stat[y]->combo, i) != -1)
                    has_alt = log_add_exp(has_alt, prhap->data[seti]);
                else
                    not_alt = log_add_exp(not_alt, prhap->data[seti]);
            }
            variant_print(&output, var_set, i, seen, rcount, acount, total, has_alt, not_alt);
        }
    }

    if (verbose)
    {
        for (readi = 0; readi < read_list->len; readi++)
        {
            if (read_data[readi]->prgu == 0 && read_data[readi]->prgv == 0 && read_data[readi]->pout == 0)
                continue; // unprocessed read
            flockfile(stderr);
            fprintf(stderr, "%s\t%s\t%d\t", read_data[readi]->name, read_data[readi]->chr, read_data[readi]->pos);
            fprintf(stderr, "%f\t%f\t%f\t", read_data[readi]->prgu, read_data[readi]->prgv, read_data[readi]->pout);
            for (i = 0; i < read_data[readi]->n_cigar; i++)
                fprintf(stderr, "%d%c", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
            fprintf(stderr, "\t");
            if (read_data[readi]->multimapXA != NULL)
                fprintf(stderr, "%s\t", read_data[readi]->multimapXA);
            else
                fprintf(stderr, "%d\t", read_data[readi]->multimapNH);
            if (read_data[readi]->flag != NULL)
                fprintf(stderr, "%s\t", read_data[readi]->flag);
            else
                fprintf(stderr, "NONE\t");
            fprintf(stderr, "[");
            for (i = 0; i < stat[read_data[readi]->index]->combo->len; i++)
            {
                variant_t *v = var_data[stat[read_data[readi]->index]->combo->data[i]];
                fprintf(stderr, "%s,%d,%s,%s;", v->chr, v->pos, v->ref, v->alt);
            }
            fprintf(stderr, "]\n");
            funlockfile(stderr);
        }
    }

    for (i = 0; i < combo->len; i++)
        vector_int_free(combo->data[i]);
    vector_free(combo); //not destroyed because previously vector_int_free all elements
    vector_int_free(haplotypes);
    vector_double_free(prhap);
    vector_destroy(read_list);
    free(read_list);
    read_list = NULL;
    vector_destroy(stats);
    free(stats);
    stats = NULL;
    return output;
}