#include<stdio.h>
#include<stdlib.h>

int main(){
    return 0;
}

static void new_calc_likelihood(stats_t *stat, vector_t *var_set, 
                            const char *refseq, const int refseq_length, 
                            read_t **read_data, const int nreads, int readi, int *seqnt_map) {
    size_t i;//readi
    stat->ref = 0;
    stat->alt = 0;
    stat->het = 0;
    stat->ref_count = 0;
    stat->alt_count = 0;
    stat->seen = 0;

    variant_t **var_data = (variant_t **)var_set->data;
    double log_nv = log((double)var_set->len);

    int has_indel = 0;

    /* Alternative sequence */
    int altseq_length = 0;
    char *altseq = NULL;
    if (has_indel || dp) altseq = construct_altseq(refseq, refseq_length, stat->combo, var_data, &altseq_length); 

    /* Aligned reads */
    //for (readi = 0; readi < nreads; readi++) {
    /*---------------------------------------------------------------------------------------------------------*/
    /*要做一段選取variant, from start to (end+read_length)*/
    /*---------------------------------------------------------------------------------------------------------*/
//for each position
    /*每回合先重製prob_matrix*/
    double is_match[read_data[readi]->length], no_match[read_data[readi]->length];
    for (i = 0; i < read_data[readi]->length; i++) {
        is_match[i] = p_match[read_data[readi]->qual[i]];
        no_match[i] = p_mismatch[read_data[readi]->qual[i]];
        if (dp) {
            double n = is_match[i] - 1;
            is_match[i] += 1 - n;
            no_match[i] += 1 - n;
        }
    }
    /* Read probability matrix */
    double readprobmatrix[NT_CODES * read_data[readi]->length];
    set_prob_matrix(readprobmatrix, read_data[readi], is_match, no_match, seqnt_map, bisulfite);

        /* Outside Paralog Exact Formuation: Probability that read is from an outside the reference paralogous "elsewhere", f in F.  Approximate the bulk of probability distribution P(r|f):
           a) perfect match = prod[ (1-e) ]
           b) hamming/edit distance 1 = prod[ (1-e) ] * sum[ (e/3) / (1-e) ]
           c) lengthfactor = alpha ^ (read length - expected read length). Length distribution, for reads with different lengths (hard clipped), where longer reads should have a relatively lower P(r|f):
        P(r|f) = (perfect + hamming_1) / lengthfactor */
    double delta[read_data[readi]->length];
    for (i = 0; i < read_data[readi]->length; i++) delta[i] = no_match[i] - is_match[i];
    double a = sum_d(is_match, read_data[readi]->length);
    double elsewhere = log_add_exp(a, a + log_sum_exp(delta, read_data[readi]->length)) - (LGALPHA * (read_data[readi]->length - read_data[readi]->inferred_length));

    double prgu, prgv;
        //for (i =0; i < stat->combo->len; i++) { variant_t *v = var_data[stat->combo->data[i]]; printf("%d;%s;%s;", v->pos, v->ref, v->alt); }
        //printf("\t%s\t%d\t%d\t%s\n", read_data[readi]->name, read_data[readi]->pos, read_data[readi]->length, read_data[readi]->qseq);
    //依據起始位置選取variant set
    int i=0;
    while(i<var_set->len){
        if(var_set[i]==*(read_data[readi]+i)){
            picked_vector.push_back(i);
            改當前reference genome;
        }//match
        ++i;
    }
    //改readprobmatrix

    //計算score, new_calc_prob_dp下層要改只跑一次就好, 因為這層改為每個position跑一回
    /*---------------------------------------------------------------------------------------------------------*/
    prgu = new_calc_prob_dp(readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, gap_op, gap_ex, seqnt_map);
    prgv = new_calc_prob_dp(readprobmatrix, read_data[readi]->length, altseq, altseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, gap_op, gap_ex, seqnt_map);
    /*---------------------------------------------------------------------------------------------------------*/

    if (dp) {
        prgu = calc_prob_dp(readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, gap_op, gap_ex, seqnt_map);
        prgv = calc_prob_dp(readprobmatrix, read_data[readi]->length, altseq, altseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, gap_op, gap_ex, seqnt_map);
    }
    else if (has_indel) {
        prgu = calc_prob(readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
        prgv = calc_prob(readprobmatrix, read_data[readi]->length, altseq, altseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
    }
    else {
        calc_prob_snps(&prgu, &prgv, stat->combo, var_data, readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
    }//prg=probability genome, (u,v) stand for different strand
        //printf("%f\t%f\n\n", prgv, prgu);
    double pout = elsewhere;

    /* Multi-map alignments from XA tags: chr8,+42860367,97M3S,3;chr9,-44165038,100M,4; */
    //這邊不確定要怎麼改 因為這個好像跟下層dp計算偏移的理念不同
    if (read_data[readi]->multimapXA != NULL) {
        int xa_pos, n;
        char *s, xa_chr[strlen(read_data[readi]->multimapXA) + 1];
        for (s = read_data[readi]->multimapXA; sscanf(s, "%[^,],%d,%*[^;]%n", xa_chr, &xa_pos, &n) == 2; s += n + 1) {
            pout = log_add_exp(pout, elsewhere); // the more multi-mapped, the more likely it is the read is from elsewhere (paralogous), hence it scales (multiplied) with the number of multi-mapped locations
            if (strcmp(xa_chr, read_data[readi]->chr) != 0 && abs(xa_pos - read_data[readi]->pos) < read_data[readi]->length) { // if secondary alignment does not overlap primary aligment
                fasta_t *f = refseq_fetch(xa_chr, fa_file);
                if (f == NULL) continue;
                char *xa_refseq = f->seq;
                int xa_refseq_length = f->seq_length;

                double *p_readprobmatrix = readprobmatrix;
                double *newreadprobmatrix = NULL;
                if ((xa_pos < 0 && !read_data[readi]->is_reverse) || (xa_pos > 0 && read_data[readi]->is_reverse)) { // opposite of primary alignment strand
                    newreadprobmatrix = reverse(readprobmatrix, read_data[readi]->length * NT_CODES);
                    p_readprobmatrix = newreadprobmatrix;
                }

                xa_pos = abs(xa_pos);
                double readprobability = calc_prob(p_readprobmatrix, read_data[readi]->length, xa_refseq, xa_refseq_length, xa_pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
                prgu = log_add_exp(prgu, readprobability);
                prgv = log_add_exp(prgv, readprobability);
                free(newreadprobmatrix); newreadprobmatrix = NULL;
            }
            if (*(s + n) != ';') break;
        }
    }
    else if (read_data[readi]->multimapNH > 1) { // scale by the number of multimap positions
        double n = log(read_data[readi]->multimapNH - 1);
        double readprobability = prgu + n;
        pout = log_add_exp(pout, elsewhere + n);
        prgu = log_add_exp(prgu, readprobability);
        prgv = log_add_exp(prgv, readprobability);
    }

    /* Mixture model: probability that the read is from elsewhere, outside paralogous source */
    pout += lgomega;
    prgu = log_add_exp(pout, prgu);
    prgv = log_add_exp(pout, prgv);

    /* Track combination with highest variant likelihood */
    //這邊要改為紀錄最大值
    /*if (prgv > read_data[readi]->prgv) {
        read_data[readi]->index = seti;
        read_data[readi]->prgu = (float)prgu;
        read_data[readi]->prgv = (float)prgv;
        read_data[readi]->pout = (float)pout;
    }*/

    /* Mixture model: heterozygosity or heterogeneity as explicit allele frequency mu such that P(r|GuGv) = (mu)(P(r|Gv)) + (1-mu)(P(r|Gu)) */
    double phet   = log_add_exp(LOG50 + prgv, LOG50 + prgu);
    double phet10 = log_add_exp(LOG10 + prgv, LOG90 + prgu);
    double phet90 = log_add_exp(LOG90 + prgv, LOG10 + prgu);
    if (phet10 > phet) phet = phet10;
    if (phet90 > phet) phet = phet90;

    /* Priors */
    prgu += ref_prior;
    prgv += alt_prior - log_nv;
    phet += het_prior - log_nv;
    stat->ref += prgu;
    stat->alt += prgv;
    stat->het += phet;

    vector_double_add(stat->read_prgv, log_add_exp(prgv, phet));

    /* Read count incremented only when the difference in probability is not ambiguous, > ~log(2) difference and more likely than pout */
    if (prgv > prgu && prgv - prgu > 0.69 && prgv - pout > 0.69) stat->alt_count += 1;
    else if (prgu > prgv && prgu - prgv > 0.69 && prgu - pout > 0.69) stat->ref_count += 1;
    //拿掉一段debug
    //}, related with for (readi = 0; readi < nreads; readi++) {
    stat->mut = log_add_exp(stat->alt, stat->het);
    free(altseq); altseq = NULL;
}