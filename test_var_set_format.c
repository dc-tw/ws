/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "vector.h"
#include "util.h"
#include "calc.h"
#include "heap.h"

/* Constants */
#define VERSION "1.1.1"
#define ALPHA 1.3 // Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source

/* Precalculated log values */
#define M_1_LOG10E (1.0 / M_LOG10E)
#define M_1_LN10 (1.0 / M_LN10)
#define LOG50 (log(0.5))
#define LOG10 (log(0.1))
#define LOG90 (log(0.9))
#define LGALPHA (log(ALPHA))

/* Command line arguments */
static int debug;
static char *vcf_file;
static char *bam_file;
static char *fa_file;
static char *out_file;
static int nthread;
static int sharedr;
static int distlim;
static int maxdist;
static int maxh;
static int mvh;
static int pao;
static int isc;
static int nodup;
static int splice;
static int verbose;
static int lowmem;
static int phred64;
static int bisulfite;
static int const_qual;
static double hetbias;
static double omega, lgomega;
static int dp, gap_op, gap_ex;
static int rc;
static double ref_prior, alt_prior, het_prior;

/* Time info */
static time_t now;
static struct tm *time_info;
#define print_status(M, ...)     \
    time(&now);                  \
    time_info = localtime(&now); \
    fprintf(stderr, M, ##__VA_ARGS__);

KHASH_MAP_INIT_STR(rsh, fasta_t *) // hashmap: string key, vector value
static khash_t(rsh) * refseq_hash; // pointer to hashmap
static pthread_mutex_t refseq_lock;

static int bam_fetch_last(const char *bam_file, const char *chr, const int pos1, const int pos2)
{
    /* Reads in variant j = i + 1 region coordinates */
    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL)
    {
        exit_err("failed to open BAM file %s\n", bam_file);
    }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0)
    {
        exit_err("bad header %s\n", bam_file);
    }
    hts_idx_t *bam_idx = sam_index_load(sam_in, bam_file); // bam index
    if (bam_idx == NULL)
    {
        exit_err("failed to open BAM index %s\n", bam_file);
    }

    int last = -1;
    int tid = bam_name2id(bam_header, chr);
    hts_itr_t *iter = sam_itr_queryi(bam_idx, tid, pos1 - 1, pos2); // read iterator
    if (iter != NULL)
    {
        bam1_t *aln = bam_init1(); // initialize an alignment
        while (sam_itr_next(sam_in, iter, aln) >= 0)
        {
            if (aln->core.tid < 0)
                continue; // not mapped
            last = aln->core.pos + aln->core.l_qseq;
        }
        bam_destroy1(aln);
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    return (last);
}

typedef struct
{
    vector_t *queue, *results;
    pthread_mutex_t q_lock;
    pthread_mutex_t r_lock;
    size_t len;
} work_t;

static void *pool(void *work)
{
    work_t *w = (work_t *)work;
    /*
    //int size = var_set->queue->size;
    print_status("len=%d\n", work->len);
    vector_t *queue = work->queue;
    int size = work->queue->size;
    while(!queue){
        //char *chr, int pos, char *ref, char *alt, from variant_create
        printf("chr=%c pos=%d ref=%c alt=%c\n",queue, queue+ sizeof(char), 
                                                queue+ sizeof(char)+ sizeof(int), 
                                                queue+ 2*sizeof(char)+ sizeof(int));
        queue += work->queue->size;
        //queue += 3*sizeof(char)+ sizeof(int);
    }
    */
    size_t n = w->len / 10;
    while (1)
    { //pthread_t ptid = pthread_self(); uint64_t threadid = 0; memcpy(&threadid, &ptid, min(sizeof (threadid), sizeof (ptid)));
        pthread_mutex_lock(&w->q_lock);
        vector_t *var_set = (vector_t *)vector_pop(w->queue);
        pthread_mutex_unlock(&w->q_lock);
        if (var_set == NULL)
            break;
        
        //char *outstr = evaluate(var_set);
        if (outstr != NULL)
        {
            pthread_mutex_lock(&w->r_lock);
            vector_add(w->results, outstr);
            pthread_mutex_unlock(&w->r_lock);
        }
        vector_free(var_set); //variants in var_list so don't destroy
    }
    return NULL;
}

static void process(const vector_t *var_list, FILE *out_fh)
{
    size_t i, j;

    variant_t **var_data = (variant_t **)var_list->data;

    i = 0;
    vector_t *var_set = vector_create(var_list->len, VOID_T);
    if (sharedr == 1)
    { /* Variants that share a read: shared with a given first variant */
        while (i < var_list->len)
        {
            vector_t *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);

            /* Reads in variant i region coordinates */
            int i_last = bam_fetch_last(bam_file, var_data[i]->chr, 
                                        var_data[i]->pos, var_data[i]->pos);

            j = i + 1;
            while (j < var_list->len && strcmp(var_data[i]->chr, var_data[j]->chr) == 0)
            { // while last read in i will reach j
                if (var_data[j]->pos > i_last)
                    break;
                vector_add(curr, var_data[j++]);
            }
            i = j;
            vector_add(var_set, curr);
        }//依據bam_fetch_last將variants分段放入var_set<--------------------------------------------------------------
    }
    else if (sharedr == 2)
    { /* Variants that share a read: shared with any neighboring variant */
        while (i < var_list->len)
        {
            vector_t *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);

            j = i + 1;
            while (j < var_list->len && strcmp(var_data[i]->chr, var_data[j]->chr) == 0)
            { // while last read in i will reach j
                /* Reads in variant i region coordinates */
                int i_last = bam_fetch_last(bam_file, var_data[i]->chr, var_data[i]->pos, var_data[i]->pos);
                if (var_data[j]->pos > i_last)
                    break;
                vector_add(curr, var_data[j]);
                i++;
                j++;
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    else
    { /* Variants that are close together as sets */
        while (i < var_list->len)
        {
            vector_t *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);
            j = i + 1;
            while (distlim > 0 && j < var_list->len && strcmp(var_data[j]->chr, var_data[j - 1]->chr) == 0 && abs(var_data[j]->pos - var_data[j - 1]->pos) <= distlim)
            {
                if (maxdist > 0 && abs(var_data[j]->pos - var_data[i]->pos) > maxdist)
                    break;
                vector_add(curr, var_data[j++]);
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    /* Heterozygous non-reference variants as separate entries */
    int flag_add = 1;
    while (flag_add)
    {
        flag_add = 0;
        for (i = 0; i < var_set->len; i++)
        {
            vector_t *curr_set = (vector_t *)var_set->data[i];
            if (curr_set->len == 1)
                continue;

            int flag_nonset = 1;
            for (j = 0; j < curr_set->len - 1; j++)
            { // check if all entries have the same position
                variant_t *curr = (variant_t *)curr_set->data[j];
                variant_t *next = (variant_t *)curr_set->data[j + 1];
                if (curr->pos == next->pos && strcmp(curr->chr, next->chr) == 0 && strcmp(curr->ref, next->ref) == 0 && strcmp(curr->alt, next->alt) == 0)
                    vector_del(curr_set, j + 1); // delete duplicate entries
                else if (curr->pos != next->pos)
                    flag_nonset = 0;
            }
            if (flag_nonset)
            { // only 1 entry, with multiple heterozygous non-reference variants
                while (curr_set->len > 1)
                {
                    variant_t *curr = (variant_t *)vector_pop(curr_set);
                    vector_t *dup = vector_create(8, VARIANT_T);
                    vector_add(dup, curr);
                    vector_add(var_set, dup);
                }
            }
            else
            { // multiple entries comprising a set
                for (j = 0; j < curr_set->len - 1; j++)
                {
                    variant_t *curr = (variant_t *)curr_set->data[j];
                    variant_t *next = (variant_t *)curr_set->data[j + 1];
                    if (curr->pos == next->pos)
                    {
                        flag_add = 1;
                        vector_t *dup = vector_dup(curr_set);
                        vector_del(curr_set, j);
                        vector_del(dup, j + 1);
                        vector_add(var_set, dup);
                    }
                }
            }
        }
    }
    if (sharedr == 1)
    {
        print_status("# Variants with shared reads to first in set: %i entries\t%s", (int)var_set->len, asctime(time_info));
    }
    else if (sharedr == 2)
    {
        print_status("# Variants with shared reads to any in set: %i entries\t%s", (int)var_set->len, asctime(time_info));
    }
    else
    {
        print_status("# Variants within %d (max window: %d) bp: %i entries\t%s", distlim, maxdist, (int)var_set->len, asctime(time_info));
    }

    print_status("# Options: maxh=%d mvh=%d pao=%d isc=%d nodup=%d splice=%d bs=%d lowmem=%d phred64=%d\n", maxh, mvh, pao, isc, nodup, splice, bisulfite, lowmem, phred64);
    print_status("#          dp=%d gap_op=%d gap_ex=%d\n", dp, gap_op, gap_ex);
    print_status("#          hetbias=%g omega=%g cq=%d\n", hetbias, omega, const_qual);
    print_status("#          verbose=%d\n", verbose);
    print_status("# Start: %d threads \t%s\t%s", nthread, bam_file, asctime(time_info));

    vector_t *queue = vector_create(var_set->len, VOID_T);
    vector_t *results = vector_create(var_set->len, VOID_T);
    for (i = 0; i < var_set->len; i++)
        vector_add(queue, var_set->data[i]);

    work_t *w = malloc(sizeof(work_t));
    w->queue = queue;
    w->results = results;
    w->len = var_set->len;

    pthread_mutex_init(&w->q_lock, NULL);
    pthread_mutex_init(&w->r_lock, NULL);

    pthread_t tid[nthread];
    for (i = 0; i < nthread; i++)
        pthread_create(&tid[i], NULL, pool, w);
    for (i = 0; i < nthread; i++)
        pthread_join(tid[i], NULL);

    pthread_mutex_destroy(&w->q_lock);
    pthread_mutex_destroy(&w->r_lock);

    free(w);
    w = NULL;
    vector_free(var_set); //variants in var_list so don't destroy

    qsort(results->data, results->len, sizeof(void *), nat_sort_vector);
    fprintf(out_fh, "# SEQ\tPOS\tREF\tALT\tReads\tRefReads\tAltReads\tProb\tOdds\tSet\n");
    for (i = 0; i < results->len; i++)
        fprintf(out_fh, "%s", (char *)results->data[i]);
    vector_destroy(queue);
    free(queue);
    queue = NULL;
    vector_destroy(results);
    free(results);
    results = NULL;
    print_status("# Done:\t%s\t%s", bam_file, asctime(time_info));
}

static void print_usage()
{
    printf("\nUsage: eagle [options] -v variants.vcf -a alignment.bam -r reference.fasta\n\n");
    printf("Required:\n");
    printf("  -v --vcf      FILE   Variants VCF file. [stdin]\n");
    printf("  -a --bam      FILE   Alignment data bam files, ref-coord sorted with bai index file.\n");
    printf("  -r --ref      FILE   Reference sequence, fasta file with fai index file.\n");
    printf("Options:\n");
    printf("  -o --out      FILE   Output file. [stdout]\n");
    printf("  -t --nthread  INT    Number of threads. [1]\n");
    printf("  -s --sharedr  INT    Group nearby variants that share a read, 0:distance based/off, 1:shared with first, 2:shared with any. [0]\n");
    printf("  -n --distlim  INT    Group nearby variants within n bases, 0:off. [10]\n");
    printf("  -w --maxdist  INT    Maximum number of bases between any two variants in a set of hypotheses, 0:off. [0]\n");
    printf("  -m --maxh     INT    Maximum number of combinations in the set of hypotheses, instead of all 2^n. [1024]\n");
    printf("     --mvh             Output the maximum likelihood hypothesis in the set instead of marginal probabilities.\n");
    printf("     --pao             Primary alignments only.\n");
    printf("     --isc             Ignore soft-clipped bases.\n");
    printf("     --nodup           Ignore marked duplicate reads (based on SAM flag).\n");
    printf("     --splice          RNA-seq spliced reads.\n");
    printf("     --bs       INT    Bisulfite treated reads. 0: off, 1: top/forward strand, 2: bottom/reverse strand, 3: both. [0]\n");
    printf("     --dp              Use dynamic programming to calculate likelihood instead of the basic model.\n");
    printf("     --gap_op   INT    DP gap open penalty. [6]. Recommend 2 for long reads with indel errors.\n");
    printf("     --gap_ex   INT    DP gap extend penalty. [1].\n");
    printf("     --verbose         Verbose mode, output likelihoods for each read seen for each hypothesis to stderr.\n");
    printf("     --lowmem          Low memory usage mode, the default mode for snps, this may be slightly slower for indels but uses less memory.\n");
    printf("     --phred64         Read quality scores are in phred64.\n");
    printf("     --hetbias  FLOAT  Prior probability bias towards non-homozygous mutations, between [0,1]. [0.5]\n");
    printf("     --omega    FLOAT  Prior probability of originating from outside paralogous source, between [0,1]. [1e-6]\n");
    printf("     --cq       INT    Constant quality as a phred score, ignoring the quality field in SAM. [0 is off]\n");
    printf("     --rc              Wrapper for read classification settings: --omega=1.0e-40 --isc --mvh --verbose --lowmem.\n");
    printf("     --version         Display version.\n");
}

int main(int argc, char **argv)
{
    /* Command line parameters defaults */
    debug = 0;
    vcf_file = NULL;
    bam_file = NULL;
    fa_file = NULL;
    out_file = NULL;
    nthread = 1;
    sharedr = 0;
    distlim = 10;
    maxdist = 0;
    maxh = 1024;
    mvh = 0;
    pao = 0;
    isc = 0;
    nodup = 0;
    splice = 0;
    bisulfite = 0;
    verbose = 0;
    lowmem = 0;
    phred64 = 0;
    dp = 0;
    gap_op = 6;
    gap_ex = 1;
    hetbias = 0.5;
    omega = 1.0e-6;
    const_qual = 0;
    rc = 0;

    static struct option long_options[] = {
        {"debug", optional_argument, NULL, 'd'},
        {"vcf", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"ref", required_argument, NULL, 'r'},
        {"out", optional_argument, NULL, 'o'},
        {"nthread", optional_argument, NULL, 't'},
        {"sharedr", optional_argument, NULL, 's'},
        {"distlim", optional_argument, NULL, 'n'},
        {"maxdist", optional_argument, NULL, 'w'},
        {"maxh", optional_argument, NULL, 'm'},
        {"maxh", optional_argument, NULL, 'm'},
        {"mvh", no_argument, &mvh, 1},
        {"pao", no_argument, &pao, 1},
        {"isc", no_argument, &isc, 1},
        {"nodup", no_argument, &nodup, 1},
        {"splice", no_argument, &splice, 1},
        {"verbose", no_argument, &verbose, 1},
        {"phred64", no_argument, &phred64, 1},
        {"lowmem", no_argument, &lowmem, 1},
        {"dp", no_argument, &dp, 1},
        {"gap_op", optional_argument, NULL, 981},
        {"gap_ex", optional_argument, NULL, 982},
        {"hetbias", optional_argument, NULL, 990},
        {"omega", optional_argument, NULL, 991},
        {"bs", optional_argument, NULL, 992},
        {"cq", optional_argument, NULL, 993},
        {"rc", no_argument, &rc, 1},
        {"version", optional_argument, NULL, 999},
        {0, 0, 0, 0}};

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "d:v:a:r:o:t:s:n:w:m:", long_options, &opt)) != -1)
    {
        switch (opt)
        {
        case 0:
            //if (long_options[option_index].flag != 0) break;
            break;
        case 'd':
            debug = parse_int(optarg);
            break;
        case 'v':
            vcf_file = optarg;
            break;
        case 'a':
            bam_file = optarg;
            break;
        case 'r':
            fa_file = optarg;
            break;
        case 'o':
            out_file = optarg;
            break;
        case 't':
            nthread = parse_int(optarg);
            break;
        case 's':
            sharedr = parse_int(optarg);
            break;
        case 'n':
            distlim = parse_int(optarg);
            break;
        case 'w':
            maxdist = parse_int(optarg);
            break;
        case 'm':
            maxh = parse_int(optarg);
            break;
        case 981:
            gap_op = parse_int(optarg);
            break;
        case 982:
            gap_ex = parse_int(optarg);
            break;
        case 990:
            hetbias = parse_double(optarg);
            break;
        case 991:
            omega = parse_double(optarg);
            break;
        case 992:
            bisulfite = parse_int(optarg);
            break;
        case 993:
            const_qual = parse_int(optarg);
            break;
        case 999:
            printf("EAGLE %s\n", VERSION);
            exit(0);
        default:
            exit_usage("Bad options");
        }
    }
    if (optind > argc)
    {
        exit_usage("Bad program call");
    }

    FILE *vcf_fh = stdin;
    if (vcf_file != NULL)
    { // default vcf file handle is stdin unless a vcf file option is used
        vcf_fh = fopen(vcf_file, "r");
        if (vcf_fh == NULL)
        {
            exit_err("failed to open VCF file %s\n", vcf_file);
        }
    }
    else
    {
        vcf_file = "stdin";
    }
    if (bam_file == NULL)
    {
        exit_usage("Missing alignments given as BAM file!");
    }
    if (fa_file == NULL)
    {
        exit_usage("Missing reference genome given as Fasta file!");
    }
    if (nthread < 1)
        nthread = 1;
    if (sharedr < 0 || sharedr > 2)
        sharedr = 0;
    if (distlim < 0)
        distlim = 10;
    if (maxdist < 0)
        maxdist = 0;
    if (maxh < 0)
        maxh = 0;
    if (gap_op <= 0)
        gap_op = 6;
    if (gap_ex <= 0)
        gap_ex = 1;
    if (hetbias < 0 || hetbias > 1)
        hetbias = 0.5;
    if (omega < 0 || omega > 1)
        omega = 1e-6;
    if (rc)
    {
        omega = 1e-40;
        isc = 1;
        mvh = 1;
        verbose = 1;
        lowmem = 1;
    }
    lgomega = (log(omega) - log(1.0 - omega));

    ref_prior = log(0.5);
    alt_prior = log(0.5 * (1 - hetbias));
    het_prior = log(0.5 * hetbias);

    FILE *out_fh = stdout;
    if (out_file != NULL)
        out_fh = fopen(out_file, "w"); // default output file handle is stdout unless output file option is used

    init_seqnt_map(seqnt_map);
    init_q2p_table(p_match, p_mismatch, 50);

    /* Start processing data */
    clock_t tic = clock();
    vector_t *var_list = vcf_read(vcf_fh);
    print_status("# Read VCF: %s\t%i entries\t%s", vcf_file, (int)var_list->len, asctime(time_info));
    /*---------*/
    //C read from ref-genome, input C to A
    int position=0, cur=0;
    while (position<read->length)
    {
        if(var_list->data[cur])
        ++position;
    }
    /*---------*/


    refseq_hash = kh_init(rsh);

    pthread_mutex_init(&refseq_lock, NULL);
    process(var_list, out_fh);
    if (out_file != NULL)
        fclose(out_fh);
    else
        fflush(stdout);
    pthread_mutex_destroy(&refseq_lock);

    khiter_t k;
    for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); k++)
    {
        if (kh_exist(refseq_hash, k))
        {
            fasta_destroy(kh_val(refseq_hash, k));
            free(kh_val(refseq_hash, k));
            kh_val(refseq_hash, k) = NULL;
        }
    }
    kh_destroy(rsh, refseq_hash);
    vector_destroy(var_list);
    free(var_list);
    var_list = NULL;

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
