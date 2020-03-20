static vector_t *ref_genome_read(vector_t ref_genome, int ref_genome_length, vector_t *var_list)
{
    //vector_t *var_list = vector_create(8, VARIANT_T);
    int pos = 0;

    while (pos < ref_genome_length)
    {
        if(ref_genome[pos] == 'C' || ref_genome[pos] == 'G')
        {
            char alt;
            if(ref_genome[pos] == 'C')
            {
                alt = 'A';
            }
            else if (ref_genome[pos] == 'G')
            {
                alt = 'T';
            }
            
            variant_t *v = variant_create(var_list->data[0], pos, ref_genome[pos], alt);
            vector_add(var_list, v);
        }
        ++pos;
    }
    qsort(var_list->data, var_list->len, sizeof(void *), nat_sort_variant);//sort it
    return var_list;
}