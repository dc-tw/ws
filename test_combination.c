#include<stdio.h>
#include<stdlib.h>
typedef struct {
    size_t len, size;
    int  **data;
    //enum type type;
} vector_t;
typedef struct {
    size_t len, size;
    int *data;
} vector_int_t;

void vector_int_init(vector_int_t *a, int initial_size) {
    a->len = 0;
    a->size = initial_size;
    int tmp = initial_size*sizeof(int);
    a->data = malloc(tmp);
}
vector_int_t *vector_int_create(int initial_size) {
    vector_int_t *a = malloc(sizeof (vector_int_t));
    //vector_int_t *a;
    vector_int_init(a, initial_size);
    return a;
}
void vector_int_add(vector_int_t *a, int entry) {
    if (a->len >= a->size) {
        a->size *= 2;
        //int *p;
        int *p = realloc(a->data, a->size * sizeof (int));
        if (p == NULL) { //exit_err("failed to realloc in vector_add\n"); 
		}
        else { a->data = p; }
    }
    a->data[a->len++] = entry;
}

void vector_add(vector_t *a, void *entry) {
	if (a->len >= a->size) {
        a->size *= 2;
        //void **p;
        void **p = realloc(a->data, a->size * sizeof (void *));
        if (p == NULL) { //exit_err("failed to realloc in vector_add\n"); 
		}
        else { a->data = *p; }
    }
    printf("123");
    a->data[a->len++]= entry;
}

void combinations(vector_t *combo, int k, int n) {
    int i, c[k];
    for (i = 0; i < k; i++) c[i] = i; // first combination
    while (1) { // while (next_comb(c, k, n)) {
        // record the combination
        vector_int_t *v = vector_int_create(k);
        for (i = 0; i < k; i++) vector_int_add(v, c[i]);
        vector_add(combo, v);
        i = k - 1;
        c[i]++;
        while ((i >= 0 && i < k) && (c[i] >= n - k + 1 + i)) {
            i--;
            c[i]++;
        }
        /* Combination (n-k, n-k+1, ..., n) reached. No more combinations can be generated */
        if (c[0] > n - k) break; // return 0;
        /* c now looks like (..., x, n, n, n, ..., n), turn it into (..., x, x + 1, x + 2, ...) */
        for (i = i + 1; i < k; i++) c[i] = c[i - 1] + 1;
        // return 1;
    }
}


int main(void){
	vector_t *combo;// = vector_create(5, 0);//n=4
	printf("123");
    combinations(combo, 1, 2);
    int i;
    printf("123");
    /*for(i=0; i<2; i++){
    	printf("%d", combo->data[i]);
	}*/
	return 0;
}


