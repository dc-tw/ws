#include <iostream>
#include <vector>

using namespace std;
void derive_combo(vector<vector<int> > *combo, vector<int> *prev, int n) { // Derive the combinations in k+1 that contain the previous elements
    if (prev->size() + 1 >= n) return;

    int k = prev->size() + 1;
    int i, c[k];

	vector<int> tmp = *prev;
    for (i = 0; i < prev->size(); i++)c[i] = tmp[i]; // first combination

    c[prev->size()] = c[prev->size() - 1] + 1;

    while (c[prev->size()] < n) { // generate and record combinations
        //for (i = 0; i < k; i++) { fprintf(stderr, "%d;", c[i]); } fprintf(stderr, "\n");
        vector<int > v(k,0);
        for (i = 0; i < k; i++) v[i] += c[i];
        combo->push_back(v);
        c[prev->size()]++;
    }
    //int ii, jj; for (ii = 0; ii < combo->size; ii++) { vector_int_t **c = (vector_int_t **)combo->data; fprintf(stderr, "%d\t", (int)ii); for (jj = 0; jj < c[ii]->size; jj++) { fprintf(stderr, "%d;", c[ii]->data[jj]); } fprintf(stderr, "\n"); } fprintf(stderr, "\n");
}


int main() {
    vector<vector<int> > c;//line 530 default as 8­Ó 0s
    vector<int> z(1, 0);
    for(int i=0; i<8; ++i)c.push_back(z);
    // c is the new input for calc_likelihood
    vector<int> prev;
    for(int i=0; i<6;++i)prev.push_back(i);
    //vector_t *c = vector_create(8, VOID_T);
	derive_combo(&c, &prev, 6);
    int i = c.size(), j = c[0].size(), m, n;
    for(m=0; m<i; m++){
    	for(n=0; n<j; n++){
    		cout<<c[m][n]<<"   ";
		}
		cout<<endl;
	}
    return 0;
}
