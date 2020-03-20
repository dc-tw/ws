#include <iostream>
#include <vector>

using namespace std;


void combinations(vector<vector<int> > *combo, int k, int n) {
    int i, c[k];
    for (i = 0; i < k; i++) c[i] = i; // first combination
    while (1) { // while (next_comb(c, k, n)) {
        // record the combination
        vector<int> v(k, 0);
        for (i = 0; i < k; i++) v[i] = c[i];
        combo->push_back(v);
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

int main() {
    vector<vector<int> > combo;
    combinations(&combo, 1, 5);
    int i = combo.size(), j = combo[0].size(), m, n ;
    for(m=0; m<i; m++){
    	for(n=0; n<j; n++){
    		cout<<combo[m][n]<<"   ";
		}
		cout<<endl;
	}
    return 0;
}

/*-------------------------*/




