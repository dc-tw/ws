#include<iostream>
#include<stdlib.h>
#include<algorithm>
#include<vector>

using namespace std;

int main(){
	int variants[3][3];//{{1,3}{5,2},{8,2}};
	variants[0][0] = 1;
	variants[0][1] = 3;
	variants[0][2] = 10;
	variants[1][0] = 5;
	variants[1][1] = 2;
	variants[1][2] = 5;
	variants[2][0] = 8;
	variants[2][1] = 1;
	variants[2][2] = 3;
	for(int i=0; i<3; i++){
		cout<<variants[i][0]<<" "<<variants[i][1]<<endl;
	}
	int pick[3];
	int arr[10] = {1,2,3,4,1,2,3,4,1,2};
	for(int i=0; i<3; i++){
		if(arr[variants[i][0]] == variants[i][1])pick[i] = 1;
		else pick[i] = 0;
	}
	cout<<"------------"<<endl;
	cout<<pick[0]<<" "<<pick[1]<<" "<<pick[2]<<" ";
	return 0;
}
