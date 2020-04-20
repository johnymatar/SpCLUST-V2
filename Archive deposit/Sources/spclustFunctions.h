#include <cstring>
#include <algorithm>
#include "GaussianMixture.h"
#include "eigen_jacobi.h"
using namespace std;

bool profilingOn=false; //Set to true to print the timestamp at the starting and ending of each computation phase
bool internalValidationOn=false; //Set to true to print the internal validation indices of the resulting clustering

void printTimestamp(string msg){
	if(profilingOn){
		// declaring argument of time() 
		time_t my_time = time(NULL); 
  
		// ctime() used to give the present time 
		cout<<msg<<" "<<ctime(&my_time); 
	}
}

char tolowerchar(char in){
	if(in>='A' && in <= 'Z')
		in = in + 'a' - 'A';
	return in;
}

string tolowercstr(char *in){
	string lowres="";
	for(int i=0; i<strlen(in); i++)
		lowres+=tolowerchar(in[i]);
	return lowres;
}

string tolowerstr(string in){
	string lowres="";
	for(int i=0; i<in.length(); i++)
		lowres+=tolowerchar(in[i]);
	return lowres;
}

double distanceS(string mat, char *ali1, char *ali2, double gapOpen = -10, double gapExtend = -0.5) {
	/*
	This function returns the distance between two aligned sequences. 
	This function is a subfunction of "matriceDistances".
	*/
	//Validating the choice of distance matrix
	/*if(mat!="EDNAFULL" && mat!="BLOSUM62" && mat!="PAM250"){
		cout<<"Invalid distance matrix. Available matrices are: EDNAFULL, BLOSUM62, and PAM250.\n";
		return 999;
	}*/
	
	int l1 = strlen(ali1), l2 = strlen(ali2), k, scoreMax, dist;
	double score = 0.0;
	
	//Nucleotides and proteins comparison scores matrices
	int EDNAFULL[225][3] = { { 'A', 'A', 5 },{ 'A', 'C', -4 },{ 'A', 'B', -4 },{ 'A', 'D', -1 },{ 'A', 'G', -4 },{ 'A', 'H', -1 },{ 'A', 'K', -4 },{ 'A', 'M', 1 },{ 'A', 'N', -2 },{ 'A', 'S', -4 },{ 'A', 'R', 1 },{ 'A', 'T', -4 },{ 'A', 'W', 1 },{ 'A', 'V', -1 },{ 'A', 'Y', -4 },{ 'C', 'A', -4 },{ 'C', 'C', 5 },{ 'C', 'B', -1 },{ 'C', 'D', -4 },{ 'C', 'G', -4 },{ 'C', 'H', -1 },{ 'C', 'K', -4 },{ 'C', 'M', 1 },{ 'C', 'N', -2 },{ 'C', 'S', 1 },{ 'C', 'R', -4 },{ 'C', 'T', -4 },{ 'C', 'W', -4 },{ 'C', 'V', -1 },{ 'C', 'Y', 1 },{ 'B', 'A', -4 },{ 'B', 'C', -1 },{ 'B', 'B', -1 },{ 'B', 'D', -2 },{ 'B', 'G', -1 },{ 'B', 'H', -2 },{ 'B', 'K', -1 },{ 'B', 'M', -3 },{ 'B', 'N', -1 },{ 'B', 'S', -1 },{ 'B', 'R', -3 },{ 'B', 'T', -1 },{ 'B', 'W', -3 },{ 'B', 'V', -2 },{ 'B', 'Y', -1 },{ 'D', 'A', -1 },{ 'D', 'C', -4 },{ 'D', 'B', -2 },{ 'D', 'D', -1 },{ 'D', 'G', -1 },{ 'D', 'H', -2 },{ 'D', 'K', -1 },{ 'D', 'M', -3 },{ 'D', 'N', -1 },{ 'D', 'S', -3 },{ 'D', 'R', -1 },{ 'D', 'T', -1 },{ 'D', 'W', -1 },{ 'D', 'V', -2 },{ 'D', 'Y', -3 },{ 'G', 'A', -4 },{ 'G', 'C', -4 },{ 'G', 'B', -1 },{ 'G', 'D', -1 },{ 'G', 'G', 5 },{ 'G', 'H', -4 },{ 'G', 'K', 1 },{ 'G', 'M', -4 },{ 'G', 'N', -2 },{ 'G', 'S', 1 },{ 'G', 'R', 1 },{ 'G', 'T', -4 },{ 'G', 'W', -4 },{ 'G', 'V', -1 },{ 'G', 'Y', -4 },{ 'H', 'A', -1 },{ 'H', 'C', -1 },{ 'H', 'B', -2 },{ 'H', 'D', -2 },{ 'H', 'G', -4 },{ 'H', 'H', -1 },{ 'H', 'K', -3 },{ 'H', 'M', -1 },{ 'H', 'N', -1 },{ 'H', 'S', -3 },{ 'H', 'R', -3 },{ 'H', 'T', -1 },{ 'H', 'W', -1 },{ 'H', 'V', -2 },{ 'H', 'Y', -1 },{ 'K', 'A', -4 },{ 'K', 'C', -4 },{ 'K', 'B', -1 },{ 'K', 'D', -1 },{ 'K', 'G', 1 },{ 'K', 'H', -3 },{ 'K', 'K', -1 },{ 'K', 'M', -4 },{ 'K', 'N', -1 },{ 'K', 'S', -2 },{ 'K', 'R', -2 },{ 'K', 'T', 1 },{ 'K', 'W', -2 },{ 'K', 'V', -3 },{ 'K', 'Y', -2 },{ 'M', 'A', 1 },{ 'M', 'C', 1 },{ 'M', 'B', -3 },{ 'M', 'D', -3 },{ 'M', 'G', -4 },{ 'M', 'H', -1 },{ 'M', 'K', -4 },{ 'M', 'M', -1 },{ 'M', 'N', -1 },{ 'M', 'S', -2 },{ 'M', 'R', -2 },{ 'M', 'T', -4 },{ 'M', 'W', -2 },{ 'M', 'V', -1 },{ 'M', 'Y', -2 },{ 'N', 'A', -2 },{ 'N', 'C', -2 },{ 'N', 'B', -1 },{ 'N', 'D', -1 },{ 'N', 'G', -2 },{ 'N', 'H', -1 },{ 'N', 'K', -1 },{ 'N', 'M', -1 },{ 'N', 'N', -1 },{ 'N', 'S', -1 },{ 'N', 'R', -1 },{ 'N', 'T', -2 },{ 'N', 'W', -1 },{ 'N', 'V', -1 },{ 'N', 'Y', -1 },{ 'S', 'A', -4 },{ 'S', 'C', 1 },{ 'S', 'B', -1 },{ 'S', 'D', -3 },{ 'S', 'G', 1 },{ 'S', 'H', -3 },{ 'S', 'K', -2 },{ 'S', 'M', -2 },{ 'S', 'N', -1 },{ 'S', 'S', -1 },{ 'S', 'R', -2 },{ 'S', 'T', -4 },{ 'S', 'W', -4 },{ 'S', 'V', -1 },{ 'S', 'Y', -2 },{ 'R', 'A', 1 },{ 'R', 'C', -4 },{ 'R', 'B', -3 },{ 'R', 'D', -1 },{ 'R', 'G', 1 },{ 'R', 'H', -3 },{ 'R', 'K', -2 },{ 'R', 'M', -2 },{ 'R', 'N', -1 },{ 'R', 'S', -2 },{ 'R', 'R', -1 },{ 'R', 'T', -4 },{ 'R', 'W', -2 },{ 'R', 'V', -1 },{ 'R', 'Y', -4 },{ 'T', 'A', -4 },{ 'T', 'C', -4 },{ 'T', 'B', -1 },{ 'T', 'D', -1 },{ 'T', 'G', -4 },{ 'T', 'H', -1 },{ 'T', 'K', 1 },{ 'T', 'M', -4 },{ 'T', 'N', -2 },{ 'T', 'S', -4 },{ 'T', 'R', -4 },{ 'T', 'T', 5 },{ 'T', 'W', 1 },{ 'T', 'V', -4 },{ 'T', 'Y', 1 },{ 'W', 'A', 1 },{ 'W', 'C', -4 },{ 'W', 'B', -3 },{ 'W', 'D', -1 },{ 'W', 'G', -4 },{ 'W', 'H', -1 },{ 'W', 'K', -2 },{ 'W', 'M', -2 },{ 'W', 'N', -1 },{ 'W', 'S', -4 },{ 'W', 'R', -2 },{ 'W', 'T', 1 },{ 'W', 'W', -1 },{ 'W', 'V', -3 },{ 'W', 'Y', -2 },{ 'V', 'A', -1 },{ 'V', 'C', -1 },{ 'V', 'B', -2 },{ 'V', 'D', -2 },{ 'V', 'G', -1 },{ 'V', 'H', -2 },{ 'V', 'K', -3 },{ 'V', 'M', -1 },{ 'V', 'N', -1 },{ 'V', 'S', -1 },{ 'V', 'R', -1 },{ 'V', 'T', -4 },{ 'V', 'W', -3 },{ 'V', 'V', -1 },{ 'V', 'Y', -3 },{ 'Y', 'A', -4 },{ 'Y', 'C', 1 },{ 'Y', 'B', -1 },{ 'Y', 'D', -3 },{ 'Y', 'G', -4 },{ 'Y', 'H', -1 },{ 'Y', 'K', -2 },{ 'Y', 'M', -2 },{ 'Y', 'N', -1 },{ 'Y', 'S', -2 },{ 'Y', 'R', -4 },{ 'Y', 'T', 1 },{ 'Y', 'W', -2 },{ 'Y', 'V', -3 },{ 'Y', 'Y', -1 } };
	int PAM250[529][3]={{'A','A',2},{'R','A',-2},{'N','A',0},{'D','A',0},{'C','A',-2},{'Q','A',0},{'E','A',0},{'G','A',1},{'H','A',-1},{'I','A',-1},{'L','A',-2},{'K','A',-1},{'M','A',-1},{'F','A',-3},{'P','A',1},{'S','A',1},{'T','A',1},{'W','A',-6},{'Y','A',-3},{'V','A',0},{'B','A',0},{'Z','A',0},{'X','A',0},{'A','R',-2},{'R','R',6},{'N','R',0},{'D','R',-1},{'C','R',-4},{'Q','R',1},{'E','R',-1},{'G','R',-3},{'H','R',2},{'I','R',-2},{'L','R',-3},{'K','R',3},{'M','R',0},{'F','R',-4},{'P','R',0},{'S','R',0},{'T','R',-1},{'W','R',2},{'Y','R',-4},{'V','R',-2},{'B','R',-1},{'Z','R',0},{'X','R',-1},{'A','N',0},{'R','N',0},{'N','N',2},{'D','N',2},{'C','N',-4},{'Q','N',1},{'E','N',1},{'G','N',0},{'H','N',2},{'I','N',-2},{'L','N',-3},{'K','N',1},{'M','N',-2},{'F','N',-3},{'P','N',0},{'S','N',1},{'T','N',0},{'W','N',-4},{'Y','N',-2},{'V','N',-2},{'B','N',2},{'Z','N',1},{'X','N',0},{'A','D',0},{'R','D',-1},{'N','D',2},{'D','D',4},{'C','D',-5},{'Q','D',2},{'E','D',3},{'G','D',1},{'H','D',1},{'I','D',-2},{'L','D',-4},{'K','D',0},{'M','D',-3},{'F','D',-6},{'P','D',-1},{'S','D',0},{'T','D',0},{'W','D',-7},{'Y','D',-4},{'V','D',-2},{'B','D',3},{'Z','D',3},{'X','D',-1},{'A','C',-2},{'R','C',-4},{'N','C',-4},{'D','C',-5},{'C','C',12},{'Q','C',-5},{'E','C',-5},{'G','C',-3},{'H','C',-3},{'I','C',-2},{'L','C',-6},{'K','C',-5},{'M','C',-5},{'F','C',-4},{'P','C',-3},{'S','C',0},{'T','C',-2},{'W','C',-8},{'Y','C',0},{'V','C',-2},{'B','C',-4},{'Z','C',-5},{'X','C',-3},{'A','Q',0},{'R','Q',1},{'N','Q',1},{'D','Q',2},{'C','Q',-5},{'Q','Q',4},{'E','Q',2},{'G','Q',-1},{'H','Q',3},{'I','Q',-2},{'L','Q',-2},{'K','Q',1},{'M','Q',-1},{'F','Q',-5},{'P','Q',0},{'S','Q',-1},{'T','Q',-1},{'W','Q',-5},{'Y','Q',-4},{'V','Q',-2},{'B','Q',1},{'Z','Q',3},{'X','Q',-1},{'A','E',0},{'R','E',-1},{'N','E',1},{'D','E',3},{'C','E',-5},{'Q','E',2},{'E','E',4},{'G','E',0},{'H','E',1},{'I','E',-2},{'L','E',-3},{'K','E',0},{'M','E',-2},{'F','E',-5},{'P','E',-1},{'S','E',0},{'T','E',0},{'W','E',-7},{'Y','E',-4},{'V','E',-2},{'B','E',3},{'Z','E',3},{'X','E',-1},{'A','G',1},{'R','G',-3},{'N','G',0},{'D','G',1},{'C','G',-3},{'Q','G',-1},{'E','G',0},{'G','G',5},{'H','G',-2},{'I','G',-3},{'L','G',-4},{'K','G',-2},{'M','G',-3},{'F','G',-5},{'P','G',0},{'S','G',1},{'T','G',0},{'W','G',-7},{'Y','G',-5},{'V','G',-1},{'B','G',0},{'Z','G',0},{'X','G',-1},{'A','H',-1},{'R','H',2},{'N','H',2},{'D','H',1},{'C','H',-3},{'Q','H',3},{'E','H',1},{'G','H',-2},{'H','H',6},{'I','H',-2},{'L','H',-2},{'K','H',0},{'M','H',-2},{'F','H',-2},{'P','H',0},{'S','H',-1},{'T','H',-1},{'W','H',-3},{'Y','H',0},{'V','H',-2},{'B','H',1},{'Z','H',2},{'X','H',-1},{'A','I',-1},{'R','I',-2},{'N','I',-2},{'D','I',-2},{'C','I',-2},{'Q','I',-2},{'E','I',-2},{'G','I',-3},{'H','I',-2},{'I','I',5},{'L','I',2},{'K','I',-2},{'M','I',2},{'F','I',1},{'P','I',-2},{'S','I',-1},{'T','I',0},{'W','I',-5},{'Y','I',-1},{'V','I',4},{'B','I',-2},{'Z','I',-2},{'X','I',-1},{'A','L',-2},{'R','L',-3},{'N','L',-3},{'D','L',-4},{'C','L',-6},{'Q','L',-2},{'E','L',-3},{'G','L',-4},{'H','L',-2},{'I','L',2},{'L','L',6},{'K','L',-3},{'M','L',4},{'F','L',2},{'P','L',-3},{'S','L',-3},{'T','L',-2},{'W','L',-2},{'Y','L',-1},{'V','L',2},{'B','L',-3},{'Z','L',-3},{'X','L',-1},{'A','K',-1},{'R','K',3},{'N','K',1},{'D','K',0},{'C','K',-5},{'Q','K',1},{'E','K',0},{'G','K',-2},{'H','K',0},{'I','K',-2},{'L','K',-3},{'K','K',5},{'M','K',0},{'F','K',-5},{'P','K',-1},{'S','K',0},{'T','K',0},{'W','K',-3},{'Y','K',-4},{'V','K',-2},{'B','K',1},{'Z','K',0},{'X','K',-1},{'A','M',-1},{'R','M',0},{'N','M',-2},{'D','M',-3},{'C','M',-5},{'Q','M',-1},{'E','M',-2},{'G','M',-3},{'H','M',-2},{'I','M',2},{'L','M',4},{'K','M',0},{'M','M',6},{'F','M',0},{'P','M',-2},{'S','M',-2},{'T','M',-1},{'W','M',-4},{'Y','M',-2},{'V','M',2},{'B','M',-2},{'Z','M',-2},{'X','M',-1},{'A','F',-3},{'R','F',-4},{'N','F',-3},{'D','F',-6},{'C','F',-4},{'Q','F',-5},{'E','F',-5},{'G','F',-5},{'H','F',-2},{'I','F',1},{'L','F',2},{'K','F',-5},{'M','F',0},{'F','F',9},{'P','F',-5},{'S','F',-3},{'T','F',-3},{'W','F',0},{'Y','F',7},{'V','F',-1},{'B','F',-4},{'Z','F',-5},{'X','F',-2},{'A','P',1},{'R','P',0},{'N','P',0},{'D','P',-1},{'C','P',-3},{'Q','P',0},{'E','P',-1},{'G','P',0},{'H','P',0},{'I','P',-2},{'L','P',-3},{'K','P',-1},{'M','P',-2},{'F','P',-5},{'P','P',6},{'S','P',1},{'T','P',0},{'W','P',-6},{'Y','P',-5},{'V','P',-1},{'B','P',-1},{'Z','P',0},{'X','P',-1},{'A','S',1},{'R','S',0},{'N','S',1},{'D','S',0},{'C','S',0},{'Q','S',-1},{'E','S',0},{'G','S',1},{'H','S',-1},{'I','S',-1},{'L','S',-3},{'K','S',0},{'M','S',-2},{'F','S',-3},{'P','S',1},{'S','S',2},{'T','S',1},{'W','S',-2},{'Y','S',-3},{'V','S',-1},{'B','S',0},{'Z','S',0},{'X','S',0},{'A','T',1},{'R','T',-1},{'N','T',0},{'D','T',0},{'C','T',-2},{'Q','T',-1},{'E','T',0},{'G','T',0},{'H','T',-1},{'I','T',0},{'L','T',-2},{'K','T',0},{'M','T',-1},{'F','T',-3},{'P','T',0},{'S','T',1},{'T','T',3},{'W','T',-5},{'Y','T',-3},{'V','T',0},{'B','T',0},{'Z','T',-1},{'X','T',0},{'A','W',-6},{'R','W',2},{'N','W',-4},{'D','W',-7},{'C','W',-8},{'Q','W',-5},{'E','W',-7},{'G','W',-7},{'H','W',-3},{'I','W',-5},{'L','W',-2},{'K','W',-3},{'M','W',-4},{'F','W',0},{'P','W',-6},{'S','W',-2},{'T','W',-5},{'W','W',17},{'Y','W',0},{'V','W',-6},{'B','W',-5},{'Z','W',-6},{'X','W',-4},{'A','Y',-3},{'R','Y',-4},{'N','Y',-2},{'D','Y',-4},{'C','Y',0},{'Q','Y',-4},{'E','Y',-4},{'G','Y',-5},{'H','Y',0},{'I','Y',-1},{'L','Y',-1},{'K','Y',-4},{'M','Y',-2},{'F','Y',7},{'P','Y',-5},{'S','Y',-3},{'T','Y',-3},{'W','Y',0},{'Y','Y',10},{'V','Y',-2},{'B','Y',-3},{'Z','Y',-4},{'X','Y',-2},{'A','V',0},{'R','V',-2},{'N','V',-2},{'D','V',-2},{'C','V',-2},{'Q','V',-2},{'E','V',-2},{'G','V',-1},{'H','V',-2},{'I','V',4},{'L','V',2},{'K','V',-2},{'M','V',2},{'F','V',-1},{'P','V',-1},{'S','V',-1},{'T','V',0},{'W','V',-6},{'Y','V',-2},{'V','V',4},{'B','V',-2},{'Z','V',-2},{'X','V',-1},{'A','B',0},{'R','B',-1},{'N','B',2},{'D','B',3},{'C','B',-4},{'Q','B',1},{'E','B',3},{'G','B',0},{'H','B',1},{'I','B',-2},{'L','B',-3},{'K','B',1},{'M','B',-2},{'F','B',-4},{'P','B',-1},{'S','B',0},{'T','B',0},{'W','B',-5},{'Y','B',-3},{'V','B',-2},{'B','B',3},{'Z','B',2},{'X','B',-1},{'A','Z',0},{'R','Z',0},{'N','Z',1},{'D','Z',3},{'C','Z',-5},{'Q','Z',3},{'E','Z',3},{'G','Z',0},{'H','Z',2},{'I','Z',-2},{'L','Z',-3},{'K','Z',0},{'M','Z',-2},{'F','Z',-5},{'P','Z',0},{'S','Z',0},{'T','Z',-1},{'W','Z',-6},{'Y','Z',-4},{'V','Z',-2},{'B','Z',2},{'Z','Z',3},{'X','Z',-1},{'A','X',0},{'R','X',-1},{'N','X',0},{'D','X',-1},{'C','X',-3},{'Q','X',-1},{'E','X',-1},{'G','X',-1},{'H','X',-1},{'I','X',-1},{'L','X',-1},{'K','X',-1},{'M','X',-1},{'F','X',-2},{'P','X',-1},{'S','X',0},{'T','X',0},{'W','X',-4},{'Y','X',-2},{'V','X',-1},{'B','X',-1},{'Z','X',-1},{'X','X',-1}};
	int BLOSUM62[400][3]={{'C','C',9},{'S','C',-1},{'T','C',-1},{'P','C',-3},{'A','C',0},{'G','C',-3},{'N','C',-3},{'D','C',-3},{'E','C',-4},{'Q','C',-3},{'H','C',-3},{'R','C',-3},{'K','C',-3},{'M','C',-1},{'I','C',-1},{'L','C',-1},{'V','C',-1},{'F','C',-2},{'Y','C',-2},{'W','C',-2},{'C','S',-1},{'S','S',4},{'T','S',1},{'P','S',-1},{'A','S',1},{'G','S',0},{'N','S',1},{'D','S',0},{'E','S',0},{'Q','S',0},{'H','S',-1},{'R','S',-1},{'K','S',0},{'M','S',-1},{'I','S',-2},{'L','S',-2},{'V','S',-2},{'F','S',-2},{'Y','S',-2},{'W','S',-3},{'C','T',-1},{'S','T',1},{'T','T',4},{'P','T',1},{'A','T',-1},{'G','T',1},{'N','T',0},{'D','T',1},{'E','T',0},{'Q','T',0},{'H','T',0},{'R','T',-1},{'K','T',0},{'M','T',-1},{'I','T',-2},{'L','T',-2},{'V','T',-2},{'F','T',-2},{'Y','T',-2},{'W','T',-3},{'C','P',-3},{'S','P',-1},{'T','P',1},{'P','P',7},{'A','P',-1},{'G','P',-2},{'N','P',-1},{'D','P',-1},{'E','P',-1},{'Q','P',-1},{'H','P',-2},{'R','P',-2},{'K','P',-1},{'M','P',-2},{'I','P',-3},{'L','P',-3},{'V','P',-2},{'F','P',-4},{'Y','P',-3},{'W','P',-4},{'C','A',0},{'S','A',1},{'T','A',-1},{'P','A',-1},{'A','A',4},{'G','A',0},{'N','A',-1},{'D','A',-2},{'E','A',-1},{'Q','A',-1},{'H','A',-2},{'R','A',-1},{'K','A',-1},{'M','A',-1},{'I','A',-1},{'L','A',-1},{'V','A',-2},{'F','A',-2},{'Y','A',-2},{'W','A',-3},{'C','G',-3},{'S','G',0},{'T','G',1},{'P','G',-2},{'A','G',0},{'G','G',6},{'N','G',-2},{'D','G',-1},{'E','G',-2},{'Q','G',-2},{'H','G',-2},{'R','G',-2},{'K','G',-2},{'M','G',-3},{'I','G',-4},{'L','G',-4},{'V','G',0},{'F','G',-3},{'Y','G',-3},{'W','G',-2},{'C','N',-3},{'S','N',1},{'T','N',0},{'P','N',-2},{'A','N',-2},{'G','N',0},{'N','N',6},{'D','N',1},{'E','N',0},{'Q','N',0},{'H','N',-1},{'R','N',0},{'K','N',0},{'M','N',-2},{'I','N',-3},{'L','N',-3},{'V','N',-3},{'F','N',-3},{'Y','N',-2},{'W','N',-4},{'C','D',-3},{'S','D',0},{'T','D',1},{'P','D',-1},{'A','D',-2},{'G','D',-1},{'N','D',1},{'D','D',6},{'E','D',2},{'Q','D',0},{'H','D',-1},{'R','D',-2},{'K','D',-1},{'M','D',-3},{'I','D',-3},{'L','D',-4},{'V','D',-3},{'F','D',-3},{'Y','D',-3},{'W','D',-4},{'C','E',-4},{'S','E',0},{'T','E',0},{'P','E',-1},{'A','E',-1},{'G','E',-2},{'N','E',0},{'D','E',2},{'E','E',5},{'Q','E',2},{'H','E',0},{'R','E',0},{'K','E',1},{'M','E',-2},{'I','E',-3},{'L','E',-3},{'V','E',-3},{'F','E',-3},{'Y','E',-2},{'W','E',-3},{'C','Q',-3},{'S','Q',0},{'T','Q',0},{'P','Q',-1},{'A','Q',-1},{'G','Q',-2},{'N','Q',0},{'D','Q',0},{'E','Q',2},{'Q','Q',5},{'H','Q',0},{'R','Q',1},{'K','Q',1},{'M','Q',0},{'I','Q',-3},{'L','Q',-2},{'V','Q',-2},{'F','Q',-3},{'Y','Q',-1},{'W','Q',-2},{'C','H',-3},{'S','H',-1},{'T','H',0},{'P','H',-2},{'A','H',-2},{'G','H',-2},{'N','H',1},{'D','H',1},{'E','H',0},{'Q','H',0},{'H','H',8},{'R','H',0},{'K','H',-1},{'M','H',-2},{'I','H',-3},{'L','H',-3},{'V','H',-2},{'F','H',-1},{'Y','H',2},{'W','H',-2},{'C','R',-3},{'S','R',-1},{'T','R',-1},{'P','R',-2},{'A','R',-1},{'G','R',-2},{'N','R',0},{'D','R',-2},{'E','R',0},{'Q','R',1},{'H','R',0},{'R','R',5},{'K','R',2},{'M','R',-1},{'I','R',-3},{'L','R',-2},{'V','R',-3},{'F','R',-3},{'Y','R',-2},{'W','R',-3},{'C','K',-3},{'S','K',0},{'T','K',0},{'P','K',-1},{'A','K',-1},{'G','K',-2},{'N','K',0},{'D','K',-1},{'E','K',1},{'Q','K',1},{'H','K',-1},{'R','K',2},{'K','K',5},{'M','K',-1},{'I','K',-3},{'L','K',-2},{'V','K',-3},{'F','K',-3},{'Y','K',-2},{'W','K',-3},{'C','M',-1},{'S','M',-1},{'T','M',-1},{'P','M',-2},{'A','M',-1},{'G','M',-3},{'N','M',-2},{'D','M',-3},{'E','M',-2},{'Q','M',0},{'H','M',-2},{'R','M',-1},{'K','M',-1},{'M','M',5},{'I','M',1},{'L','M',2},{'V','M',-2},{'F','M',0},{'Y','M',-1},{'W','M',-1},{'C','I',-1},{'S','I',-2},{'T','I',-2},{'P','I',-3},{'A','I',-1},{'G','I',-4},{'N','I',-3},{'D','I',-3},{'E','I',-3},{'Q','I',-3},{'H','I',-3},{'R','I',-3},{'K','I',-3},{'M','I',1},{'I','I',4},{'L','I',2},{'V','I',1},{'F','I',0},{'Y','I',-1},{'W','I',-3},{'C','L',-1},{'S','L',-2},{'T','L',-2},{'P','L',-3},{'A','L',-1},{'G','L',-4},{'N','L',-3},{'D','L',-4},{'E','L',-3},{'Q','L',-2},{'H','L',-3},{'R','L',-2},{'K','L',-2},{'M','L',2},{'I','L',2},{'L','L',4},{'V','L',3},{'F','L',0},{'Y','L',-1},{'W','L',-2},{'C','V',-1},{'S','V',-2},{'T','V',-2},{'P','V',-2},{'A','V',0},{'G','V',-3},{'N','V',-3},{'D','V',-3},{'E','V',-2},{'Q','V',-2},{'H','V',-3},{'R','V',-3},{'K','V',-2},{'M','V',1},{'I','V',3},{'L','V',1},{'V','V',4},{'F','V',-1},{'Y','V',-1},{'W','V',-3},{'C','F',-2},{'S','F',-2},{'T','F',-2},{'P','F',-4},{'A','F',-2},{'G','F',-3},{'N','F',-3},{'D','F',-3},{'E','F',-3},{'Q','F',-3},{'H','F',-1},{'R','F',-3},{'K','F',-3},{'M','F',0},{'I','F',0},{'L','F',0},{'V','F',-1},{'F','F',6},{'Y','F',3},{'W','F',1},{'C','Y',-2},{'S','Y',-2},{'T','Y',-2},{'P','Y',-3},{'A','Y',-2},{'G','Y',-3},{'N','Y',-2},{'D','Y',-3},{'E','Y',-2},{'Q','Y',-1},{'H','Y',2},{'R','Y',-2},{'K','Y',-2},{'M','Y',-1},{'I','Y',-1},{'L','Y',-1},{'V','Y',-1},{'F','Y',3},{'Y','Y',7},{'W','Y',2},{'C','W',-2},{'S','W',-3},{'T','W',-3},{'P','W',-4},{'A','W',-3},{'G','W',-2},{'N','W',-4},{'D','W',-4},{'E','W',-3},{'Q','W',-2},{'H','W',-2},{'R','W',-3},{'K','W',-3},{'M','W',-1},{'I','W',-3},{'L','W',-2},{'V','W',-3},{'F','W',1},{'Y','W',2},{'W','W',11}};
	
	k = 0;
	while (k<l1) {
		if (ali1[k] == '-' && ali2[k] == '-') {
			l1--;
			l2--;
			for (int i = k; i<l1; i++)
				ali1[i] = ali1[i + 1];
			ali1[l1] = '\0';
			for (int i = k; i<l2; i++)
				ali2[i] = ali2[i + 1];
			ali2[l2] = '\0';
		}
		else { k++; }
	}
	scoreMax = 5 * l1;
	for (k = 160; k>0; k--) {
		int l = 0;
		while (l<l1 - k + 1) {
			bool m1 = true, m2 = true;
			for (int i = l; i<l + k; i++) {
				if (ali1[i] != '-') {
					m1 = false;
					break;
				}
			}
			for (int i = l; i<l + k; i++) {
				if (ali2[i] != '-') {
					m2 = false;
					break;
				}
			}
			if (m1 || m2) {
				score += gapOpen + gapExtend*k;
				for (int i = l; i<l1 - k; i++) {
					ali1[i] = ali1[i + k];
				}
				l1 = l1 - k;
				ali1[l1] = '\0';
				for (int i = l; i<l2 - k; i++) {
					ali2[i] = ali2[i + k];
				}
				l2 = l2 - k;
				ali2[l2] = '\0';
			}
			else { l++; }
		}
	}
	
	if(tolowerstr(mat)=="ednafull"){
		for (k = 0; k<l1; k++) {
			dist=0;
			for (int i = 0; i<225; i++)
				if (EDNAFULL[i][0] == ali1[k] && EDNAFULL[i][1] == ali2[k]) {
					dist = EDNAFULL[i][2];
					break;
				}
			score+=dist;
		}
	}else if(tolowerstr(mat)=="pam250"){
		for (k = 0; k<l1; k++) {
			dist=-8;
			for (int i = 0; i<529; i++)
				if (PAM250[i][0] == ali1[k] && PAM250[i][1] == ali2[k]) {
					dist = PAM250[i][2];
					break;
				}
			score+=dist;
		}
	}else{ //mat=="BLOSUM62"
		for (k = 0; k<l1; k++) {
			dist=-4;
			for (int i = 0; i<400; i++)
				if (BLOSUM62[i][0] == ali1[k] && BLOSUM62[i][1] == ali2[k]) {
					dist = BLOSUM62[i][2];
					break;
				}
			score+=dist;
		}
	}
	
	return double(scoreMax - score) / scoreMax;
}

std::string exec(const char* cmd) {
	/*
	This function executes an external shell command and retrieves its output in a string.
	Input:
	-cmd : The command to execute.
	Output:
	-result : A string holding the command output.
	*/
	char buffer[128];
	std::string result = "";
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	FILE* pipe = popen(cmd, "r");
#else
	FILE* pipe = _popen(cmd, "r");
#endif
	if (!pipe) throw std::runtime_error("popen() failed!");
	try {
		while (!feof(pipe)) {
			if (fgets(buffer, 128, pipe) != NULL)
				result += buffer;
		}
	}
	catch (...) {
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		pclose(pipe);
#else			
		_pclose(pipe);
#endif
		throw;
	}
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	pclose(pipe);
#else			
	_pclose(pipe);
#endif
	return result;
}

string getStringFromDouble(double num){
	ostringstream streamObj;
	streamObj << std::fixed;
	streamObj << std::setprecision(18);
	streamObj << num;
	string strObj = streamObj.str();
	return strObj;
}

void toAffinityMat(double *mat, int n, string type="rwnl"){
	if(tolowerstr(type)=="rwnl"){ //Random Walk Normalized Laplacian
		double deg;
		for(int l=0; l<n; l++){
			deg=0.0;
			for(int m=0; m<n; m++)
				deg+=mat[l*n+m];
			for(int m=0; m<n; m++)
				if(l==m)
					mat[l*n+m]=(deg-mat[l*n+m])/deg;
				else
					mat[l*n+m]=(0.0-mat[l*n+m])/deg;
		}
	} else if(tolowerstr(type)=="ul"){ //Unnormalized laplacian
		double deg;
		for(int l=0; l<n; l++){
			deg=0.0;
			for(int m=0; m<n; m++)
				deg+=mat[l*n+m];
			for(int m=0; m<n; m++)
				if(l==m)
					mat[l*n+m]=deg-mat[l*n+m];
				else
					mat[l*n+m]=0.0-mat[l*n+m];
		}
	} else if(tolowerstr(type)=="mod"){ //Modularity
		double K=0.0,*k;
		k=new double [n];
		for(int l=0; l<n; l++){ //Calculating the degrees vector k and the total degree K
			k[l]=0.0;
			for(int m=0; m<n; m++){
				k[l]+=mat[l*n+m];
				K+=mat[l*n+m];
			}
		}
		for(int l=0; l<n; l++){//M=(A-kkT/K)/K
			for(int m=0; m<n; m++){
				mat[l*n+m]=(mat[l*n+m]-(k[l]*k[m])/K)/K;
			}
		}
		delete k;
	} else if(tolowerstr(type)=="bh"){ //Bethe Hessian
		//r= (sqrt(c)) where c is the average degreee of the graph
		//In our case, all the nodes are linked so the degree of each node is the number of nodes - 1
		//therefore r = sqrt(c) = sqrt(n-1) and r^2 = n-1
		double r=sqrt(n-1), *k;
		k=new double [n]; // 1 + SUM i!=j ((Wij^2)/(r^2 - Wij^2)) this constitutes the first part of the equation where i=j. k=all the nodes except i because all the nodes are neighbors to i
		for(int l=0; l<n; l++){ //Calculating the degrees vector k and the total degree K
			k[l]=1.0;
			for(int m=0; m<n; m++){
				if(m!=l){
					k[l]+=(pow(mat[l*n+m],2))/((n-1)-(pow(mat[l*n+m],2)));
				}
			}
		}
		for(int l=0; l<n; l++){//M=(A-kkT/K)/K
			for(int m=0; m<n; m++){
				mat[l*n+m]=-(r*mat[l*n+m])/((n-1)-(pow(mat[l*n+m],2)));
				if(l==m){
					mat[l*n+m]=mat[l*n+m]+k[l];
				}
			}
		}
		delete[] k;
	} else { //Non-Backtracking
	}
}

bool is_number(const string& s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void eigenSolver(int nbSequences, double *affMat, int it_max, double *eigVects, double *eigVals, int& it_num, int& rot_num){
	
}

int featuresCalc(double *affMat, double **featuresMat, int& nbFeatures, int& nbClusters, int nbSamples, string EVMethod="delta", float EVCutOff=0.01, string affMatType="rwnl"){
	
	//Calculate the Eigenvalues and Eigenvectors
	double *eigVals, *eigVects; //a holds the matrix, d holds the eigenvalues, v holds the eigenvectors
	int it_max=100; //maximum number of iterations
	int it_num; //used number of iterations
	int rot_num; //used number of rotations
		
	//Allocate space for eigenvalues
	if (!(eigVals = new double[nbSamples])){
		cout << "Allocation for 'eigVals' failed. \n";
		return 1;
	} 
	//Allocate space for eigenvectors
	if (!(eigVects = new double[nbSamples*nbSamples])){
		cout << "Allocation for 'eigVects' failed. \n";
		delete []eigVals;
		return 1;
	}
		
	jacobi_eigenvalue ( nbSamples, affMat, it_max, eigVects, eigVals, it_num, rot_num );
	//Print the eingenvalues
	/*for(int i=0; i<nbSamples; i++)
		cout<<eigVals[i]<<"  ";*/

	//Calculating the number of eigenvectors to use based on the chosen methos and EVCutOff
	if(EVMethod=="delta"){
		nbFeatures=1;
		while (eigVals[nbFeatures+1] > eigVals[nbFeatures] + EVCutOff){
            nbFeatures+=1;
            if (nbFeatures == (nbSamples-2))
               break;
		}
	} else if(EVMethod=="log"){
		nbFeatures=int(log(nbSamples));
	} else if(EVMethod=="bh"){//Bethe Hessian: count the number of negative eigenvectors as nbClusters and use their eigenvalues
		nbFeatures=1;
		while (eigVals[nbFeatures]<0){
			nbFeatures+=1;
            if (nbFeatures == (nbSamples-1))
               break;
		}
		nbClusters=nbFeatures;
	}
	
	//Allocating and initializing the features matrix
	for(int i=0; i<nbSamples; i++)
		if (!(featuresMat[i] = new double[nbFeatures])){
		cout << "Allocation for 'featuresMat' failed. \n";
		return 1;
	} 
	for(int i=0; i<nbSamples; i++){
		for(int j=0; j<nbFeatures; j++){
			if(tolowerstr(affMatType)!="mod"){ //If NOT using the modularity matrix then take the eigenvectors of the smallest eigenvalues
				featuresMat[i][j]=eigVects[j*nbSamples+i];
			} else { //If using the modularity matrix then take the eigenvectors of the highest eigenvalues
				featuresMat[i][j]=eigVects[(nbSamples-1-j)*nbSamples+i];
			}
		}
	}

	return 0;
}

int getBestNbClusters(double **data, int number_data, int number_features, int rndseed, string criterion = "BIC", string type_covariance="full", int max_iterations=1000){
	int *p = new int[number_data];
	int nbClusters;
	GaussianMixture *gmm;
	double crit=INFTY;
	if (criterion=="BIC"){
       for(int j = 1;j<number_data+1;j++){
			gmm = new GaussianMixture(j,number_data, number_features, type_covariance,max_iterations,0.01,0.001,rndseed);
			(*gmm).train(data,p,number_data,number_features);
			if(crit>(*gmm).bic()){
				crit=(*gmm).bic();
				nbClusters = j;
			}else{
				delete gmm;
				break;
			}
			delete gmm;
			//cout<<"For nb = "<<j<<" BIC = "<<crit<<endl;
		}
	}else if (criterion=="AIC"){
        for(int j = 1;j<number_data+1;j++){
			gmm = new GaussianMixture(j,number_data, number_features, type_covariance,max_iterations,0.01,0.001,rndseed);
			(*gmm).train(data,p,number_data,number_features);
			if(crit>(*gmm).aic()){
				crit=(*gmm).aic();
				nbClusters = j;
			}else{
				delete gmm;
				break;
			} 
			delete gmm;
		}    
	}
	delete[] p;
	return nbClusters;
}

double silhouette(double **coords, int *clusterNbs, int nbPts, int dim){
	double sumSi=0.0, Si, ai, bi, dij, sdij;
	int nbClstrs=0, Ci, Ck;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	for(int i=0; i<nbPts; i++){
		Ci=0;
		sdij=0.0;
		for(int j=0; j<nbPts; j++)
			if(clusterNbs[j]==clusterNbs[i])
				Ci++;
		if(Ci!=1){
			for(int j=0; j<nbPts; j++){
				if(i!=j && clusterNbs[j]==clusterNbs[i]){
					dij=0.0;
					for(int l=0; l<dim; l++){
						dij=dij+pow((coords[i][l]-coords[j][l]),2);
					}
					dij=sqrt(dij);
					sdij=sdij+dij;
				}
			}
			ai=(1.0/double(Ci-1))*sdij;
			
			bi=INFTY;
			for(int k=0; k<=nbClstrs; k++){
				if(k!=clusterNbs[i]){
					Ck=0;
					sdij=0.0;
					for(int j=0; j<nbPts; j++)
						if(clusterNbs[j]==k)
							Ck++;
					if(Ck>0){
						for(int j=0; j<nbPts; j++){
							if(clusterNbs[j]==k){
								dij=0.0;
								for(int l=0; l<dim; l++){
									dij=dij+pow((coords[i][l]-coords[j][l]),2);
								}
								dij=sqrt(dij);
								sdij=sdij+dij;
							}
						}
						if(bi>((1.0/double(Ck))*sdij))
							bi=(1.0/double(Ck))*sdij;
					}
				}
			}
			if(ai<bi)
				Si=1.0-(ai/bi);
			else
				Si=(bi/ai)-1.0;
			sumSi=sumSi+Si;
		}
	}
	return sumSi/double(nbPts);
}

int getMedoidIndex(int clstrNb, double **coords, int *clusterNbs, int nbPts, int dim){ //Returns the medoid index of cluster clstrNb
	double sdyxi = INFTY;
	int mIndex = -1;
	double dij, sdij;
	for(int i=0; i<nbPts; i++){
		if(clusterNbs[i]==clstrNb){
			sdij=0.0;
			for(int j=0; j<nbPts; j++){
				if(i!=j && clusterNbs[j]==clusterNbs[i]){
					dij=0.0;
					for(int l=0; l<dim; l++){
						dij=dij+pow((coords[i][l]-coords[j][l]),2);
					}
					dij=sqrt(dij);
					sdij=sdij+dij;
				}
			}
			if(sdij<sdyxi){
				sdyxi=sdij;
				mIndex=i;
			}
		}
	}
	return mIndex;
}

int getDsCentroidIndex(double **coords, int nbPts, int dim){ //Returns the centoid index of the data set
	double sdyxi = INFTY;
	int cIndex = -1;
	double dij, sdij;
	for(int i=0; i<nbPts; i++){
		sdij=0.0;
		for(int j=0; j<nbPts; j++){
			if(i!=j){
				dij=0.0;
				for(int l=0; l<dim; l++){
					dij=dij+pow((coords[i][l]-coords[j][l]),2);
				}
				dij=sqrt(dij);
				sdij=sdij+dij;
			}
		}
		if(sdij<sdyxi){
			sdyxi=sdij;
			cIndex=i;
		}
	}
	return cIndex;
}

double daviesbouldin(double **coords, int *clusterNbs, int nbPts, int dim){
	double *S, **R, DB, dmij, sdmij;
	int nbClstrs=0, realNbClstrs=0, clstrSize;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	S=new double[nbClstrs+1]; //nbClstrs+1 because nbClstrs hold the highest cluser number and the numbers could be starting from 0.
	R=new double*[nbClstrs+1];
	for(int i=0; i<=nbClstrs; i++)
		R[i]=new double[nbClstrs+1];
	for(int i=0; i<=nbClstrs; i++){
		int mi = getMedoidIndex(i,coords,clusterNbs,nbPts,dim);
		clstrSize=0;
		sdmij=0.0;
		if(mi!=(-1)){//Cluster i not empty;
			realNbClstrs++;//The number of non-empty clusters to take into consideration
			for(int j=0; j<nbPts; j++){
				if(clusterNbs[j]==i){
					clstrSize++;
					dmij=0.0;
					for(int l=0; l<dim; l++){
						dmij=dmij+pow((coords[mi][l]-coords[j][l]),2);
					}
					dmij=sqrt(dmij);
					sdmij=sdmij+dmij;
				}
			}
			S[i]=sdmij/double(clstrSize);
		} else {S[i]=-1;}//Cluster i is empty
	}
	for(int i=0; i<=nbClstrs; i++){
		for(int j=i; j<=nbClstrs; j++){
			if(j==i){
				R[i][j]=-INFTY; //Will never be max because we have max i!=j in the DB sum
			} else{
				int mi = getMedoidIndex(i,coords,clusterNbs,nbPts,dim);
				int mj = getMedoidIndex(j,coords,clusterNbs,nbPts,dim);
				if(mi!=-1 && mj!=-1){//Both clusters are not empty
					double dmimj=0.0;
					for(int l=0; l<dim; l++){
						dmimj=dmimj+pow((coords[mi][l]-coords[mj][l]),2);
					}
					dmimj=sqrt(dmimj);
					R[i][j]=(S[i]+S[j])/dmimj;
					R[j][i]=R[i][j];
				} else {R[i][j]=-INFTY; R[j][i]=-INFTY;} //Will never be taken into consideration when chosing max
			}
		}
	}
	double sMaxRij=0.0;
	double maxRij;
	for(int i=0; i<=nbClstrs; i++){
		maxRij=-INFTY;
		for(int j=0; j<=nbClstrs; j++)
			if(R[i][j]>maxRij)
				maxRij=R[i][j];
		if(maxRij!=-INFTY) //Omit empty clusters
			sMaxRij = sMaxRij + maxRij;
	}
	DB=sMaxRij/double(realNbClstrs);

	for(int i=0; i<=nbClstrs; i++)
		delete [] R[i];
	delete [] R;
	delete [] S;
	return DB;
}

double calinskiharabasz(double **coords, int *clusterNbs, int nbPts, int dim){
	double CH, sqrdmij, ssqrdc, SSw=0.0, tss=0.0, SSb;
	int nbClstrs=0, realNbClstrs=0, ciMedInd, dsCInd;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	for(int i=0; i<=nbClstrs; i++){
		ssqrdc=0.0;
		ciMedInd=getMedoidIndex(i,coords,clusterNbs,nbPts,dim);
		if(ciMedInd!=-1){ //Cluster exists (not empty)
			realNbClstrs++;
			for(int j=0; j<nbPts; j++){
				sqrdmij=0.0;
				if(clusterNbs[j]==i){
					for(int l=0; l<dim; l++){
						sqrdmij=sqrdmij+pow((coords[ciMedInd][l]-coords[j][l]),2);
					}
				}
				ssqrdc = ssqrdc + sqrdmij;
			}
		}
		SSw = SSw + ssqrdc;
	}

	dsCInd=getDsCentroidIndex(coords,nbPts,dim);
	if(dsCInd!=-1){ //Data set not empty
		for(int j=0; j<nbPts; j++){
			sqrdmij=0.0;
			for(int l=0; l<dim; l++){
				sqrdmij=sqrdmij+pow((coords[dsCInd][l]-coords[j][l]),2);
			}
			tss = tss + sqrdmij;
		}
	}
	SSb=tss-SSw;

	CH=(SSb/SSw)*(double(nbPts-realNbClstrs)/double(realNbClstrs-1));

	return CH;
}

double dunn(double **coords, int *clusterNbs, int nbPts, int dim){
	double D, minSep = INFTY, maxDmtr = -INFTY, dij;
	int nbClstrs=0;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	for(int i=0; i<=nbPts; i++){
		for(int j=i+1; j<nbPts; j++){
			if(clusterNbs[i]!=clusterNbs[j]){
				dij=0.0;
				for(int l=0; l<dim; l++){
					dij=dij+pow((coords[i][l]-coords[j][l]),2);
				}
				dij = sqrt(dij);
				if(minSep>dij)
					minSep=dij;
			}
		}
	}

	for(int i=0; i<=nbPts; i++){
		for(int j=i+1; j<nbPts; j++){
			if(clusterNbs[i]==clusterNbs[j]){
				dij=0.0;
				for(int l=0; l<dim; l++){
					dij=dij+pow((coords[i][l]-coords[j][l]),2);
				}
				dij = sqrt(dij);
				if(maxDmtr<dij)
					maxDmtr=dij;
			}
		}
	}

	D=minSep/maxDmtr;

	return D;
}

void printInternalValidationIndices(double **coords, int *clusterNbs, int nbPts, int dim){	
	if(internalValidationOn){
		cout<<"The mean Silhouette is: "<<silhouette(coords,clusterNbs,nbPts,dim)<<endl;
		cout<<"Davies-Bouldin index is: "<<daviesbouldin(coords,clusterNbs,nbPts,dim)<<endl;
		cout<<"Calinski-Harabasz index is: "<<calinskiharabasz(coords,clusterNbs,nbPts,dim)<<endl;
		cout<<"Dunn index is: "<<dunn(coords,clusterNbs,nbPts,dim)<<endl;
	}
}