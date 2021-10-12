//#include<amconf.h>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<mpi.h>
#include<time.h>
#include"spclustFunctions.h"
using namespace std;

bool installed=false;
int Aligned_seq_length; //Alligned sequence length
char progPath[255]=""; //Path of the program directory
string mDist; //Distance matrix used for distance calculation
double gapOpenI, gapExtendI; //Parameters used for distance calculation
string *nc; //Array that will hold the input genomes' names
double *MatSimil; //MatSimil will hold the similarity matrix
int nbSequences; //The number of sequences in the input fasta file initialized in similarity() function
char *buffChar1,*buffChar2; //Buffer for sequences distance calculation
int *dCount; //Nunber of distances to be calculated by each process
//MPI communication buffers and variables
int myid;
char *aliSeqs; //Aligned sequences buffer
int numprocs; //Nunber of processes
MPI_Status status;
double *buffDouble; //Calculated distances buffer
int progEnd=0; //Exit flag


void sendDataToWorkers(){ //Send the required data to workors processes
	for(int i=1; i<numprocs; i++){
		MPI_Send(&progEnd,1,MPI_INT,i,99,MPI_COMM_WORLD);
		MPI_Send(&nbSequences,1,MPI_INT,i,98,MPI_COMM_WORLD);
		MPI_Send(&Aligned_seq_length,1,MPI_INT,i,97,MPI_COMM_WORLD);
		MPI_Send(aliSeqs,nbSequences*(Aligned_seq_length+1)+1,MPI_CHAR,i,96,MPI_COMM_WORLD);
	}
}

void slaveDistCalc(int n=nbSequences){
	//Determine the number of distances to be calculated by each slave
	int nCalcs=((n*n)+n)/2; //Number of calculations
	for(int i=0; i<numprocs; i++){
		if(nCalcs%numprocs>i)
			dCount[i]=ceil(float(nCalcs)/float(numprocs));
		else
			dCount[i]=floor(float(nCalcs)/float(numprocs));
	}
	
	//Allocate the distance vector of this slave
	buffDouble=new double[dCount[myid]];
	
	//Fill this slave's distances vector
	int sInd=0; //Starting calculation item's index
	for(int i=0; i<myid; i++)
		sInd+=dCount[i];
	int cInd=0; //Calculation item's index
	for (int i = 0; i<n; i++)
		for (int j = i; j<n; j++,cInd++){
			if(cInd>=sInd && cInd<sInd+dCount[myid]){
				strncpy(buffChar1,aliSeqs+(i*(Aligned_seq_length+1)),Aligned_seq_length); //Copy the i th sequence
				buffChar1[Aligned_seq_length]='\0';
				strncpy(buffChar2,aliSeqs+(j*(Aligned_seq_length+1)),Aligned_seq_length); //Copy the j th sequence
				buffChar2[Aligned_seq_length]='\0';
				buffDouble[cInd-sInd]=distanceS(mDist,buffChar1,buffChar2,gapOpenI,gapExtendI);
			}
		}
}

void matriceDistances(string **dicoMuscle, double *matDist, int n = nbSequences) {
	/*
	This function is a subfunction of "similarity".
	Input:
	-dicoMuscle : A dictionary of sequences in which the keys are the names of the sequences.
	-nc : The names of the sequences.
	Output:
	-Matrix : A matrix of distance.
	*/
	
	Aligned_seq_length=dicoMuscle[0][1].length();
	buffChar1=new char[Aligned_seq_length+1];
	buffChar2=new char[Aligned_seq_length+1];
	if(numprocs == 1){
		for (int i = 0; i<n; i++)
			for (int j = i; j<n; j++){
				strcpy(buffChar1,dicoMuscle[i][1].c_str());
				strcpy(buffChar2,dicoMuscle[j][1].c_str());
				matDist[i*n+j] = distanceS(mDist,buffChar1,buffChar2,gapOpenI,gapExtendI);
			}
	} else {
		// Vectorize the aligned sequences
		aliSeqs=new char[nbSequences*(Aligned_seq_length+1)+1];
		strcpy(aliSeqs,dicoMuscle[0][1].c_str());
		strcat(aliSeqs,"+");
		for(int i=1; i<nbSequences; i++){
			strcat(aliSeqs,dicoMuscle[i][1].c_str());
			strcat(aliSeqs,"+");
		}
		strcat(aliSeqs,"\0");

		//Send the required data to slaves
		sendDataToWorkers();

		//Do the required calculations from this slave
		slaveDistCalc();

		//Collect the results from slaves and fill the lower part of the distance matrix
		int curProc=0; //Current process' results
		int toInd=dCount[curProc];
		int curInd=0, vectInd=0;
		for (int i = 0; i<n; i++)
			for (int j = i; j<n; j++,curInd++,vectInd++){
				if(curInd==toInd){ //Get and continue with the next process' results
					curProc++;
					toInd+=dCount[curProc];
					vectInd=0;
					MPI_Recv(buffDouble,dCount[curProc],MPI_DOUBLE,curProc,curProc*100+11,MPI_COMM_WORLD,&status);
				}
				matDist[i*n+j] = buffDouble[vectInd];
			}
	}

	//Fill the upper part of the symmetric distance matrix
	for (int i = 0; i<n; i++) 
		for (int j = 0; j<i; j++)
			matDist[i*n+j] = matDist[j*n+i];

	delete aliSeqs; delete buffDouble; delete buffChar1; delete buffChar2; //Free memory from dynamic allocations
}

void listGenomes(string alig){
	/*ifstream inFile(fListe.c_str());
	if(!inFile)
	{
		cout<<"Couldn't open input fasta file"<<endl;
		exit(1);
	}*/
	istringstream inFile(alig);
	string line;
	int i=0;
	while(getline(inFile, line)){
		if(line=="")
			continue;
        if(line.at(0)=='>'||line.at(0)==';'){
			line.erase(0,1);
			#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
				if(line[strlen(line.c_str())-1]=='\r')
					line[strlen(line.c_str())-1]='\0';
			#endif
			nc[i]=line.c_str();
			i++;
		}
	}
}

int similarity(string fListe, string alignMode="maxPrecision") {
	/*
	This function returns the similarity matrix.
	Input:
	-fListe : Name of the fasta file holding a list of sequences with their names.
	Output:
	-MatSimil : The similarity matrix (global variable).
	-nc : Sequential number with the names of the sequences (global variable).
	Output:
	-0 = success or 1 = alignment failed
	*/

	char cmd[300] = "";
	if(!installed){
	#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Under Windows OS
		strcat(cmd,"\"");
	#endif
		strcat(cmd,"\"");
		strcat(cmd,progPath);
		strcat(cmd,"muscle\" ");
	}else
		strcat(cmd,"muscle ");
	if(tolowerstr(alignMode)=="maxprecision")
		strcat(cmd,"-quiet -in ");
	else if(tolowerstr(alignMode)=="moderate")
		strcat(cmd,"-quiet -maxiters 4 -in ");
	else if(tolowerstr(alignMode)=="fast")
		strcat(cmd,"-quiet -maxiters 2 -in ");

	string liste="";

	if(tolowerstr(alignMode)!="none"){
		//Perform the alignment
		strcat(cmd,"\"");
		strcat(cmd, fListe.c_str());
		strcat(cmd,"\"");
		//cout<<cmd<<endl;
		liste = exec(cmd);
		// Save the result
		/*std::ofstream out("aligned.txt");
		out << liste;
		out.close();*/
	} else {
		// Read the aligned sequences fron input file
		ifstream ifs(fListe.c_str());
		string inliste( (std::istreambuf_iterator<char>(ifs) ),
						   (std::istreambuf_iterator<char>()    ) );
		liste=inliste;
	}
	if(liste==""){
		return 1;
	}
	int index;
	string **dicoMuscle,tempGen;
	double max = 0.0;
	nbSequences=0;
	for (int i = 0; i<liste.size(); i++)
		if ((liste[i] == '>' && liste[i+1] != '>') || (liste[i] == ';' && liste[i+1] != ';')) //Taking into consideration the error of duplication by checking liste[i+1]
			nbSequences++;
	//Allocate space for pointers to rows of a matrix
	if (!(MatSimil = new double[nbSequences*nbSequences])){
		cout << "Allocation for 'MatSimil' failed. \n";
		return 3;
	}
	if (!(dicoMuscle = new string*[nbSequences])){
		cout << "Allocation for 'dicoMuscle' failed. \n";
		return 3;
	}
	if (!(nc = new string[nbSequences])){
		cout << "Allocation for 'nc' failed. \n";
		return 3;
	}

	//initialize the list of sequences references as in input file order
	listGenomes(liste);
	
	//Shuffle the sequences order to minimize the ordering interference on GMM
	customRND myrnd;
	random_shuffle(nc,nc+nbSequences-1,[&](int i){return (myrnd.rnd()%i);});
	
	for (int i = 0; i<nbSequences; i++) {
		dicoMuscle[i] = new string[2];
		dicoMuscle[i][0]="";//Initialize to empty string in order to detect if the current row is used in the case of genes with duplicate name
		dicoMuscle[i][1]="";
	}
	for (int i = 0; i<liste.size();) {
		if (liste[i] == '>'||liste[i]==';') {
			tempGen="";
			index=0;
			i++;
			do {
				tempGen += liste[i];
				i++;
			} while (liste[i] != '\n'&&liste[i]!='\r');
			while(tempGen.compare(nc[index])!=0 || dicoMuscle[index][0]!=""){//strcmp(nc[index].c_str(),tempGen.c_str())!=0){
				index++;
			}
			dicoMuscle[index][0]=tempGen;
		}
		else {
			do {
				if (liste[i] != '\n' && liste[i] != '\r') {
					dicoMuscle[index][1] += liste[i];
				}
				i++;
			} while (liste[i] != '>'&&liste[i] != ';'&&i<liste.size());
		}
	}

	//Validate that the aligned sequences all have the same size
	int aligseqsize = dicoMuscle[0][1].length();
	
	for (int i = 1; i<nbSequences; i++){
		if(dicoMuscle[i][1].length()!=aligseqsize){
			//Free memory from dynamically allocated matrixes
			for (int i = 0; i<nbSequences; i++) {
				delete[] dicoMuscle[i];
			}
			delete[] dicoMuscle;
			delete[] nc;
			return 2;
		}
	}
	
	matriceDistances(dicoMuscle, MatSimil, nbSequences); //Calculates the distance matrix
	
	for (int i = 0; i<nbSequences; i++)
		for (int j = 0; j<nbSequences; j++)
			if (MatSimil[i*nbSequences+j]>max)
				max = MatSimil[i*nbSequences+j];
	for (int i = 0; i<nbSequences; i++)
		for (int j = 0; j<nbSequences; j++)
			MatSimil[i*nbSequences+j] /= max;
	for (int i = 0; i<nbSequences; i++)
		for (int j = 0; j<nbSequences; j++)
			MatSimil[i*nbSequences+j] = 1 - MatSimil[i*nbSequences+j];

	//Free memory from dynamically allocated matrixes
	for (int i = 0; i<nbSequences; i++) {
		delete[] dicoMuscle[i];
	}
	delete[] dicoMuscle;
	return 0;
}

void killWorkers(){ //Send a message to the worker processes that the work is done
	progEnd=1;
	for(int i=1; i<numprocs; i++)
		MPI_Send(&progEnd,1,MPI_INT,i,99,MPI_COMM_WORLD);
}

string getClusteringString(int number_data, int nbClusters, string *refs, int *p){
	/*for (int i = 0; i < nbClusters; i++){
		clustersString +="[";
		for (int j = 0; j < number_data; j++)
			if(p[j]==i)
				clustersString+="'" + refs[j] + "', ";
		clustersString = clustersString.substr(0, clustersString.size()-2);
		clustersString +="]\r\n";
	}*/
	vector<vector<string>> clusters(nbClusters);
	string clustersString="";

	for(int i = 0;i<number_data;i++)
		clusters[int(p[i])].push_back(refs[i]);

	for (int i = 0; i < nbClusters; i++){
		clustersString +="[";
		for (int j = 0; j < clusters[i].size(); j++)
			clustersString += "'" + clusters[i][j] + "', ";
		clustersString = clustersString.substr(0, clustersString.size()-2);
		clustersString +="]\r\n\r\n";
	}
	return clustersString;
}

string GMM_Clustering(double **data, string *refs, int number_data, int number_features, string ccCriterion = "bestBIC", int nbRuns = 500, int neStop = 50, string criterion = "BIC", string type_covariance="full", int max_iterations=1000, int nbClusters=-1){

	int *p = new int[number_data];
	bool calculateBestNbClusters=(nbClusters==-1); //if nbClusters != -1 then it has not been calculated and provided by other methods. It will be calculated here based on the best BIC.
	GaussianMixture *gmm;
	string clustersString="";

	if(tolowerstr(ccCriterion)=="fast"){

		int rndseed=320; // Fixed seed chosen randomly in order to avoid getting different clustering at each run
		if(calculateBestNbClusters)
			nbClusters = getBestNbClusters(data, number_data, number_features, rndseed, criterion, type_covariance, max_iterations);
		gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, rndseed);    
		(*gmm).train(data,p,number_data,number_features);
		clustersString = getClusteringString(number_data, nbClusters, refs, p);
		delete gmm;

	} else if(tolowerstr(ccCriterion)=="bestbic"){

		int neTries=0;
		double bestBIC=INFTY;
		int bestSeed;
		//int totalTries=0;
		for(int i=0; i<nbRuns && neTries<neStop; i++){ // "i" will be used as the random seed for each loop
			//totalTries++;
			if(calculateBestNbClusters)
				nbClusters = getBestNbClusters(data, number_data, number_features, i, criterion, type_covariance, max_iterations);
			gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, i);    
			(*gmm).train(data,p,number_data,number_features);
			if(bestBIC>(*gmm).bic()){
				bestBIC=(*gmm).bic();
				neTries=1;
				bestSeed=i;
			} else {
				neTries++;
			}
			delete gmm;
		}
		//cout<<"Stopped after "<<totalTries<<" tries\n";
		if(calculateBestNbClusters)
			nbClusters = getBestNbClusters(data, number_data, number_features, bestSeed, criterion, type_covariance, max_iterations);
		gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, bestSeed);    
		(*gmm).train(data,p,number_data,number_features);
		clustersString = getClusteringString(number_data, nbClusters, refs, p);
		delete gmm;

	} else{ // The choice is mostFreq

		string *cStrings = new string[nbRuns];
		int *occurences = new int[nbRuns];
		double *bestBIC = new double[nbRuns];
		for(int i=0;i<nbRuns;i++){
			cStrings[i]="";
		}
		int highestFreq=0, bestCluteringIndex;

		for(int i=0; i<nbRuns; i++){ // i will be used as the random seed for each loop
			if(calculateBestNbClusters)
				nbClusters = getBestNbClusters(data, number_data, number_features, i, criterion, type_covariance, max_iterations);
			gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, i);    
			(*gmm).train(data,p,number_data,number_features);
			clustersString = getClusteringString(number_data, nbClusters, refs, p);

			for(int j=0; j<nbRuns; j++){
				if(clustersString==cStrings[j]){
					occurences[j]++;
					if(bestBIC[j]>(*gmm).bic())
						bestBIC[j]=(*gmm).bic();
					break;
				} else if(cStrings[j]==""){
					cStrings[j]=clustersString;
					occurences[j]=1;
					bestBIC[j]=(*gmm).bic();
					break;
				}
			}
			delete gmm;
		}

		for(int i=0; i<nbRuns; i++){ // Getting the index of the best clustering: the most frequently occurent one, and if many we take the one with the best BIC
			if(cStrings[i]==""){//cout<<"stopped at "<<i<<endl;
				break;
			}else if(occurences[i]>highestFreq){
				highestFreq=occurences[i];
				bestCluteringIndex=i;
			} else if(occurences[i]==highestFreq){
				if(bestBIC[i]<bestBIC[bestCluteringIndex])
					bestCluteringIndex=i;
			}
		}//cout<<"\nThe best index is "<<bestCluteringIndex<<" and the highest freq is "<<highestFreq<<endl;

		clustersString=cStrings[bestCluteringIndex];
		delete[] cStrings; delete[] occurences; delete[] bestBIC;

	}

	printInternalValidationIndices(data,p,number_data,number_features);

	delete[] p;

	return clustersString;	
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	
	dCount=new int[numprocs];
	//Setting default value and reading the input arguments needed by all processes
	mDist="EDNAFULL";
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-mdist")==0){
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				MPI_Finalize();
				return 0;
			} else
				mDist=argv[i+1];
		}
	gapOpenI=-10;
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-gapopen")==0)
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				MPI_Finalize();
				return 0;
			} else
				gapOpenI=atof(argv[i+1]);
	
	
	gapExtendI=-0.5;
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-gapextend")==0)
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				MPI_Finalize();
				return 0;
			} else
				gapExtendI=atof(argv[i+1]);

	if(myid==0){ //This is the main process
		string fListe,fGroupes;
		string groupes;
		string alignMode;
		string matType;
		string ccCriterion;
		string nbRunsStr;
		int nbRuns;
		string neStopStr;
		int neStop;

		//Displaying usage message in case called without arguments
		if(argc==1){
			//cout<<"You are running "<<PACKAGE_STRING<<endl;
			cout<<"SpClust performs biological sequences clustering using GMM.\n\n"
			<<"usage: mpispclust -in [input fasta file] -out [output clustering file] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop] -matType [affinity matrix type]\n"
			<<"   or: mpiexec -n [number of slave processes] spclust -in [input fasta file] -out [output clustering file] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop] -matType [affinity matrix type]\n\n"
			<<"Available scoring matrices are EDNAFULL, BLOSUM62, and PAM250. Defaults to EDNAFULL if not specified.\n"
			<<"Available alignment modes are: none, fast, moderate, and maxPrecision. none considers the input sequences aligned and does not perforn any alignment. fast and moderate limit the number of iterations for the alignment to 2 and 4 respectively (using MUSCLE). Defaults to maxPrecision if not specified.\n"
			<<"Available clustering choice criteria are: bestBIC, mostFreq, and fast. Defaults to bestBIC.\n"
			<<"   nbRuns and neStop are both used for the clustering choice criterion bestBIC to indicate the stop condition for the best BIC choice.\n"
			<<"   neStop is ignored for the clustering choice criterion mostFreq. GMM will be run nbRuns times and the most occurent clustering will be chosen.\n"
			<<"   nbRuns and neStop are ignored for the clustering choice criterion fast. GMM will be run only once.\n"
			<<"   If not specified, nbRuns defaults to 500 and neStop defaults to 50.\n"
			<<"Available affinity matrices types are UL (Unnormalized Laplacian), RWNL (Random Walk Normalized Laplacian), MOD (Modularity), and BH (Bethe Hessian). Defaults to RWNL if not specified.\n\n"
			/*<<"Note: parameters are case sensitive, e.g. using -alignmode instead of -alignMode will cause this parameter to be disregarded, and using blosum62 instead of BLOSUM62 will be mentioned as an error.\n\n"*/;
			killWorkers();
			MPI_Finalize();
			return 0;
		}

		//Checking if our program is installed
	#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		char chkcmd[20];
		string spclustres, /*spclustGMMres, */muscleres;
		strcpy(chkcmd,"whereis spclust");
		spclustres = exec(chkcmd);
		/*strcpy(chkcmd,"whereis spclustGMM");
		spclustGMMres = exec(chkcmd);*/
		strcpy(chkcmd,"whereis muscle");
		muscleres = exec(chkcmd);
		if(spclustres.length()>9 && /*spclustGMMres.length()>12 && */muscleres.length()>9){
			installed=true;
			//Setting the working directory
			strcpy(progPath,exec("pwd").c_str());
			progPath[strlen(progPath)-1]='/';
		}
	#endif
	
		if(!installed){
		//Setting the running path of our executables
		#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
			realpath(argv[0], progPath);
		#else			
			_fullpath(progPath, argv[0], sizeof(progPath));
		#endif
			progPath[strlen(progPath)-7]='\0';
		}

		//Setting the default values for the arguments
		fListe = string(progPath);
		fListe += "sequences.fasta";
		fGroupes = string(progPath);
		fGroupes += "Clustering.txt";
		alignMode = "maxPrecision";
		ccCriterion = "bestBIC";
		nbRunsStr = "500";
		neStopStr = "50";
		matType = "RWNL";

		//Reading the input arguments specific to the main process only
		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-in")==0)
				if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				killWorkers();
				MPI_Finalize();
				return 0;
			} else
				fListe=argv[i+1];

		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-out")==0)
				if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				killWorkers();
				MPI_Finalize();
				return 0;
			} else
				fGroupes=argv[i+1];

		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-alignmode")==0)
				if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				killWorkers();
				MPI_Finalize();
				return 0;
			} else
				alignMode=argv[i+1];

		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-cccriterion")==0)
				if(i+1>=argc){
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
					killWorkers();
					MPI_Finalize();
					return 0;
				} else
					ccCriterion=argv[i+1];

		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-nbruns")==0)
				if(i+1>=argc){
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
					killWorkers();
					MPI_Finalize();
					return 0;
				} else
					nbRunsStr=argv[i+1];

		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-nestop")==0)
				if(i+1>=argc){
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
					killWorkers();
					MPI_Finalize();
					return 0;
				} else
					neStopStr=argv[i+1];

		for(int i=1; i<argc; i+=2)
			if(strcmp(tolowercstr(argv[i]).c_str(),"-mattype")==0)
				if(i+1>=argc){
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
					killWorkers();
					MPI_Finalize();
					return 0;
				} else
					matType=argv[i+1];

		//Validate the input arguments
		if(tolowerstr(mDist)!="ednafull" && tolowerstr(mDist)!="blosum62" && tolowerstr(mDist)!="pam250"){
			cout<<"Error: Invalid distance matrix. Available matrices are: EDNAFULL, BLOSUM62, and PAM250.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		
		if(tolowerstr(alignMode)!="fast" && tolowerstr(alignMode)!="moderate" && tolowerstr(alignMode)!="maxprecision" && tolowerstr(alignMode)!="none"){
			cout<<"Error: Invalid alignmnet mode. Available modes are: fast, moderate, maxPrecision, and none.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		
		if (fListe.substr(fListe.find_last_of(".") + 1) != "fasta" && fListe.substr(fListe.find_last_of(".") + 1) != "dat"){
			cout << "Error: Invalid input filename! Please input a fasta file containing the genomes using the argument -in.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		/*// Removed restriction on the output file extension
		if (fGroupes.substr(fGroupes.find_last_of(".") + 1) != "txt" && fGroupes.substr(fGroupes.find_last_of(".") + 1) != "dat"){
			cout << "Error: Invalid output filename! Please name a txt file for the results using the argument -out.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		*/
		
		if(tolowerstr(ccCriterion)!="fast" && tolowerstr(ccCriterion)!="bestbic" && tolowerstr(ccCriterion)!="mostfreq"){
			cout<<"Error: Invalid clustering choice criterion. Available criteria are: bestBIC, mostFreq, and fast.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		
		if(!is_number(nbRunsStr)){
			cout<<"Error: The input number of runs must be a positive integer.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		nbRuns = atoi(nbRunsStr.c_str());
		
		if(!is_number(neStopStr)){
			cout<<"Error: The input no-emprovement stop parameter must be a positive integer.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		neStop = atoi(neStopStr.c_str());
		
		if(tolowerstr(matType)!="ul" && tolowerstr(matType)!="rwnl" && tolowerstr(matType)!="mod" && tolowerstr(matType)!="bh"){
			cout<<"Error: Invalid affinity matrix choice. Available choices are: LU, RWNL, MOD, and BH.\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
	
		//Check if the input fasta file is accessible and the required modules are present
		std::ifstream infile(fListe.c_str());
		if(!infile.good()){
			cout<<"Error: The input fasta file is missing or not accessible.\n\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
		char modsPath[300]="";//For storing the paths of the required files in order to check their existance		
	#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		//Check muscle
		strcpy(modsPath,progPath);
		strcat(modsPath,"muscle");
		std::ifstream infile2(modsPath);
		if(!infile2.good() && !installed){
			cout<<"Error: muscle executable is missing or not accessible.\n\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}

	#else
		//Check muscle
		strcpy(modsPath,progPath);
		strcat(modsPath,"muscle.exe");
		std::ifstream infile2(modsPath);
		if(!infile2.good()){
			cout<<"Error: muscle executable is missing or not accessible.\n\n";
			killWorkers();
			MPI_Finalize();
			return 0;
		}
	#endif

		//Calculate the similarity matrix
		printTimestamp("Alignment and similarity matrix calculation started at: ");
		int errres = similarity(fListe,alignMode); //Capture the return value to check for errors

		//Check for alignment errors
		if(errres==1){
			killWorkers();
			cout<<"Sequences alignment failed. Job aborted!\n";
			MPI_Finalize();
			return 0;
		} else if(errres==2){
			killWorkers();
			cout<<"Input aligned sequences must have the same size. Job aborted!\n";
			MPI_Finalize();
			return 0;
		} else if(errres==3){
			killWorkers();
			cout<<"Memory allocation error. Job aborted!\n";
			MPI_Finalize();
			return 0;
		}
		printTimestamp("Alignment and similarity matrix calculation ended at: ");
/*
		string matriceSimilitude="";
		string refsGenomes="";
	
		for(int i=0; i<nbSequences; i++){
			for(int j=0; j<nbSequences; j++){
				matriceSimilitude+=getStringFromDouble(MatSimil[i*nbSequences+j]);
				matriceSimilitude+=" ";
			}
			matriceSimilitude+="\r\n";
		
			refsGenomes+=nc[i];
			refsGenomes+=",";
		}
	
		std::ofstream out1("matSimil.txt");
		out1 << matriceSimilitude;
		out1.close();
		std::ofstream out2("refs.txt");
		out2 << refsGenomes;
		out2.close();
*/
		
		printTimestamp("Affinity matrix calculation started at: ");
		//Calculate the affinity matrix
		toAffinityMat(MatSimil, nbSequences,matType);
		
		printTimestamp("Affinity matrix calculation ended and features matrix calculation started at: ");
		
		int chNbEigenVects, nbClusters=-1; //if nbClusters is provided as -1 to GMM_Clustering() then it will be calculated within GMM_Clustering() based on the best BIC.
		double **vecPropT;
		if (!(vecPropT = new double*[nbSequences])){ //The features matrix has number of rows equal to the number of nodes (sequences). The number of columns are allocated in the next function.
			killWorkers();
			cout<<"Allocation for 'vecPropT' failed. \n";
			MPI_Finalize();
			return 0;
		}

		string EVMethod="delta";
		if(tolowerstr(matType)=="mod")
			EVMethod="log";

		errres=featuresCalc(MatSimil,vecPropT,chNbEigenVects,nbClusters,nbSequences,EVMethod,0.01,matType);
		//Check for errors
		if(errres==1){
			killWorkers();
			cout<<"Memory allocation error. Job aborted!\n";
			MPI_Finalize();
			return 0;
		}
		

		//Free the memory used by MatSimil to leave enough memory for the GMM
		if(MatSimil!=NULL)
			delete [] MatSimil;

		printTimestamp("Features calculation ended and clustering started at: ");
		string clustersString = GMM_Clustering(vecPropT, nc, nbSequences, chNbEigenVects, ccCriterion, nbRuns, neStop, "BIC", "full", 1000, nbClusters);
		printTimestamp("Clustering ended at: ");
		ofstream outfile;
		outfile.open(fGroupes.c_str());
		outfile<<clustersString;
		outfile.close();
		//cout<<"done";system("pause");

		//Free danamically allocated memory
		for(int i=0; i<nbSequences; i++)
			if(vecPropT[i]!=NULL)
				delete []vecPropT[i];
		delete []vecPropT;
		if(nc!=NULL)
			delete[] nc;

	} else { //This is a calculation worker process

		MPI_Recv(&progEnd, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status); //Check if there is no work to do
		if(progEnd==1){
			MPI_Finalize();
			return 0; //If the is no work then exit (end process)
		}

		MPI_Recv(&nbSequences, 1, MPI_INT, 0, 98, MPI_COMM_WORLD, &status); //Receive the aligned sequence size
		MPI_Recv(&Aligned_seq_length, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status); //Receive the aligned sequence size
		aliSeqs=new char[nbSequences*(Aligned_seq_length+1)+1]; //Allocate aliSeqs to receive the aligned sequences
		MPI_Recv(aliSeqs, nbSequences*(Aligned_seq_length+1)+1, MPI_CHAR, 0, 96, MPI_COMM_WORLD, &status); //Receive the aligned sequences
		
		buffChar1=new char[Aligned_seq_length+1];
		buffChar2=new char[Aligned_seq_length+1];
		
		//Do the required calculations from this slave
		slaveDistCalc();

		//Send the results to the master
		MPI_Send(buffDouble,dCount[myid],MPI_DOUBLE,0,myid*100+11,MPI_COMM_WORLD);

		delete aliSeqs; delete buffDouble; delete buffChar1; delete buffChar2; //Free dynamic allocations
	}
	delete dCount;
	MPI_Finalize();
	return 0;
}
