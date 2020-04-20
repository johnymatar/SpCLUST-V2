# SpCLUST-V2
SpCLUST-V2 is a package for divergent biological sequences clustering. Conversely from traditional clustering methods that focuses on the speed of clustering highly similar sequences, SpCLUST uses a Machine Learning Gaussian Mixture Model and targets the clustering accuracy of divergent sequences with the best possible speed.
SpCLUST-V2 includes the option of using Edgar, R.C.'s MUSCLE module (www.drive5.com) for sequences alignment.

# Prerequisite
SpCLUST-V2 uses MPI for parallel computation and the executable building for the installation package. Below are some basic instructions for installing MPI on your system.
- For Linux users:
  •	install mpich
  •	install openmpi
  •	install openmpi-devel
  •	echo "export PATH=$PATH:/usr/lib64/openmpi/bin" >> ~/.bashrc
- For Windows users:
  •	Download MS-MPI SDK and Redist installers from Microsoft's website: https://msdn.microsoft.com/en-us/library/bb524831.aspx
  •	Install the downloaded packages

# Installation on Linux
- Get the installation package from the "Linux" folder in our repository: "wget https://github.com/johnymatar/SpCLUST-V2/raw/master/Linux/install.tar.xz"
- Extract the package: "tar -xvf install.tar.xz"
- Run the following commands: "cd SpCLUST-V2", "./configure", "make"
- Run the following command as a sudoer: "make install"
- You can now call the executables with the desired arguments
- For serial computation use "spclust" with the desired arguments
- For parallel computation use "mpispclust" with the desired arguments
- To use the graphical interface, install mono (run "apt install mono-complete" as a sudoer) and then call "guispclust"

# Usage without installation on Linux
- Get the standalone package from the "Linux" folder in our repository: "wget https://github.com/johnymatar/SpCLUST-V2/raw/master/Linux/standalone.tar.xz"
- Extract the package: "tar -xvf standalone.tar.xz"
- Keep the extracted files together in a same directory and, for each use, browse to that directory from the console: e.g. "cd ~/SpCLUST-V2"
- For serial computation use "./spclust" with the desired arguments
- For parallel computation use "./mpispclust" with the desired arguments
- To use the graphical interface, install mono (run "apt install mono-complete" as a sudoer) and then call "./guispclust" or "mono ./guispclust"

# Usage without installation on M.S. Windows
- Download and extract SpCLUST-V2.zip from the "Windows" folder in our repository
- For using the command line interface, browse to the extracted folder from your open CMD
- Call "spclust" with the desired arguments (from the open CMD)
- For parallel computation call "mpispclust" with the desired arguments (from the open CMD)
- For using the graphical interface, run "guispclust"

# Integration into Galaxy
Follow the instructions in the README file in the Galaxy folder

# Current version features
- Cross-platform: tested on Linux and Windows. Although not tested, the source files should also compile and run on MAC OS
- Portable (can be used without installation)
- Parallel computation for the distance matrix using MPI, thus enabling its use on omputation clusters.
- Compatible for Galaxy integration
- usage: mpispclust -in [input fasta file] -out [output clustering file] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop]  -matType [affinity matrix type]
     or: mpiexec -n [number of slave processes] spclust -in [input fasta file] -out [output clustering file] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop]  -matType [affinity matrix type]

  Available scoring matrices are EDNAFULL, BLOSUM62, and PAM250. Defaults to EDNAFULL if not specified.
  Available alignment modes are: none, fast, moderate, and maxPrecision. none considers the input sequences aligned and does not perforn any alignment. fast and moderate limit the number of iterations for the alignment to 2 and 4 respectively (using MUSCLE). Defaults to maxPrecision if not specified.
  Available clustering choice criteria are: bestBIC, mostFreq, and fast. Defaults to bestBIC.
    nbRuns and neStop are both used for the clustering choice criterion bestBIC to indicate the stop condition for the best BIC choice.
    neStop is ignored for the clustering choice criterion mostFreq. GMM will be run nbRuns times and the most occurent clustering will be chosen.
    nbRuns and neStop are ignored for the clustering choice criterion fast. GMM will be run only once.
    If not specified, nbRuns defaults to 500 and neStop defaults to 50.
  Available affinity matrices types are UL (Unnormalized Laplacian), RWNL (Random Walk Normalized Laplacian), MOD (Modularity), and BH (Bethe Hessian). Defaults to RWNL if not specified.
