- Download the Linux standalone version and extract it to the home directory (SpCLUST-V2's executables should be reachable in ~/SpCLUST-V2/ path)
  In case "spclust" executable in the standalone package does not run on your machine, it should be rebuilt from the source files (mpic++ spclust.cpp -o spclust)
- Edit the spclust-v2.xml file and fix the full path of the executable in the command tag at line 4 (e.g. /home/[your user directory]/SpCLUST-V2/spclust-v2.sh $input $output $...)
- Create ~/galaxy/tools/SpCLUST-V2/ directory and place the spclust-v2.xml file in it
- Edit ~/galaxy/config/tool_conf.xml to add SpCLUST-V2 tool (a sample entry is present at the end of the posted tool_conf.xml file)
