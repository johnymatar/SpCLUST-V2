n=$(nproc --all)
BASEDIR=$(dirname "$0")
arg1=$1
arg2=$2
arg3=$3
arg4=$4
arg5=$5
arg6=$6
arg7=$7
arg8=$8
mpirun -np $n --allow-run-as-root $BASEDIR/spclust -in $arg1 -out $arg2 -alignmode $arg3 -mdist $arg4 -cccriterion $arg5 -nbruns $arg6 -nestop $arg7 -mattype $arg8 
