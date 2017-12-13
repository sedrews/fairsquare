PARAMS="-mi -mf -z -e 0.15"
FILELIST=$1
RESDIR=$2
TIMEOUT=$3
mkdir $RESDIR
if [ $? -ne 0 ]; then
    echo Specify a non-existing directory to create and store output files
    exit
fi
cd ../src/
for F in $(cat ../oopsla/$FILELIST); do
    OUT=$(echo $F | sed 's/^.*\///' | sed 's/\.fr$//')
    CMD="python fairProve.py -f ../oopsla/$F $PARAMS -o ../oopsla/$RESDIR/$OUT"
    if [ -z $TIMEOUT ]; then
        $CMD
    else
        trap 'kill -INT -$pid' INT
        timeout $TIMEOUT $CMD &
        pid=$!
        wait $pid
    fi
done
