n=${1}
for((i=0;i<${n};++i)) ; do
	qsub submit_simple.sh
done
