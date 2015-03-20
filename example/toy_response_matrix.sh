clean=${1}
exe=example/toy_response_matrix

ofile=example/response_matrix.dat
if [ "${clean}" == "clean" ] ; then rm ${ofile}; fi
if [ ! -e ${ofile} ] ; then
	${exe} -n 100000000 -o ${ofile}
fi

# ofile=response_matrix.root
# if [ ! -e ${ofile} ] ; then
	# ${exe} -n 100000000 -o ${ofile}
# fi
 
# ofile=response_matrix_gun.root
# if [ ! -e ${ofile} ] ; then
	# ${exe} -n 100000000 -gun -o ${ofile}
# fi

