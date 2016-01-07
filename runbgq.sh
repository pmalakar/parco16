#!/bin/bash -x
#COBALT --disable_preboot
 
#export L1P_POLICY=std
#export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

#Free bootable blocks
boot-block --reboot
 
PROG=iotree
NODES=$1

#for PROG in iotree #iotree_2coll iotree_1coll iotree_gather
#do
for iter in 1 # 2 
do
for THRD in 1 #2 4 8 16 
do
for ppn in 16 #8 #4 #2 1
do
for MSG in 512 #512 #1024 #2048 4096 8192
#for MSG in 1 2 4 512 1024 #2048 4096 8192
do 
 for coalesced in 1 0 
 do
 for streams in 1 0 #2
 do
 for blocking in 0 #1
 do
 rm -f dummy*
 for type in 0 1 #2 
 do
	if [ $type -gt 0 ] && [ $coalesced -eq 1 ] 
	then
			continue;
	fi
	if [ $type -ge 0 ] && [ $streams -gt 0 ] && [ $coalesced -eq 0 ] 
	then
			continue;
	fi
	RANKS=`echo "$NODES*$ppn"|bc`
	OUTPUT=${PROG}_N${NODES}_R${ppn}_${MSG}_${coalesced}_${blocking}_${type}_${streams}
	rm -f ${OUTPUT}.cobaltlog ${OUTPUT}.output ${OUTPUT}.error
	echo 
	echo "* * * * *"
	echo 
	echo "Starting $OUTPUT with numthreads=$THRD ppn=$ppn args=${MSG} ${coalesced} ${blocking} ${type} ${streams}"
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=4 MUSPI_NUMRECFIFOS=4 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_4
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=8 MUSPI_NUMRECFIFOS=8 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_8
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=4194304" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_4M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=8388608" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_8M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=16777216" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_16M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=4194304" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_4M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=8388608" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_8M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=16777216" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_16M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}

	continue;

  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV_LOCAL=4M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_local_4M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV_LOCAL=8M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_local_8M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=4M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_remote_4M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=8M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_remote_8M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs PAMID_RZV=16M PAMID_STATISTICS=1 PAMID_VERBOSE=1 : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_remote_16M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=4M" --envs "PAMID_RZV_LOCAL=4M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_4M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=8M" --envs "PAMID_RZV_LOCAL=8M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_8M
  	#rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=4M" --envs "PAMID_RZV_LOCAL=4M" --envs "PAMID_THREAD_MULTIPLE=1" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_ptm
  	#rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_THREAD_MULTIPLE=1" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_ptm
	echo 
	echo "* * * * *"
	echo
done 
done 
done 
done 
done 
done 
done 
done

exit

