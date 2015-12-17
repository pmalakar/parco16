#!/bin/bash -x
#COBALT --disable_preboot
 
#export L1P_POLICY=std
#export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

#Free bootable blocks
boot-block --reboot
 
PROG=iotree
NODES=$1

for iter in 1 # 2 
do
for THRD in 1 #2 4 8 16 
do
for ppn in 8 # #4 #2 1
do
for MSG in 16 #32 64 128 #256 #1024 # 8 16 32 64 # 128 256 512 #384 512 768 1024 2048 4096 8192
do 
 for coalesced in 0 1 
 do
 for blocking in 0 #1
 do
 rm -f dummy*
 for type in 2 1 0
 do
	if [ $type -gt 0 ] && [ $coalesced -eq 1 ] 
	then
			continue;
	fi
	RANKS=`echo "$NODES*$ppn"|bc`
	OUTPUT=${PROG}_N${NODES}_R${ppn}_${MSG}_${coalesced}_${blocking}_${type}
	rm -f ${OUTPUT}.cobaltlog ${OUTPUT}.output ${OUTPUT}.error
	echo 
	echo "* * * * *"
	echo 
	echo "Starting $OUTPUT with numthreads=$THRD ppn=$ppn args=${MSG} ${coalesced} ${blocking} ${type}"
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} > ${OUTPUT}
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

exit

