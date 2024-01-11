#!/bin/bash

usage(){
	echo ""
	echo "Usage: $0 all [pid] | PID1 PID2 ... PIDn"
	echo ""
	echo "	all	: lists all the AF2 folders for all the current GPU running PIDs"
	echo "	pid	: adds the PIDs for each AF2 folder"
	echo "	PIDn	: a list of PIDs to get AF2 folder locations"
	echo ""
	exit 1
}

[ $# -lt 1 ] && usage

# Pass all as first argument to get all PIDs
if [ "$1" == all  ]; then
	pids_list=($(nvidia-smi --query-compute-apps=pid --format=csv,noheader))
else
	pids_list="$@"
fi


# Pass PIDs as positional arguments
for arg in ${pids_list[@]}; do
	[ "$2" == "pid" ] && echo "PID $arg:"
	ps aux | grep $arg | grep -o '\./AF2/[^ ]*'
	[ $? -ne 0 ] && echo "No process found"
done
