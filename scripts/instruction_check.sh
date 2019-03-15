#!/bin/bash

# In the instructions.c, you can uncomment the PRINT_INSTRUCTIONS macro.
# This will print out the instructions in the output file during runs.
# This script can be used to compare the used instructions between different runs.
# The script expects two arguments (i.e. file 1 and file 2)

begin1=($(grep -n "#START INSTR" $1 | sed 's/:.*$//'))
end1=($(grep -n "#END INSTR" $1 | sed 's/:.*$//'))
begin2=($(grep -n "#START INSTR" $2 | sed 's/:.*$//'))
end2=($(grep -n "#END INSTR" $2 | sed 's/:.*$//'))

if [ ${#begin1[@]} -ne ${#end1[@]} ]; then
    >&2 echo "Unfinished instructions in file $1, skipping the incomplete instructions"
fi
if [ ${#begin2[@]} -ne ${#end2[@]} ]; then
    >&2 echo "Unfinished instructions in file $2, skipping the incomplete instructions"
fi

if [ ${#begin1[@]} -ne ${#begin2[@]} ]; then
    >&2 echo "File $1 has ${#begin1[@]} while file $2 has ${#begin2[@]}."
    >&2 echo "Will only compare overlapping set of instructions."
    NR_INSTR=$((${#begin1[@]}<${#begin2[@]}?${#begin1[@]}:${#begin2[@]}))
else
    NR_INSTR=${#begin1[@]}
fi

echo "Comparing $NR_INSTR sets of instructions"
differences=0
for i in $(seq 0 $(($NR_INSTR-1))); do
	difference1=$((${end1[$i]} - ${begin1[$i]}))
	difference2=$((${end2[$i]} - ${begin2[$i]}))

	echo "Instruction set $i"
	if [ $difference1 -ne $difference2 ]; then
	    >&2 echo "Set has different amount of instructions ($difference1 neq $difference2)"
	    differences=$(($differences+1))
	else
	    # This is faster than first head then tail.
	    # While tail is running, head can already start reading. 
	    # Once head has read enough lines, it will send a signal to 
	    # terminate tail, and sort can be started.
	    #
	    # If doing first head, then tail, this is not possible.
	    instr1=$(tail -n+${begin1[$i]} $1 | head -n$difference1 | sort)
	    instr2=$(tail -n+${begin2[$i]} $2 | head -n$difference2 | sort)
	    diff <(echo "$instr1") <(echo "$instr2")
	    if [ $? -ne 0 ]; then
		differences=$(($differences+1))
	    fi
	fi
done

echo "$differences different sets were found."
exit $differences
