#!/bin/bash
# Run script example: ./script.sh genera/sulfate_reducing > sulfate_reducing_orders
while read thing; do 
	taxonomy=$(grep "$thing" fullnamelineage.dmp | 
	awk 'BEGIN {FS="\t"}; {print $5}' | 
	awk 'NR==3');
	#echo $taxonomy
	numgroup=$(echo $taxonomy |awk 'BEGIN {FS="; "}; {print $1 $2 $3 $4 $5}' | grep -o "group"| wc -l)
	#echo $numgroup
	if [[ "$numgroup" == 2 ]]; then
		order=$(echo $taxonomy | awk 'BEGIN {FS="; "}; {print $7}')
		#echo "two groups found"
	elif [[ "$numgroup" == 1 ]]; then
		#echo "one group found"
		order=$(echo $taxonomy | awk 'BEGIN {FS="; "}; {print $6}');
	elif echo $taxonomy | grep -q ";"; then
		#echo "no groups found"
		order=$(echo $taxonomy | awk 'BEGIN {FS="; "}; {print $5}');
	else
                order="Unknown Order"
		#echo "no taxonomy found"
	fi;
	order=$(echo $order | sed 's/;//')
	echo -e $order "\t" $thing;
done < $1
