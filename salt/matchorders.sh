#!/bin/bash
while read line; do
	#echo $line
	order=$(echo $line | awk 'BEGIN {FS=" "}; {print $2}')
	virus=$(echo $line | awk 'BEGIN {FS=" "}; {print $1}')
	#echo $order
	type=""
	sr=$(grep -c "$order" orders/sulfredlist); if [ $sr -gt 0 ]; then type="Sulfate Reducing"; fi
	so=$(grep -c "$order" orders/sulfoxlist); if [ $so -gt 0 ]; then type="Sulfur Oxidizing"; fi
	n=$(grep -c "$order" orders/nitrlist); if [ $n -gt 0 ]; then type="Nitrifying"; fi
	io=$(grep -c "$order" orders/ironoxlist); if [ $io -gt 0 ]; then type="Iron Oxidizing"; fi
	echo -e "$virus\t$order\t$type"
done < saltorders
