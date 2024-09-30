#!/bin/bash

# get orders for known genera
# Run script example: ./script.sh genera/sulfate_reducing > sulfate_reducing_orders
while read thing; do 
        taxonomy=$(grep "$thing" fullnamelineage.dmp | 
        awk 'BEGIN {FS="\t"}; {print $5}' | 
        awk 'NR==3');
        numgroup=$(echo $taxonomy |awk 'BEGIN {FS="; "}; {print $1 $2 $3 $4 $5}' | grep -o "group"| wc -l)
        if [[ "$numgroup" == 2 ]]; then
                order=$(echo $taxonomy | awk 'BEGIN {FS="; "}; {print $7}')
        elif [[ "$numgroup" == 1 ]]; then
                order=$(echo $taxonomy | awk 'BEGIN {FS="; "}; {print $6}');
        elif echo $taxonomy | grep -q ";"; then
                order=$(echo $taxonomy | awk 'BEGIN {FS="; "}; {print $5}');
        else
                order="Unknown Order"
        fi;
        order=$(echo $order | sed 's/;//')
        echo -e $order "\t" $thing;
done < $1

# next, had to manually search for unknown ones and add them to the previous output files

# get clean lists of the different types of known orders
awk 'BEGIN {FS="\t"}; (NR>1) {print $1}' ordergenera/sulfur_oxidizing_orders | sort -u > sulfoxlist
awk 'BEGIN {FS="\t"}; (NR>1) {print $1}' ordergenera/sulfate_reducing_orders | sort -u > sulfredlist
awk 'BEGIN {FS="\t"}; (NR>1) {print $1}' ordergenera/nitrifying_orders | sort -u > nitrlist
awk 'BEGIN {FS="\t"}; (NR>1) {print $1}' ordergenera/iron_oxidizing_orders | sort -u > ironoxlist


# get list of orders from hostpredictions
awk -v OFS='\t' 'BEGIN {FS=","}; (NR>1) {print $1,$6}' hostprediction.csv | awk '(NR>1) {print $0}' | awk '$2!=""' > saltorders

# search for matches to known orders
while read line; do
        order=$(echo $line | awk 'BEGIN {FS=" "}; {print $2}')
        virus=$(echo $line | awk 'BEGIN {FS=" "}; {print $1}')
        type=""
        sr=$(grep -c "$order" orders/sulfredlist); if [ $sr -gt 0 ]; then type="Sulfate Reducing"; fi
        so=$(grep -c "$order" orders/sulfoxlist); if [ $so -gt 0 ]; then type="Sulfur Oxidizing"; fi
        n=$(grep -c "$order" orders/nitrlist); if [ $n -gt 0 ]; then type="Nitrifying"; fi
        io=$(grep -c "$order" orders/ironoxlist); if [ $io -gt 0 ]; then type="Iron Oxidizing"; fi
        echo -e "$virus\t$order\t$type"
done < saltorders

# get counts of each host order predicted(includes those that did not match to known order)
awk -F '\t' '{print $2}' hosts | sort | uniq -c > hostordercounts

