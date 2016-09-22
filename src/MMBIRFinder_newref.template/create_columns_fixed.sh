#!/bin/sh
ll=`grep -n '@PG' trial2/unaligned_2.sam | awk 'BEGIN { FS = ":"} ; { print $1 }'`
less trial2/unaligned_2.sam | sed "1,${ll}d" | awk '{print $1 "\t" $10}' > columns.txt

