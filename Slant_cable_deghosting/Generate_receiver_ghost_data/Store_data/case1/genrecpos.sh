#!/bin/bash

for i in $(seq 0 1 100) 
do
	x=$(($i*8+1500))
        z=$(($i*1+10))
        echo  "$x	-$z"
done

