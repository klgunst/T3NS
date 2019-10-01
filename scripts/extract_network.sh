#!/bin/bash

file=$1
output=$(awk '/Site to orbital: /,/#/' $file | sed '/#/,$d')
sites=$(echo "$output" | sed -n 's~Site to orbital: ~~p')
nsites=$(echo "$sites" | wc -w)
bsites=$(echo "$sites" | tr -dc '*' | wc -c)
psites=$(($nsites - $bsites))
bonds=$(echo "$output" | tr -dc '>' | wc -c)
echo "NR_SITES = $nsites"
echo "NR_PHYS_SITES = $psites"
echo "NR_BONDS = $bonds"

echo "$output" | sed 's~Site to orbital: ~/\n~;s/*$/* /;s/-> //;s~Bonds :~/~'
