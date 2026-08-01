#!/bin/bash

awk 'condition {action}'
awk --> is a loop
in python --> 
for row in data:
    if layout == "Paired":
    print (row)
    if layout == "Paired"  and size <500:
    print (row)

in awk -->
awk '$2 == "Paired" {print $1}'
awk '$2 == "Paired" && $3 <500 {print $1}'
for finding columns -->
awk -F'\t' '{print $2}'  --> tab separated file
awk -F',' '{print $2}' --> comma separated file
NR == 1 --> first line
