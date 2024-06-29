#!/bin/bash

END=80
for((i=40;i<=END;i++))
do
cd analysis"$i"/multisite
python ../../matlab.py
cd ../../
done
