#!/bin/sh
FileName=$1
newName=${FileName::-3}.pdf
pandoc -s $FileName -o $newName;
