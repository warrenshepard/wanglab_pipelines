#!/bin/bash
# This converts mango files to WashU format
# This also unzips .txt files for WashU viewing

# RUN FROM Hichipper FOLDER

mkdir -p WashU

for f in *"all.mango"; do
awk '{if($8 < 0.01) print $0}' ${f} > ${f}_FDR_0.01.mango
awk '{if($8 < 0.05) print $0}' ${f} > ${f}_FDR_0.05.mango
done

for f in *".mango"; do
awk 'NR>1{print $1,$2,$3,$4":"$5"-"$6","$7,$8}' OFS="\t" ${f} > WashU/${f}.txt
done

for f in *".interaction.txt.gz"; do
cp $f WashU/$f
cd WashU
gunzip $f
cd ..
done