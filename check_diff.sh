!#/bin/bash

for file in *; do
    diff -y --suppress-common-lines ${file} ../mcdif_test_new_rta/${file} > ${file}.txt
    echo $file
    read czekaj
    mcedit ${file}.txt

done