#!/bin/bash

echo "w" >> diff.txt
ed - $1 < diff.txt
diff -y --suppress-common-lines $1 ../mcdif_test_new_rta/$1

