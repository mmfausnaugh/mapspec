#!/usr/bin/bash

#run the test in this directory
if [ -f mapspec.params ]; then
   rm mapspec.params
fi
do_map -r ref.smooth.txt -w oiii.window   -s speclist_use -o mapspec.params

for i in $(ls scale*); do
    diff "$i" output_checks/"$i" > /dev/null
    [[ 0 -eq $? ]] && echo "File $i passes test." || \
    echo "Error: file $i differs from what is in output_checks"
done

diff mapspec.params output_checks/mapspec.params > /dev/null
[[ 0 -eq $? ]] && echo "File mapspec.params  passes test." || echo "Error: file mapspec.params differs from what is in output_checks"

diff chains/ output_checks/chains/ > /dev/null
[[ 0 -eq $? ]] && echo "Directory 'chains' passes test."|| echo "Error: directory chains differs from what is in output_checks"
