#!/bin/bash
if [[ $# -ne 1 || ! -d $1 ]]; then 
    echo USAGE: $0 /path/to/datastore/folder
    exit 1
fi

path=$1
pushd $path
#get the "key" component
key=`basename $path | sed 's/.*\.//'`
md5_name=CHECKSUM.$key.md5

if [[ -e $md5_name ]]; then rm $md5_name; fi
#FIXME: this will work on LIS using BSD. if anyone else ever uses it, will need to test OS and adjust accordingly I think
md5 -r * > $md5_name
popd
