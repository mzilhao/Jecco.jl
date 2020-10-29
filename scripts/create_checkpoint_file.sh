#!/usr/bin/env bash

# creates a checkpoint file by combining the files with the bulk, boundary and
# gauge data in a single hdf5 file

if [ "$#" = 0 ]; then
    echo "Usage: $(basename $0)  <iteration>"
    exit 0
fi

if ! [ -x "$(command -v h5copy)" ]; then
    echo "h5copy could not be found"
    exit 1
fi

it="$1"
padit="$(printf %08d $it)"

bdry_file="boundary_${padit}.h5"
gauge_file="gauge_${padit}.h5"
bulk_file="bulk_${padit}.h5"

# check if all needed files exist in this folder

if ! [ -f "${bdry_file}" ]; then
    echo "${bdry_file} does not exit"
    exit 1
fi

if ! [ -f "${gauge_file}" ]; then
    echo "${gauge_file} does not exit"
    exit 1
fi

if ! [ -f "${bulk_file}" ]; then
    echo "${bulk_file} does not exit"
    exit 1
fi

# copy data over to the (new) checkpoint file

checkpoint_file="checkpoint_it${padit}.h5"

cp -v ${bulk_file} ${checkpoint_file}

h5copy -v -i ${bdry_file} -o ${checkpoint_file} -s /data/${it}/fields/a4 -d /data/${it}/fields/a4
h5copy -v -i ${bdry_file} -o ${checkpoint_file} -s /data/${it}/fields/fx2 -d /data/${it}/fields/fx2
h5copy -v -i ${bdry_file} -o ${checkpoint_file} -s /data/${it}/fields/fy2 -d /data/${it}/fields/fy2
h5copy -v -i ${gauge_file} -o ${checkpoint_file} -s /data/${it}/fields/xi -d /data/${it}/fields/xi
