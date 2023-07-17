#! /usr/bin/env bash

set -eo pipefail

# wrapper script that mounts local directories into container and forwards arguments

# EXAMPLE: ./run-angler.sh \
#  --input_dir /home/my-user/AnglerVersion4/app/static/mrna_fasta \
#  --output_dir /home/my-user/angler-results-tmp \
#  --debug \
#  --log_dir /home/my-user/AnglerVersion4/app/log
#

args=()
input_dir=
output_dir=
log_dir=

while [ -n "$1" ]; do
    if [[ "$1" == '--input_dir' ]]; then
        shift
        input_dir=$1
        shift

    elif [[ "$1" == '--output_dir' ]]; then
        shift
        output_dir=$1
        shift

     elif [[ "$1" == '--log_dir' ]]; then
        shift
        log_dir=$1
        shift

    else
        args+=("$1")
        shift
    fi
done

docker run --rm \
    -v "${input_dir}":/tmp/input:ro \
    -v "${output_dir}":/tmp/output \
    -v "${log_dir}":"/angler4/log" \
    angler4:latest \
        --input_dir /tmp/input \
        --output_dir /tmp/output \
        --log_dir /angler4/log \
        "${args[@]}"
