#!/bin/bash

tar -cvzf $SCRATCH/alphafast/alphafast_hits_$(date +%Y%m%d%H%M).tar.gz \
    --absolute-names \
    $(find $SCRATCH/alphafast/*/output_hits/ -mindepth 1 -maxdepth 1 -type d)
