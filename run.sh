#!/bin/bash

# Figure out where generep.nf is, assuming it's in the same directory as
# this script. Otherwise, set the GR_HOME variable manually.
prog=$(readlink -f $0)
GR_HOME=$(dirname $prog)
export PYTHONPATH=$PYTHONPATH:${GR_HOME}/bin

CONF=$1

if [ -f "$CONF" ];
then
  nextflow -c $CONF run $GR_HOME/generep.nf -resume
else
  echo Usage: run.sh configfile
fi

