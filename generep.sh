#!/bin/bash

# Figure out where generep.nf is, assuming it's in the same directory as
# this script. Otherwise, set the GR_HOME variable manually.
prog=$(readlink -f $0)
GR_HOME=$(dirname $prog)
export PYTHONPATH=$PYTHONPATH:${GR_HOME}/bin

CONF=$1

if [ ! -f "$CONF" ];
then
  echo Usage: run.sh configfile
  exit 1
fi

NF_CONF=$($GR_HOME/bin/make_conf.py $CONF)

nextflow -c $NF_CONF run $GR_HOME/generep.nf -resume

