# Generep
Genome-wide network analysis

## Introduction
Generep is a complete pipeline for gene network reconstruction from 
gene expression data. It takes as input a matrix of gene expression
values (genes on rows, samples on columns) and outputs a network of
gene-gene interactions in different formats (see Output below).

## Prerequisites
* Nextflow. Generep is implemented as a [NextFlow](https://www.nextflow.io/) script, it therefore requires
a recent version of nextflow to run. 

* Aracne2. Generep relies on the Aracne2 program to compute Mutual Information (MI)
values for gene-gene interactions, and it therefore expects the aracne2 executable
to be found in $PATH. Aracne2 can be installed from http://califano.c2b2.columbia.edu/aracne.

* Python 2.7, with the numpy and scipy libraries. Support for Python3 is coming soon.

## Configuration file

Generep requires a configuration file describing the details of the desired
analysis. An example configuration file is provided as `sample.conf`. Here is
a summary of what the configuration file should contain:

```
[General]
# Prefix for output files (will be placed in results/)
title = final

# Gene expressions table (genes on rows, patients on columns)
expressions = data/gene-expressions.txt

# File containing IDs of transcription factors, one per line
tfs = data/hs-tf.txt

# File mapping IDs to gene names, for final output
translation = data/ids.csv

# Number of columns in expressions table before the start of data.
# First column is assumed to contain gene identifiers. If there are
# additional non-data columns, set this to a number greater than 1.
genecols = 1

# Number of shuffled datasets
nshuffled = 3

[Bootstrap]
# Number of bootstrap rounds
nrounds = 100

[Aracne]
# DPI threshold (see Aracne documentation)
dpi = 0.1

[CX]
# Attributes for generation of CX files.
(see sample.conf for the full list)
```
