#!/bin/bash

sed 's, ,_,g' -i spikeprot*.fasta

mafft --thread 8 --op 10 spikeprot*.fasta > spike_GISAID_aligned.fasta