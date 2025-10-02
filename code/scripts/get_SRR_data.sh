#!/bin/bash
echo $1
which fastq-dump
fastq-dump --split-3 --gzip $1
