#!/usr/bin/bash

if [ ! -d "course_data" ]; then
  if [ ! -f "course_data.tar.gz" ]; then
    wget https://single-cell-transcriptomics.s3.eu-central-1.amazonaws.com/course_data.tar.gz
  fi
  tar -xvf course_data.tar.gz
fi
