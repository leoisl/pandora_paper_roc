#!/usr/bin/env bash

wget https://www.dropbox.com/s/7804udh9c6th7v9/data_and_metadata.zip?dl=1 -O data_and_metadata.zip
cd .. && unzip scripts/data_and_metadata.zip
rm scripts/data_and_metadata.zip