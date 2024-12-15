#!/bin/bash

# extract data
mkdir GSE250061 GSE230375 GSE216784
tar -xvf GSE250061_RAW.tar -C GSE250061
tar -xvf GSE230375_RAW.tar -C GSE230375
tar -xvf GSE216784_RAW.tar -C GSE216784

echo "Done!"
