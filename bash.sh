#!/bin/bash

python scripts/gen.py -c X1 -t train 
python scripts/gen.py -c X1 -t val
python scripts/gen.py -c X2 -t train
python scripts/gen.py -c X2 -t val