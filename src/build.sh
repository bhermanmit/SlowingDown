#!/bin/bash

cmake -H. -Bbuild -Doptimize=on
make -s -C build
