#!/bin/bash

mpb plots

mv dashboard $1
mv $1/report.html $1/index.html

scp -r $1 wrtabb@t3.unl.edu:/home/hep/wrtabb/public_html/DrellYan/$1
rm -r $1
