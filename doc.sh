#!/bin/sh

# This is a small shell script that automates the documentation update process.

./waf configure --ntl /usr/local --build-docs
./waf build
rm -rf doc/html
git clone https://github.com/umontreal-simul/latticetester.git doc/html
cd doc/html
git checkout -b gh-pages origin/gh-pages
cd ../..
cp -r build/doc/html/* doc/html
cd doc/html
read
git add .
git commit -m 'Update to documentation'
git push
