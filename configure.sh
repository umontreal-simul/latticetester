#!/bin/sh

# This is a small script that can be used to configure the build when a shell
# scripting langage is available on the system.
# This is meant to allow the user to not write a really long command every time
# he wants to change it's configuration.

read -p "Installation prefix (default: '/usr/local'):" PREFIX
read -p "Build output directory (default: './build'):" OUT
read -p "NTL installation prefix (default: '/usr/local'):" NTL
read -p "BOOST installation prefix (default: empty)" BOOST
read -p "GMP installation prefix (default: empty)" GMP
read -p "Do you want to link the program staticaly (y/N):" STATIC
read -p "Do you want to build the documentation locally (y/N):" DOCS

if [ -z $PREFIX ]; then
  PREFIX=" "
else
  PREFIX="--prefix "$PREFIX
fi
if [ -z $OUT ]; then
  OUT=" "
else
  OUT="--out "$OUT
fi
if [ -z $NTL ]; then
  NTL="--ntl /usr/local/"
else
  NTL="--ntl "$NTL
fi
if [ -z $BOOST ]; then
  BOOST=" "
else
  BOOST="--boost "$BOOST
fi
if [ -z $GMP ]; then
  GMP=" "
else
  GMP="--gmp "$GMP
fi
if [ -z $STATIC ]; then
  STATIC=" "
elif [ $STATIC = "y" -o $STATIC = "Y" ]; then
  STATIC="--link-static "
fi
if [ -z $DOCS ]; then
  DOCS=" "
elif [ $DOCS = "y" -o $DOCS = "Y" ]; then
  DOCS="--build-docs "
fi

./waf configure$PREFIX$OUT$NTL$BOOST$GMP$STATIC$DOCS
