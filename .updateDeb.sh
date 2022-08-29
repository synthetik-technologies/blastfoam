#! /bin/bash

version="$1"
versionString="Version: $1" 
sed -i "s/Version:.*/$versionString/g" debian/blastfoam/DEBIAN/control
dch -v $version "Automated Release"
sed -i -E "s/(blastfoam_)(.*)(_.*)/\1$version\3/g" debian/files

echo "Debug"
cat debian/blastfoam/DEBIAN/control
cat debian/files
cat debian/changelog
