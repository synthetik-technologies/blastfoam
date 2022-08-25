#! /bin/bash

version="$1"
versionString="Version: $1" 
sed -i "s/Version:.*/$versionString/g" debian/blastfoam/DEBIAN/control
dch -v $version " "
# sed -i -E "s/(blastfoam_)(.*)(_.*)/\1$version\3/g" debian/files
