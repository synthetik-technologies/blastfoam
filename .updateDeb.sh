#! /bin/bash

version="Version: $1"
sed -i "s/Version:.*/$version/g" debian/blastfoam/DEBIAN/control