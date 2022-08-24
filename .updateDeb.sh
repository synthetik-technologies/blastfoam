#! /bin/bash

sed -i "s/Version:.*/Version: $1/ debian/blastfoam/DEBIAN/control"