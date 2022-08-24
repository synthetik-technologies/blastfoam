#! /bin/bash

sed -i "s/Version:.*/Version: $1/g" debian/blastfoam/DEBIAN/control