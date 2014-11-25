#!/bin/sh
rm Doxyfile
sed "s/__PROJECT_NUMBER__/\"`git log --pretty=format:"%ai" -1`\"/g" Doxyfile.in > Doxyfile
doxygen


