#!/bin/sh
rm Doxyfile
sed "s/__PROJECT_NUMBER__/\"`git log --pretty=format:"%ai" -1`\"/g" Doxyfile.in > Doxyfile
doxygen
echo "... compressing"
cd doxygen
tar czf ../API_EG.tar.gz html
cd ..
echo "... copying to ftp site"
scp API_EG.tar.gz ftp.mech.kuleuven.be:.
echo "... executing server side script"
ssh ftp.mech.kuleuven.be 'bash -s' < extract_API_EG.sh


