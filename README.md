# FACS


## install pfunit 
currently not as part of the repo but globally installed
-get it from:
https://github.com/Goddard-Fortran-Ecosystem/pFUnit

tested version: 932130e787fcee0b8070515ede6912b35e0db6a9 (HEAD -> master, tag: v4.1.5)
-install commands:
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
cd pfUnit
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr ..
make
make tests
sudo make install
