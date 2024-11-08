echo "NTK"
g++ -std=c++20 -O3 -o dynamont_NTK dynamont_NTK.cpp utils.cpp
echo "NT_banded"
g++ -std=c++20 -O3 -o dynamont_NT dynamont_NT.cpp utils.cpp
echo "NT"
g++ -std=c++20 -O3 -o dynamont_NT_banded dynamont_NT_banded.cpp utils.cpp
