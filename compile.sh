echo "NTK"
g++ -I include/ -O3 -Wall -o src/dynamont/dynamont_NTK src/dynamont/dynamont_NTK.cpp src/dynamont/utils.cpp
echo "NT"
g++ -I include/ -O3 -Wall -o src/dynamont/dynamont_NT src/dynamont/dynamont_NT.cpp src/dynamont/utils.cpp
