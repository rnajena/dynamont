# run cpp tests

# TODO: change these paths to build proper package
gtest=/home/yi98suv/miniconda3/envs/science/pkgs/gtest-1.15.2-h434a139_0/include/
library=/home/yi98suv/miniconda3/envs/science/pkgs/gtest-1.15.2-h434a139_0/lib/

dynamont=$(dirname $0)/../src/dynamont/
testscript=$(dirname $0)/test.cpp
testexe=$(dirname $0)/test
utilsscript=$(dirname $0)/../src/dynamont/utils.cpp

g++ -std=c++20 -DUNIT_TESTING -I $gtest -I $dynamont $testscript $utilsscript -L $library -lgtest -lgtest_main -pthread -o $testexe
export LD_LIBRARY_PATH=$library:$LD_LIBRARY_PATH
$testexe
