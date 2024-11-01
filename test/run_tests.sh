source ~/.bashrc
conda activate science

echo "Start CPP tests..."
echo ""

# run cpp tests
g++ -DUNIT_TESTING -I /home/yi98suv/miniconda3/envs/science/pkgs/gtest-1.15.2-h434a139_0/include/ -I ../src/dynamont/ test.cpp ../src/dynamont/utils.cpp -L /home/yi98suv/miniconda3/envs/science/pkgs/gtest-1.15.2-h434a139_0/lib/ -lgtest -lgtest_main -pthread -o test
export LD_LIBRARY_PATH=/home/yi98suv/miniconda3/envs/science/pkgs/gtest-1.15.2-h434a139_0/lib/:$LD_LIBRARY_PATH
./test

echo ""
echo "Start Python tests..."
echo ""

# run python tests
cd ..
python -m pytest
