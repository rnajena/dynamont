# source ~/.bashrc
# conda activate science

echo "Start CPP tests..."
echo ""
$(dirname $0)/run_cpp_test.sh

echo ""
echo "Start Python tests..."
echo ""
$(dirname $0)/run_python_test.sh