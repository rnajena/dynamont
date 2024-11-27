import subprocess
import sys
from os.path import join
from src.python.segmentation.__init__ import __build__

def main():
    # Path to the C++ executable
    exe_path = join(__build__, 'dynamont_NT')

    # Execute the C++ executable with command-line arguments
    result = subprocess.run([str(exe_path)] + sys.argv[1:], capture_output=True, text=True)

    # Print the output of the C++ executable
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)

    # Exit with the return code of the C++ executable
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()