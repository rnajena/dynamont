from os.path import join, dirname
from subprocess import call
from sys import exit, argv

def main():
    exe_path = join(dirname(__file__), "../bin/dynamont-NTC")
    exit(call([exe_path] + argv[1:]))

if __name__ == "__main__":
    main()