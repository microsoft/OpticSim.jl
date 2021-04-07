import os
import sys
import multiprocessing

# set up the project (should already be installed in base through docker so should be quick)
ret = os.system(f"xvfb-run julia --project -e \"using Pkg; Pkg.resolve()\"")
if os.WEXITSTATUS(ret) != 0:
    sys.exit(os.WEXITSTATUS(ret))

# run the script
ret = os.system(f"JULIA_NUM_THREADS={multiprocessing.cpu_count()} xvfb-run julia --project {' '.join(sys.argv[1:])}")
sys.exit(os.WEXITSTATUS(ret))
