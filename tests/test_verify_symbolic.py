import sys
import os
import runpy

def test_verify_symbolic():
    # Execute verify_symbolic.py as a script
    script_path = os.path.join(os.path.dirname(__file__), "..", "src", "verify_symbolic.py")
    try:
        runpy.run_path(script_path, run_name="__main__")
    except SystemExit as e:
        # verify_symbolic.py calls sys.exit(1) on failure, sys.exit(0) on success
        assert e.code == 0 or e.code is None
