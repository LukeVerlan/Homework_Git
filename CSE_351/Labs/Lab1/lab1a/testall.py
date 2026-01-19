#!/usr/bin/env python3
"""
All in one test script runner for Lab 1a.

By Samantha Dreussi <sdreussi@cs.washington.edu>
"""

import subprocess
import re

compilation = 0
test = 0
style = 0

compilation_run = subprocess.run(["make"], capture_output=True, text=True)
compilation_err = compilation_run.stderr
if compilation_err:
    print("----- COMPILATION RUN -----")
    print(compilation_err)
    compilation = 1
else:
    test_harness = subprocess.run(["./ptest"], capture_output=True, text=True)
    test_err = re.search("fail", test_harness.stdout)
    if test_err:
        print("----- TEST HARNESS -----")
        print(test_harness.stdout)
        test = 1


    dlc = subprocess.run(["./dlc.py"], capture_output=True, text=True)
    dlc_err = re.match(r"\[.+\]", dlc.stdout)
    if dlc_err:
        print("----- STYLE CHECK -----")
        print(dlc.stdout)
        style = 1



print("----- RESULTS -----")
if compilation: print("Compilation failed! See results above for compilation error details.")
else:
    if test: print("ptest failed! See results above for test details.")
    if style: print("dlc style check failed! See results above for test details.")

if (not (compilation or test or style)):
    print("All checks passed! Ready to submit to Gradescope")
