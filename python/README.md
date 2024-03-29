# Testing the code here

## Running everything

To run all tests (unit tests, system tests and pylint), simply run `make
all` from this directory.

## Python environment

use the "npl" conda environment on cheyenne

## Unit and system tests

Unit and system tests can be run in one of two ways; these do the same
thing, but support different options:

1. via `make test`

   You can specify a few arguments to this:
   
   - python version: `make python=python3.9 test` (defaults to `python3`; you should expect errors if trying to run with python2)
   - verbose: `make verbose=true test`
   - debug: `make debug=true test`

   Note that unit tests and system tests can be run separately with
   `make utest` or `make stest`, or they can all be run with `make
   test`.

2. via `./run_slim_py_tests`

   You can specify various arguments to this; run `./run_slim_py_tests
   -h` for details

## pylint

You can run pylint on everything in the slim package with `make lint`.

## black

You can run black on everything in the slim package with `make black`.

