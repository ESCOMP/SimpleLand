# Makefile for running tests on the python code here

# These variables can be overridden from the command-line
python = not-set
verbose = not-set
debug = not-set

ifneq ($(python), not-set)
    PYTHON=$(python)
else
    PYTHON=python3
endif

ifneq ($(debug), not-set)
    TEST_ARGS+=--debug
endif
ifneq ($(verbose), not-set)
    TEST_ARGS+=--verbose
endif

PYLINT=pylint
PYLINT_ARGS=-j 4 --rcfile=slim/.pylintrc
DIR_SRC = \
	slim \
	../cime_config/slim_cime_py

all: lint black test
test: utest stest

.PHONY: utest
utest: FORCE
	$(PYTHON) ./run_slim_py_tests $(TEST_ARGS) --unit

.PHONY: stest
stest: FORCE
	$(PYTHON) ./run_slim_py_tests $(TEST_ARGS) --sys

.PHONY: lint
lint: FORCE
	$(PYLINT) $(PYLINT_ARGS) $(DIR_SRC)

.PHONY: black
# Run black check on all of the python files here and undeneath.
# Should be identical to run_black option below, other than the --check option added in
# Use the black configure file to explicitly set a few things and specifiy the exact files.
black: FORCE
	black --check --config pyproject.toml $(DIR_SRC)

.PHONY: run_black
# Run black on all of the python files here and undeneath -- CHANGING THE FILES!
# Should be identical to above, other than the "--check" option is removed.
# Use the black configure file to explicitly set a few things and specifiy the exact files.
run_black: FORCE
	black --config pyproject.toml $(DIR_SRC)


.PHONY: clean
clean: FORCE
	find . -name '*.pyc' -exec rm {} \;

FORCE:

