#!/bin/bash

# Test with micromamba
eval "$(micromamba shell hook --shell bash)"

# The channel_priority needs to be set to flexible

TRAVIS_PYTHON_VERSION=3.12
for dir in $(find . -name "__pycache__"); do
    rm -r $dir
done
if [ -e build ]; then
    rm -r build
fi
micromamba create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda python=$TRAVIS_PYTHON_VERSION
micromamba activate pgt_test_${TRAVIS_PYTHON_VERSION}
requirement_file=requirements_CI.txt
micromamba install --yes -c bioconda --file ${requirement_file}
python -m pip install .
coverage run -m py.test
coverage html
coverage-badge -f -o docs/coverage.svg
micromamba deactivate


TRAVIS_PYTHON_VERSION=3.14
for dir in $(find . -name "__pycache__"); do
    rm -r $dir
done
if [ -e build ]; then
    rm -r build
fi
micromamba create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda python=$TRAVIS_PYTHON_VERSION
micromamba activate pgt_test_${TRAVIS_PYTHON_VERSION}
micromamba install --yes -c bioconda bedtools
requirement_file=requirements_CI.txt
python -m pip install -r ${requirement_file}
python -m pip install .
py.test pygenometracks --doctest-modules
micromamba deactivate
