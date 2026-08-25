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

# TRAVIS_PYTHON_VERSION=3.11
# conda create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda -c conda-forge python=$TRAVIS_PYTHON_VERSION
# conda activate pgt_test_${TRAVIS_PYTHON_VERSION}
# conda install --yes -c conda-forge -c bioconda bedtools
# pip install -r requirements_CI.txt
# python setup.py install
# py.test pygenometracks --doctest-modules
# conda deactivate
