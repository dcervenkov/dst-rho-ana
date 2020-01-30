#!/bin/bash

DIRS="../DSRhoBackground ../DSRhoCPFit ../DSRhoEfficiency ../DSRhoLifetime ../DSRhoPeek ../DSRhoYield"

BOLD=$(tput bold)
NORMAL=$(tput sgr0)
RETURN_CODE=0

for DIR in ${DIRS}; do
    cd ${DIR}
    echo "${BOLD}$(basename ${DIR})${NORMAL}"
    echo "Building..."
    make -j4 &> /dev/null

    if [ $? -eq 0 ]; then
        if [[ -f tests/run_tests.py ]]; then
            cd tests
            echo "Running tests..."
            ./run_tests.py 2>/dev/null | grep --after-context=100 -e "^Summary$"

            if [ $? -ne 0 ]; then
                echo "ERROR: $(basename ${DIR}) tests failed!"
                RETURN_CODE=1
            fi

            cd ..
        fi
    else
        echo "ERROR while building $(basename ${DIR}), skipping..."
        RETURN_CODE=1
    fi
    echo
done

exit $RETURN_CODE