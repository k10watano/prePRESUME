#!/bin/bash
# -*- cofing: utf-8 -*-
THIS_PATH=$(cd $(dirname $0); pwd)
PRESUME="${THIS_PATH}/../PRESUME.py"
OUTDIR=${THIS_PATH}/test
TREE="${THIS_PATH}/../example/example_1/PRESUMEout/PRESUMEout.nwk"
# run PRESUME
mkdir ${OUTDIR}
python3 ${PRESUME} --output ${OUTDIR} --cassiopeir -n 100
echo "Done!"