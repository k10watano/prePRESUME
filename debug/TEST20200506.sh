#!/bin/bash
# -*- cofing: utf-8 -*-
THIS_PATH=$(cd $(dirname $0); pwd)
PRESUME="${THIS_PATH}/../PRESUME.py"
OUTDIR=${THIS_PATH}/test
TREE="${THIS_PATH}/../example/example_1/PRESUMEout/PRESUMEout.nwk"
# run PRESUME
mkdir ${OUTDIR}
python3 ${PRESUME} --output ${OUTDIR} --indel  ${TREE} --ram_I 0.1 --ram_D 0.1 --prop_q 0.5 --int_r 1.0
cat ${TREE} | sed -e 's/[0-9]//g' > ${OUTDIR}/reference1
cat ${OUTDIR}/PRESUMEout_from_tree/test.nwk | sed -e 's/[0-9]//g' > ${OUTDIR}/reference2
diff -s ${OUTDIR}/reference1 ${OUTDIR}/reference2
echo "Done!"