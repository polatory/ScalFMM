#!/bin/bash
mkdir platypus 
cd platypus
git clone git+ssh://mkhannou@scm.gforge.inria.fr/gitroot/scalfmm/scalfmm.git >/dev/null 2>&1
cd scalfmm
git checkout mpi_implicit >/dev/null 2>&1
cd ..
rm -rf scalfmm/.git > /dev/null
tar czf scalfmm.tar.gz scalfmm
scp scalfmm.tar.gz scm.gforge.inria.fr:/home/groups/scalfmm/htdocs/orgmode/implicit >/dev/null
ssh scm.gforge.inria.fr "cd /home/groups/scalfmm/htdocs/orgmode/; chmod og+r implicit -R;"
cd ..
rm -rf platypus >/dev/null
