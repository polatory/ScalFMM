#!/bin/bash
cd /home/mkhannou/scalfmm/Doc/noDist/implicit
emacs implicit.org --batch -f org-html-export-to-html --kill
ssh scm.gforge.inria.fr "cd /home/groups/scalfmm/htdocs/orgmode/; rm -rf implicit"
cd ..
rsync -e ssh -avz --delete-after implicit scm.gforge.inria.fr:/home/groups/scalfmm/htdocs/orgmode/
ssh scm.gforge.inria.fr "cd /home/groups/scalfmm/htdocs/orgmode/; chmod og+r implicit -R;"

#Put an up-to-date tarbal
cd /home/mkhannou/scalfmm
./export_tarbal.sh &
