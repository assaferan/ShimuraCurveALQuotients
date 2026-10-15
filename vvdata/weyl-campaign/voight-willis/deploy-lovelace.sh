#!/bin/sh
# Copy this checkout to lovelace under a directory of its own and rewrite the absolute paths the
# scripts carry.  A directory of its own because the shared checkout there has other sessions'
# jobs running from it, and updating a tree under a running job mixes two versions of the code.
# The paths have to be rewritten because Magma's `load` takes a string literal, so the library
# paths cannot be passed in as an argument.
#
#   sh deploy-lovelace.sh [remote-directory]
set -e
SRC=/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign
DST=${1:-/home/assaferan/vwrun/campaign}
HOST=lovelace

ssh $HOST "mkdir -p $DST"
rsync -az --delete \
    --exclude '.git' --exclude 'worktrees' --exclude 'polymake' \
    --exclude 'vvdata/weyl-campaign/voight-willis/out' \
    "$SRC/" "$HOST:$DST/"
ssh $HOST "cd $DST && mkdir -p vvdata/weyl-campaign/voight-willis/out && \
    grep -rl '$SRC' vvdata/weyl-campaign/voight-willis | xargs sed -i 's#$SRC#$DST#g' && \
    echo rewritten: \$(grep -rl '$DST' vvdata/weyl-campaign/voight-willis | wc -l) files"
echo "deployed to $HOST:$DST"
