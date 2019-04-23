#!/bin/sh
# compute dependencies for the directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

if test $# = 0
then
    dirs=" src" 
          
else
    dirs=$*
fi

for DIR_ in $dirs
do
    DIR=`echo $DIR_ | sed 's?/??' `
    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    DEPENDS=""
    case $DIR in 
        src )
                  DEPENDS="$DEPENDS "            ;;
    esac

    # generate dependencies file
    if test -d $TOPDIR/../$DIR
    then
        cd $TOPDIR/../$DIR
       
        $TOPDIR/moduledep.sh $DEPENDS > make.depend
        $TOPDIR/includedep.sh $DEPENDS >> make.depend
    fi

    rm -f make.depend.tmp

    # check for missing dependencies
    if grep @ make.depend
    then
        notfound=1
        echo WARNING: dependencies not found in directory $DIR
    else
        echo directory $DIR : ok
    fi
done

if test "$notfound" = ""
then
    echo all dependencies updated successfully
else
    echo not all dependencies were found
fi
