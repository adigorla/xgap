#!/bin/sh

usage(){
    echo -e "\n\tUSAGE: ./`basename $0` /path/to/xgap/root/dir /path/to/python3.X/interpreter\n"
    exit
}

if [ $# -ne 2 ]
then
    echo -e "\n\tERROR: Please refer to the proper usage shown below and provide 2 input args\n" 
    usage
fi

if [ ! -d "$1" ]
then
    echo -e "\n\tERROR: INVALID PATH to xgap's directory or directory dosen't exist; please provide valid path to xgap's ROOT directory\n"
    exit
elif [ ! -f "$2" ]
then
    echo -e "\n\tERROR: INVALID PATH to python3.X interpreter or file dosen't exist\n"
    exit
fi

XGAPPATH=${1}

PYPATH="#!${2}"

INITF=__init__.py

CHKPTF=checkpoint.py

SETUF=setup_utils.py 

if [ -f $XGAPPATH/bin/xgap ]
then
    hashbang=$(awk 'NR==1{print;exit}' $XGAPPATH/bin/xgap)
    echo -e "\n\tReplacing $hashbang with $PYPATH in file ${XGAPPATH}/bin/xgap\n"
    sed -i'' -e "1 s|$hashbang|$PYPATH|" $XGAPPATH/bin/xgap
else
    echo -e "\n\tERROR: the path to xgap's directory is NOT the path path to xgap's ROOT directory! Please provide path to xgap's ROOT directory\n"
fi

for path in $(find ${XGAPPATH}/xgap/ -type f -name '*.py'); do
    if [[ "$INITF" != "`basename ${path}`" ]] && [[ "$CHKPTF" != "`basename ${path}`" ]] && [[ "$SETUF" != "`basename ${path}`" ]] && [[ "$(find ${XGAPPATH}/xgap/scheduler/ -type f -name '*.py')" != *"${path}"* ]]
    then
	hashbang=$(awk 'NR==1{print;exit}' $path)
	echo -e "\tReplacing $hashbang with $PYPATH in file $path\n"
	sed -i'' -e "1 s|$hashbang|$PYPATH|" $path
    fi
done

echo -e "\n\tSuccessfully Completed!\n"
