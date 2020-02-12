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
    echo -e "\n\tERROR: INVALID PATH to xgap's directory or directory dosen't exist; please provide valid path to winston's ROOT directory\n"
    exit
elif [ ! -f "$2" ]
then
    echo -e "\n\tERROR: INVALID PATH to python3.X interpreter or file dosen't exist\n"
    exit
fi

WINSPATH=${1}

PYPATH="#!${2}"

INITF=__init__.py

CHKPTF=checkpoint.py

SETUF=setup_utils.py 

if [ -f $WINSPATH/bin/xgap ]
then
    hashbang=$(awk 'NR==1{print;exit}' $WINSPATH/bin/xgap)
    echo -e "\nReplacing $hashbang with $PYPATH in file ${WINSPATH}/bin/xgap"
    sed -i'' -e "1 s|$hashbang|$PYPATH|" $WINSPATH/bin/xgap
else
    echo -e "\n\tERROR: the path to xgap's directory is NOT the path path to winston's ROOT directory! Please provide path to winston's ROOT directory\n"
fi

for path in $(find ${WINSPATH}/xgap/ -type f -name '*.py'); do
    if [[ "$INITF" != "`basename ${path}`" ]] && [[ "$CHKPTF" != "`basename ${path}`" ]] && [[ "$SETUF" != "`basename ${path}`" ]] && [[ "$(find ${WINSPATH}/xgap/scheduler/ -type f -name '*.py')" != *"${path}"* ]]
    then
	hashbang=$(awk 'NR==1{print;exit}' $path)
	echo -e "\nReplacing $hashbang with $PYPATH in file $path"
	sed -i'' -e "1 s|$hashbang|$PYPATH|" $path
    fi
done
