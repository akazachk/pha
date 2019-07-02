#!/bin/bash
# Takes as an (optional) argument the directory where you wish to install the COIN-OR software

if [ "$SHELL" != "/bin/bash" ]
  then echo "This script assumes a bash shell. Comment out the exit below and proceed on your own. (For example, with sh, the >& redirection will not work below.)"
  exit
fi

## User needs to define
#PHA_DIR="~/repos/pha"
COIN_STUB="Cgl"
COIN_VERSION="0.59"
#COIN_REVISION=
COIN_URL="https://projects.coin-or.org/svn/${COIN_STUB}/stable/${COIN_VERSION}"
COIN_DIR_NAME="${COIN_STUB}-${COIN_VERSION}"

## Check arguments
if [ -z "$PHA_DIR" ]
  then echo "Need to define PHA_DIR. Exiting."
  exit
fi
if [ -z "$1" ]
then
	COIN_DIR="${PHA_DIR}/lib/${COIN_DIR_NAME}"
else
	COIN_DIR="${1}/${COIN_DIR_NAME}"
fi

## Some extra arguments (can be commented out)
EXTRA_ARGS=
if [ ! -z ${ENV_LAPACK_LIB} ]
then
  EXTRA_ARGS+=" with-lapack=-L${ENV_LAPACK_LIB} -llapack"
fi
if [ ! -z ${ENV_BLAS_LIB} ]
then
  EXTRA_ARGS+=" with-blas=-L${ENV_BLAS_LIB} -lblas"
fi

## Ignore below unless you wish to use OsiCpxSolverInterface
#UNAME=`uname`
#if [ "$UNAME" = "Darwin" ]
#then
#	CPLEX_ARCH="x86-64_osx"
#	CPLEX_DIR="/Applications/CPLEX_Studio128"
#else
#	CPLEX_ARCH=x86-64_linux
#	CPLEX_DIR=/home/ibm/cplex-studio/12.8.0.0
#fi
#CPLEX_INC="$CPLEX_DIR/cplex/include/ilcplex"
#CPLEX_LIB_DIR="$CPLEX_DIR/cplex/lib/$CPLEX_ARCH/static_pic"
#CPLEX_LIB="-L$CPLEX_LIB_DIR -lcplex -lm -ldl -lpthread"
#EXTRA_ARGS+=" --with-cplex-incdir=$CPLEX_INC --with-cplex-lib=$CPLEX_LIB"

echo "    Installing ${COIN_STUB} into ${COIN_DIR}"
mkdir -p $COIN_DIR
if [ -z ${COIN_REVISION} ]
then
  svn co $COIN_URL $COIN_DIR
else
  svn co -r ${COIN_REVISION} $COIN_URL $COIN_DIR
fi
cd $COIN_DIR

echo "    Configuring ${COIN_DIR}/build"
mkdir -p build
cd build
../configure -C ${EXTRA_ARGS} >& last_config.txt
make
make install
cd ..

## Uncomment below if you wish to compile in debug mode
echo "    Configuring ${COIN_DIR}/buildg (debug mode)"
ADD_CXXFLAGS='-fkeep-inline-functions'
EXTRA_ARGS+=" --enable-debug=yes"
mkdir -p buildg
cd buildg
../configure -C ${EXTRA_ARGS} >& last_config.txt
make
make install

echo ""
echo "Done with installing COIN-OR files into $COIN_DIR"
