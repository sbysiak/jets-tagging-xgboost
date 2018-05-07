#!/bin/bash
DATASET=$1
RANDOMIZE=$2
RECOMPILE=$3

# Compile aliroot
if [ "$RECOMPILE" = 1 ]; then
  cd /opt/alice/ali-master/sw/BUILD/AliPhysics-latest/AliPhysics
  make -j 4 PWGHFjetsHF && make -j 4 install PWGHFjetsHF

  if [ $? -eq 0 ]; then
    cd -
  else
    cd -
    exit 10
  fi
fi


if [ ! -e "env.sh" ]
then
  echo "ERROR: env.sh not found."
  exit 3
fi
source env.sh

if [ "$DATASET" = 1 ]; 
  then
  export TEST_DIR='LHC17f8g'
  export DATASET_FILE=$TEST_DIR"_root_archive_AliAOD.txt"
elif [[ $DATASET == LHC* ]];
  then
  export TEST_DIR=$DATASET
  export DATASET_FILE=$TEST_DIR"_root_archive_AliAOD.txt"
else
  echo "Option missing: 1 - LHC17f8g"
  exit 2
fi

if [ ! -e "generate.C" ]
then
  echo "ERROR: generate.C not found."
  exit 3
fi

# Copy train files to seperate folder
rm -r ./train/*
rmdir ./train
mkdir ./train
cp ./env.sh ./train/
cp ./generate.C ./train/
cp ./handlers.C ./train/
cp ./globalvariables.C ./train/
cp ./MLTrainDefinition.cfg ./train/

if [ "$RANDOMIZE" = 1 ]
then
  python -c "from random import shuffle; f_in=open('$DATASET_FILE'); lines=f_in.readlines(); shuffle(lines); f_in.close(); f_out=open('train/'+'$DATASET_FILE','w'); f_out.writelines(lines); f_out.close()"
else
  cp ./$DATASET_FILE ./train
fi
cd ./train

aliroot -b -q generate.C\(\)

if [ ! -e "lego_train.sh" ]
then
  echo "ERROR: lego_train.sh not found."
  exit 3
fi
bash ./lego_train.sh

if [ ! -e "lego_train_validation.sh" ]
then
  echo "ERROR: lego_train_validation.sh not found."
  exit 3
fi
bash ./lego_train_validation.sh


cd -
