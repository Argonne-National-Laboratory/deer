#!/bin/bash

NCPU="10"
CHECKPOINT_FILE_NAME="testfile.cpr"
EXEXCUTIONER_TYPE="dbg"

## Step 2-3: break and split the mesh this seems safe
$EXECUTIONER -i RVE_test_break.i --split-mesh 1 --split-file $CHECKPOINT_FILE_NAME

##Step 4: runnign the simulation for the first time
mpirun -np $NCPU $EXECUTIONER -i RVE_test_run0_split.i  Mesh/msh/file=$CHECKPOINT_FILE_NAME Outputs/nemesis=true Outputs/checkpoint=true

#Step 5: restarting the simulation using the msae input file, notice the problem paramter to relaod stateful var
RESTART_CHECKPOINT_FILE_LOCATION="RVE_test_run0_split_out_cp"
CHECKPOINT_ID="0002"
mpirun -np $NCPU $EXECUTIONER -i RVE_test_run0_split.i Mesh/msh/file=$RESTART_CHECKPOINT_FILE_LOCATION/${CHECKPOINT_ID}_mesh.cpr  Problem/restart_file_base=$RESTART_CHECKPOINT_FILE_LOCATION/$CHECKPOINT_ID Outputs/nemesis=true Outputs/checkpoint=true Outputs/file_base="restart"
