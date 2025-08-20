#!/bin/bash

# Exit on any error
set -e

# Define variables
EXPNO="12"
TARGET_DIR="Experiments/$EXPNO"
SOURCE_FILE="./src/ice.f90"   # Change this to your source file name
EXECUTABLE="zoupa"
POST_FILE="./output_post_files/output_vp_test_$EXPNO"

# Step 1: Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Step 2: Copy the source file to the directory
cp "$SOURCE_FILE" "$TARGET_DIR/"


# Step 4: Compile the model (Fortran example using gfortran)
make

# Step 5: Run the model
nohup ./"$EXECUTABLE" > "$POST_FILE" &