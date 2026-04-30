#!/bin/bash

# Exit on any error
set -e

# Define variables
EXPNO="19"
TARGET_DIR="Experiments/$EXPNO"
SOURCE_FILE="namelistSIM"   # Change this to your source file name
EXECUTABLE="zoupa"
POST_FILE="./output_post_files/output_wp_$EXPNO"

# Step 0: Clean the previous build
make clean

# Step 1: Compile the model
make

# Step 2: Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Step 3: Copy the source file to the directory
cp "$SOURCE_FILE" "$TARGET_DIR/"

# Step 4: Run the model
nohup ./"$EXECUTABLE" < input_1dsim > "$POST_FILE" &

wait &

grep 'Number of viscous and plastic grid cells' output_post_files/output_mu_test_"$EXPNO" > ./simplotting1d/num_visc_plas_"$EXPNO".out