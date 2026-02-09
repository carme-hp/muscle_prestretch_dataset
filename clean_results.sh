#!/bin/bash

TARGET_DIR="results"

# Check if the directory exists
if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: '$TARGET_DIR' does not exist"
  exit 1
fi

# Delete all files but keep directories
find "$TARGET_DIR" -type f -delete

echo "All files deleted from '$TARGET_DIR', directories preserved."