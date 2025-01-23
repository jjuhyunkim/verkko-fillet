#!/bin/bash

# Ensure the directory argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

# Set the directory to work with
directory=$(realpath "$1")

# Get the current user's username
user=$(whoami)

# Process each file in the directory
find "$directory" -type f | while read -r file; do
    # Check if the file exists and is accessible
    if [ ! -e "$file" ]; then
        echo "File $file does not exist or is not accessible"
        continue
    fi

    # Get the file owner
    owner=$(stat -c "%U" "$file")

    # If the file is owned by the current user, remove write permissions
    if [ "$owner" == "$user" ]; then
        if chmod a-w "$file"; then
            echo "Removed write permission from $file"
        else
            echo "Failed to remove write permission from $file"
        fi
    else
        echo "Skipped $file (not owned by $user)"
    fi
done

