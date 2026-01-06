#!/bin/bash

# Script to auto-pull, build, and run HPC assignment

set -e  # Exit on any error

# Configuration
REPO_DIR="."
SSH_KEY="$HOME/.ssh/hpc_deploy_key"
PROGRAM="./main"

# Check if ssh-agent is already running, if not start it
if [ -z "$SSH_AUTH_SOCK" ]; then
    eval $(ssh-agent -s)
    trap "ssh-agent -k" EXIT  # Kill agent when script exits
    ssh-add "$SSH_KEY"
elif ! ssh-add -l &>/dev/null; then
    # Agent is running but no keys loaded
    ssh-add "$SSH_KEY"
fi

# Navigate to repo
cd "$REPO_DIR"

# Pull latest changes
echo "Pulling latest changes..."
git pull origin main  # or whatever your branch is

# Build
echo "Building project..."
make clean
make all BLAS=atlas

# Run
echo "Running program..."
$PROGRAM

echo "Done!"