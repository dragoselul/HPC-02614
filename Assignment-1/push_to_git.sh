#!/bin/bash

SSH_KEY="$HOME/.ssh/hpc_deploy_key"
# Check if ssh-agent is already running, if not start it
if [ -z "$SSH_AUTH_SOCK" ]; then
    eval $(ssh-agent -s)
    trap "ssh-agent -k" EXIT  # Kill agent when script exits
    ssh-add "$SSH_KEY"
elif ! ssh-add -l &>/dev/null; then
    # Agent is running but no keys loaded
    ssh-add "$SSH_KEY"
fi

git push