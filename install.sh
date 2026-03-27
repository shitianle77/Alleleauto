#!/bin/bash
set -e

SUBMODULE_DIR="genetribe"

echo "=========================================="
echo "Installing Alleleauto dependencies"
echo "=========================================="

# If submodule is missing, try to initialize it automatically
if [ ! -f "$SUBMODULE_DIR/install.sh" ]; then
    echo "Submodule not found. Initializing..."
    git submodule sync --recursive
    git submodule update --init --recursive
    if [ ! -f "$SUBMODULE_DIR/install.sh" ]; then
        echo "ERROR: Failed to initialize submodule. Please check your network and try again."
        exit 1
    fi
    echo "✓ Submodule initialized."
else
    echo "✓ Submodule found."
fi

# Enter submodule and run its own installation script
cd "$SUBMODULE_DIR"
./install.sh
cd ../..

echo "✓ GeneTribe installed."

# Other optional operations...

echo "=========================================="
echo "Installation complete."
echo "=========================================="
