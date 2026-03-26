#!/bin/bash
set -e

# Submodule path
SUBMODULE_DIR="genetribe"

echo "=========================================="
echo "Installing Alleleauto dependencies"
echo "=========================================="

# Check if submodule exists (using install.sh as indicator)
if [ ! -f "$SUBMODULE_DIR/install.sh" ]; then
    echo "ERROR: Submodule $SUBMODULE_DIR is missing."
    echo "Please clone with --recurse-submodules, or run:"
    echo "  git submodule update --init --recursive"
    exit 1
fi

echo "✓ Submodule found."

# Enter submodule and run its own install script
cd "$SUBMODULE_DIR"
./install.sh
cd ../..

echo "✓ GeneTribe installed."

# Optional: link genetribe executable to main repo's bin directory for convenience
# mkdir -p bin
# ln -sf "$(pwd)/$SUBMODULE_DIR/genetribe" bin/genetribe
# echo "GeneTribe linked to bin/genetribe."

echo "=========================================="
echo "Installation complete."
echo "=========================================="
