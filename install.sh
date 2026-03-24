#!/bin/bash
set -e

# 子模块路径
SUBMODULE_DIR="genetribe"

echo "=========================================="
echo "Installing Alleleauto dependencies"
echo "=========================================="

# 检查子模块是否存在（通过 install.sh 判断）
if [ ! -f "$SUBMODULE_DIR/install.sh" ]; then
    echo "ERROR: Submodule $SUBMODULE_DIR is missing."
    echo "Please clone with --recurse-submodules, or run:"
    echo "  git submodule update --init --recursive"
    exit 1
fi

echo "✓ Submodule found."

# 进入子模块，执行其自身的安装脚本
cd "$SUBMODULE_DIR"
./install.sh
cd ../..

echo "✓ GeneTribe installed."

# 可选：将 genetribe 可执行文件链接到主仓库的 bin 目录，方便使用
# mkdir -p bin
# ln -sf "$(pwd)/$SUBMODULE_DIR/genetribe" bin/genetribe
# echo "GeneTribe linked to bin/genetribe."

echo "=========================================="
echo "Installation complete."
echo "=========================================="
