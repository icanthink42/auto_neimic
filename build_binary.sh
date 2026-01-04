#!/usr/bin/env bash
set -euo pipefail

# Build a standalone Windows GUI binary with PyInstaller.
# Note: PyInstaller does not cross-compile; run this on Windows (or Wine) with Python and PyInstaller installed.
if ! command -v pyinstaller >/dev/null 2>&1; then
  echo "PyInstaller not found. Install with: pip install pyinstaller"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

pyinstaller --noconfirm --onefile --windowed --name auto_niemiec main.py

echo "Windows .exe created at dist/auto_niemiec.exe"

