#!/bin/bash
set -e

cd "$(dirname "$0")"

echo "🔧 Preparing gmac_rs for crates.io publish..."

TEMP_FILES=(
  "gmac_rs/README.md"
  "gmac_rs/LICENSE"
)
TEMP_DIRS=(
  "gmac_rs/docs"
)

cp README.md gmac_rs/
cp LICENSE gmac_rs/
cp -r docs gmac_rs/

# Preview what will be published
cd gmac_rs
echo "📦 Previewing crate contents:"
cargo package --list

# Confirm publish
# echo
# read -p "✅ Publish to crates.io? (y/N): " confirm
# if [[ "$confirm" =~ ^[Yy]$ ]]; then
#   cargo publish --dry-run
# else
#   echo "❌ Aborted."
# fi

# Clean up copied files
echo "🧹 Cleaning up temporary files..."
cd ..
for f in "${TEMP_FILES[@]}"; do
  rm -f "$f"
done
for d in "${TEMP_DIRS[@]}"; do
  rm -rf "$d"
done

echo "✅ Done."
