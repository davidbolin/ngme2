#!/bin/bash

set -e

# 1. Build the pkgdown site
echo "Building pkgdown site..."
Rscript -e 'pkgdown::build_site()'

# 2. Save current branch name
BRANCH=$(git rev-parse --abbrev-ref HEAD)

# 3. Fetch and switch to gh-pages branch
echo "Switching to gh-pages branch..."
git fetch origin
git checkout gh-pages

# 4. Copy site files from docs/ to root of gh-pages
echo "Copying site files..."
rsync -av --delete docs/ .  # safer than rm/cp, preserves .git and deletes removed files

# 5. Commit and push changes
echo "Committing and pushing..."
git add .
git commit -m "Update pkgdown site [skip ci]" || echo "No changes to commit"
git push origin gh-pages

# 6. Switch back to original branch
echo "Switching back to $BRANCH..."
git checkout "$BRANCH"

echo "Done! Your pkgdown site is deployed."