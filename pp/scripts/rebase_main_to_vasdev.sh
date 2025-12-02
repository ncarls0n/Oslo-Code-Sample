#!/usr/bin/env bash
git push origin vasdev
git checkout main
git rebase vasdev
git push origin main
git checkout vasdev
