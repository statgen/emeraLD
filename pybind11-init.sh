#!/bin/bash

# Git strategy for merging in pybind11 for the first time. 
# You could change TARGET_BRANCH below to merge pybind11 into
# a branch other than master. 
TARGET_BRANCH="master"

git remote add -f -t master --no-tags pybind11 https://github.com/pybind/pybind11.git && \
  git checkout pybind11/master && \
  git subtree split -P include/pybind11 -b pybind-merge-branch && \
  git checkout ${TARGET_BRANCH} && \ 
  git subtree add --squash -P src/pybind11 pybind-merge-branch && \
  git branch -D pybind-merge-branch

