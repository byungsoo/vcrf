#!/bin/bash
# VCRF installer

cd vendor

echo "== X.Ren et al [cvpr12] =="
echo "Downloading..."
wget http://homes.cs.washington.edu/~xren/research/cvpr2012_scene/scene_labeling_cvpr2012_v1.zip
echo "Extracting..."
unzip scene_labeling_cvpr2012_v1.zip
rm scene_labeling_cvpr2012_v1.zip
echo "Done."

cd ..
