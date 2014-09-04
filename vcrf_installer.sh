#!/bin/bash
# VCRF installer

# check to see if we have the right XQuartz....

echo "== X.Ren et al [cvpr12] =="
echo "Downloading..."
wget http://homes.cs.washington.edu/~xren/research/cvpr2012_scene/scene_labeling_cvpr2012_v1.zip
mv scene_labeling_cvpr2012_v1.zip ./vendor 
echo "Extracting..."
unzip ./vendor/scene_labeling_cvpr2012_v1.zip
rm ./vendor/scene_labeling_cvpr2012_v1.zip
echo "Done."

