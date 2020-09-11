# CountmRNA.jl
[![DOI](https://zenodo.org/badge/255479900.svg)](https://zenodo.org/badge/latestdoi/255479900)

Our gene noise modulation project captured a series of 3D live cell fluorescence
time-lapse image using spin-disk confocal. We aim to count mRNA number inside
each nuclei during gene expression. Thanks to many other researches' work, we apply
and adjust algorithm to split each trajectory, extract nuclei, then recognize mRNA. 
This algorithm features on morphology processing.

There are our algorithm in julia code. You are free to take anything to help your study!
Even you don't care algorithm, `julia2ims.jl` and `tiffxml.jl` might save your
day if you work with imaris/imagej.

Red fluorescence prefixed with NLS are expressed to label nuclei, while mRNA are
more brilliant. We use Andor spin-disk microscope, 60x TIRF objective,
Hamamatsus Ocra Fusion 4(roi 1900x1300). Each position capture 20 z slices and
about 110 time points(~ 18 hours).

## Files

| NAME              | DESCRIPTION |
|-------------------|----------------------------------------------------------
| splitcell.jl      | Use Laplace of Gaussian(LoG) filter extract nuclei from raw 3d image |
| lineage.jl        | Use connected component to find each track and use watershed to split |
| segmentation3d.jl | Use Otsu's method to threshold 3d nuclei |
| normalization3d.jl| Normalize minimal and mean intensity of nuclei |
| julia2ims.jl      | Useful functions to load and save imaris 5 file |
| tiffxml.jl        | Useful functions to save tiff with OME-TIFF info |
| CountmRNA.jl         | Completed pipline to extract and track nuclei  |
| deprecations.jl   | Deprecated functons collections |
| ims2info.ipynb    | Extract Imaris data to matlab .mat files |
| notebook/xxx.ipynp| Debug files respond to each julia function, you can ignore them |


## Algorithm
Algorithm is special to our data feature and analysis target, which try to low
compute complexity but not for general case. In our strategy, we first project 4D
image(XYZT) into 3D image(XYT) to simplify problem, then search connected components
in projected image to locate each cell/nuclei trajectory and border. At the end,
benefiting from known trajectory and border, segmentation can be limited in
local. We apply common Otsu's method at each marked space in original data,
which separates reliable 3D nuclei out of cytoplasm background in spite of cell
intensity variation and photobleach.

Firstly, extract binary nuclei mask from raw image. We apply medium filter(size=5) and
LoG filter(gamma=40) on z-projected image. 

Secondly, we search connected component as cell trajectories in 3D image(XYT).
Because cells move slow(large overlap between frame) and occupy large area(under 60X objective), using
connected component do work correctly. More, comparing finding closest neighbour
by cell centroid, its direct and completed logic help to handle cell come/leave
from image border. Each detected connected components are assigned unique id. By
filtering connected component length in time dimension, short trajectories are
removed. Then possible collisions are detected by scanning connected component
number among each time slice for each trajectory, then they are split using a
distance transform and resign new unique id. These trajectories are the seeds to
divide each cell in whole image.

Thirdly, we use above trajectories as seeds/markers to perform watershed at each
z-projected image. As a result, each detected cell occupy unique area without
overlap in whole image, which called land mask. (Note: border mark don't equal
to nuclei edge)

Fourthly, backing to original 3D image(XYZ), we extend each 2D land mask to 3D
and apply Otsu's method to separate nuclei. Otsu's method requires enough
background and object information in image, this is why we keep both nuclei and
cytoplasm inside land mask. We use a 512x512x20 box to mask land mask again,
then Otsu's method is applied to original 3D image defined by final mask. Due to
cell intensity variation and photobleach, threshold is calculated for each nuclei 
and time point independently.

Fifthly, normalize nuclei intensity for imaris analysis. Because imaris is not
flexible to set variable parameter for search RNA spot, we normalize minimal
intensity and mean intensity of nuclei to fixed value.


## Count mRNA spot with imaris
Before we have these code, we tried to use Imaris for all these work, but Imaris allocated
huge memory for calculation and failed to handle photobleach and cell intensity
variaction. So we decide to simplify question: extract cell using own code and 
just let imaris count mRNA spot. More, my limited experience and time also
suggest me to use imaris at this time. Before load into imaris, I still need to
remove some failed segmentation result by hand.

Just use Imaris' spot model, the creation Parameters are:
```
[Algorithm]
  Enable Region Of Interest = false
  Enable Region Growing = false
  Enable Tracking = false
  Enable Region Growing = false
  Enable Shortest Distance = false
[Source Channel]
  Source Channel Index = 1
  Estimated XY Diameter = 0.540 um
  Estimated Z Diameter = 1.50 um
  Background Subtraction = false
[Classify Spots]
  # NOTE: The value are slight different among expriments
  "Intensity Center Ch=1 Img=1" above 2000 
  or
  "Intensity Mean Ch=1 Img=1" above 1850 or 1950 or 1800
```

## Citation
MSB-20-9991, Wei Huang, Quantitative Control of Noise in Mammalian Gene Expression by Dynamic Histone Regulations


## Known Issue and Todo
1. fail to separate some collisions: may require fully watershed instead of
just distance transform
2. some dark line occur in middle of nuclei: may some biology feature
3. fill hole inside neclei


## Licence
MIT License

Copyright (c) 2020 H.F.
