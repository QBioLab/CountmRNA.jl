# CountmRNA.jl
Our gene modulation project captured a series of 3D live cell fluorescence
time-lapse image using spin-disk confocal. We aim to count mRNA number inside
each nuclei during gene expression. Thanks to many other researches' work, we apply
and modify algorithm to split each trajectory, extract nuclei, then recognize mRNA. 
This algorithm mainly uses morphology processing.

Here is julia code for our algorithm. You can take anything to help your study!
Even you don't care algorithm, `julia2ims.jl` and `tiffxml.jl` might save your
day if you work with imaris/imagej.


Red fluorescence prefixed with NLS are expressed to label nuclei, while mRNA are
more brilliant. We use Andor spin-disk microscope, 60x TIRF objective,
Hamamatsus Ocra Fusion 4(roi 1900x1300). Each position capture 20 z slices and
about 130 time points(~ 20 hours).


## Files

| NAME              | DESCRIPTION |
|-------------------|----------------------------------------------------------
| splitcell.jl      | Use Laplace of Gaussian(LoG) filter extract nuclei from raw 3d image |
| lineage.jl        | Use connected component to find each track and use watershed to split |
| segmentation3d.jl | Use Otsu's method to threshold 3d nuclei |
| normalization3d.jl| Normalize minimal and mean intensity of nuclei |
| julia2ims.jl      | Useful functions to load and save imaris 5 file |
| tiffxml.jl        | Useful functions to save tiff with OME-TIFF info |
| xxx.ipynp         | Debug file respond to each julia function, you can ignore them |

## Algorithm
Algorithm is special to our data feature and analysis target, which try to low
compute complexity but not for general case. In our strategy, we project 4D
image(XYZT) into 3D image(XYT) to simplify problem, then search connected components
in projected image to locate each cell/nuclei trajectory and border. At the end,
benefiting from known trajectory and border, segmentation can be limited in
local. We apply common Otsu's method at each marked space in original data,
which separates reliable 3D nuclei out of cytoplasm background in spite of cell
intensity variation and photobleach.

Firstly, extract binary nuclei mask from raw image. We apply medium filter(5) and
LoG filter(40) on z projected image. 

Secondly, we search connected component as cell trajectories in 3D image(XYT).
Because cells move slow and occupy large area(under 60X objective), using
connected component do work correctly. More, comparing finding closest neighbour
by cell centroid, its direct and completed logic help to handle cell come/leave
from image border. Each detected connected components are assigned unique id. By
filtering connected component length in time dimension, short trajectories are
remove. Then possible collisions are detected by scanning connected component
number among each time slice for each trajectory, then they are split using a
distance transform and resign new unique id. These trajectories are seeds to
divide each cell in whole image.

Thirdly, we use above trajectories as seeds/markers to performance watershed at
each z-projection image. As a result, each detected cell occupy unique area in
whole image without overlap.

Fourthly, pick single cell from watershed split image and apply otsu-threshold.
Because otsu-threshold requires background and object information in image, we
apply above watershed segmentation to obtain individual cell image contained
nuclei, cytoplasm. We use watershed segmentation result at XY slice to mask 3D
individual cell image from raw data, then use a 512x512x20 box to mask it again.
Otsu-threshold is apply to individual cell to extract nuclei. Due to intensity
variation and photobleach, otsu-threshold is calculated for each cell at each time
independently. Final nuclei image are intensity larger than otsu-threshold,
otherwise set to zero.

Fifthly, normalize nuclei intensity for imaris analysis. Because imaris is not
flexible to set variable parameter for search RNA spot, we set minimal intensity
and mean intensity of nuclei to certain value.

## Count mRNA spot with imaris
In fact, at the begin, we try to use imaris do all the work, but imaris allocate
huge memory during calculation and fail to handle photobleach and cell intensity
variation. So we decide to simplify question: extract cell using own code, then
just let imaris count mRNA spot. My limited experience and time also suggest me
to use imaris at this time.


Creation Parameters:
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


## Licence
MIT License

Copyright (c) 2020 H.F.
