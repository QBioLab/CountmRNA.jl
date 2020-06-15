# CountmRNA.jl

A julia script collection to tracking single cell and count mRNA inside each cell.

For our live cell fluorescence image by spin-disk confocal microscope, we only
track cell in own algorithm and recognize mRNA with imaris.

_You can take anything if can help your study! julia2ims.jl and tiffxml.jl
might save your time if you work with imaris/imagej, even you are not
intersting at cell track algorithm_

## Files

| NAME              | FUNCTION
| splitcell.jl      | Use Laplace of Gaussian(LoG) filter extract cell from raw 3d image
| lineage.jl        | Use connected component to find each track and use watershed to split
| segmentation3d.jl | Use Otsu-threshold to extract 3d nuclei 
| normalization3d.jl| Normalize minimal and mean intensity of nuclei
| julia2ims.jl      | Useful functions to load and save imaris 5 file
| tiffxml.jl        | Useful functions to save tiff with OME-TIFF info


## Algorithm
In general, I use watershed to split cell/nuclei, and scan 3d connected component
find track. To extract 3d nuclei, we apply otsu-threshold. Finally, extracted 
nuclei are normalized and loaded into imaris.

These algorithm are special to our data feature and analysis target, which low
compute cost efficiently but are not for general case. Although we are processing 4D
image stack(XYZT), we deicide to located, split, track cell/nuclei from
z-project 3D images(XYT). Finally, extra 3D nuclei from known location. Because
not much cell inside single slice(about 10 cells) and cells take large area of
view( we use 60x objective). 

Firstly, extra binary cell mask from raw image. We apply medium filter(5x5) and
LoG filter(40) on each time point of z-project image. 

Secondly, we search connected component to figure out cell trajectories in 3D
images(XYT). Because our cells move slow between adjacent frame and short time
point(compare to pixel number of X/Y dimension), we are able to
use connected component. But more meaningful point is that it use much more
direct logic and clear figure to handle cell come/leave from camera edge
instead of linking closest neighbour by cell centroid. Each detected connected
components are assigned unique id. By filtering connected component length in
t-dimension, short trajectories are remove. Then we detect possible collision
in each long trajectories, which scan connected component number for each time
slice of single trajectory. Collision are split by a distance transform and
resign new unique id. Now, we are able to get rough cell/nuclei mask. But it is
not enough fine to extract 3D nuclei.

Thirdly, we use above nuclei mask as marker and performance watershed segmentation
at each time point of z-projection image. It split whole image for different cell.  

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
Recognizing 3D object is beyond my knowledge and experiment, so I use imaris
(Test on Imaris 9.5).

Imaris Parameters:


## Licence
MIT License
Copyright (c) 2020 H.F.
