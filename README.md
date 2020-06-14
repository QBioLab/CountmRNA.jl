# CountmRNA.jl

A julia collection to tracking single cell and count mRNA inside each cell.

For our live cell fluorescence image by spin-disk confocal microscope, we only
track cell in own algorithm and recognize mRNA with imaris.

(_You can take anything if can help your study!_)

## Files
| NAME              | FUNCTION
| splitcell.jl      | Use Laplace of Gaussian(LoG) filter extract cell from raw 3d image
| lineage.jl        | Use connected component to find each track and use watershed to split
| segmentation3d.jl | Use Otsu-threshold to extract 3d nuclei 
| normalization3d.jl| Normalize minimal and mean intensity of nuclei
| julia2ims.jl      | Useful functions to load and save imaris 5 file
| tiffxml.jl        | Useful functions to save tiff with OME-TIFF info

## Algorithm

