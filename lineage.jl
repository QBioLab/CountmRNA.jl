using FileIO
using Plots
using Images
using ImageSegmentation

"""
Version Comment
0.1		initial
#TODO: add continue state space
"""


# Get living time
livingtime(stack) = [ sum(stack[:, :, i])>0 for i in 1:size(stack)[3] ]

"""
Find long-lived trajactory
"""
function find_time_line(markers_t)
    # search connected components in 3D space 
    shortest_t = 90
    time_line = label_components( markers_t.>0 )
    line_amount = maximum(time_line)
    # More advanced and fine punch and merge could be done
    # But we just select long living trajactory, remove short-lived one
    living_time =[ livingtime(time_line.==line) for line in 1:line_amount ]
    living_length = [sum(living_time[line]) for line in 1:line_amount]
    # label 0 mean background
    shortlived = living_length .< shortest_t
    longlived_label = (1:line_amount)[.~shortlived]
    #longlived_time = living_time[]
    for line in (1:line_amount)[shortlived]
        time_line .*= (time_line.≠line)
    end
    time_line, longlived_label, live_time[longlived_label]
end





"""
Choose cloest neighborhood
"""
function find_nearest(root, next_level, maxdistance)
    nucleus_distance = [sum(abs.(root .- next_level[i])) for i in 1:length(next_level)]
    mindist, mindist_index = findmin(nucleus_distance);
    if mindist < maxdistance
        return mindist_index
    else
        return nothing
    end
end

"""
Find every child with given root
"""
function buildtree(root::Int, position)
    tree = []
    height = length(position)
    level = 2
    child = find_nearest(position[1][root], position[level], 100);
    while( child ≠ nothing && level<height)
        push!(tree, child)
        child = find_nearest(position[level][child], position[level+=1], 100)
        #print("$level ")
    end
    tree
end

"""
Extract long-life path
"""
function generate_path(pos_map)
    all_path = []
    for root in 1:length(pos_map[1])
        tree = buildtree(root, pos_map);
        if length(tree) > 70
            path = [pos_map[i+1][tree[i]] for i in 1:length(tree) ];
            push!(all_path, path)
        end
    end
    all_path;
end

"""
Using given mask to export roi of cell 
"""
function export_cell(track, imgs, edge_mask)
    height, width = size(imgs[:,:,1]);
    crop_img = zeros(700, 700, length(track)*20);
    for t in 1:length(track)
        x, y = Int.(floor.(track[t]));
        xmin, xmax, ymin, ymax = x-349 , x+350,  y-349, y+350
        if xmin < 1
            xmax = xmax - xmin + 1
            xmin = 1
        end
        if xmax > 1900
            xmin = xmin - (xmax-1900)
            xmax = 1900
        end
        if ymin < 1
            ymax = ymax - ymin +1
            ymin = 1
        end
        if ymax > 1300
            ymin = ymin - (ymax-1300)
            ymax = 1300
        end
        masked_cell= (edge_mask[:,:,t] .== edge_mask[x, y, t] )
        for i in 1:20
            crop_img[:,:,20*(t-1)+i] =(masked_cell.*imgs[:,:,20*(t-1)+i])[xmin:xmax, ymin:ymax]
        end
    end
    crop_img;
end
