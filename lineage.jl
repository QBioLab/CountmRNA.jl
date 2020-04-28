using FileIO
using Plots
using Images
using ImageSegmentation

"""
Version Comment
0.1		initial
0.2		extract time-lines and split contacted branches. hf@0427
"""


# Get living time
livingtime(stack) = [ sum(stack[:, :, i])>0 for i in 1:size(stack)[3] ]
area_t(stack) = [ sum(stack[:, :, i]) for i in 1:size(stack)[3] ]

"""
Find long-lived trajactory
"""
function find_time_line(markers_t)
    # search connected components in 3D space 
    shortest_t = 90
	print("Finding connected component")
    @time time_line = label_components( markers_t.>0 )
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
Split contacted cell
"""
function split_contacted_cell!(old_time_line::Array{Int64,3},
        old_longlived_labels::Array{Int64,1}, old_living_time::Array{Array{Bool,1},1})
    t_len = size(old_time_line)[3]
    println("Detecting contacted branch")
    conn_z_number = []; # mark connected components more than 1
    for label in old_longlived_labels
        branches = time_line_2 .== label
        push!(conn_z_number, [maximum(label_components(branches[:,:, i])) .> 1 for i in 1:t_len])
    end
    contacted_labels = old_longlived_labels[sum.(conn_z_number) .> 10]
    print("Found contacted branch: ")
    println(contacted_labels)
    # if two more connected component touch to at same z slice more than 10 times
    # select area longer than 5e4
    #old_label = longlived_label_2[mean.(area_longlived) .> 5e4 ]

    # split 3d branch by split 2d cell slice by slice
    for contacted_label in contacted_labels
        dist_const = 30 # distance constant
        contacted_branch = old_time_line .==  contacted_label
        dist  = zeros(size( contacted_branch ))
        local_markers = zeros(Bool, size( contacted_branch ))
        println("Splitting branch $contacted_label now")
        for  t in 1:t_len
            dist[:,:, t] = distance_transform(feature_transform(.~contacted_branch[:,:,t]))
            local_markers[:,:,t] = dist[:,:,t] .> dist_const
            #split_water = watershed( .- dist[:,:,t], label_components(markers[:,:,t]))
            # 直接用 dist 可能会错误打断完整轨迹，需额外膨胀
            # 用分水岭 可能会误连，需额外腐蚀
        end
        local_time_line, local_longlived_labels, local_living_time = find_time_line( local_markers );
        println("Reassigning contacted branch $contacted_label ")
        # Assign new label to labels set and update mask

        if length(local_longlived_labels) > 0
            #TODO: what if all lines are shorter than 90 after breaking
            # remove old branch, old label, old living time.
            old_time_line .-= (( old_time_line .== contacted_label ).*contacted_label)
            index2remove =  old_longlived_labels .== contacted_label
            deleteat!( old_longlived_labels, index2remove )
            deleteat!( old_living_time, index2remove)
            for i in 1:length(local_longlived_labels)
                global_label = get_unique_label(old_longlived_labels)
                old_time_line .+= ((local_time_line .== local_longlived_labels[i]) .* global_label)
                push!( old_longlived_labels, global_label )
                push!( old_living_time, local_living_time[i] )
                println("Branch $contacted_label -> $global_label")
            end
        end
    end
    old_time_line, old_longlived_labels, old_time_line
end

"""
Return new unique label which don't exist in orignal set
set it larger than 100 
"""
function get_unique_label(_old_labels)
    i = 100 # start from 100 to distiguish with oldest labels
    while(true)
        if i in _old_labels
            i+=1
        else
            break
        end
    end
    i
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
