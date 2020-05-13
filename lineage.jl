using Images
using ImageSegmentation

"""
Version Comment
0.1		initial
0.2		extract time-lines and split contacted branches. hf@0427
0.2.1	add error handle to avoid no longlived_cell hf@0512
"""

# Get living time
cal_livingtime(stack) = [ sum(stack[:, :, i])>0 for i in 1:size(stack)[3] ]
area_t(stack) = [ sum(stack[:, :, i]) for i in 1:size(stack)[3] ]
"""
function component_lengths(img::AbstractArray{Int})
    n = zeros(Int64,maximum(img)+1)
    for i=1:length(img)
        n[img[i]+1]+=1
    end
    n
end
"""
function get_unique_label(_old_labels, number::Int=1)
    i = 100 # start from 100 to distiguish with oldest labels
	labels = zeros(Integer, number)
	for count in 1:number
	    while( i ∈ _old_labels)
        	i+=1
		end
		_old_labels = union(_old_labels, i)
		labels[count] = i
    end
    labels
end

"Find long-lived track by searching connected components in 3D"
function find_time_line(markers_t)
    local shortest_t = 90
    println("Finding connected component")
    local time_line = label_components( markers_t.>0 )
    local time_line_whole = copy(time_line)
    local line_amount = maximum(time_line)
    # More advanced and fine punch and merge could be done
    # But we just select long living trajactory, remove short-lived one
    # Calculate living length
    x_len, y_len, t_len = size(time_line)
    local n = zeros(Int64, t_len, line_amount+1)
    # label 0 mean background
    @inbounds Threads.@threads for t in 1:t_len
        @inbounds for x in 1:x_len
            @inbounds for y in 1:y_len
            	n[t, time_line[x, y, t]+1] += 1
            end
        end
    end
    living_time =[ n[:, line].>0 for line in 2:line_amount+1 ]
    longlived = [sum(living_time[line]) for line in 1:line_amount] .> shortest_t
    local longlived_label = (1:line_amount)[longlived]
    # Remove shortlived branches 
    @inbounds for I in CartesianIndices(size(time_line))
        if time_line[I] ∉ longlived_label 
            time_line[I] = 0
        end
    end
    time_line, longlived_label, living_time[longlived_label], time_line_whole
end

function find_index4label(_labels)
	index = 1
	indexs = zeros(Int, maximum(_labels))
	for label in _labels
		indexs[label]=index
		index+=1
	end
	indexs
end

"Split contacted cell by watershed"
function split_contacted_cell!(old_time_line::Array{Int64,3}, 
        old_longlived_labels::Array{Int64,1}, old_living_time::AbstractArray,
        old_time_line_whole::Array{Int64, 3} )
    t_len = size(old_time_line)[3]
	if length(old_longlived_labels) == 0
		println("No longlived cell is found, stopping")
	end
    println("Detecting contacted branch")
	local conn_z_t = zeros(Int, length(old_longlived_labels),t_len)
	local label2index = find_index4label(old_longlived_labels)
	@time for t in 1:t_len
		Threads.@threads for label in old_longlived_labels
			local branches = old_time_line[:, :, t] .== label
			conn_z_t[label2index[label], t] = maximum(label_components(branches))
		end
	end
	contacted_labels = old_longlived_labels[sum(conn_z_t.>1, dims=2)[:] .> 3]
    #If two more connected component touch to at same z slice more than 10 times
    print("Found contacted branch: ")
    println(contacted_labels)
    # select area longer than 5e4
    #old_label = longlived_label_2[mean.(area_longlived) .> 5e4 ]
    
    # split 3d branch by split 2d cell slice by slice
    for contacted_label in contacted_labels
        dist_const = 30 # distance constant 
        local contacted_branch = old_time_line .==  contacted_label
        local dist  = zeros(size( contacted_branch ))
        local local_markers = zeros(Bool, size( contacted_branch ))
        println("Splitting branch $contacted_label now")
        Threads.@threads for  t in 1:t_len # TODO: only split when connected component descrease
            dist[:,:, t] = distance_transform(feature_transform(.~contacted_branch[:,:,t]))
            local_markers[:,:,t] = dist[:,:,t] .> dist_const
            #split_water = watershed( .- dist[:,:,t], label_components(markers[:,:,t]))
            # 直接用 dist 可能会错误打断完整轨迹，需额外膨胀
            # 用分水岭 可能会误连，需额外腐蚀
        end
        local_time_line, local_longlived_labels, local_living_time, local_time_line_whole = find_time_line( local_markers );
        
        println("Reassigning contacted branch $contacted_label ")
        # Assign new label to labels set and update mask
        if length(local_longlived_labels) > 0
            #TODO: what if all lines are shorter than 90 after breaking
            # remove old branch, old label, old living time.
            @inbounds for point in eachindex(old_time_line)
                if old_time_line[point] == contacted_label 
                    old_time_line[point] = 0
                    old_time_line_whole[point] = 0
                end
            end
            index2remove =  old_longlived_labels .== contacted_label
            deleteat!( old_longlived_labels, index2remove )
            deleteat!( old_living_time, index2remove)
            for i in 1:length(local_longlived_labels)
				global_label = get_unique_label(old_longlived_labels)[1]
				# update time_line index 
                @inbounds Threads.@threads for point in eachindex(local_time_line)
                    if local_time_line[point] == local_longlived_labels[i]
                        old_time_line[point] = global_label
                        old_time_line_whole[point] = global_label
					end
                end
                push!( old_longlived_labels, global_label )
                push!( old_living_time, local_living_time[i] )
                println("Branch $contacted_label -> $global_label")
            end

            # add remaining markers to time_line_whole
			local_label_max = maximum(local_time_line_whole)
			global_label_max = maximum(old_time_line_whole)
			tmp = get_unique_label(1:global_label_max, local_label_max - length(local_longlived_labels))
			tmp_index = 1 
			global_labels = zeros(Integer, local_label_max+1)
            for i in 1:local_label_max
                if i ∉ local_longlived_labels
					# preserve 1th index for 0, so start at 2
					global_labels[i+1] = tmp[tmp_index]
					tmp_index += 1
				end
			end
            @inbounds Threads.@threads for point in eachindex(local_time_line)
				γ = global_labels[local_time_line_whole[point]+1] 
				if γ ≠ 0
                   old_time_line_whole[point] = γ
                end
			end
        end
    end
    #old_time_line, old_longlived_labels, old_living_time, old_time_line_whole
	nothing
end


"extract walking pathway, [xy, time, labels]"
function walking(_time_line, _longlived_labels, _livingtime)
    _tracks = zeros( 2, length(_livingtime[1]), length(_longlived_labels))
    _max_index = maximum(_time_line)
    @inbounds for t in 1: length(_livingtime[1])
        tmp = component_centroids_lables(_time_line[:, :, t], _longlived_labels, _max_index)
        @inbounds for i in 1:length(_longlived_labels)
            _tracks[:, t, i] = [ tmp[i][1]  tmp[i][2] ]
        end
    end
    _tracks #[pos, t, branch]
    #TODO: soomth pathway to make sence in biology
end


" 分封领地"
function grant_domain( _raw_imgs, _time_line, _longlived_labels, _livingtime, _time_line_whole)
    z_depth = 20
    t_len = length(_livingtime[1])
    h, w = size(_raw_imgs[:,:,1]);
    watershed_maps = zeros(Integer, h, w, t_len)
    println("Grant domain for each detected cell by watershed")
    @time @inbounds Threads.@threads for t in 1:t_len 
        #print("-")
        watershed_maps[:, :, t] = labels_map(watershed(
                .-maximum(_raw_imgs[:,:,z_depth*(t-1)+1:z_depth*t],dims=3)[:,:,1], _time_line_whole[:,:,t]) )
    end
    
	println("")
	println("Drawing longlived_maps")
    longlived_maps = zeros(Integer, h, w, t_len)
	@inbounds Threads.@threads for I in eachindex(longlived_maps)
        if watershed_maps[I] ∈ _longlived_labels
			longlived_maps[I] = watershed_maps[I]
        end
    end

	println(_longlived_labels)
    longlived_maps, watershed_maps
    #longlived_maps
end

" Using given mask to export roi of cell "
function pick_cell(_raw_imgs::Array{Gray{Normed{UInt16,16}},3}, _longlived_maps::Array{Integer,3},
					cell_label::Int64, cell_tracks::Array{Float64,2}, cell_livingtime::BitArray{1})
    z_depth = 20
    t_len = length(cell_livingtime)
    h, w = size(_raw_imgs[:, :, 1])
    cell_img = zeros(Gray{Normed{UInt16,16}}, 512, 512, t_len*z_depth)
    #cell_img = zeros(512, 512, t_len*z_depth)
    
    @inbounds Threads.@threads for t in (1:t_len)[cell_livingtime] 
    #for t in (1:t_len)[cell_livingtime]
		#print("$t ")
		# only choose frame when object exist
        bounder = box(cell_tracks[:, t], h, w)
        #cell_img[:, :, (t-1)*z_depth+1:t*z_depth] = 
		#	_raw_imgs[bounder[1], bounder[2], (t-1)*z_depth+1:t*z_depth] .*
        #    (_longlived_maps[bounder[1], bounder[2], t] .== _longlived_labels[_index])
		d₁ref = (bounder[1])[1] -1
		d₂ref = (bounder[2])[1] -1
		z = (t-1)*z_depth+1 : t*z_depth
        @inbounds for d₁ in bounder[1] 
            @inbounds for  d₂ in bounder[2]
                if _longlived_maps[d₁, d₂, t] == cell_label
					#cell_img[d₁-d₁ref, d₂-d₂ref, z] = copy(_raw_imgs[d₁, d₂, z])
					cell_img[d₁-d₁ref, d₂-d₂ref, z] = _raw_imgs[d₁, d₂, z]
                end
            end
        end
    end
    cell_img
end

"`component_centroids(labeled_array)` -> an array of centroids for selected lables, excluding the background label 0"
function component_centroids_lables(img::AbstractArray{Int,N}, labels::AbstractArray{Int, 1}, max_index::Integer) where N
    len = length(0:maximum(max_index))
    #len = length(0:maximum(labels))
    n = fill(zero(CartesianIndex{N}), len)
    counts = fill(0, len)
    @inbounds for I in CartesianIndices(size(img))
        v = img[I] + 1
        n[v] += I
        counts[v] += 1
    end
    map(v -> n[v].I ./ counts[v], labels.+1)
end

" Return box with fixed height and width"
function box(cell_center, _d₁max, _d₂max)
    α = 256 
    Α = 2*α 
    d₁, d₂ = Int.(floor.(cell_center))
    d₁min, d₁max, d₂min, d₂max = d₁-α+1, d₁+α,  d₂-α+1, d₂+α
    if d₁min < 1
        d₁max = Α
        d₁min = 1
    end
    if d₁max > _d₁max
        d₁min = _d₁max - Α + 1
        d₁max = _d₁max 
    end
    if d₂min < 1
        d₂max = Α
        d₂min = 1
    end
    if d₂max > _d₂max
        d₂min = _d₂max - Α + 1
        d₂max = _d₂max
    end
    d₁min:d₁max, d₂min:d₂max
end
