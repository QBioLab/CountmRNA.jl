"""
Normalize brightness for each time point
Version Comment
0.1     hf, initial version
"""

using Images
using Statistics
using LsqFit
using Clustering

"""
Extract feature from image
"""
function extract_feature(imgs)
    z_depth = 20
    t_len = size(imgs)[3] ÷ z_depth
    low_line = zeros(t_len)
    mean_line  = zeros(t_len)
    for t in 1:t_len
        #nozeros_pixel = [real(i) for i in imgs[:,:, (t-1)*20+1:20*t] if i ≠ 0]
		nozeros_pixel = imgs[:,:, (t-1)*20+1:20*t][imgs[:,:, (t-1)*20+1:20*t].≠0]
		if length(nozeros_pixel) == 0 
			println("no cell in time point $t")
			continue # skip empty stack
		end
        low_line[t] = minimum(nozeros_pixel)
        mean_line[t] = mean(nozeros_pixel)
    end
    low_line, mean_line
end

medianfilter1(v, ws) = [median(v[i:(i+ws-1)]) for i=1:(length(v)-ws+1) ]
medianfilter(v, ws) = medianfilter1(vcat(0,v,0),ws)


"""
Use k-means to remove outliers
Input raw data series, return index without outliers
"""
function kmeans_rm_outliers( raw_line )
	# TODO: if only one cluster?
	t_len = length(raw_line)
    X = zeros(2, t_len)
    X[1, :] = fill(0.02, t_len) # just choose 0.02 in random.
    X[2, :] = raw_line; # forgive me just add one more dimension for k-means
    result = kmeans(X, 2; init=:kmpp)
    if result.counts[1]> result.counts[2] # chose the one have more points
        index = 1
    else
        index = 2
    end
    t_nooutliers = [ i  for i in 1:t_len  if result.assignments[i] == index ]
end


"""
Removing outliers and smoothing feature
"""
function correct_feature( low_line, mean_line )
    t_len = length(low_line)
    ws = 3 # use small window to remove spark
    low_line = medianfilter(low_line, ws)
    mean_line = medianfilter(mean_line, ws)
   
	# TODO: if only one cluster? Add additional judge by std(mean_line[1:end-1].-mean_line[2:end])
	t_nooutliers = kmeans_rm_outliers(mean_line)
    mean_line_nooutliers = mean_line[t_nooutliers]
	# assume mean_line and low_line share same shape
    low_line_nooutliers = low_line[t_nooutliers]
    @. model(x, p) = p[1]*exp(-x*p[2]) # fit with expontional function
    fit_mean = curve_fit(model, t_nooutliers, mean_line_nooutliers, [0.003, 0.0003])
    fit_low = curve_fit(model, t_nooutliers, low_line_nooutliers, [0.002, 0.0003])
    corrected_mean_line = model(1:t_len, fit_mean.param)
    corrected_low_line = model(1:t_len, fit_low.param) 
	# TODO: what if mean.(fit, no_outliers) ?
    
    corrected_low_line, corrected_mean_line
end

function normalize(imgs)
    z_depth = 20
    t_len = size(imgs)[3] ÷ z_depth
    d1_len, d2_len  = size(imgs)[1],  size(imgs)[2]
    low_line, mean_line = extract_feature(imgs)
    corrected_low_line, corrected_mean_line = correct_feature(low_line, mean_line)
    imgs_norm = zeros(size(imgs))
    
    for t in 1:t_len
        imgs_norm[:,:,(t-1)*20+1:20*t] = ((imgs[:,:,(t-1)*20+1:20*t] .- corrected_low_line[t]) ./ (corrected_mean_line[t]-corrected_low_line[t]).*0.005 .+ 0.02) .* ( imgs[:, :, (t-1)*20+1:20*t] .>0)
        #print("$t ")
    end
    
	imgs_norm, Dict("low_line"=>low_line, "mean_line"=>mean_line, 
					"corrected_low_line"=>corrected_low_line, "corrected_mean_line"=>corrected_mean_line)
end

#@time s14c1 = load(File(format"TIFF", "$data_dir/$(cell_name[14])"));
#@time imgs_norm, low_line, mean_line = normalize(s14c1);
#plot(low_line)
#plot!(low_line_raw, marker=(:star))
#plot!(mean_line)
#plot!(mean_line_raw, marker=(:dot))
#plot!(ylim=(0,0.004))
