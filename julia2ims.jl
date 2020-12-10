using HDF5
using Dates 
using Images
#https://portal.hdfgroup.org/display/HDF5/Introduction+to+HDF5
# https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.format-Tuple{TimeType,AbstractString} 
# https://cpb-us-w2.wpmucdn.com/sites.wustl.edu/dist/f/1861/files/2019/02/02-imaris-image-properties-and-edit-2jj1avj.pdf
# https://github.com/tlambert03/imarispy/blob/master/imarispy/imaris.py

"""
Version Commit
0.1		hf, first verion
TODO: genrate thumb and fix timestamps
"""

"Save 4D UInt16 Array(x,y,z,t) to .ims file"
function  save2ims(array, fname::String="myfile.ims";
              subsamp=((1, 1, 1)),
              chunks=((16, 128, 128)),
              compression=2,
              thumbsize=256,
              dx=0.108, dz=0.5, dt_min=10)

    df = DateFormat("Y-mm-dd HH:MM:SS.sss")
    RecordingDate = DateTime(2020,5,7,14,19,22,0)

    nx, ny, nz, nt = size(array)
    nc = 1
    nr = 1     #nr = length(subsamp)

    GROUPS = [
        "DataSetInfo",
        "Thumbnail",
        "DataSetTimes",
        "DataSetInfo/Imaris",
        "DataSetInfo/Image",
        "DataSetInfo/TimeInfo"
    ]

    ATTRS = [
        ("/", "ImarisDataSet", "ImarisDataSet"),
        ("/", "ImarisVersion", "5.5.0"),
        ("/", "DataSetInfoDirectoryName", "DataSetInfo"),
        ("/", "ThumbnailDirectoryName", "Thumbnail"),
        ("/", "DataSetDirectoryName", "DataSet"),
        ("DataSetInfo/Imaris", "Version", "9.5"),
        ("DataSetInfo/Imaris", "ThumbnailMode", "thumbnailMIP"),
        ("DataSetInfo/Imaris", "ThumbnailSize", thumbsize),
        ("DataSetInfo/Image", "X", nx),
        ("DataSetInfo/Image", "Y", ny),
        ("DataSetInfo/Image", "Z", nz),
        ("DataSetInfo/Image", "NumberOfChannels", nc),
        ("DataSetInfo/Image", "Noc", nc),
        ("DataSetInfo/Image", "Unit", "um"),
        ("DataSetInfo/Image", "Description", "Hi"),
        ("DataSetInfo/Image", "MicroscopeModality", "Inverted Microscope"),
        ("DataSetInfo/Image", "RecordingDate", Dates.format(RecordingDate, df)),
        ("DataSetInfo/Image", "Name", "qblab"),
        ("DataSetInfo/Image", "ExtMin0", "0"),
        ("DataSetInfo/Image", "ExtMin1", "0"),
        ("DataSetInfo/Image", "ExtMin2", "0"),
        ("DataSetInfo/Image", "ExtMax0", nx * dx),
        ("DataSetInfo/Image", "ExtMax1", ny * dx),
        ("DataSetInfo/Image", "ExtMax2", nz * dz),
        ("DataSetInfo/Image", "LensPower", "40x"),
        ("DataSetInfo/TimeInfo", "DatasetTimePoints", nt),
        ("DataSetInfo/TimeInfo", "FileTimePoints", nt)
    ];

    COLORS = ("1 1 1", "1 0 1", "1 0 0", "0 0 1")
    for c in 0:nc-1
        grp = "DataSetInfo/Channel $c"
        push!(GROUPS, grp)
        push!(ATTRS, (grp, "ColorOpacity", 1))
        push!(ATTRS, (grp, "ColorMode", "BaseColor"))
        push!(ATTRS, (grp, "Color", COLORS[3]))
        push!(ATTRS, (grp, "GammaCorrection", 1))
        push!(ATTRS, (grp, "ColorRange", "0 255"))
        push!(ATTRS, (grp, "Name", "Channel $c"))
    end

    for t in 0:nt-1
        strr = Dates.format(RecordingDate + Dates.Minute(dt_min*t), df)
        #strr = "2020-05-07 {:02d}:{:02d}:{:02d}.000".format(h, m, s)
        push!(ATTRS, ("DataSetInfo/TimeInfo", "TimePoint$t", strr))
        #TODO: Imaris don't use these timestamps as time axis
    end

    h5open(fname, "w") do file
         for grp in GROUPS
            g_create(file, grp)
        end

        for info in ATTRS
            grp, key = info[1:2]
            value = "$(info[3])"
            h5a_write_S1(file[grp], key, value)
        end

        #file["Thumbnail/Data"] = UInt8.(maximum(array[:,:, 1:20,1], dims=3).>>8)
        file["Thumbnail/Data"] = UInt8.(maximum(array[:,:, 1:end,1], dims=3).>>8)

        # add data
        for t in 0:nt-1
            for c in 0:nc-1
                data = array[:, :, :, t+1]
                for r in 0:nr-1
                    edges, count = build_histogram(data, 256)
                    grp = g_create(file, "/DataSet/ResolutionLevel $r/TimePoint $t/Channel $c/")
                    grp["Histogram"] = UInt64.(count[1:end])
                    h5a_write_S1(grp, "HistogramMin", "1000")
                    h5a_write_S1(grp, "HistogramMax", "$(edges[end])")
                    grp["Data", "chunk", (256,256,4), "compress", compression] = data
                    #grp["Data", "chunk", (256,256,4)] = data
                    h5a_write_S1(grp, "ImageSizeX", "$(size(data)[1])")
                    h5a_write_S1(grp, "ImageSizeY", "$(size(data)[2])")
                    h5a_write_S1(grp, "ImageSizeZ", "$(size(data)[3])")
                end
            end
        end
    end
	fname
end

"Encode string with ASCII S1"
function h5a_write_S1(parent, key, value)
    dspace = dataspace((length(value),))
    try
        attr = a_create(parent, key,  HDF5Datatype(HDF5.H5T_C_S1), dspace)
        HDF5.writearray(attr, HDF5.H5T_C_S1, value)
    finally
        close(dspace)
    end
    nothing
end

function str2int(str)
    local tmp =""
    for i in str
        tmp=tmp*i
    end
    parse(Int, tmp)
end

function loadims(imsname::String)
    local attr_t_len = h5readattr(imsname, "DataSetInfo/TimeInfo")["DatasetTimePoints"];
    local t_len  = str2int(attr_t_len)
    local z_depth = 20
    local img = zeros(Gray{N0f16}, 1300, 1900, t_len*z_depth);
    h5open(imsname) do imsfile
        for i in 0:t_len-1
			data = imsfile["DataSet/ResolutionLevel 0/TimePoint $i/Channel 0/Data"]
            img[:, :, i*20+1:20*(i+1)] = reinterpret(N0f16,data[1:1300,1:1900,1:20]);
        end
    end
    img
end

function loadims2(imsname::String)
    local attr_t_len = h5readattr(imsname, "DataSetInfo/TimeInfo")["DatasetTimePoints"];
    local t_len  = str2int(attr_t_len)
    local z_depth = 20
    local img = zeros(Gray{N0f16}, 512, 512, t_len*z_depth);
    h5open(imsname) do imsfile
        for i in 0:t_len-1
			data = imsfile["DataSet/ResolutionLevel 0/TimePoint $i/Channel 0/Data"]
            img[:, :, i*20+1:20*(i+1)] = reinterpret(N0f16,data[1:512,1:512,1:20]);
        end
    end
    img
end
