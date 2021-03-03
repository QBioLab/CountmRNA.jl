using LightXML
using JLD
using HDF5 
"""
Usage:
    julia ims_cp.jl imaris_root_dir exp_folder_name imaris_creation_name destination

Version  Commit
0.1      first version
TODO: combine scene and image file
"""

"copy and rename imaris batch result from imaris database with Folder and Creation name"
function ims_cp(ims_root::String, ims_array::String, ims_creation::String, output::String)
    id2name = Dict()
    id2path= Dict()
    # load ims TreeData.xml
    TREEDATE = parse_file("$ims_root/TreeData.xml")
    # enter root node
    xroot = root(TREEDATE)
    
    # local -> mRNA -> ims_array    
    for folder in child_elements(xroot) 
        if (name(folder) =="Folder") && (attribute(folder, "name") == "local")
            for folder2 in child_elements(folder)
                if (name(folder2) == "Folder") && (attribute(folder2, "name") == "mRNA")
                    for folder3 in child_elements(folder2)
                        # enter single experiment folder
                        #println(name(folder3)," ", attribute(folder3, "name"))
                        if (name(folder3) == "Folder") && (attribute(folder3, "name") == ims_array)
                            # enter "Image" node for each folder
                            for e in child_elements( folder3 )
                                # travel all image
                                if name(e) == "Image"
                                    image_name = attribute(e, "name")
                                    image_id = attribute(e, "id")
                                    image_path = replace( attribute(e, "path"), 
                                                         "file:///C%3A/Imaris_Data_qblab/"=>"" )
                                    push!(id2name, image_id=>image_name)
                                    push!(id2path, image_id=>image_path)
                                elseif name(e) == "CreationParameter"
                                    # get CreationParameter id
                                    creation_parameters = attribute(e, "CreationParameters")
                                    open("$output/$ims_creation.txt", "w") do io
                                        for line in split(creation_parameters, "> <")
                                            println(io, line)
                                        end
                                    end
                                end
                            end
                            # cp batch result
                            for e in child_elements( folder3 )
                                if (name(e) == "BatchRun") && (attribute(e, "name") == "$ims_array $ims_creation")
                                    #println(e["batchrunresults"]["result"])
                                    batchrunresults = e["batchrunresults"][1]["result"]
                                    for result in batchrunresults
                                        if attribute(result, "status") == "FINISHED"
                                            id = attribute(result, "input")
                                            outputpath = attribute(result, "outputpath")
                                            src_scene = replace(outputpath, "file:/"=>"$ims_root")
                                            src_ims = "$ims_root/$(id2path[id])"
                                            out_ims = "$output/$(id2name[id]).ims";
                                            print(id2name[id], " ")
                                            cp(src_ims, out_ims, force=true)
                                            ims_f = h5open(out_ims, "r+")
                                            scene_f = h5open(src_scene, "r")
                                            if haskey(ims_f, "Scene")
                                                HDF5.delete_object(ims_f, "Scene")
                                            end
                                            if haskey(ims_f, "Scene8")
                                                HDF5.delete_object(ims_f, "Scene8")
                                            end
                                            HDF5.copy_object(scene_f, "Scene", ims_f, "Scene")
                                            HDF5.copy_object(scene_f, "Scene8", ims_f, "Scene8")
                                            close(ims_f)
                                            close(scene_f)
                                        end
                                    end
                                else
                                    #println("$ims_creation no found!")
                                end
                            end
                        else
                            #println("$ims_array no found!")
                        end
                    end
                end
            end
        end
    end
end

if length(ARGS) != 4
    println("Usage: julia ims_cp.jl imaris_root_dir exp_folder_name imaris_creation_name destination")
else
    ims_cp(ARGS[1], ARGS[2], ARGS[3], ARGS[4])
end

