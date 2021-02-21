using LightXML
using JLD

"""
Return the mapping between image name and image id of imaris file from "TreeData.xml"
Version Comment
0.1     first version
"""


"Feed TreeData.xml and return mapping dictionary"
function get_ims_mapping(file_name::String)
    # Dictionary: Dict(date=> Dict(image_name=>image_id or image_id=>image_name)) 
    mapping = Dict() 
    xdoc = parse_file(file_name)
    # enter root node
    xroot = root(xdoc)
    for c in child_nodes(xroot) 
        if is_elementnode(c)
            e = XMLElement(c)
            # enter first "Folder" node
            for c2 in child_nodes(c)
                e2 = XMLElement(c2)
                # enter second "Folder" node: mRNA
                if (attribute(e2, "name") == "mRNA")
                    for e3 in child_nodes( e2 )
                        exp = XMLElement( e3 )
                        # enter third "Folder" node: for each experiment
                        if name(exp) == "Folder"
                            exp_name = attribute(exp, "name")
                            mapping_cell = Dict()
                            println(exp_name)
                            # enter "Image" node for each view
                            for e4 in child_nodes( exp )
                                cell = XMLElement( e4 )
                                # travel all image
                                if name(cell) == "Image"
                                    image_name = attribute(cell, "name")
                                    image_id = replace(attribute(cell, "id"), "http://localhost:8048/data_service/"=>"")
                                    push!(mapping_cell, image_name=>image_id)
                                    push!(mapping_cell, image_id=>image_name)
                                    #println(image_name, " ", image_id)
                                end
                            end
                            push!(mapping, exp_name=>mapping_cell)
                        end
                    end
                    #println(attribute(e2, "name"))
                end
            end
        end
    end
    return mapping
end

mapping_dict = get_ims_mapping( ARGS[1] );

save("$(dirname(ARGS[1]))/id_mapping.jld", mapping_dict)
