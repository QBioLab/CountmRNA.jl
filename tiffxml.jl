using LightXML
"""
20200420 hf first version
"""

"""
Create and return XML string with xyzt
"""
function generatexml(SizeX::Integer, SizeY::Integer, SizeZ::Integer, SizeT::Integer)
    tiffXML = XMLDocument();
    OME = create_root(tiffXML, "OME");
    xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" 
    xmlns_xsi="http://www.w3.org/2001/XMLSchema-instance"
    Creator="OME Bio-Formats 5.2.2" 
    UUID="urn:uuid:a5ae8c1b-ac04-4544-97c1-bbdd0bdf8629" # use fixed value now
    xsi_schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd"
    OME_attributes = Dict("xmlns"=>xmlns, "xmlns:xsi"=>xmlns_xsi, "Creator"=>Creator, "UUID"=>UUID, "xsi:schemaLocation"=>xsi_schemaLocation)
    set_attributes(OME, OME_attributes)
    
    name = "test.tiff"
    Image= new_child(OME, "Image")
    set_attributes(Image, Dict("ID"=>"Image:0", "Name"=>name))
    
    Pixels = new_child(Image, "Pixels");
    #SizeC=1; SizeT=10;  SizeX=1300; SizeY=1900; SizeZ=2;
    PhysicalSizeX=PhysicalSizeY=0.108; PhysicalSizeZ=0.5;
    TimeIncrement=600;
    Pixels_attributes = Dict("ID"=>"Pixels:0","Type"=>"uint16", "BigEndian"=>"false",
    "DimensionOrder"=>"XYZCT", "SizeC"=>SizeC,"SizeX"=>SizeX, "SizeY"=>SizeY, "SizeZ"=>SizeZ, "SizeT"=>SizeT,
    "PhysicalSizeX"=>PhysicalSizeX, "PhysicalSizeY"=>PhysicalSizeY, "PhysicalSizeZ"=>PhysicalSizeZ,
    "PhysicalSizeXUnit"=>"µm", "PhysicalSizeYUnit"=>"µm", "PhysicalSizeZUnit"=>"µm",
    "TimeIncrement"=>TimeIncrement, "TimeIncrementUnit"=>"s");
    set_attributes(Pixels, Pixels_attributes)
    
    Channel_ = new_child(Pixels, "Channel");
    set_attributes(Channel_, Dict("ID"=>"Channel:0:0", "SamplesPerPixel"=>"1"));
    LightPath =new_child(Channel_, "LightPath");
    
    for i in 1:SizeC*SizeZ*SizeT
        new_child(Pixels, "TiffData")
    end
    
    tiffXML;
end

"""
Embed OME-XML to tiff with XYZT informance
"""
function embedxml(SizeX::Integer, SizeY::Integer, SizeZ::Integer, SizeT::Integer, img_name::String)
    tiffxml = generatexml(SizeX, SizeY, SizeZ, SizeT);
    print("Embedding OME-XML")
    run(`/home/hf/Bin/bftools/tiffcomment -set $tiffxml $img_name`);
    nothing;
end
