


#################### retrieve_ISA

# This rule is to retrieve the Spanish ISA rasters

rule retrieve_isa_onwind:
    message:
        "... Retrieving ISA index from MITECO for onwind carrier"    
    output:
        tiff_file="data/ISA/Clas_ISA_eol_pb.tiff",
    run:
        from pathlib import Path
        from zipfile import ZipFile
        import urllib.request
        import io

        url = "https://www.miteco.gob.es/content/dam/miteco/es/calidad-y-evaluacion-ambiental/temas/evaluacion-ambiental-de-planes-programas-y-proyectos/Zonificacion_EOL_clasificada_2024.zip"

        Path(output.tiff_file).parent.mkdir(parents=True, exist_ok=True)

        with urllib.request.urlopen(url) as response:
            with ZipFile(io.BytesIO(response.read())) as zf:
                with zf.open("Clas_ISA_eol_pb.tiff") as src, open(output.tiff_file, "wb") as dst:
                    dst.write(src.read())        


rule retrieve_isa_solar:
    message:
        "... Retrieving ISA index from MITECO for solar carrier"    
    output:
        tiff_file="data/ISA/Clas_ISA_ftv_pb.tiff",
    run:
        from pathlib import Path
        from zipfile import ZipFile
        import urllib.request
        import io

        url = "https://www.miteco.gob.es/content/dam/miteco/es/calidad-y-evaluacion-ambiental/temas/evaluacion-ambiental-de-planes-programas-y-proyectos/Zonificacion_FTV_clasificada_2024.zip"

        Path(output.tiff_file).parent.mkdir(parents=True, exist_ok=True)

        with urllib.request.urlopen(url) as response:
            with ZipFile(io.BytesIO(response.read())) as zf:
                with zf.open("Clas_ISA_ftv_pb.tiff") as src, open(output.tiff_file, "wb") as dst:
                    dst.write(src.read())        
