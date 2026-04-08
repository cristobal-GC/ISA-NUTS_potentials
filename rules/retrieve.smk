


#################### retrieve_ISA

# This rule is to retrieve the Spanish ISA rasters

rule retrieve_isa_onwind:
    message:
        "... [retrieve_isa_onwind] Retrieving ISA index from MITECO for onwind carrier"    
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
        "... [retrieve_isa_solar] Retrieving ISA index from MITECO for solar carrier"    
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


#################### retrieve_GEBCO

GEBCO_CFG = config.get("gebco", {})

rule retrieve_gebco:
    message:
        "... [retrieve_gebco] Retrieving GEBCO bathymetry data"
    params:
        url=GEBCO_CFG["url"],
    output:
        gebco=f"{GEBCO_CFG['folder']}/{GEBCO_CFG['file_name']}",
    run:
        from pathlib import Path
        import urllib.request

        output_folder = Path(output["gebco"]).parent
        output_folder.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(params.url, output["gebco"])
