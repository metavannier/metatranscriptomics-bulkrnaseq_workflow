# ----------------------------------------------
# MicrobeAnnotator: database building
# ----------------------------------------------

rule MicrobeAnnotator_build:
  output:
    microbeannotatordb = REF + config["microbeannotator"]["folder"] + "/microbeannotator.db",
    conversiondb = REF + config["microbeannotator"]["folder"] + "/conversion.db"

  params:
    foldermicrobeannotator = REF + config["microbeannotator"]["folder"],
    method = config["microbeannotator"]["method"],
    threads = config["microbeannotator"]["threads"],
    lightversion = config["microbeannotator"]["lightversion"]

  conda: 
    CONTAINER + "microbeannotator.yaml"

  shell:
    """
    python --version
    lightversion={params.lightversion}
    if [ ${{lightversion}} == 'light' ]
    then
        microbeannotator_db_builder -d {params.foldermicrobeannotator} -m {params.method} -t {params.threads} --light;
    else
        microbeannotator_db_builder -d {params.foldermicrobeannotator} -m {params.method} -t {params.threads};
    fi
    """

# ----------------------------------------------
# MicrobeAnnotator: annotation
# ----------------------------------------------

rule MicrobeAnnotator_annotation:
  output:
    MicrobeAnnotatorannotation = config["microbeannotator"]["annotationfolder"] + "/mockfile.txt",

  params:
    protfiles = REF + config["microbeannotator"]["protfiles"],
    annotationfolder = REF + config["microbeannotator"]["annotationfolder"],
    method = config["microbeannotator"]["method"],
    threads = config["microbeannotator"]["threads"],
    processes = config["microbeannotator"]["processes"],
    foldermicrobeannotator = REF + config["microbeannotator"]["folder"],
    
  conda: 
    CONTAINER + "microbeannotator.yaml"

  shell:
    """
    pip install hmmer==0.1.0
    microbeannotator -i {params.protfiles} -d {params.foldermicrobeannotator} -o {params.annotationfolder} -m {params.method} -p {params.processes} -t {params.threads}
    touch {params.annotationfolder}/mockfile.txt
    """