# ----------------------------------------------
# MicrobeAnnotator: database building
# ----------------------------------------------

rule MicrobeAnnotator_build:
  output:
    microbeannotatordb = config["microbeannotator"]["folder"] + "/microbeannotator.db",
    conversiondb = config["microbeannotator"]["folder"] + "/conversion.db"

  params:
    foldermicrobeannotator = config["microbeannotator"]["folder"],
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
    protfiles = config["microbeannotator"]["protfiles"],
    annotationfolder = config["microbeannotator"]["annotationfolder"],
    method = config["microbeannotator"]["method"],
    threads = config["microbeannotator"]["threads"],
    processes = config["microbeannotator"]["processes"],
    foldermicrobeannotator = config["microbeannotator"]["folder"],
    
  conda: 
    CONTAINER + "microbeannotator.yaml"

  shell:
    """
    pip install hmmer==0.1.0
    microbeannotator -i {params.protfiles} -d {params.foldermicrobeannotator} -o {params.annotationfolder} -m {params.method} -p {params.processes} -t {params.threads}
    touch {params.annotationfolder}/mockfile.txt
    """