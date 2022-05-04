# ##########################################################################
# This script is used to launch the compilation of the RMarkdown report
# independently of Rstudio interface
# ##########################################################################

WORKING_DIR = getwd()

STEP = "09_differential_expression"

SCRIPT_DIR = file.path( WORKING_DIR, "03_Script")
OUTPUT_DIR = file.path( WORKING_DIR, "05_Output")

RDS_FILE = snakemake@input[["rds"]]
RDS = file.path( WORKING_DIR, RDS_FILE)

rmarkdown::render( input = file.path( SCRIPT_DIR, "diffexp.Rmd"),
                   output_dir = file.path( OUTPUT_DIR, STEP),
                   output_file  = "diffexp.html",
                   quiet = FALSE)
