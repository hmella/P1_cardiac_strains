import argparse

# Initiate the parser
parser = argparse.ArgumentParser()

# add long and short argument
help_str = "Choose if paramters are generated or not"
parser.add_argument("--generate-parameters", "-gp", help=help_str,
    action='store_true')

help_str = "Choose if CSPAMM images are generated or not"
parser.add_argument("--cspamm", "-c", help=help_str,
    action='store_true')

help_str = "Choose if DENSE images are generated or not"
parser.add_argument("--dense", "-d", help=help_str,
    action='store_true')

help_str = "Choose if REF images are generated or not"
parser.add_argument("--reference", "-r", help=help_str,
    action='store_true')

help_str = "Number of spins samples"
parser.add_argument("--nb-samples", "-Ns", help=help_str, type=int)

help_str = "Choose initial data"
parser.add_argument("--initial-data", "-i", help=help_str, type=int)

help_str = "Choose final data"
parser.add_argument("--final-data", "-f", help=help_str, type=int)

# read arguments from the command line
args = parser.parse_args()