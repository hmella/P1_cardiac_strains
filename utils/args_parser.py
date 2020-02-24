import argparse

# Initiate the parser
parser = argparse.ArgumentParser()

# add long and short argument
parser.add_argument("--generate-params", "-gp", help="Choose if paramters are generated or not", type=bool)
parser.add_argument("--nb-data", "-Nd", help="Number of data", type=int)
parser.add_argument("--nb-samples", "-Ns", help="Number of spins samples", type=int)

# read arguments from the command line
args = parser.parse_args()
