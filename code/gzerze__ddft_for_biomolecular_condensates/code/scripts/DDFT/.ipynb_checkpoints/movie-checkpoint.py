import argparse
from Protein_RPA.utils.graphics import *




def create_movie(args):
    output_dir=args.o
    input_dir=args.i
    save_movie(input_dir,output_dir)


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description='Take output filename to produce Movie')
     parser.add_argument('--i',help="Directory of PNG files", required = True);

     parser.add_argument('--o',help="Directory of Movie files", required = True);
     #parser.add_argument('--s',help="Name of state folder", required = True);
     args = parser.parse_args();

     create_movie(args);
