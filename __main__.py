#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 01:18:34 2022

@author: dozeduck
"""
import getopt
import sys
from PDB_reader_writer import PDBfile

x = PDBfile()

args=sys.argv[1:]
try:
   opts, args = getopt.getopt(args,"h:i:o:",["help",
                                            "input_pdb=",
                                            "output_name="])
except getopt.GetoptError:
   print ('PDBreader.py -i <inputprotein> -o <output_name>')
   sys.exit(2)                                                                      # Exiting the program raises a SystemExit exception, 0 means normal exit, and others are abnormal exits.
 
for opt, arg in opts:                                                               # Generate several pairs of value, e.g: opr,arg = -i,PLDXPAL
   if opt == '-h':
      print ('PDBreader.py -i <inputprotein> -o <output_name>')
      sys.exit()
   elif opt in ("-i", "--input_pdb"):
      f = str(arg)
   elif opt in ("-o", "--outputfile"):
      out_put_name = str(arg)
      
x.PDBreader(f)
x.Define_Chains()
x.PDBwriter(out_put_name)