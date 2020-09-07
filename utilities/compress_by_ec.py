#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

import re
import sys, os
from optparse import OptionParser, OptionGroup
from glob import glob

script_name = sys.argv[0]
usage= script_name + """ --pin <input pf file> --pout <output pf file>  \n --rin <reduced input file> 
                         --rout <reduced output file>
                          This script takes the 0.pf and reduced.txt files to create a new set of 
                          0.pf and  reduced.txt
                     """
parser = OptionParser(usage)

parser.add_option( "--pin", dest="pfin", default =None, 
                  help='the input 0.pf file')

parser.add_option( "--pout", dest="pfout", default = None, 
                  help='the output .pf file')

parser.add_option( "--rin", dest="redin", default = None, 
                  help='reduced input file')

parser.add_option( "--rout", dest="redout", default = None, 
                  help='reduced out file')


parser.add_option( "--level", dest="level", default = '2',  choices=['0', '1', '2'],
                  help='two orfs are equivalent  \n' + \
                        '0: if ECs the similar, \n ' +  \
                        '1: if product/function strings are the similar\n' 
                        '2: if product strings prefix similar\n' 
                  )

FIELDS  = ['ID', 'NAME', 'STARTBASE', 'ENDBASE', 'FUNCTION', 'PRODUCT-TYPE', 'EC']
def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)


def check_arguments(opts, args):
    if opts.pfin == None:
         return False
    if opts.redin == None:
         return False

    return True

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


def  fix_pgdb_input_files(pgdb_folder, pgdbs = []):
     pgdb_list = glob(pgdb_folder + '*cyc/*/input/organism.dat')     

     for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
def read_pf_file(filel):
     patternID =  re.compile(r'^ID\t.*')

     annotations = []
     ecs= {}
     try:
         orgfile = open(filel,'r')
     except IOError:
         print("ERROR : Cannot open organism file" + str(filel))
         return 
     lines = orgfile.readlines()
     orfids = []
     annotation={}
     products ={}

          
     orfinfo = {}
     for line in lines:
         fields = line.strip().split('\t')

         if len(fields)==1 and re.search('^//', line.strip()):
             if 'ID' in orfinfo \
                 and 'NAME' in orfinfo \
                 and 'STARTBASE' in orfinfo \
                 and 'ENDBASE' in orfinfo \
                 and 'FUNCTION' in orfinfo \
                 and 'PRODUCT-TYPE' in orfinfo:

                 orfid = orfinfo['ID']
                 products[orfid] = re.sub('  ', ' ', orfinfo['FUNCTION'].lower())

                 orfids.append(orfid)
                 annotation[orfid] =[]
                 for field in FIELDS:
                     if field in orfinfo:
                         if fields[0]=="FUNCTION":
                             products[orfid] = re.sub('  ', ' ', fields[1].lower())
                         annotation[orfid].append(field + '\t' + orfinfo[field])
             orfinfo = {}
                     
         if len(fields)==2:
             orfinfo[fields[0]] = fields[1]

     orgfile.close()


     orfids.sort(key = lambda x: products[x] )

     #for orfid in orfids:
     #      print(annotation[orfid], products[orfid])

     return orfids, products, annotation, ecs
     

def read_reduced_file(filel):

     try:
         orgfile = open(filel,'r')
     except IOError:
         print("ERROR : Cannot open organism file" + str(filel))
         return 
     lines = orgfile.readlines()
     equivalent={}

     for line in lines:
         fields = line.strip().split('\t')
         if len(fields)!=2:
            continue
         equivalent[fields[0]] =fields[1]

     orgfile.close()

     return equivalent


def write_new_file(lines, output_file):
    
    print("Fixing file " + output_file )
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print("ERROR :Cannot open output file "  + output_file)
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()


def  modify_files(filerenames, source, target):

    for file, newname in filerenames.iteritems():
        print("<<<<" , file)
        with open(file) as f:
            lines = f.readlines()
            outputfile = open(file + '.tmp','w')
            for line in lines:
                newline = line.replace(source, target)
                fprintf(outputfile, "%s", newline)
            outputfile.close();
            os.remove(file)
            os.rename(file+'.tmp', newname)


def get_files_to_modify(folder):
   all_files = []
   for root, directories, filenames in os.walk( folder +'/'):
     for filename in filenames: 
        all_files.append(os.path.join(root,filename))

   return all_files
   

def  get_new_file_names(files, source, target):
     sourceRe = re.compile(source)
     filerename ={}

     for file in files:
         if sourceRe.search(file):
            newname = file.replace(source+'.', target + '.', 2)
            filerename[file] = newname 
         else:
            filerename[file] = file 
     return filerename


def cleanFunctionStrings(_products):
     funarray = []
     for orf, prod in _products.items():
        #prod = re.sub('-', ' ', prod)
        prod = re.sub(',', '', prod)
        prod = re.sub('\/', ' ', prod)
        prod = re.sub(':', ' ', prod)
        ncharwords=0
        
        for word in prod.split():
            if len(word)==1:    
               ncharwords += 1

        if ncharwords > 1 or len(prod.split()) > 5  :
            prod = "hypothetical protein"
        funarray.append( (orf, prod) )

     funarray.sort(key = lambda x: x[1] )
     
     prevfun = None
     products = {}
     for orf, fun in funarray:
       if prevfun != None and  fun.startswith(prevfun):
          products[orf] = prevfun
       else:
          products[orf] = fun
          prevfun = fun

       #print(fun, '\t\t', prevfun)
     return products 

# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print(usage)
       sys.exit(0)

    orfids, _products, annotations, ecs = read_pf_file(opts.pfin)

    # read the reduced file
    equivalent  = read_reduced_file(opts.redin)

    # orf-> function dictionary
    products = cleanFunctionStrings(_products)

    seenecs ={}
    seenproducts ={}
    neworfs =[]
    prodstring = None 
    for orfid in orfids:

       if opts.level == '0':
  
          if orfid not in ecs:
             continue
           
          if ecs[orfid] in seenecs: 
             if not orfid in equivalent:
                equivalent[orfid] = seenecs[ecs[orfid]]
          else:
              neworfs.append(orfid)
              seenecs[ecs[orfid]] = orfid
 

       if opts.level == '1':
          if orfid not in products:
             continue

          if products[orfid] in seenproducts:
              if not orfid in equivalent:
                 equivalent[orfid] = seenproducts[products[orfid]]
                 if orfid in ["O_83746_0", "O_63745_0"]:
                    print("not equivalent 1", orfid, products[orfid])
          else:
               seenproducts[products[orfid]] = orfid
               neworfs.append(orfid)
               if orfid in ["O_83746_0", "O_63745_0"]:
                 print("equivalent 1", orfid, _products[orfid])

       if opts.level == '2':
           if prodstring == None or not  products[orfid].startswith(prodstring):
              prodstring = products[orfid]
              prodorfid = orfid
              neworfs.append(orfid)
           else:
              equivalent[orfid] = prodorfid


    write_new_pf_file(neworfs, annotations,  opts.pfout)

#   print  len(orfids), len(neworfs), len(equivalent)
    print("Original # ORFs   : {}".format(len(orfids)))
    print("Compressed # ORFs : {}".format(len(neworfs)))
    print("Equiv # ORFS      : {}".format(len(equivalent)))
    print("Annotations #     : {}".format(len(annotations)))

    write_new_reduced_file(equivalent, opts.redout)


def write_new_reduced_file(equivalent, filename):
    outputfile = open(filename,'w')
    for orfid in equivalent: 
       if orfid !=  equivalent[orfid]:
          fprintf(outputfile, "%s\t%s\n",orfid, equivalent[orfid])
    outputfile.close()

def write_new_pf_file(neworfs, annotations, filename):
    with open(filename,'w') as outputfile:
        for orfid in neworfs: 
           for line in annotations[orfid]:
               fprintf(outputfile, "%s\n",line)
           fprintf(outputfile, "//\n")
        outputfile.close()

    #print len(equivalent.keys()), len(orfids), len(ecs)
    #for ec in ecs:
    #   print ecs[ec]


# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

