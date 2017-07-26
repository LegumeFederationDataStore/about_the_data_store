#!/usr/local/bin/python

#UNDER DEVELOPMENT  

'''
Extracts from a structures Excel file, information about an expression dataset
and writes to markdown file.

Sudhansu Dash: May 2017

'''

#TO DO:
# output file name as arg??
#USAGE:  

#-----------------------------------------------------------------------------

import sys
import re
#For reading Excel file
import openpyxl as opx  


#----------------------

#Excel File to read
xlfile = sys.argv[1]

# *****  PROVIDE:  OUTPUT FILE NAME   ******
outputfile = open('cicar-SRP017394-atlas-on-ICC4958-dataset.md', 'w')
#----------------

print "Excel file: ",xlfile

#Read xlfile
wb = opx.load_workbook(xlfile, read_only=True)
print wb.get_sheet_names()  #prints all the worksheets in workbook

# Get the worksheets from workbook
ws_dataset = wb.get_sheet_by_name('dataset')
ws_datasource = wb.get_sheet_by_name('datasource')
ws_sample = wb.get_sheet_by_name('sample')
ws_method = wb.get_sheet_by_name('method')
ws_expdesign = wb.get_sheet_by_name('expdesign')


# Dataset:
#---------

dict_dataset = {}  ## to append to this dict of attrib=value

print "\n" + "##  " + "DATASET" + "\n"
outputfile.write("\n" + "##  " + "DATASET" + "\n")

for row in tuple(ws_dataset.rows):   # creates a tuple of all rows in ws_dataset
    
    #print row[0].value #debug
    if (row[0].value and re.match('^#', row[0].value)):   # skip commented line, 1st col/cel is '#'
          ##row[0].value must be a string for re.match; Caution about 'NoneObject'
        continue
    
    else:
        k = row[1].value
        v = row[2].value
        
        if (k):
            #k = k.rstrip()  #key k in 2nd cell
            print '####  ' + k
            outputfile.write('####  ' + k + "\n")
        
        if (v):
            #v = v.rstrip()  #value v in 3rd cell
            print v
            outputfile.write(v + "\n")
      
    print "\n"  # \n after every sheet-row (item)
    dict_dataset[k] = v
#end for
#i


# Datasource:
#------------

dict_dsource = {}  ## to append to this dict of attrib=value

print "\n" + "##  " + "DATASOURCE" + "\n"
outputfile.write("\n" + "##  " + "DATASOURCE" + "\n")

for row in tuple(ws_datasource.rows):   # creates a tuple of all rows in ws_dataset
    if (row[0].value and re.match('^#', row[0].value)):   # skip commented line, 1st col/cel is '#'
        continue
    else:
        k = row[1].value
        v = row[2].value
        
        if (k):
            #k = k.rstrip()  #key k in 2nd cell
            print '####  ' + k
            outputfile.write('####  ' + k + "\n")
        
        if (v):
            #v = v.rstrip()  #value v in 3rd cell
            print v
            outputfile.write(v + "\n")
            
    print "\n"  # \n after every sheet-row (item)
    dict_dsource[k] = v
#end for
#


print "\n" + "##  " + "METHOD" + "\n"
outputfile.write("\n" + "##  " + "METHOD" + "\n")


# Method:
#--------

dict_method = {}

for row in tuple(ws_method.rows):   # creates a tuple of all rows in ws_dataset
    if (row[0].value and re.match('^#', row[0].value)):   # skip commented line, 1st col/cel is '#'
        continue
    else:
        k = row[1].value
        v = row[2].value
        
        if (k):
            #k = k.rstrip()  #key k in 2nd cell
            print '####  ' + k
            outputfile.write('####  ' + k + "\n")
        
        if (v):
            #v = v.rstrip()  #value v in 3rd cell
            print v
            outputfile.write(v + "\n")
            
    print "\n"  # \n after every sheet-row (item)
    dict_method[k] = v
#end for
#



##SAMPLES INTO A LIST-OF-DICTIONARY
#----------------------------------


## cells in row starting with 'sample_name'
attribute_names = []  #a list, to append to
for row in tuple(ws_sample.rows):   # creates a tuple of all rows in ws_sample
    if (row[0].value == 'sample_name'):   # r=row, if r 1st col is ?sample_name?
        for cell in row:    #start loop through that row
            #print cell.value  # print cell value
            attribute_names.append(cell.value)


attribute_names_str = "\t".join(attribute_names)    #list to string with '\t' as sep
print "Sample attribute names: \n", attribute_names


data_list_dict = []  ## To append to a list of dicts(each row is a dict)

for row in tuple(ws_sample.rows):   # creates a tuple of all rows in ws_sample
    if (row[0].value and re.match('^#', row[0].value)):   # skip commented line, 1st col/cel is '#'
        continue
    
    rowvalues = []    #list
    
    for cell in row:    #start loop through that row
        rowvalues.append(cell.value)
    #End- inner for
    
    dict_col_value = dict(zip(attribute_names, rowvalues))  ##Creates a dict col-names and this row's values
    data_list_dict.append(dict_col_value)  ## Appends this dict to a list (A list of dict at the end)
#End-outer for    

print "\n" + "##  " + "SAMPLES" + "\n"
outputfile.write("\n" + "##  " + "SAMPLES" + "\n")

for samp in data_list_dict:
    row = samp['sample_name'] + "\t" + samp['sample_uniquename'] + "\t" + samp['description'] \
    + "\t" + samp['treatment'] + "\t" + samp['tissue'] + "\t" + samp['dev_stage']  \
    + "\t" + str(samp['age']) + "\t" \
    + samp['organism'] + "\t" + samp['infraspecies'] + "\t" + samp['cultivar'] \
    + "\t" + samp['sra_run'] + "\t" + samp['biosample_accession'] + "\t" + samp['sra_accession'] \
    + "\t" + samp['bioproject_accession'] + "\t" + samp['sra_study'] 
         # + samp[''] + samp['']
    print row + "\n"
    outputfile.write(row + "\n")



outputfile.write("\n"+"\n"+"\n")

outputfile.close()


##-------------scratch pad------------
'''
for row in tuple(ws_dataset.rows):   # creates a tuple of all rows in ws_dataset
    if (re.match('^#', row[1].value)):   # skip commented line, 1st col/cel is '#'
        #print re.match('^#', row[0].value)
        #print row[0].value
        continue
    else:
        #print row[0].value,":",row[1].value,":",row[2].value,"\n"
        print "cell2: ", row[2].value
#
#
'''
