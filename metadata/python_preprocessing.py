# Script to preprocess metadata for scikit learn
import csv
import sys

sub=583



clin_out_fname='pat_FR_'+str(sub)+'_clinical_szrs_detailed.csv'
subclin_out_fname='pat_FR_'+str(sub)+'_subclinical_szrs_detailed.csv'

# generate mat filenames for Python script to load
all_files = range(1,150)
clin_files = []
subclin_files = []


f = open(clin_out_fname, 'rt')
try:
    reader = csv.reader(f)
    for row in reader:
        clin_files.append(int(row[0]))
finally:
    f.close()

f = open(subclin_out_fname, 'rt')
try:
    reader = csv.reader(f)
    for row in reader:
        subclin_files.append(int(row[0]))

finally:
    f.close()


ictal_files = clin_files + subclin_files 

print "ALL SEIZURE FILES: "
for fnum in sorted(set(ictal_files)):
    print '"EDMSE_ideal_pat_FR_' + str(sub) + "_" + str(fnum) + '.mat",'


clin_files = set(clin_files)
subclin_files  = set(subclin_files)



print "CLINICAL SEIZURE FILES: "
for fnum in sorted(clin_files):
    print '"EDMSE_ideal_pat_FR_' + str(sub) + "_" + str(fnum) + '.mat",'

print "SUBCLINICAL SEIZURE FILES: "
for fnum in sorted(subclin_files):
    print '"EDMSE_ideal_pat_FR_' + str(sub) + "_" + str(fnum) + '.mat",'    


unique_ictal_files = set(ictal_files)
dupl_ictal_files = set([x for x in ictal_files if ictal_files.count(x) > 1])
print "Multi-ictal Files: ", len(dupl_ictal_files), dupl_ictal_files
print "Unique Ictal Files: ", len(unique_ictal_files), unique_ictal_files


print "INTERICTAL FILES: "
interictal_files = set(all_files)- set(clin_files) #- set(subclin_files)
for fnum in sorted(interictal_files):
    print '"EDMSE_ideal_pat_FR_' + str(sub) + "_" + str(fnum) + '.mat",'



print set(ictal_files)
print set(all_files)

print (set(all_files) - set(ictal_files))
