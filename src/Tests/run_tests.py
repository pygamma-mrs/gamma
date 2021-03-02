#run_tests.py
from __future__ import division
from __future__ import print_function

import os
import glob
import sys
import subprocess
import numpy
import scipy
import scipy.io as sio
import scipy.linalg.interpolative as sli
import filecmp as fc
import time
from optparse import OptionParser



# According to the python 2.6.4 online documentation,
# The following names have currently been registered:
# 'posix', 'nt', 'mac', 'os2', 'ce', 'java', 'riscos'
os_id = os.name


def is_numeric(x):
    if(len(x) == 0):
        return False
    else:
        count = 0
        has_decimal = False
        for n in x:
            if count == 0 and (n=='+' or n=='-'):
                    continue
            elif n=='.':
                if has_decimal == False:
                    has_decimal = True
                    continue
                else:
                    return False
            elif n.isdigit():
                continue
            else:
                return False
            
    return True            
     
# define comparison function to be used by run_tests.py
def cmp_nn(file1, file2):

    mrv =  0.000002
    dncp = 0.0000001

    # Compare two files 
    # ignoring the newlines (lf vs cr-lf)

    f1 = open(file1)
    f2 = open(file2)

    L1 = f1.readline()
    L2 = f2.readline()

    result = True
    while L1 or L2:
        if L1 == None:
            L2 = L2.strip("\r\n\f")
            L2 = L2.strip()
            if len(L2) > 0:
                result == False
                break
            else:
                continue
        if L2 == None:
            L1 = L1.strip("\r\n\f")
            L1 = L1.strip()
            if len(L1) > 0:
                result == False
                break
            else:
                continue
        L1 = L1.rstrip("\r\n\f")
        L2 = L2.rstrip("\r\n\f")
        if L1 != L2:
            # Tokenize the string and compare them - element by element.
            tokens1 = L1.split()
            tokens2 = L2.split()
            if(len(tokens1) != len(tokens2)):
                result = False
                break
            
            for t1, t2 in zip(tokens1, tokens2):
                s1 = t1.strip()
                s2 = t2.strip()                

                if s1 == s2:
                    continue
                else:
                    # print s1
                    # print s2
                    if is_numeric(s1) and is_numeric(s2):
                        float_1 = float(s1)
                        float_2 = float(s2) 

                        if(compare_values(float_1, float_2, mrv, dncp) == False):
                            result = False
                            break
                    else:
                        result = False
                        break

        if(result == False):
            break
        L1 = f1.readline()
        L2 = f2.readline()   

    f1.close()
    f2.close()

    if is_verbose and result == False:
        print(file1 + ": '" + L1 + "'")
        print(file2 + ": '" + L2 + "'")

    return result


# define matlab comparison function to be used by run_tests.py
def cmp_mat(file1, file2):
    # Compare two matlab files 
    # ignoring the newlines (lf vs cr-lf)    

    contents1 = sio.loadmat(file1)
    contents2 = sio.loadmat(file2)

    # Warning: Cannot guarantee that the data will be the 0th element in this array!!

    keylist = list(contents1.keys())

    for key in keylist:
        if key.startswith("test_"):        
            
            data1 = contents1[ key ]
            data2 = contents2[ key ]

            #if scipy.rank(data1) != scipy.rank(data2):
            if sli.estimate_rank(data1, eps=1e-10) != sli.estimate_rank(data2, eps=1e-10):
                print("Data sets in the two files are not of the same rank")
                return False

            if data1.shape != data2.shape:
                print("Data sets in the two files do not have the same shape")
                return False

            #print "About to run compare_data()"
            return compare_data(data1, data2)

    print("No data to compare!")
    return True


def compare_data(data_1, data_2):

    max_relative_diff =  0.000002
    dontcare_point     = 0.0000001
    # See below for the analogous values for complex numbers
    # I.e.: component_relative_diff AND component_dontcare_point

    l = len(data_1.shape)

    if l==0:
        #print "l==0"
        return compare_values(data_1, data_2, max_relative_diff, dontcare_point)

    if l == 1:
        #print "l==1"
        return compare_lists(data_1, data_2, max_relative_diff, dontcare_point)        
    else:
        #print "l>1"
        i = 0
        for d1 in data_1:
            d2 = data_2[i]
            result = compare_data(d1, d2)
            i += 1
            if result == True:
                continue
            else:
                return False

    return True


def compare_lists(sources, targets, max_relative_difference, do_not_care_point):

    for source, target in zip(sources, targets):
        if compare_values(source, target, max_relative_difference, do_not_care_point, True) == True:
            continue
        else:
            return False

    #print "returning true"
    return True



def compare_values(test_value, expected_value, max_relative_difference, do_not_care_point, do_print=False):
   
    # Compares two numbers (usually floats) for fuzzy equality.

    tvt = type(test_value)

    if tvt == float or tvt == numpy.float32 or tvt == numpy.float64:
        return compare_floats(test_value, expected_value, max_relative_difference, do_not_care_point, do_print)
    elif tvt == complex or tvt == numpy.complex64 or tvt == numpy.complex128:
        return compare_complex(test_value, expected_value, max_relative_difference, do_not_care_point, do_print)
    else:
        if test_value != expected_value:
            if do_print == True:
                print("Other:")
                print("test_value = " + str(test_value) + " and expected_value = " + str(expected_value))
            return False
        else:
            return True
                
def compare_floats(test_value, expected_value, max_rel_diff, do_not_care_point, do_print = False):

    abs_difference = abs(test_value - expected_value)

    # See if the difference is out of tolerance, and also meaningful
    # i.e. near-zero values sometimes show up as noise, like 7.382239e-304)

    if abs_difference > do_not_care_point:
        if abs_difference > ( abs(expected_value) * max_rel_diff ):
            if do_print == True:
                print("float:")
                print("test_value = %.13f and expected_value = %.13f" % (test_value, expected_value))
            return False
    
    return True

def compare_complex(test_value, expected_value, max_rel_diff, dncp, do_print = False):

    abs_difference = abs(test_value - expected_value)

    tR = test_value.real
    eR = expected_value.real

    tI = test_value.imag
    eI = expected_value.imag

    if abs_difference > dncp:
        if abs_difference > ( abs(expected_value) * max_rel_diff ):
            if do_print == True:
                print("complex-A:")
                print("I values are %.29f, %.29f" % (tI, eI))
                print("R values are %.29f, %.29f" % (tR, eR))
                print("test_value = (%.11f + %.11fj)" % (tR, tI))
                print("expected_value = (%.11f + %.11fj)" % (eR, eI))
            return False

    component_relative_diff =  0.000002
    component_dontcare_point = 0.00000005

    compareR = compare_floats(tR, eR, component_relative_diff, component_dontcare_point)
    compareI = compare_floats(tI, eI, component_relative_diff, component_dontcare_point)

    if compareR == False or compareI == False:
        if do_print == True:
            print("complex-B:")
            print("compareR == " + str(compareR) + "  compareI == " + str(compareI))
            print("I values are %.29f, %.29f" % (tI, eI))
            print("R values are %.29f, %.29f" % (tR, eR))
            print("test_value = (%.11f + %.11fj)" % (tR, tI))
            print("expected_value = (%.11f + %.11fj)" % (eR, eI))
        return False
    else:
        return True

# parse the command line.
parser = OptionParser()
parser.add_option("-p", "--path", dest="path", action="store", default="",
                  help="path to executable test progams", metavar="PATH")
parser.add_option("-v", "--verbose", dest="verbose", 
                    action="store_true", default=False,
                    help="print status messages to stdout")

(options, args) = parser.parse_args()


filelist = glob.glob('*.suite')
total_failures = 0
total_comparisons = 0
is_verbose = options.verbose

print("")

for filename in filelist:
    bad_thing_happened = False
    current_suite = ""
    f = open(filename)

    execlist = []
    delfiles = []

    # read through and get a list of files to delete
    # and organize the executable lines into a shorter list.
    for line in f.readlines():

        line = line.strip()
                
        if len(line) > 0:
            if line[0] == '#':  # ignore lines that begin with a '#'
                continue
        else:
            continue  # ignore totally blank lines

        #split into command and detailed information.
        substrings = line.split(':')

        # validate that this data line has a valid format.
        if len(substrings) != 2:
            continue

        s = substrings[0].strip()

        # FIXME: 
        # We should add a check to make sure
        # the file does not try to compare files before
        # running the command, etc.

        # Is this a command we understand?
        if s == "compare" or s == "suite" or s == "command":

            execlist.append(line)

            # make a list of files to delete.
            if s == "compare":
                st = substrings[1].strip()
                filenames = st.split(",")
                if len(filenames) >= 2: # make sure there are 2 files to compare.
                    # then we should delete the first one.
                    delfiles.append( filenames[0].strip() )
        else:
            print("")
            print("*** Warning: Do not understand instruction, " + s + " !  ***")
            print("")         
            continue


    print("Deleting old output files: ")
    for dfile in delfiles:
        if dfile.startswith("golden") or dfile.startswith("./golden") or dfile.startswith(".\golden"):
            print("")
            print("    *** WARNING: Do not want to delete your golden files, " + dfile)
            print("    Double check how you wrote your suite, " + filename)
            print("")
            continue
        if os.path.exists(dfile) == True:
            print("    " + dfile)
            os.remove(dfile)

    print("")

    # Now execute all the executable lines...

    for line in execlist:

        # A little duplication here, but it's faster to code...for now.

        substrings = line.split(':')

        # double check that there are two items in this list.
        if len(substrings) != 2:
            continue

        s = substrings[0].strip()

        if s == "suite":
            current_suite = substrings[1].strip()
            print("Running Suite: " + current_suite)
        if s == "command":
            command = ""
            command += options.path
            command += substrings[1].lstrip()
            print("")
            print("running command: " + command)
            sys.stdout.flush()



            t1 = time.time()
            retcode = subprocess.call(command, shell=True)
            t2 = time.time()
            if retcode < 0:
                s1 = "Attempt to call " + substrings[1] + " was terminated by signal " + str(retcode)
                print(s1, file=sys.stderr)
            else:
                print('\'%s\' took %0.3f seconds' %(command, (t2-t1)))

        if s == "compare":
            total_comparisons += 1
            st = substrings[1]
            st = st.lstrip()
            st = st.rstrip()

            filenames = st.split(",")
            if len(filenames) < 2:
                print("Error while attempting to compare files")
                print("Filecount was too small")
                total_failures += 1
                continue

            file1 = filenames[0]
            file2 = filenames[1].lstrip()
            if is_verbose == True:
                print("Comparing: " + file1 + " " + file2)
                sys.stdout.flush()

            filebad = False
            fileout = ""
            if os.path.exists(file1) == False:
                fileout = file1
                filebad = True
            elif os.path.exists(file2) == False:
                fileout = file2
                filebad = True
            if filebad:
                sys.stdout.write(" *** TEST FAILED ***\n")
                print("")
                print("Warning: Could not find file, " + fileout + " required for comparison")
                print("")
                total_failures += 1
                continue

            if file1.endswith(".mat"):
                if file2.endswith(".mat"):
                    result = cmp_mat(file1, file2)
                else:
                    print("")
                    print("Both files need to be of type matlab (i.e. 'xxx.mat').")
                    print("")
                    result = False
            else:
                result = cmp_nn(file1, file2)

            if result == False:
                if is_verbose == True:
                    print("*** TEST FAILED ***")
                    print("")
                    sys.stdout.flush()
                total_failures += 1
                continue

    f.close()

    if bad_thing_happened == True:
        print("")
        print("*** Error while Running the \'" + current_suite + "\' Suite, Tests not complete!")

print("\nSummary:")
print("\tRan a total of " + str(total_comparisons) + " comparisons")
print("\tThere were " + str(total_failures) + " failure(s)")
print("")

