import os
import subprocess
from subprocess import call
#os.system("ls -l")

#call(["ls", "-l"])
#abc = "leapr.cpp"
#call(["vim", abc])

def get_name_of_all_files():
    x = subprocess.check_output(['ls']).split()
    l = []
    for entry in x:
        entry = entry.decode('utf-8')
        l.append(entry) 
    return l

def get_all_dirs():
    l = get_name_of_all_files()
    dirs = []
    for entry in l:
        if '.' not in entry:
            dirs.append(entry)
    return dirs
    
def get_all_test_files():
    l = get_name_of_all_files()
    tests = []
    for entry in l:
        if '.test.cpp' in entry:
            tests.append(entry)
    return tests

def run_test(test):
    my_str = "g++ -std=c++14 " + test 
    my_str = "g++ " + test 
    os.system(my_str)
    try: 
        output = subprocess.check_output(["./a.out"])
        if 'All tests passed' in output.decode('utf-8'):
            return True
    except:
        return False
    
    


def testing_func(failed=[]):
    dirs = get_all_dirs()
    tests = get_all_test_files()
    for test in tests:
        print(" - ",test )
        if not run_test(test):
            print("Oh no! "+test+" failed :( ")
            failed.append(test)
    for directory in dirs:
        os.chdir("./"+directory)
        print('Checking',directory)
        failed = testing_func(failed)
        os.chdir("../")
    return failed

output = testing_func()
if len(output) == 0:
    print( "\n:)\n" )
else:
    print("\nThe following tests failed: ")
    print(output)
    print( ":(" )






