import sys

try:
    with open(sys.argv[1], 'r') as f:
        for lines in f.readlines():
            line = lines.split()
            try:
                if line[0] == 'Parameter:':
                    if line[1] == 'CHRG':
                        #print(line[-1])
                        #print(line[-1].split('"'))
                        print(abs(round(float(line[-1].split('"')[1]))))
            except IndexError:
                None
except FileNotFoundError:
    None

