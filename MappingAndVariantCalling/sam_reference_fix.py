import sys
#import os
#rom stat import S_IS

def fix(line):
    if "@SQ\t" in line:
        split_line = line.split("\t")
        length = split_line[2]
        length = length.split(":")
        length[1] = int(length[1])
        length[1] -= 1
        length[1] = str(length[1])
        split_line[2] = ":".join(length)
        return "\t".join(split_line)+"\n"
    return line

if __name__=="__main__":
    for line in sys.stdin:
        print(fix(line), end="")
