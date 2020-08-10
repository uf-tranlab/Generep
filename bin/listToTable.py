#!/usr/bin/env python

import sys

def main(l):
    ls = [p.strip() for p in l.strip("[]").split(",")]
    for i in range(0, len(ls), 3):
        sys.stdout.write("{}\t{}\t{}\n".format(ls[i], ls[i+1], ls[i+2]))

if __name__ == "__main__":
    main(sys.argv[1])
