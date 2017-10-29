#! /usr/bin/env python
import csv
import sys

csv1 = csv.reader(open(sys.argv[1], 'rt'))
csv2 = csv.reader(open(sys.argv[2], 'rt'))

d1 = {}
for row in csv1:
    d1[row[0]] = row[1:]

d2 = {}
for row in csv2:
    d2[row[0]] = row[1:]

kk = set(d1.keys()).intersection(d2.keys())

for k in sorted(kk):
    if d1[k] != d2[k]:
        print(k)
        print('\t', d1[k])
        print('\t', d2[k])
        print('--')
