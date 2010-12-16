#!/usr/bin/env python
from numpy import *
import sys
ri=map(eval,sys.argv[1:])
nr=len(ri)
print """
vertices=
{"""
print " { 0 , 0 },"
for r in ri[:-1]:
    print " { %f, 0 },"  % r
    print " { 0, %f }," % r
print " { %f, 0 }, " %  ri[-1]
print " { 0, %f } " %  ri[-1],
print"""
}

elements=
{"""
id=0
print " { 0, 1, 2, %d }," % id
for i in range(1,2*nr-3,2):
    id=id+1
    print " { %d, %d, %d, %d, %d }," % (i,i+2,i+3,i+1,id)
id=id+1
i=2*nr-3
print " { %d, %d, %d, %d, %d }" % (i,i+2,i+3,i+1,id),
print"""
}

boundaries=
{
"""
id=1
print " { 0, 1, 1},"
print " { 0, 2, 2},"
id=3
for i in range(1,2*nr-1,2):
    print " { %d, %d, %d }," % (i,i+2,id)
    id=id+1
    print " { %d, %d, %d }," % (i+1,i+3,id)
    id=id+1
i=2*nr-1
print " { %d, %d, %d }" % (i,i+1,id)

print """
}
curves=
{"""
for i in range(1,2*nr-1,2):
    print " { %d, %d, 90 }," % (i,i+1)
i=2*nr-1
print " { %d, %d, 90 }" % (i,i+1)
print "}"

