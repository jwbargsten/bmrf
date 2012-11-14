import string
import sys

filterfile=sys.argv[1]
infile=sys.argv[2]

def read_filter(file):
 f=open(file,'r')
 data=f.readlines()
 f.close()
 list=[]
 for i in data:
  s=string.split(i)
  p,q=s[0],s[1]
  if list.count(p)==0:
   list.append(p)
  if list.count(q)==0:
   list.append(q)
 return list

def read_data(file,filt):
 f=open(file,'r')
 data=f.readlines()
 f.close()
 for i in data:
  s=string.split(i)
  p=s[0]
  if filt.count(p)>0:
   print i[:-1]

flist=read_filter(filterfile)
read_data(infile,flist)
