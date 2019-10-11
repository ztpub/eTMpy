import os

import sys
fa_name = sys.argv[1]
mir_name = sys.argv[2]
P_score = sys.argv[3]
print(fa_name)
print(mir_name)
P_score = float(P_score)
print(P_score)

if P_score<0:
    P_score = 2
    print "P_score should be non-negative value, and this parameter is reset as default 2"


f = open(fa_name)
lines = f.readlines()
f.close()

Hash = {}
title = ''
start = 0
record = ''
for i in range(0,len(lines)):
    line = lines[i].replace("\n","").replace("\r","").replace("\S","").replace("U","T")
    if line.startswith(">"):
        if len(title)>0 and len(record)>0:
            Hash[title] = record
            record = ''
        
        title = line
    else:
        record = record + line

    if i==len(lines)-1:
        Hash[title] = record

sf = open("my_py_positive.fasta","w")
for title in Hash.keys():
    line = Hash[title]
    #print line
    if len(line)>200:
        plus = ''
        for k in range(0,len(line)):
            if line[k]=='A':
                plus = plus + 'T'
            if line[k]=='T':
                plus = plus + 'A'
            if line[k]=='C':
                plus = plus + 'G'
            if line[k]=='G':
                plus = plus + 'C'
        sf.write(title+"\n")
        sf.write(plus+"\n")
sf.close()


sf = open("my_py_negative.fasta","w")
for title in Hash.keys():
    line = Hash[title]
    #print line
    if len(line)>200:
        minus = line[::-1]
        sf.write(title+"\n")
        sf.write(minus+"\n")
sf.close()

########################################################        

def buildseq_plus(sequence):
   newseq = {};

   L = len(sequence)
   sequence = sequence[::-1]
   newseq[0] = sequence[0:(L-9)] + "---" + sequence[(L-9):L]
   newseq[1] = sequence[0:(L-10)] + "---" + sequence[(L-10):L]
   #newseq[0] = sequence[0:(L-11)] + "---" + sequence[(L-11):L]
   newseq[2] = sequence[0:(L-10)] + "--" + sequence[L-10] + "-" + sequence[(L-9):L]
   newseq[3] = sequence[0:(L-10)] + "-" + sequence[L-10] + "--" + sequence[(L-9):L]

   newseq[4] = sequence[0:(L-11)] + "--" + sequence[(L-11):(L-9)] + "-" + sequence[(L-9):L]
   newseq[5] = sequence[0:(L-11)] + "-" + sequence[(L-11):(L-9)] + "--" + sequence[(L-9):L]
   newseq[6] = sequence[0:(L-11)] + "--" + sequence[(L-11)] + "-" + sequence[(L-10):L]
   newseq[7] = sequence[0:(L-11)] + "-" + sequence[(L-11)] + "--" + sequence[(L-10):L]
   newseq[8] = sequence[0:(L-11)] + "-" + sequence[(L-11)] + "-" + sequence[(L-10)] + "-"+ sequence[(L-9):L]

   #for k in range(0,9):
   #    tmp = newseq[k]
   #    newseq[k] = tmp[::-1]
   #print "----"

   return newseq

#############################################################
def buildseq_minus(sequence):
   newseq = {};

   L = len(sequence)
   sequence = sequence[::-1]
   newseq[0] = sequence[0:(L-9)] + "---" + sequence[(L-9):L]
   newseq[1] = sequence[0:(L-10)] + "---" + sequence[(L-10):L]
   newseq[2] = sequence[0:(L-11)] + "---" + sequence[(L-11):L]  ################################### minus version
   newseq[3] = sequence[0:(L-10)] + "--" + sequence[L-10] + "-" + sequence[(L-9):L]
   newseq[4] = sequence[0:(L-10)] + "-" + sequence[L-10] + "--" + sequence[(L-9):L]

   newseq[5] = sequence[0:(L-11)] + "--" + sequence[(L-11):(L-9)] + "-" + sequence[(L-9):L]
   newseq[6] = sequence[0:(L-11)] + "-" + sequence[(L-11):(L-9)] + "--" + sequence[(L-9):L]
   newseq[7] = sequence[0:(L-11)] + "--" + sequence[(L-11)] + "-" + sequence[(L-10):L]
   newseq[8] = sequence[0:(L-11)] + "-" + sequence[(L-11)] + "--" + sequence[(L-10):L]
   newseq[9] = sequence[0:(L-11)] + "-" + sequence[(L-11)] + "-" + sequence[(L-10)] + "-"+ sequence[(L-9):L]

   #for k in range(0,9):
   #    tmp = newseq[k]
   #    newseq[k] = tmp[::-1]
   #print "----"

   return newseq

#################################################################

def sub(string,p,c):
        new = []
        for s in string:
            new.append(s)
        new[p] = c
        return ''.join(new)
    

def buildGUpairs(text):
   mat = {};

   L = len(text)
   mat[0] = text[::-1]
   #print(text+"\n"+mat[0]+"\n")

   i = 1;
   pos = -1;
   while(pos>=-1):
       tt = text
       pos = tt.find('G', pos+1)
       if pos==-1:
           break
       Index = pos
       tt = sub(tt, pos, "A")
       mat[str(i)+" G"] = tt[::-1]

       if(mat[str(i)+" G"]):
           Gnext = pos
           Gtext = tt[:]
           while(Gnext>=-1):
               Gtext = tt
               Gnext = Gtext.find('G', Gnext+1)
               if Gnext==-1:
                   break
               Gtext = sub(Gtext, Gnext, "A")
               mat[str(i)+" "+ str(Gnext) + " 2G"] = Gtext[::-1]

       if(mat[str(i)+" G"]):
           Gnext = -1
           Gtext = tt[:]
           while(Gnext>=-1):
               Gtext = tt
               Gnext = Gtext.find('T',Gnext+1)
               if Gnext==-1:
                   break
               Gtext = sub(Gtext, Gnext, "C")
               mat[str(i)+" "+ str(Gnext) + " GT"] = Gtext[::-1]

       i = i + 1
       #print i

    
   i = 1
   pos = -1
   while(pos>=-1):
       tt = text
       pos = tt.find('T',pos+1)
       if pos==-1:
           break
       Index = pos
       tt = sub(tt, pos, "C")
       mat[str(i)+" T"] = tt[::-1]

       if(mat[str(i)+" T"]):
           Gnext = pos
           Gtext = tt[:]
           while(Gnext>-1):
               Gtext = tt
               Gnext = Gtext.find('T',Gnext+1)
               if Gnext==-1:
                   break
               Gtext = sub(Gtext, Gnext, "C")
               mat[str(i)+" "+ str(Gnext) + " 2T"] = Gtext[::-1]


       i = i + 1

   #print(mat)
   return mat

###########################################################################
###########################################################################
sf = open("eTM_py_list_all.txt", "w")
sf1 = open("eTM_py_list.txt", "w")

f = open("my_py_positive.fasta")
lines = f.readlines()
f.close()

EGC = {}
title = ''
start = 0
record = ''
for i in range(0,len(lines)):
    line = lines[i].replace("\n","").replace("\r","").replace("\S","").replace("U","T")
    if line.startswith(">"):
        if len(title)>0 and len(record)>0:
            EGC[title] = record
            record = ''
        
        title = line
    else:
        record = record + line

    if i==len(lines)-1:
        EGC[title] = record


        
f = open(mir_name)
lines = f.readlines()
f.close()

MIR = {}
title = ''
start = 0
record = ''
for i in range(0,len(lines)):
    line = lines[i].replace("\n","").replace("\r","").replace("\S","").replace("U","T")
    if line.startswith(">"):
        if len(title)>0 and len(record)>0:
            MIR[title] = record.replace("U","T")
            record = ''
        
        title = line
    else:
        record = record + line

    if i==len(lines)-1:
        MIR[title] = record.replace("U","T")



for e in MIR.keys():
    #print e, MIR[e]
    #print("@@@@@@@@@@@@@@@@@"+e)
    
    
    myseq = buildseq_plus(MIR[e])

    #print(myseq)
    
    match = MIR[e]
    match = match[1:8]
    #print match
    #print "---"
    match = buildGUpairs(match)

    #for xe in match.keys():
    #    print xe, match[xe]

    #ontinue

    L = len(MIR[e])
    L = L + 3
    for ee1 in match.keys():
        for ee2 in EGC.keys():
            #print EGC[ee2]
            #print match[ee1]
            pos = -1
            while pos >= -1:
               #print EGC[ee2]
               #print match[ee1]
               
               pos = EGC[ee2].find(match[ee1],pos+1)
               #print pos
               if pos+7 > len(EGC[ee2]):
                   break
               if pos==-1:
                   break
               index_match = pos - L + 8
               if index_match<0:
                   break

               mm = EGC[ee2]
               mm = mm[index_match:(index_match+L)]
              
               #print "??????????"
               #print match[ee1]
               #print mm

               ##
               for se in myseq.keys():
                   seq = myseq[se]
                   #print se
                   #print seq
                   #print mm
                   
                   sa = 0
                   saa = 0
                   mark = 0

                   while sa <= len(mm):
                       mn1 = mm[sa:(sa+1)]
                       mn2 = seq[sa:(sa+1)]
                       #print mn1, mn2
                       if not mn1==mn2:
                           saa = saa + 1
                       sa = sa + 1
                   #print "saa=", saa
                   if saa>=3 and saa<=6:
                       #print "???",seq

                       mark = 1
                       mm1 = mm

                       plus = ''
                       for k in range(0,len(mm1)):
                           if mm1[k]=='A':
                               plus = plus + 'T'
                           if mm1[k]=='T':
                               plus = plus + 'A'
                           if mm1[k]=='C':
                               plus = plus + 'G'
                           if mm1[k]=='G':
                               plus = plus + 'C'
                       mm1 = plus

                       start = index_match + 1
                       end = index_match + L

	               #@ss=split(/\s/,$key);
	               #@newtit=split(/\-/,$ss[1]);
	               #$s1=$newtit[0]+$start-1;
	               #$e1=$newtit[0]+$end-1;

                       s1 = start - 1
                       e1 = end - 1


                       mirs = seq[::-1]
                       tars = mm1[::-1]

                       penalty_score = 0

                       #print mirs[9:12]
                       if mirs[9]=='-' and mirs[10]=='-' and mirs[11]=='-':
                           penalty_score = penalty_score + 0
                       else:
                           if mirs[9]=='-' and mirs[10]=='-':
                               penalty_score = penalty_score + 0
                           else:
                               if mirs[10]=='-' and mirs[11]=='-':
                                   penalty_score = penalty_score + 0
                               else:
                                   if mirs[9]=='-':
                                       penalty_score = penalty_score + 0.5
                                   if mirs[10]=='-':
                                       penalty_score = penalty_score + 0.5
                                   if mirs[11]=='-':
                                       penalty_score = penalty_score + 0.5

                       for pk in range(0,len(mirs)):
                           if pk==9 or pk==10 or pk==11:
                               continue
                           if mirs[pk]=='-':
                               penalty_score = penalty_score + 0.5

                       for pk in range(0,len(mirs)):
                           if pk==9 or pk==10 or pk==11:
                               continue
                           if mirs[pk]=='T' and tars[pk]=='G':
                               penalty_score = penalty_score + 0.5
                           if mirs[pk]=='G' and tars[pk]=='T':
                               penalty_score = penalty_score + 0.5

                           if mirs[pk]=='A' and tars[pk]=='C':
                               penalty_score = penalty_score + 1
                           if mirs[pk]=='C' and tars[pk]=='A':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='A' and tars[pk]=='G':
                               penalty_score = penalty_score + 1
                           if mirs[pk]=='G' and tars[pk]=='A':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='T' and tars[pk]=='C':
                               penalty_score = penalty_score + 1
                           if mirs[pk]=='C' and tars[pk]=='T':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='A' and tars[pk]=='A':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='T' and tars[pk]=='T':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='C' and tars[pk]=='C':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='G' and tars[pk]=='G':
                               penalty_score = penalty_score + 1

                               
                       #print(ee2+str(s1)+'..'+str(e1)+"(+)"+mm1+e+seq)
                       #print ee1
                       sf.write(ee2+"\t"+str(s1)+'..'+str(e1)+"(+)"+"\t5'"+mm1+"3'\t"+e+"\t3'"+seq+"5'\n")
                       if penalty_score<=P_score:
                           sf1.write(ee2+"\t"+str(s1)+'..'+str(e1)+"(+)"+"\t5'"+mm1+"3'\t"+e+"\t3'"+seq+"5'\n")
                           
                   if mark==1:
                      break

f = open("my_py_negative.fasta")
lines = f.readlines()
f.close()

EGC = {}
title = ''
start = 0
record = ''
for i in range(0,len(lines)):
    line = lines[i].replace("\n","").replace("\r","").replace("\S","").replace("U","T")
    if line.startswith(">"):
        if len(title)>0 and len(record)>0:
            EGC[title] = record
            record = ''
        
        title = line
    else:
        record = record + line

    if i==len(lines)-1:
        EGC[title] = record


        
f = open(mir_name)
lines = f.readlines()
f.close()

MIR = {}
title = ''
start = 0
record = ''
for i in range(0,len(lines)):
    line = lines[i].replace("\n","").replace("\r","").replace("\S","").replace("U","T")
    if line.startswith(">"):
        if len(title)>0 and len(record)>0:
            MIR[title] = record.replace("U","T")
            record = ''
        
        title = line
    else:
        record = record + line

    if i==len(lines)-1:
        MIR[title] = record.replace("U","T")
        
for e in MIR.keys():
    #print e, MIR[e]
    #print("@@@@@@@@@@@@@@@@@"+e)
    
    
    myseq = buildseq_minus(MIR[e])

    #print(myseq)
    
    match = MIR[e]
    match = match[1:8]
    #print match
    #print "---"
    match = buildGUpairs(match)

    #for xe in match.keys():
    #    print xe, match[xe]

    #ontinue

    L = len(MIR[e])
    L = L + 3
    for ee1 in match.keys():
        for ee2 in EGC.keys():
            #print ee2
            #print EGC[ee2]
            #print match[ee1]
            pos = -1
            while pos >= -1:
               #print EGC[ee2]
               #print match[ee1]
               
               pos = EGC[ee2].find(match[ee1],pos+1)
               #print pos
               if pos+7 > len(EGC[ee2]):
                   break
               if pos==-1:
                   break
               index_match = pos - L + 8
               if index_match<0:
                   break

               mm = EGC[ee2]
               mm = mm[index_match:(index_match+L)]
              
               #print "??????????"
               #print match[ee1]
               #print mm

               ##
               for se in myseq.keys():
                   seq = myseq[se]
                   #print se
                   #print seq
                   #print mm
                   
                   sa = 0
                   saa = 0
                   mark = 0

                   while sa <= len(mm):
                       mn1 = mm[sa:(sa+1)]
                       mn2 = seq[sa:(sa+1)]
                       #print mn1, mn2
                       if not mn1==mn2:
                           saa = saa + 1
                       sa = sa + 1
                   #print "saa=", saa
                   if saa>=3 and saa<=6:
                       #print "???",seq

                       mark = 1
                       mm1 = mm

                       plus = ''
                       for k in range(0,len(mm1)):
                           if mm1[k]=='A':
                               plus = plus + 'T'
                           if mm1[k]=='T':
                               plus = plus + 'A'
                           if mm1[k]=='C':
                               plus = plus + 'G'
                           if mm1[k]=='G':
                               plus = plus + 'C'
                       mm1 = plus

                       
                       start = index_match + 1
                       end = index_match + L

                       ############################################################################ minus version
                       LL = len(EGC[ee2])
                       new_start = LL - end
                       new_end = LL - start
                       s1 = new_start
                       e1 = new_end
                       
	               #@ss=split(/\s/,$key);
	               #@newtit=split(/\-/,$ss[1]);
	               #$s1=$newtit[0]+$start-1;
	               #$e1=$newtit[0]+$end-1;

                       mirs = seq[::-1]
                       tars = mm1[::-1]

                       penalty_score = 0

                       #print mirs[9:12]
                       if mirs[9]=='-' and mirs[10]=='-' and mirs[11]=='-':
                           penalty_score = penalty_score + 0
                       else:
                           if mirs[9]=='-' and mirs[10]=='-':
                               penalty_score = penalty_score + 0
                           else:
                               if mirs[10]=='-' and mirs[11]=='-':
                                   penalty_score = penalty_score + 0
                               else:
                                   if mirs[9]=='-':
                                       penalty_score = penalty_score + 0.5
                                   if mirs[10]=='-':
                                       penalty_score = penalty_score + 0.5
                                   if mirs[11]=='-':
                                       penalty_score = penalty_score + 0.5

                       for pk in range(0,len(mirs)):
                           if pk==9 or pk==10 or pk==11:
                               continue
                           if mirs[pk]=='-':
                               penalty_score = penalty_score + 0.5

                       for pk in range(0,len(mirs)):
                           if pk==9 or pk==10 or pk==11:
                               continue
                           if mirs[pk]=='T' and tars[pk]=='G':
                               penalty_score = penalty_score + 0.5
                           if mirs[pk]=='G' and tars[pk]=='T':
                               penalty_score = penalty_score + 0.5

                           if mirs[pk]=='A' and tars[pk]=='C':
                               penalty_score = penalty_score + 1
                           if mirs[pk]=='C' and tars[pk]=='A':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='A' and tars[pk]=='G':
                               penalty_score = penalty_score + 1
                           if mirs[pk]=='G' and tars[pk]=='A':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='T' and tars[pk]=='C':
                               penalty_score = penalty_score + 1
                           if mirs[pk]=='C' and tars[pk]=='T':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='A' and tars[pk]=='A':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='T' and tars[pk]=='T':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='C' and tars[pk]=='C':
                               penalty_score = penalty_score + 1

                           if mirs[pk]=='G' and tars[pk]=='G':
                               penalty_score = penalty_score + 1

                       #print(ee2+str(s1)+'..'+str(e1)+"(+)"+mm1+e+seq)
                       sf.write(ee2+"\t"+str(s1)+'..'+str(e1)+"(-)"+"\t5'"+mm1+"3'\t"+e+"\t3'"+seq+"5'\n")
                       if penalty_score<=P_score:
                           sf1.write(ee2+"\t"+str(s1)+'..'+str(e1)+"(+)"+"\t5'"+mm1+"3'\t"+e+"\t3'"+seq+"5'\n")
                           
                   if mark==1:
                      break
                    
sf.close()
                  
    


