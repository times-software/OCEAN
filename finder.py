#This is the automated one
import shutil, os, sys
from distutils.dir_util import copy_tree

pseudo = "psps"
#this will be where all the pseudofiles are


#znucl code:

a = os.path.join("/flash", "jtv1", "OCEAN", "bin", pseudo)
aalist = os.listdir(a)
alist = sorted(aalist)
numa = map(int, alist)
#lists everything in pseudo directory    

commonz = open("Common/znucl", "r")
searchz = commonz.readlines()
commonz.close()

#looks at znucl file in Common


print "Available elements (by atomic number):"
print alist


for number in numa:
        if number == 1:
                name = "H"
        elif number == 2:
                name = "He"
        elif number == 3:
                name = "Li"
        elif number == 4:
                name = "Be"
        elif number == 5:
                name = "B"
        elif number == 6:
                name = "C"
        elif number == 7:
                name = "N"
        elif number == 8:
                name = "O"
        elif number == 9:
                name = "F"
        elif number == 10:
                name = "Ne"
        elif number == 11:
                name = "Na"
        elif number == 12:
                name = "Mg"
        elif number == 13:
                name = "Al"
        elif number == 14:
                name = "Si"
        elif number == 15:
                name = "P"
        elif number == 16:
                name = "S"
        elif number == 17:
                name = "Cl"
        elif number == 18:
                name = "Ar"
        elif number == 19:
                name = "K"
        elif number == 20:
                name = "Ca"
        elif number == 21:
                name = "Sc"
	elif number == 22:
                name = "Ti"
        elif number == 23:
                name = "V"
        elif number == 24:
                name = "Cr"
        elif number == 25:
                name = "Mn"
        elif number == 26:
                name = "Fe"
        elif number == 27:
                name = "Co"
        elif number == 28:
                name = "Ni"
        elif number == 29:
                name = "Cu"
        elif number == 30:
                name = "Zn"
        #need to enter more just in case he enters a lot in the future
        else:
                name == "not inputed yet."
        #by atomic number it lists abbreviated name

        print "Atomic number " + str(number) + " is " + name
#prints what atomic number is what symbol for the user


def pathmaker():
	oblist = []
	print "in pathmaker"
	semicore = semic() #for now
	qual = fakequalitylist() #for now
	for object in zdict:
		start = os.path.join(a, zdict[object])
		print start
		
		next = os.path.join(start, semicore[object])
		print next
	
		last = os.path.join(next, qual[object])
		print last
		
		oblist.append(last)
		#need to grab the letter associated with the object of zdict for the dict semicore
#		for thing in semicore:
#			next = os.path.join(start, semicore[thing])
#			print next
#			#for now
#			for thingy in qual: #change to qualitylist()
#				last = os.path.join(next, qual[thingy])
#				print last
#				oblist.append(last)
	print oblist
        
#	for object in semicorelist:
#                ob = os.path.join(item, object)
#                oblist.append(ob)
#                return oblist
#supposed to combine paths from znucllist with each of the strings from semicorelist to make the next path for quality

def semic():
	commons = open("Common/semicore", "r")
	searchs = commons.readlines()
	commons.close()
	z = 0
	y = 0
	sdict = {}
	newletters = []
	print "Will semicore file be included?"
	ssslist = []	
	for index, line in enumerate(searchs):
		nothing = index, line
		print nothing[1] 
		letters = nothing[1].rsplit(' ', x)
		print letters
		for l in letters:
			print l
			if "\n" in l:
				l = l[:-1]
			newletters.append(l)		
		print newletters
		#should break down the semmicore text into the individual letters based on how many znucls there are
		
		for bl in newletters:
			print bl
			z = z + 1 
			y = y + 1
			sdict[y] = bl
		print "Semicore dictionary?"
		print sdict
		return sdict
		
		for bl in sdict:
			print "A letter?"
			print sdict[bl]
		
		#for semicore in sdict:
		#	if semicore in slist:
				
		#		print semicore
	    	#		if semicore.lower() in "t":	
		#			ssslist.append('T')
		#		        
		#		else:
                #			ssslist.append('F')
		#		return ssslist
		#for now i'm not going to check that it's actually one of the options, im just going to pick it			
#looks at semicore file in Common and for each znucl it sees if for it's semicore options which one it should
#pick based on the semicores that correspond to each znucl in Common 

def qualitylist():
     	q = item
       	qqlist = os.listdir(q)
       	qlist = sorted(qqlist)
       	lengthq = len(qqlist)  
       	print "\nAvailable options for quality:"
       	print qlist
	qqqlist = []

	commonq = open("Common/pp.quality", "r")
	searchq = commonq.readlines()
	commonq.close()

	for index, line in enumerate(searchq):
		everything = index, line 
		print everything[1]
		for quality in qlist:
			if quality in everything[1]:
				for i in range(len(znucllist)):
					qqqlist.append(str(quality))
					return qqqlist 
def fakequalitylist():
	t = 0
	qdict = {}
	for thing in zdict:
		print thing
		t = t + 1
		qdict[t] = ("40")
	return qdict
	print sorted(qdict)		

#makes path and then lists quality files. They're sorted and the length is measured. The sorted list is printed 

zlist = []

for index, line in enumerate(searchz):
	something = index, line
	print "\nList of everything found in znucl file (in Common):"
        print something[1]
	for znucl in alist:
        	if znucl in something[1]:
                	print znucl #not needed atm
			zlist.append(znucl)
print "The complete list of znucls"
print zlist
#looks through znucl file and grabs the second item in the list it returns (the string with all the data)
#then if any of the stuff in that list is one of the options listed in the psps directory, it prints each of the options
#that matched something

znucllist = []
for item in zlist:
	zzz = os.path.join(a, str(item))
	znucllist.append(zzz)

print "\nList of znucl paths:"
print znucllist
print len(znucllist)
#makes an empty list and then takes all the znucls from above, makes a path for them, and adds that path to a list

print "\nPath for each znucl:"
for item in znucllist:
	print item

print "\nDictionary for znucls:"
x = 0
zdict = {}
for nucl in zlist:
	x = x + 1
	zdict[x] = nucl
print sorted(zdict)
#makes the list of znucls into a dictionary so that it's easier to keep track of what's what 
	

#print "\nList of semicores:"
#semicorelist = []
#for item in znucllist:
#	s = item
#	sslist = os.listdir(s)
#	slist = sorted(sslist)
#	print "List of options in znucl directories:"
#	print slist
#	lengths = len(sslist)
#	#makes path and lists all of semicore options, sorts it, and sees how long it is
#
#	if lengths == 1:		
#        	sssc = 'F'
#		semicorelist.append(sssc)
#		pathmaker()
#        	#if there's only one option for a znucl, that option is picked
#	elif lengths == 2:
#		print "List of semicores from Common (for ti in this case):"
#	     	for sssc in semic():
#			print sssc
#			semicorelist.append(sssc)
#			pathmaker()
#  		#if there's two options for a znucl, code looks to see if it's T or F 
#	else:
#       		print "Something went wrong."
#       		sys.exit(0)
#currently not really working

#finalist = []
#for item, i in pathmaker():
#	print item
#	for thing in qualitylist():
#		fff = os.path.join(item, thing)			
#		finalist.append(fff)
		
#	print finalist
pathmaker()	
	
						
