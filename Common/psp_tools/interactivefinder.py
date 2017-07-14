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


def element():	
        type = int(input("Type which element (by number):"))
	if type in numa:
		return type
	else: 
		return 0    
#if user input matches one of the options, that option is selected. If the input does not, the program is stopped


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


nucl = element()
#nucl determines the outcome of begin() which in turn decides semicore files

if nucl != 0:
	znucl = str(nucl)
elif nucl == 0:
	print "That option is not available. Try again."
	sys.exit(0)
else:
	print "Error."
	sys.exit(0)
#znucl is used for the final path
#the if statement is just to make sure the user gave a valid answer before adding that to the final path






#semicore code:

def semic():
        word = raw_input("Type yes to include semicore file:")
        if word.lower() in "yes":
                return True
        else:
                return False
#for the semicore file



s = os.path.join(a, znucl)
sslist = os.listdir(s)
slist = sorted(sslist)
lengths = len(sslist)
#makes path and lists all of semicore options, sorts it, and sees how long it is

if lengths == 1:
	semicore = slist[0]
	#if there's only one option, that option is picked
elif lengths == 2:
	if semic() == True:
		semicore = "T"
	#if there's two options it asks the user if they want it or not
	else:	
		semicore = "F"
	#if they don't say they want it, it's not included
else:
	print "Something went wrong."
	sys.exit(0)
	#if there's more than two options something is weird and the code just ends. There should only be T and F or F
#decides if semicore is included or not





#quality code:

def qualitylist():
        q = os.path.join(a, znucl, semicore)
        qqlist = os.listdir(q)
	qlist = sorted(qqlist)
	lengthq = len(qqlist)  
        print "Available options for quality:"
        print qlist
#makes path and then lists quality files. They're sorted and the length is measured. The sorted list is printed for the user.	




def qual():
	q = os.path.join(a, znucl, semicore)
        qqlist = os.listdir(q)
        qlist = sorted(qqlist)
        lengthq = len(qqlist)
        numq = map(int, qqlist)
#repeated so i can access these variables

        if ask in numq:	
		return str(ask)
	#currently this doesn't really work
	#i need to loop it through all different options and if ask = any one of the the possible options, str(ask) is returned
	elif lengthq == 1:
                return qlist[0]	
        #will choose the only option because there is only 1 directory listed by qlist
	elif ask > 0 or ask < 150:
		if ask <= 40:
			outcome = 40
                	if outcome in numq:
				return str(outcome)
			else:
				return "0"
                	#will choose 40 if list is longer and user didn't give an exact number, as long as 40 matches the list
        	elif ask > 40 and ask <= 80:
                	outcome = 80
			if outcome in numq:
				return str(outcome)
			else:
				return "0"
                	#will choose 80 if list is longer and user didn't give an exact number, as long as 80 matches the list 
		elif ask > 80:
			outcome = 100
			if outcome in numq:
				return str(outcome)
			else:
				return "0"
			#will choose 100 if list is longer and user didn't give an exact number, as long as 100 matches the list
        	else:
                	sys.exit(0)
                	#exit
	#if the users input exists in the quality directory
	else:
		return "0"
		sys.exit(0)
		#ends code if it's not valid 


qualitylist()
#lists options for user
ask = int(input("What quality would you like:"))
checkask = ask
#asks for input. It was necessary to remove it from qual() for checker(). checkask is just to avoid problems within checker()
firstq = qual()
#takes value from qual() and adds it to the final path


def checker():
	if firstq == "0":
		print "No valid quality was returned. Try again."
		sys.exit(0)
		#if the user's input didn't end up on any of the available options, the code ends.
	elif str(checkask) != firstq:
		print "Your request, " + str(checkask) + ", was not available. " + firstq + " was picked instead."
		retry = raw_input("Type 'retry' to re-enter quality or accept value and continue:")
		if "re" in retry.lower():
			ask = int(input("What quality would you like now:"))
			qual()
			quality = qual()
		elif "quit" in retry.lower():
			print "Code ended."
			sys.exit(0)
		else:
			print firstq + " was chosen."
			quality = qual()
	else:
		quality = qual()
#checks to make sure user's input was valid and if it matched the listed options. If it doesn't match, they're notified.

checker()
quality = qual()




#for the final directory with the needed files:

finalpath = os.path.join("/flash", "jtv1", "OCEAN", "bin", pseudo, znucl, semicore, quality)
print finalpath
print "Here's the directory."


os.chdir("/flash/maz2")
#changes current working directory so that files can be copied into there

userchoice = raw_input("Where would you like to copy your files? Type quit to stop:")
if "quit" in user:
	print "Stopped."
else:
	copy_tree(finalpath, userchoice)
#copies files from final path to users choice of directory or to a newly created directory
	print "The pseudofiles have been copied for you."
