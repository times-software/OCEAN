import os, sys, shutil
from distutils.dir_util import copy_tree


pseudo = "psps"

a = os.path.join("tester", pseudo)
aalist = os.listdir(a)
alist = sorted(aalist)
numa = range(1, 121)

print alist


#znucl code:

def element():
	i = int(input("Enter which element you would like to add (by atomic number):"))
	if i in numa:
		return i
	else:
		return 0
		sys.exit(0)
	

nucl = element()
znucl = str(nucl)

p = znucl
#changes based on which thing needs to be checked. in this case it's znucl

def permission():
        if p not in alist:
                return "A"
                #for allowed
        elif p in alist:
                return "D"
                #for denied
        else:
                return "?"

allowed = permission()

s = os.path.join(a, znucl)

if allowed == "A":
        os.makedirs(s)
	sc = True
	#user creates a new directory
elif allowed == "D":
        print "That directory already exists."
        inp = raw_input("Would you like to continue without overwriting?")
        if "y" in inp[0]:
                print "Continuing in " + znucl
		sc = False
        #user continues in already existing directory
        elif "n" in inp[0]:
                overwrite = raw_input("Warning, this will overwrite an already existing directory. Type yes to overwrite:")
        	if "y" in overwrite[0]:
			shutil.rmtree(s)         
			os.makedirs(s)
			print "The directory was overwritten."
			sc = True
		else:
			print "The directory was kept as is."
			sc = False
        #Allows them to overwrite but warns them       
        else:
                print "Something went wrong."
                sys.exit(0)
else:
        print "Huh."
        sys.exit(0)







#long code for semicore:

T = os.path.join(s, "T")
F = os.path.join(s, "F")
#for later use with the user either going into False or True

def continuer():
	 inpu = raw_input("Which would you like to continue in (T or F)?")
         if "t" in inpu.lower():
         	return "1"
		q = os.path.join(s, "T")
         elif "f" in inpu.lower():
		return "2"
         	q = os.path.join(s, "F")
         else:
                print "That's not an option."
                sys.exit(0)
	#if True and False exist, user has to decide which one to enter


#determining semicore:

if sc == True:
	scinput = raw_input("Is semicore an option for this element?")
	if "y" in scinput:
		print "Adding T and F to " + znucl + " directory."
		os.makedirs(T)
		os.makedirs(F)
		print "Continuing on in 'T'. Make sure to add to 'F' later."
		q = os.path.join(s, "T")
	#the user chooses to have T and F in new znucl directory
	elif "y" not in scinput:
		print "Adding F to " + znucl + " directory."
		os.makedirs(F)
		q = os.path.join(s, "F")
	#the user chooses to have only F in new znucl directory
	else:
		print "Something went wrong."
		sys.exit(0)
#the user made a new directory for an element so they decide if there's T and F or just F as options

elif sc == False:
	slist = os.listdir(s)
	lengths = len(slist)

	if lengths == 0:
		scinput = raw_input("Would you like to add 'T' as an option?")
                if "y" in scinput:
                        print "Adding 'T' to " + znucl + " directory."
                        os.makedirs(T)
			os.makedirs(F)
			print "Both 'T' and 'F' have been added. Continuing to 'T', make sure to add to 'F' later."
                        #user adds True and False to directory
			q = os.path.join(s, "T")
			#there's just an empty False for now. This puts trust in users to not just leave an empty False
		elif "y" not in scinput:
                	print "Adding only 'F' to "  + znucl + " directory."
                        os.makedirs(F)
			q = os.path.join(s, "F")
			#user adds only False to directory
		else:
			print "Can't add nothing."
			sys.exit(0)

	elif lengths == 1:
		scinput = raw_input("Would you like to add 'T' as an option?")
		if "y" in scinput:
			print "Adding 'T' to " + znucl + " directory."
			os.makedirs(T)
			#user adds True to directory
			q = os.path.join(s, "T")
			#user is directed to True because they gotta add stuff to it
		elif "y" not in scinput:
			print "Continuing on into 'F'"
			q = os.path.join(s, "F")
		#user doesn't add anything and continues into False	
		else:
			sys.exit(0)
	#False is only currently existing option, user has to choose to add True or not

	elif lengths == 2:
		print "True and False already exist."
		if continuer() == "1":
		        q = os.path.join(s, "T")
		elif continuer() == "2":
        		q = os.path.join(s, "F")

		#directs user to method that decides if they are going in True or False 
	else:
		print "Something's wrong with the directory."
		sys.exit(0)
#the directory already existed for their choice of znucl, so now they need to decide if they will add T to semicore,
#continue on in T or F if the directory has both, or to continue on in F (or exit i guess)






#quality code:

def filecopy():
	finallist = os.listdir(finalpath)
        print "In final directory: " + finalpath
        #takes finalpath of where the user wants to add files to and lists current files in it and what it is
	#copy code from finder.py to decide if it was a valid input and such
        userchoice = raw_input("Where are the pseudofiles you are copying into here?")
        userlist = os.listdir(userchoice)
	#asks where user wants to find these files and creates a list of all the files that will be copied	

        for file in userlist:
		filepath = os.path.join(userchoice, file)
                if file not in finallist:
                        shutil.copy(filepath, finalpath)
			print file + " was copied."
                        #if file doesn't exist in finalpath, it goes ahead and copies the files given by the user
                #add part after it lists files copied where it asks user which is which and then writes info.txt
		elif file in finallist:
                        print file + "was not copied because it already exists."
                        owcopy = raw_input("Would you like to overwrite this file?")
                        if "y" in owcopy:
                                shutil.copy(filepath, finalpath)
				print file + " was copied."
                        if "n" in owcopy:
                                print file + " was not copied."
                        #warns user if file already exists in finalpath and asks if they want to copy it anyway
                else:
                        print "Something went wrong."
#just checks files to see if they already exist in finalpath and whether to overwrite them, it's used in quality code

print q

qqlist = os.listdir(q)
qlist = sorted(qqlist)
lengthq = len(qlist)
numq = map(int, qqlist)


if lengthq == 0:
	print "There are no pre-existing directories"
	userinput = "y"
else:
	print "The available options are:" 
	print qlist
	userinput = raw_input("Would you like to add a different quality (or overwrite a pre-existing one)?")
#if the list is empty, it tells the user there are none and makes them add a quality.
#if the list has stuff, it lists them and then asks the user if they would like to add another quality.


if "y" in userinput:
	numinput = int(input("What quality would you like to add?"))

	if numinput not in numq and numinput in range(1, 1001):
		finalpath = os.path.join(q, str(numinput))
		os.makedirs(finalpath)
		print str(numinput) + " was added."
		#directory is made
		filecopy()
	#if the user's choice of quality does not exist and is in range, the directory is made. files are copied	
	elif numinput in numq:
		print "That quality already exists."
		owquality = raw_input("Would you like to overwrite the pre-existing directory " + str(numinput) + ":") 
		if "y" in owquality[0]:
			finalpath = os.path.join(q, str(numinput))
			shutil.rmtree(finalpath)
                        os.makedirs(finalpath)
                        print "The directory was overwritten."
			#overwrites directory
                        filecopy()
		#the user decides to overwrite pre-existing directory and they copy their files in
                else:
                        print "The directory was kept as is."
                        finalpath = os.path.join(q, str(numinput))
			#continue in pre-existing directory		
			filecopy()
		#if the user decides not to change quality, they continue in the pre-existing directory and copy files

	else:
		print "The quality given was out of range."
		sys.exit(0)
	#if quality wasn't in the range given, the code ends
#if user said they wanted to add another quality

elif "y" not in userinput:
	ask = int(input("Which quality would you like to continue into?"))
	
#########################################################################################################
	def qual():
        
#repeat variable list if needed here

        	if ask in numq:
                	return str(ask)
        	#if the user's input existed in the list, that option is chosen
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

	checkask = ask
	#asks for input. It was necessary to remove it from qual() for checker(). checkask is just to avoid problems within checker()
	firstq = qual()

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

#########################################################################################################	
#using finder.py code to check that quality made sense

	finalpath = os.path.join(q, quality)
	filecopy()
#if user picked a valid quality, they continue into that quality and copy files 

else:
	print "why"	

