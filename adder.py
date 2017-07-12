#used to add stuff, like an opposite of finder.py

import os, sys

pseudo = "psps"

a = os.path.join("/flash", "jtv1", "OCEAN", "bin", pseudo)
aalist = os.listdir(a)
alist = sorted(aalist)


print alist

nucl = int(input("Enter which element you would like to add (by atomic number):")
znucl = str(nucl)

p = znucl
#changes based on which thing needs to be checked. in this case it's znucl

def permission():
	if p in alist:
		return "A"
		#for allowed
	elif p not in alist:
		return "D"
		#for denied
	else:
		return "?"

allowed = permission()

s = os.path.join(a, znucl)

if allowed == "A":  
	os.makedirs(s)

elif allowed == "D":
	print "That directory already exists."
	input = raw_input("Would you like to continue in existing directory?")
	if "y" in input:
		print "Continuing in" + znucl
	#user continues in already existing directory
	elif "y" not in input:
		print "You cannot overwrite the directory. Try again."
		sys.exit(0)
	#currently doesn't allow them to overwrite but can be changed to allow someone to 	
	else:
		print "Something went wrong."
		sys.exit(0)
else:
	print "Huh."
	sys.exit(0)	
