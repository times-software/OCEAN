# Psi

## the vector psi
In ocean psi is the vector of electron+hole pairs for the BSE Hamiltonian. 
The goal here is to deal with both core and valence BSE in a single set of 
functions. 

### Globals
It is acceptable for psi to have globals that pertain to the physical system. 
This is because for each given system there are certain universal variables, 
e.g., xmesh, kmesh, cell volume, etc. 

### Function list
1. init
2. new
3. delete




#### psi_init
Takes the system and MPI info and initializes all of the globals. Hopefully 
there will be no need to pass system or MPI info into any other psi routine. 
