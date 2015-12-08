# Psi

## the vector psi
In ocean psi is the vector of electron+hole pairs for the BSE Hamiltonian. 
The goal here is to deal with both core and valence BSE in a single set of 
functions. 

### Globals
It is acceptable for psi to have globals that pertain to the physical system. 
This is because for each given system there are certain universal variables, 
e.g., xmesh, kmesh, cell volume, etc. 

### Public function list
1. init
2. new
3. delete
4. load




#### psi_init
Takes the system and MPI info and initializes all of the globals. Hopefully 
there will be no need to pass system or MPI info into any other psi routine. 


#### psi_new
Optionally takes in a psi to copy over. Allocates the vector spaces and 
either sets values or zeros it out

#### psi_kill/delete

#### psi_load
psi_load is a front end, which at the moment, only will call into the legacy 
ocean (cksc/cksv for core). In the future this is where things can branch 
out to accept improvements or make more general.

