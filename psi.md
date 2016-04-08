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

### Type OCEAN_vector
The OCEAN_vector is the object that holds the core and valence coefficients. 
The data can be distributed in 3 ways. The object tracks which data holders 
are allocated, which are valid, and this module handles moving data between 
them. This allows the external interaction routines to always be able to 
ask for and recieve the expected vector full of valid data. The object also 
keeps track of the MPI communications: communicator, requests, statuses, etc.




#### Data components
1. r & i -- hold the entire vector for the core (nband, nkpt, nalpha)
2. write_r & write_i -- buffers for receiving 
3. store_r & store_i -- fully distributed vector so that it takes minimal space

#### Parameters


#### Status flags
1. storage_type -- bitwise flag for what is allocated
2. valid_storage -- bitwise flag for which parts of vector have good data
3. update -- logical for is no data *good*; requires a summation across all procs

#### Comms
1. core_comm
2. core_myid
 


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

