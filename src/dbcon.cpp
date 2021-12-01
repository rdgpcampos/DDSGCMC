/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

------------------------------------------------------------------------- */
// Concurrent MC calculation for database construction
//


#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include <dirent.h>
#include <algorithm>
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

#include <random>           // these are include files for MC
#include "mcvar2.h"

using namespace LAMMPS_NS;
int main(int narg, char **arg)
{

    if (narg != 4) {
        printf("Syntax: npar in.lammps mc.ctl rnd_seed\n");
        exit(1);
    }

    char *infile = arg[1];
    char *mc_ctl = arg[2];
	char *rnd_seed = arg[3];

    // Initial settings here
    MCVars mc;
    try
    {
        mc.read_ctl(mc_ctl);
    }
    catch (char* str)
    {
        printf("%s\n",str);
        exit(1);
    }

    MPI_Init(&narg,&arg);
    int me,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if(me==0) 
    {
        fprintf(stdout,"### Atom-swapping Monte Carlo simulator ###\n");
        fprintf(stdout,"### Temperature = %.1f K\n",mc.temp);
        fprintf(stdout,"### Concentration = %f\n",mc.c0);
        fprintf(stdout,"### Iteration = %d cycles\n",mc.niter);
        fprintf(stdout,"### Resume from %d steps\n",mc.nstart);
        fprintf(stdout,"### Output interval = %d steps\n",mc.nout);
        fprintf(stdout,"### Rejection limit = %d cycles\n",mc.rej_limit);
        fprintf(stdout,"### Boltzmann factor = %g\n",mc.beta);
        fprintf(stdout,"###########################################\n");
    }

    const int ninstance = 2;
    if (nprocs % ninstance) {
        if (me == 0)
            printf("ERROR: Total procs must be divisble by 2\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    if(me==0)
    {
        FILE *fp;
        try{
            if((fp=fopen("log.mc","w"))==NULL) throw 1;
            fprintf(fp,"#  MC_step      E1              E2              E1*             E2*            dE1             dE2           Nrej   t1 t2     flag\n");
            fclose(fp);
            fp=NULL;
        }
        catch(int e)
        {
            // exception: output file cannot open
            fprintf(stderr,"Error in output (header)\n");
            return EXIT_FAILURE;
        }
    }

    // create one communicator per instance each with P/N procs
    const int instance = me*ninstance / nprocs;
    const int npar = nprocs/ninstance;
    const int me_lammps = me%npar;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,instance,me_lammps,&comm_lammps);

//if(me==0) fprintf(stderr,"me nprocs instance npar me_lammps\n");
//fprintf(stderr,"%d %d %d %d %d\n",me,nprocs,instance,npar,me_lammps);
  
    // each instance: unique screen file, log file, temperature

    char str1[32],str2[32],str3[32];

    char **lmparg = new char*[5];
    lmparg[0] = NULL;                 // required placeholder for program name
    lmparg[1] = (char *) "-screen";
    //sprintf(str1,"screen.%d",instance);
    lmparg[2] = (char *) "none";
    lmparg[3] = (char *) "-log";
    sprintf(str2,"log.lammps.%d",instance);
    lmparg[4] = str2;

    // open N instances of LAMMPS
    // either of these methods will work

	// Each instance corresponds to a system!!!

    LAMMPS *lmp = new LAMMPS(5,lmparg,comm_lammps);
    delete [] lmparg;

    // run input script thru all instances of LAMMPS (initial settings)
    char *eachinfile = new char[32];
    // input file names = inlammps_A.lmp & inlammps_B.lmp
    if(instance==0) sprintf(eachinfile,"%s_A.lmp",infile);
    if(instance==1) sprintf(eachinfile,"%s_B.lmp",infile);
    if(me_lammps==0) fprintf(stderr,"reading %s...\n",eachinfile);
    lammps_file(lmp,eachinfile); 
    if(me_lammps==0) fprintf(stderr,"Input file %s read\n",eachinfile);
    delete [] eachinfile;
   

    // for different number of atoms in each instance
    int natoms = static_cast<int> (lmp->atom->natoms);
    int *tmp = new int[ninstance];
    int *all_tmp = new int[ninstance];
    for(int i=0;i<ninstance;i++) tmp[i] = all_tmp[i] = 0;
    if(me_lammps==0) tmp[instance] = natoms;
    MPI_Allreduce(tmp,all_tmp,ninstance,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    int natomA = all_tmp[0];
    int natomB = all_tmp[1];

    delete [] tmp;
    delete [] all_tmp;

    double *xA = new double[3*natomA];
    double *xB = new double[3*natomB];
    double *xA_start = new double[3*natomA];
    double *xB_start = new double[3*natomB];
    double *vA = new double[3*natomA];
    double *vB = new double[3*natomB];
    double *vA_start = new double[3*natomA];
    double *vB_start = new double[3*natomB];
    int    *tA = new int[natomA];
    int    *tB = new int[natomB];

    if(instance==0) lammps_gather_atoms(lmp,(char *) "type",0,1,tA); // if lammps version is 29Oct20 or above, write 1,1 instead of 0,1
    if(instance==1) lammps_gather_atoms(lmp,(char *) "type",0,1,tB);
    
    MPI_Bcast(xA,3*natomA,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(tA,natomA,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(xB,3*natomB,MPI_DOUBLE,npar,MPI_COMM_WORLD);
    MPI_Bcast(tB,natomB,MPI_INT,npar,MPI_COMM_WORLD);

//-------------------MD----------------------
// Timestep dt needs to be such that n*dt >= eq_t, where n is the number of steps in MD run
// and eq_t is the equilibration time defined in the NVT command
    char command1[64], command2[64], command3[64];
    lmp->input->one("timestep 0.001"); 
    lmp->input->one("run 0");
    lmp->input->one("thermo_style custom step vol press temp pe atoms");
    lmp->input->one("thermo 100");

//-------------------------------------- Set initial velocity of the atoms --------------------------------------------------
	if(instance==0) sprintf(command1,"velocity all create %.1f 14236 rot yes mom yes dist gaussian",mc.temp*2);
	if(instance==1) sprintf(command1,"velocity all create %.1f 54379 rot yes mom yes dist gaussian",mc.temp*2);
    lmp->input->one(command1);
//-------------------------------------- NVT ensemble (If needed, change to NPT) --------------------------------------------
    double *press1 = (double *) lammps_extract_compute(lmp,(char *) "thermo_press",0,0);
    double p0 = *press1;
    if(me==0) fprintf(stderr,"Initial pressure: %.2f\n",p0);
    press1=NULL;
	sprintf(command3,"fix nf all nvt temp %.1f %.1f 1",mc.temp,mc.temp); // Set for NVT or NPT
    lmp->input->one(command3);
//-------------------------------------- Avoid "drifting ice cube" problem --------------------------------------------------
    lmp->input->one("fix mom all momentum 50 linear 1 1 1 angular");
    if(me_lammps==0) fprintf(stderr,"Clear\n");
//-------------------------------------- Thermal equilibration --------------------------------------------------------------
    lmp->input->one("run 2000");
    lmp->input->one("unfix nf");
	char restart_command[128];
    sprintf(restart_command,"fix nf all nvt temp %.1f %.1f 0.1",mc.temp,mc.temp); // Set for NVT or NPT
    lmp->input->one(restart_command);
//-------------------------------------- Get atomic positions and velocities before MC calculation --------------------------
    if(instance==0) lammps_gather_atoms(lmp,(char *) "x",1,3,xA);
    if(instance==1) lammps_gather_atoms(lmp,(char *) "x",1,3,xB);
    if(instance==0) lammps_gather_atoms(lmp,(char *) "x",1,3,xA_start);
    if(instance==1) lammps_gather_atoms(lmp,(char *) "x",1,3,xB_start);
    if(instance==0) lammps_gather_atoms(lmp,(char *) "v",1,3,vA);
    if(instance==1) lammps_gather_atoms(lmp,(char *) "v",1,3,vB);
    if(instance==0) lammps_gather_atoms(lmp,(char *) "v",1,3,vA_start);
    if(instance==1) lammps_gather_atoms(lmp,(char *) "v",1,3,vB_start);
//--------------------//

//-------------------------------------- Calculate potential energy ---------------------------------------------------------
    lmp->input->one("compute calcPE all pe");
    double *ptr1 = (double *) lammps_extract_compute(lmp,(char *) "calcPE",0,0);

    double e0 = *ptr1;
    ptr1=NULL;

//-----------------------------------------
    // settings before MC simulation
    int nrej=0;
    int count=0;
    int id1,id2;
    int t1,t2;
    double p=0.0;
    int acc=0;
    double de=0.0;
    double r=0.0;
    double *es      = new double[ninstance];
    double *all_es  = new double[ninstance];
    double *e0s     = new double[ninstance];
    double *all_e0s = new double[ninstance];
    int *tA0 = new int[natomA]; 
    int *tB0 = new int[natomB]; 

    for(int i=0;i<ninstance;i++) e0s[i] = all_e0s[i] = 0.0;
    if(me_lammps==0) e0s[instance] = e0;
    MPI_Reduce(e0s,all_e0s,ninstance,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    // random variables settings
    //std::random_device rnd; // device dependent random variable
    //unsigned rnd = 314159265; // fixed root random seed
    unsigned rnd = static_cast<unsigned>(*rnd_seed); // fixed root random seed
    unsigned seed = rnd;
    std::mt19937 mt(seed);
    std::mt19937 mt_prob(seed+1);
    std::uniform_int_distribution<int> sidA(0,natomA-1);
    std::uniform_int_distribution<int> sidB(0,natomB-1);
    std::uniform_real_distribution<double> prob(0.0,1.0);
    int *elmA = new int[mc.ntype];
    int *elmB = new int[mc.ntype];
    for(int i=0;i<mc.ntype;i++) elmA[i] = 0;
    for(int i=0;i<mc.ntype;i++) elmB[i] = 0;
    for(int i=0;i<natoms;i++) elmA[tA[i]-1]++;
    for(int i=0;i<natoms;i++) elmB[tB[i]-1]++;
 //-------------------------------------- Starts MC loop here ---------------------------------
    if(me==0) fprintf(stderr,">>>Random seed = %u\n",seed);
    if(me==0) fprintf(stderr,">>>Start\n");

    for(int n=mc.nstart; n<mc.niter; n++)
    {

        if(me==0&&n%mc.nout==0&&nrej==0) fprintf(stderr,">>>>Iteration %d\n",n);
        if(nrej>mc.rej_limit)
        {
		    if(me==0) fprintf(stderr,"Reached rejection limit at step %d\n",n);
		    break;
        }
        
//--------------------------------------- Save atomic type data to vector ---------------------
	    memcpy(tA0,tA,natomA*sizeof(int));
	    memcpy(tB0,tB,natomB*sizeof(int));

//--------------------------------------- Continue MD run -------------------------------------
        if(instance==0) lammps_scatter_atoms(lmp,(char *) "v",1,3,vA);
        if(instance==1) lammps_scatter_atoms(lmp,(char *) "v",1,3,vB);
	    lmp->input->one("run 500"); // Number of MD steps
//--------------------------------------- Save atomic velocities to vectors -----
        if(instance==0) lammps_gather_atoms(lmp,(char *) "v",1,3,vA);
        if(instance==1) lammps_gather_atoms(lmp,(char *) "v",1,3,vB);
//--------------------------------------- Get energy immediately prior to MC cycle ------------
	    double *ptr = (double *) lammps_extract_compute(lmp,(char *) "calcPE",0,0);
        double ec = *ptr;
        ptr=NULL;
	
        for(int i=0;i<ninstance;i++) e0s[i] = all_e0s[i] = 0.0;
        if(me_lammps==0) e0s[instance] = ec;
        MPI_Reduce(e0s,all_e0s,ninstance,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        
//--------------------------------------- MC cycle (# steps = # of atoms) ---------------------
        for(int j=0; j<natoms; j++) 
        {
    	
            // id selection
            // trial move selection
	        //double trial_move = prob(mt);
	        double trial_move = 1;
	        int flag; // Flag to define which type of trial move to perform
	        if(trial_move>=0.5) 
            {
		        flag = 1; // Trial move between replicas
	        } else if (trial_move >= 0.25) 
            {
		        flag = 2; // Trial move in replica A
	        } else 
            {
		        flag = 3; // Trial move in replica B
	        }
	        if(flag==1)
            { // Exchange between replicas (Single processor)
			    t1=t2=0;	            

			    while(t1==t2) 
                {
          		    id1 = sidA(mt);
               	    id2 = sidB(mt);
               	    t1 = tA[id1];
               	    t2 = tB[id2];
	    	    }

	        } else if(flag==2) 
            { // Exchange within replica A
			    t1=t2=0;

			    while(t1==t2) 
                {
				    id1 = sidA(mt);
				    id2 = sidA(mt);
				    t1 = tA[id1];
				    t2 = tA[id2];
			    }
		
        	
	        } else if(flag==3) 
            { // Exchange within replica B
			    t1=t2=0;

			    while(t1==t2) {
				    id1 = sidB(mt);
				    id2 = sidB(mt);
				    t1 = tB[id1];
				    t2 = tB[id2];
			    }
	        }

	        // Perform exchange (according to each move type)
	        if(flag==1) 
            {
        	    tA[id1] = t2;
        	    tB[id2] = t1;
	        } else if(flag==2) 
            {
        	    tA[id1] = t2;
        	    tA[id2] = t1;	
	        } else 
            {
        	    tB[id1] = t2;
        	    tB[id2] = t1;
	        } 
            
            if(instance==0 && flag!=3) lammps_scatter_atoms(lmp,(char *) "type",0,1,tA);
            if(instance==1 && flag!=2) lammps_scatter_atoms(lmp,(char *) "type",0,1,tB);
            
            lmp->input->one("run 0"); // Update variable values

            // STRUCTURAL RELAXATION (Maybe unnecessary)
	        //lmp->input->one("min_style cg");
	        //lmp->input->one("minimize 1.0e-5 0.0 100000 100000");
            //lmp->input->one("run 0");
	
            // Sample free energy
            double *ptr = (double *) lammps_extract_compute(lmp, (char *) "calcPE",0,0);
            double ec = *ptr; 
            ptr=NULL;

            // calc acceptance ratio (energy difference)
            for(int i=0;i<ninstance;i++) es[i] = all_es[i] = 0.0;
            if(me_lammps==0) es[instance] = ec;
            MPI_Reduce(es,all_es,ninstance,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); // adding energy calculated in each processor
            if(me==0)
            {
       	        de = (all_es[0]+all_es[1])-(all_e0s[0]+all_e0s[1]);
           	    p  = de > 0 ? exp(-de*mc.beta) : 1.0;
           	    r = prob(mt_prob);
        	    acc = (p>r) ? 1 : 0;
            }
            MPI_Bcast(&acc,1,MPI_INT,0,MPI_COMM_WORLD);
            count++;

            if(me==0 && count%10==0)
            {
            	FILE *fp;
            	try{
                	if((fp=fopen("log.mc","a"))==NULL) throw 1;
               		fprintf(fp,"%8d%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e%8d   %3d%3d%8d\n"
                            ,count,all_e0s[0],all_e0s[1],all_es[0],all_es[1],all_es[0]-all_e0s[0],all_es[1]-all_e0s[1],nrej,t1,t2,flag);
                	fclose(fp);
                	fp=NULL;
            	} catch (int e) {
                	fprintf(stderr,"Error in output\n");
                	break;
            	}
            }
            
            if(!(count%mc.nout))
            {
              	//char str[32];
               	//sprintf(str,"reset_timestep %d",count);
               	char dumpstr[BUFSIZ];
               	char data[32]="mass type xs ys zs";
               //	lmp->input->one(str);
               	// for db output
               	sprintf(dumpstr,"write_dump all cfg ./dump/mc%d_%07d.cfg.gz %s %s %s",instance,count,data,mc.mod1,mc.mod2);
               	lmp->input->one(dumpstr);
               	//sprintf(str,"reset_timestep %d",n);
              	//lmp->input->one(str);
            }

            // Accept trial move    
            if(acc)
            {
            // Update number of atoms of each type
		        if (flag==1) {
            	    elmA[t1-1]--;
	    		    elmB[t1-1]++;
            	    elmA[t2-1]++;
	    		    elmB[t2-1]--;
		        }

            	for(int i=0;i<ninstance;i++) e0s[i] = all_e0s[i] = 0.0;
            	if(me_lammps==0) e0s[instance] = ec;
            	MPI_Reduce(e0s,all_e0s,ninstance,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            	nrej=0;
            } else 
            { // Trial move rejected, return system to previous state
		        if (flag==1) 
                {
            		tA[id1] = t1;
            		tB[id2] = t2;
		        } else if (flag==2) 
                {
			        tA[id1] = t1;
			        tA[id2] = t2;
		        } else 
                {
			        tB[id1] = t1;
			        tB[id2] = t2;
		        }
            	if(instance==0 && flag!=3) lammps_scatter_atoms(lmp,(char *) "type",0,1,tA);
            	if(instance==1 && flag!=2) lammps_scatter_atoms(lmp,(char *) "type",0,1,tB);
            	nrej++;
            }
        }
    } // for n

    if(me==0) fprintf(stderr,">>>MC iteration finished\n");

    delete [] xA;
    delete [] xB;
    delete [] xA_start;
    delete [] xB_start;
    delete [] vB;
    delete [] vA;
    delete [] vB_start;
    delete [] vA_start;
    delete [] tA;
    delete [] tB;
    delete [] all_es;
    delete [] es;
    delete [] all_e0s;
    delete [] e0s;

  // delete LAMMPS instances
    delete lmp;

  // close down MPI

    MPI_Comm_free(&comm_lammps);
    MPI_Finalize();

    //if(me==0) fprintf(stderr,">>>Bye!\n");

    return EXIT_SUCCESS;
}
