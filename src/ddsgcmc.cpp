#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"

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
        printf("Syntax: npar N in.lammps mc.ctl rnd_seed\n");
        exit(1);
    }

    char *infile = arg[2];
    char *mc_ctl = arg[3];
    char *rnd_seed = arg[4];

    // Initial settings here
    MCVars mc;
    try
    {
        mc.read_ctl(mc_ctl);
    }
    catch (char* str)
    {
        fprintf(stderr,"%s\n",str);
        exit(1);
    }

    // database allocation
    try
    {
        mc.alloc_db();
    }
    catch (char* str)
    {
        fprintf(stderr,"%s\n",str);
        return EXIT_FAILURE;
    }
fprintf(stderr,"Database loaded\n");

    MPI_Init(&narg,&arg);
    int me,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if(me==0)
    {
        fprintf(stdout,"### Atom-swapping Monte Carlo simulator ###\n");
        fprintf(stdout,"### Temperature = %.1f K\n",mc.temp);

        fprintf(stdout,"### Concentration = ");
		for(int i=0;i<mc.ntype;i++) fprintf(stdout,"%16.8f",mc.clist[i]);
		fprintf(stdout,"\n");

        fprintf(stdout,"### Iteration = %d cycles\n",mc.niter);
        fprintf(stdout,"### Resume from %d steps\n",mc.nstart);
        fprintf(stdout,"### Output interval = %d steps\n",mc.nout);
        fprintf(stdout,"### Rejection limit = %d cycles\n",mc.rej_limit);
        fprintf(stdout,"### Boltzmann factor = %g\n",mc.beta);
        fprintf(stdout,"###########################################\n");

        if(mc.LCOout)
        {
            mc.gen_LCOscr();
            fprintf(stdout,"# LCO.mod was generated.\n");
        }
    }


    int ninstance = 1;
    if (nprocs % ninstance) {
        if (me == 0)
            printf("ERROR: Total procs must be divisble by N\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    // header preparation
    if(me==0)
    {
        FILE *fp;
        try{
            if((fp=fopen("log.mc","w"))==NULL) throw 1;
            fprintf(fp,"#   MC_step    E               E*              dEtar       ID    #Rej      t1        t2\n");
            fclose(fp);
            fp=NULL;
        }
        catch(int e)
        {
            // exception: output file cannot open
            fprintf(stderr,"Error in output (header)\n");
            MPI_Finalize();
            return EXIT_FAILURE;
        }

        // LCO output
        if(mc.LCOout)
        {
            try{
                if((fp=fopen("log.LCO","w"))==NULL) throw 2;
                fprintf(fp,"#   MC_Step  ");
                for(int i=0;i<mc.ntype;++i)
                {
                    for(int j=i;j<mc.ntype;++j)
                    {
                        fprintf(fp,"     a%d%d    ",i+1,j+1);
                    }
                }
                fprintf(fp,"\n");
                fclose(fp);
                fp=NULL;
            }catch(int e)
            {
                fprintf(stderr,"Error in log.LCO(header)\n");
                MPI_Finalize();
                return EXIT_FAILURE;
            }
        }

    }

    // create one communicator per instance each with P/N procs
    int instance = me*ninstance / nprocs;
    int npar = nprocs/ninstance;
    int me_lammps = me%npar;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,instance,me_lammps,&comm_lammps);
  
    // each instance: unique screen file, log file, temperature

    char str1[32],str2[32],str3[32];

    char **lmparg = new char*[5];
    lmparg[0] = NULL;                 // required placeholder for program name
    lmparg[1] = (char *) "-screen";
    sprintf(str1,"screen.%d",instance);
    lmparg[2] = str1;
    lmparg[3] = (char *) "-log";
    sprintf(str2,"log.lammps.%d",instance);
    lmparg[4] = str2;

    // open N instances of LAMMPS

    LAMMPS *lmp = new LAMMPS(5,lmparg,comm_lammps);
    fprintf(stderr,"%d(%d@%d)>>>each object prepared\n",me,me_lammps,instance);

    delete [] lmparg;

    // run input script through all instances of LAMMPS (initial settings)
    lammps_file(lmp,infile);
    fprintf(stderr,"%d(%d@%d)>>>preprocess settings done\n",me,me_lammps,instance);

    if(mc.LCOout) lmp->input->one("include 'LCO.mod'");



    // settings before MC simulation
    int nrej=0;
    int win_ins,win_rank;
    int id1,t1,t2;
    double r_t = 0.0;
    double dErsv = 0.0;
    int acc = 0;
    double p = 0;
    int natoms = static_cast<int> (lmp->atom->natoms);
    double *x = new double[3*natoms];
    double *v = new double[3*natoms];
    int    *t = new int[natoms];
    double *x0 = new double[3*natoms];
    int    *t0 = new int[natoms];
    int cnt = 0;
    // random variables settings
    //std::random_device rnd; // device dependent random variable
    unsigned rnd = static_cast<unsigned>(*rnd_seed);
    unsigned seed = rnd;
    std::mt19937 mt(seed);
    std::mt19937 mt_prob(seed+1);
    std::uniform_int_distribution<int> sid(0,natoms-1);
    std::uniform_real_distribution<double> prob(0.0,1.0);
    
    lammps_gather_atoms(lmp,(char *) "x",1,3,x);
    lammps_gather_atoms(lmp,(char *) "v",1,3,v);
//------------------- Insert below values of chemical potentials when using traditional SGC ----------------
    double *tmp_mu;
    double *mu = new double[mc.ntype];
    char mu_str[64];
    for (int i=0; i<mc.ntype; ++i)
    {
        sprintf(mu_str,"mu%d",i+1);
        tmp_mu = (double *) lammps_extract_variable(lmp, (char *) mu_str, NULL);
        mu[i] = *tmp_mu;
        tmp_mu = NULL;
        if(me==0) fprintf(stdout,"Chemical potential of atom type %d: %f\n",i,tmp_mu);
    }
    delete [] tmp_mu;
 
    lammps_gather_atoms(lmp,(char *) "type",0,1,t);
    memcpy(t0,t,natoms*sizeof(int));

    // number of atoms of each element
    int *elm = new int[mc.ntype];
    for(int i=0;i<mc.ntype;i++) elm[i] = 0;
    for(int i=0;i<natoms;i++) {
	elm[t[i]-1]++;
    }
    // composition info
    if(me==0)
    {
        FILE *fp;
        try{
            if((fp=fopen("composition.dat","w"))==NULL) throw 1;
            fprintf(fp,"#   Iter   ");
            for(int i=0;i<mc.ntype;i++) fprintf(fp,"   elem%d",i+1);
            fprintf(fp,"\n");
            fprintf(fp,"%8d",0);
            for(int i=0;i<mc.ntype;i++) fprintf(fp,"%8d",elm[i]);
            fprintf(fp,"\n");
            fclose(fp);
            fp=NULL;
        }
        catch(int e)
        {
            fprintf(stderr,"Error in composition output before run.\n");
            delete [] x;
            delete [] t;
            delete [] elm;
            delete lmp;
            MPI_Comm_free(&comm_lammps);
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    }

    // for sgcmc
    std::uniform_int_distribution<int> stype(1,mc.ntype);


    char command1[256],command2[256];

    lmp->input->one("thermo_style custom step vol press temp pe");
    lmp->input->one("thermo 100");
    lmp->input->one("run 0");

    if(!isRestart)
    {
        sprintf(command1,"velocity all create %.1f %5d rot yes mom yes dist gaussian",mc.temp,46242);
        lmp->input->one(command1);

        if(mc.ensemble == 1)
        {
            double *press1 = (double *) lammps_extract_compute(lmp,(char *) "thermo_press",0,0);
            double p0 = *press1;
            if(me==0) fprintf(stderr,"Initial pressure: %.2f\n",p0);
            press1=NULL;
            sprintf(command2,"fix nf all npt temp %.1f %.1f 0.1 iso %.2f %.1f 1.0",mc.temp,mc.temp,p0,0.0);
        }else
        {
            sprintf(command2,"fix nf all nvt temp %.1f %.1f 1",mc.temp,mc.temp); // Needs to be the same setting as the one used during DB construction
        }
        lmp->input->one(command2);

        lmp->input->one("fix mom all momentum 50 linear 1 1 1 angular rescale");
        sprintf(command1,"run %d",mc.nMDinit);
        lmp->input->one(command1);
        lmp->input->one("unfix nf");
    }

    char restart_command[128];
    if(mc.ensemble == 1)
    {
        sprintf(restart_command,"fix nf all npt temp %.1f %.1f 0.1 iso 0.0 0.0 1.0",mc.temp,mc.temp);
    }else
    {
        sprintf(restart_command,"fix nf all nvt temp %.1f %.1f 0.1",mc.temp,mc.temp);
    }
    lmp->input->one(restart_command);

    lmp->input->one("compute calcPE all pe");
    lammps_gather_atoms(lmp,(char *) "x",1,3,x);
    lammps_gather_atoms(lmp,(char *) "v",1,3,v);


    // Starts MC loop here
    if(me==0) fprintf(stderr,">>>Start\n");
    for(int n=mc.nstart;n<mc.niter;n++)
    {
    //if(me==0) fprintf(stderr,">>>>Iteration %d\n",n);
        if(nrej>mc.rej_limit)
        {
            if(me==0) fprintf(stderr,"Reached rejection limit\n");
            break;
        }
//-------------------------------------- Continue MD run ----------------------------------------
	    lammps_scatter_atoms(lmp,(char *) "v",1,3,v);
        sprintf(command1,"run %d",mc.nMDmc);	
	    lmp->input->one(command1);
	    lammps_gather_atoms(lmp,(char *) "v",1,3,v);
//-----------------------------------------------------------------------------------------------
        
//------------------------ Get energy immediately prior to the MC cycle -------------------------
	    double *ptr = (double *) lammps_extract_compute(lmp,(char *) "calcPE",0,0);
        double e0 = *ptr;
        ptr=NULL;
//----------------------------- MC cycle (# of steps = # of atoms) ------------------------------
        for(int j=0;j<natoms;j++)
        {
            if(me==0) fprintf(stderr,"Checkpoint 1\n");
	        cnt++;
            // id selection
	        t1=t2=0;

       	    while(t1==t2)
       	    {
		        id1 = sid(mt);
		        t1 = t[id1];
		        t2=0;
		        r_t=prob(mt);
		        
                while(t2<mc.ntype)
		        {
			        if(mc.ccum[t2]>r_t) break;
			        t2++;
		        }
		        // convert to LAMMPS atom type
         	    t2++;		
            }
            if(me==0) fprintf(stderr,"Checkpoint 2\n");
            // atom type swapping (with reservoir)
       	    t[id1] = t2;
            if(me==0) fprintf(stderr,"%d\n",t2);
            lammps_scatter_atoms(lmp,(char *) "type",0,1,t);	 
	        lmp->input->one("run 0");
            
            // Structural relaxation (possibly unnecessary, depends on the system, and it drastically reduces speed)
	        /*
            lmp->input->one("min_style cg");
            lmp->input->one("minimize 0.0 0.01 10000 10000");
            lmp->input->one("run 0");
            */
            ptr = (double *) lammps_extract_compute(lmp,(char *) "calcPE",0,0);
            double ec = *ptr;
            ptr=NULL;

      	    // calc acceptance ratio
      	    if (me==0)
      	    {
                // target energy difference
          	    double de = ec-e0;
          	    // reservoir energy difference
          	    // BE AWARE!!!: DB gives the energy from the point of view of reservoir
          	    double rp = prob(mt_prob);
                dErsv = mc.reservoir_energy(t2-1,t1-1,rp); // DD-SGC acceptance rule
          	    //dErsv = mu[t1-1]-mu[t2-1]; // SGC acceptance rule 
                de += dErsv;
          	    p  = (de < 0.0) ? 1.0 : exp(-de*mc.beta);
          	    double rt = prob(mt_prob);
	  	        acc = (p>rt) ? 1 : 0;  
      	    }
      	    MPI_Bcast(&acc,1,MPI_INT,0,MPI_COMM_WORLD);
            if(me==0) fprintf(stderr,"Checkpoint 4\n");
        
	        if(!(cnt%mc.nout))
            {
        	    if(instance==0) 
                {
                	char dumpstr[BUFSIZ];
                	char data[32]="mass type xs ys zs";
                	sprintf(dumpstr,"write_dump all cfg ./dump/mc%07d.cfg.gz %s %s %s",cnt,data,mc.mod1,mc.mod2);
                	lmp->input->one(dumpstr);
              	}
      	    }
            
            if(mc.nrestart>0 && !(cnt%mc.nrestart))
            {
                if(instance==0)
                {
                    char dumpstr[BUFSIZ];
                    sprintf(dumpstr,"write_restart MCrestart.%d",cnt);
                    lmp->input->one(dumpstr);
                }
            }

            // if accept, update structure
            if(acc==1)
	        {
                memcpy(t0,t,natoms*sizeof(int));
	            elm[t1-1]--;
                elm[t2-1]++;
                e0 = ec;
                nrej=0;
            } else 
            {
                // if reject, revert structure
                memcpy(t,t0,natoms*sizeof(int));
                lammps_scatter_atoms(lmp,(char *) "type",0,1,t);
                nrej++;
            }
	
	        if(me==0 && cnt%10==0) 
            {
                FILE *fp;
      		    try
                {
      			    if((fp=fopen("composition.dat","a"))==NULL) throw 1;
            		fprintf(fp,"%8d",n);
            		for(int i=0;i<mc.ntype;i++) fprintf(fp,"%8d",elm[i]);
            		fprintf(fp,"\n");
             		fclose(fp);
             		fp=NULL;
     		    }
  		        catch(int e)
                {
           		    fprintf(stderr,"Error in composition output\n");
        	   	    break;
  		        }
	        }
          
            if(me==0 && cnt%10==0)
            {
                FILE *fp;
                // energy info
                try
                {
                    if((fp=fopen("log.mc","a"))==NULL) throw 1;
                    double de=ec-e0;
                    fprintf(fp,"%8d%16.8e%16.8e%16.8e%4d%8d%8d%8d\n"
                                  ,cnt,e0,ec,de,nrej,id1,t1,t2);
                    fclose(fp);
                    fp=NULL;
                }
                catch(int e)
                {
                    // exception: output file cannot open
                    fprintf(stderr,"Error in output\n");
                    break;
                }
	        }

            if(mc.LCOout && me==0 && cnt%10==0)
            {
                FILE *fp;
                try
                {
                    if((fp=fopen("log.LCO","a"))==NULL) throw 1;
                    fprintf(fp,"%8d",cnt);
                    for(int i=0;i<mc.ntype;++i)
                    {
                        for(int j=i;j<mc.ntype;++j)
                        {
                            char *alpha;
                            sprintf(alpha,"a%d%d",i+1,j+1);
                            ptr = (double *) lammps_extract_variable(lmp,alpha,0,0);
                            double a_ij = *ptr;
                            ptr = NULL;
                            fprintf(fp,"%12.4f",a_ij);
                        }
                    }
                    fprintf(fp,"\n");
                    fclose(fp);
                    fp = NULL;
                }
                catch(int e)
                {
                    fprintf(stderr,"Error in log.LCO\n");
                    break;
                }   
            }

        }	
    } // for n

    if(me==0) fprintf(stderr,">>>MC iteration finished\n");

    delete [] x;
    delete [] t;
    delete [] x0;
    delete [] t0;
    // delete LAMMPS instances
    delete lmp;
    delete mu;    
    // close down MPI
    MPI_Comm_free(&comm_lammps);
    MPI_Finalize();

    return EXIT_SUCCESS;
}
