#ifndef __MCVARS_H__
#define __MCVARS_H__

#include "stdlib.h"
#include "string.h"
#include <vector>
#include <math.h>
#define MAXARG 10

class MCVars {
    private:
        const double kb = 8.617333262e-5;
        //const int Ncum = 100; // readable value in the future
        //int Ncum;
        int split_str(char str[],char **dcp);
    public:
        MCVars();
        ~MCVars();
        double temp;
        double c0;
        int nstart,niter,nout,rej_limit,ntype;
        int nMDinit,nMDmc;
        double beta;
        char *mod1;  // label for dump command
        char *mod2;  // label for dump command
        char *str;   // useful string
        double *clist,*ccum; // concentration array

        void read_ctl(char *ctlfile);   // read controling parameters

        // for DB-GCMC
        int Ncum;
        double *Pcum; // cumulative function
        double *Ecum; // corresponding energy difference
        char *dbfile; // file name of database
        void alloc_db();       // allocate database array
        double reservoir_energy(int t1, int t2, double r);
};

MCVars::MCVars()
{
    temp=0.0;
    c0=0.0;
    nstart=0;
    nout=1;
    niter=0;
    rej_limit=1000;
    ntype=0;
    Ncum=0;
    nMDinit = 1000;
    nMDmc = 500;
    str = new char[BUFSIZ];

    dbfile = new char[BUFSIZ];
    Ecum=NULL;
    Pcum=NULL;
}
MCVars::~MCVars()
{
    if(mod1!=NULL) delete [] mod1;
    if(mod2!=NULL) delete [] mod2;
    if(str !=NULL) delete [] str;

    if(Ecum  !=NULL) delete [] Ecum;
    if(Pcum  !=NULL) delete [] Pcum;
    if(dbfile!=NULL) delete [] dbfile;
    if(clist !=NULL) delete [] clist;
    if(ccum  !=NULL) delete [] ccum;
}

void MCVars::read_ctl(char* ctlfile)
{
    int nMDinit_default = nMDinit;
    int nMDmc_default = nMDmc;

    FILE *fp;

    if((fp=fopen(ctlfile,"r"))==NULL)
    {
        sprintf(str,"%s can not open.\n",ctlfile);
        throw str;
    }

    char *cstr;
    cstr = new char[BUFSIZ];
    mod1=new char[BUFSIZ];
    mod2=new char[BUFSIZ];

    while(fgets(str,BUFSIZ,fp)!=NULL){
        if(strstr(str,"TEMP") !=NULL) sscanf(str,"%*s %*s %lf",&temp);
		// need modification for C0 reader to read multiple times
		// divide string to array 
        //if(strstr(str,"C0"  ) !=NULL) sscanf(str,"%*s %*s %lf",&c0);
        if(strstr(str,"C0"  ) !=NULL) strcpy(cstr,str);
		//

        if(strstr(str,"NINI") !=NULL) sscanf(str,"%*s %*s %d", &nstart);
        if(strstr(str,"NITR") !=NULL) sscanf(str,"%*s %*s %d", &niter);
        if(strstr(str,"NOUT") !=NULL) sscanf(str,"%*s %*s %d", &nout);
        if(strstr(str,"LREJ") !=NULL) sscanf(str,"%*s %*s %d", &rej_limit);
        if(strstr(str,"NELM") !=NULL) sscanf(str,"%*s %*s %d", &ntype);
        if(strstr(str,"NDBS") !=NULL) sscanf(str,"%*s %*s %d", &Ncum);
        if(strstr(str,"NMD0") !=NULL) sscanf(str,"%*s %*s %d", &nMDinit);
        if(strstr(str,"NMD1") !=NULL) sscanf(str,"%*s %*s %d", &nMDmc);

//        if(strstr(str,"element") !=NULL) sscanf(str,"%s",mod2);
        if(strstr(str,"element") !=NULL) memcpy(mod2,str,BUFSIZ);
    }
    fclose(fp);

    // Temperature error
    if(!(temp>0.0))
    {
        sprintf(str,"Temperature must be positive / temp=%f\n",temp);
        throw str;
    }

    // Iteration error
    if(!(niter>0))
    {
        sprintf(str,"Niter must be greater than 0\n");
        throw str;
    }

    // Element error
    if(ntype<2)
    {
        sprintf(str,"Number of atomic species must be >1\n");
        throw str;
    }

    clist = new double[ntype];

    char **dcp;
    dcp = new char*[MAXARG];
    for(int i=0;i<MAXARG;i++) dcp[i] = new char[32];

    int narg = split_str(cstr,dcp);
    
	if(narg-2!=ntype-1){
		sprintf(str,"number of types not match\n");
        delete [] dcp;
        delete [] cstr;
		throw str;
	}

    // MD steps error
    if(nMDinit<0)
    {
        fprintf(stderr,"Warning: Number of MD step must be 0 or positive. Using default NMD0 value...\n");
        nMDinit = nMDinit_default;
    }
    if(nMDmc<0)
    {
        fprintf(stderr,"Warning: Number of MD step must be 0 or positive. Using default NMD1 value...\n");
        nMDmc = nMDmc_default;
    }

    // set solute concentration
    for(int i=0;i<ntype-1;i++) sscanf(dcp[i+2],"%lf",&clist[i+1]);
    delete [] dcp;
    delete [] cstr;

    // base material
	clist[0] = 1.0;
	for(int i=1;i<ntype;i++) clist[0]-=clist[i];
	//
	ccum = new double[ntype];
	for(int i=0;i<ntype;i++) ccum[i]=0.0;
	ccum[0]=clist[0];
	for(int i=1;i<ntype;i++) ccum[i]+=(ccum[i-1]+clist[i]);

    beta = 1.0/(temp*kb);

    int npad=static_cast<int>(log10(niter))+1;
    sprintf(mod1,"modify sort id pad %d",npad);

}

void MCVars::alloc_db()
{
    // type error
    if(ntype==0)
    {
        sprintf(str,"Allocation error: number of types must be >0\n");
        throw str;
    }
    
    // Database error
    if(Ncum==0)
    {
        sprintf(str,"Number of sample in database should be >0\n");
        throw str;
    }
    // I need to make sure how to specify database files from total element number
    int dbsz = ntype*ntype*Ncum;
    Ecum = new double[dbsz];
    Pcum = new double[dbsz];

    for(int i=0;i<dbsz;i++) Ecum[i]=0.0;
    for(int i=0;i<dbsz;i++) Pcum[i]=0.0;

    double cum,de;

    for(int t1=0;t1<ntype;t1++)
    {
        for(int t2=0;t2<ntype;t2++)
        {
            if(t2==t1) continue;
            sprintf(dbfile,"dErsv_%dto%d.dat",t1+1,t2+1);
            int idx=ntype*Ncum*t1+Ncum*t2;
            int l=0;
            FILE *fp=fopen(dbfile,"r");
            if(fp==NULL)
            {
                sprintf(str,"Database %s can not open.\n",dbfile);
                throw str;
            }
            while(fgets(str,BUFSIZ,fp)!=NULL){
                sscanf(str,"%lf %*lf %lf",&de,&cum);
                Ecum[idx+l] = de;
                Pcum[idx+l] = cum;

                if(++l>Ncum)
                {
                    sprintf(str,"Database %s: No. of data not match.\n",dbfile);
                    fclose(fp);
                    throw str;
                    break; // just in case
                }
            }
            fclose(fp);
        }
    }

}

double MCVars::reservoir_energy(int t1, int t2, double r)
{
    double dErsv = 0.0;

    // here energy difference in reservoir is chosen
    int k=ntype*Ncum*t1+Ncum*t2;

    for(int i=0;i<Ncum;i++)
    {
        if(Pcum[++k]>r) break;
    }

    // back just before
    --k;

    dErsv=Ecum[k];

    return dErsv;
}

int MCVars::split_str(char str[], char **iargv){
	int tc;
	char *tp;
	char *dummy;

	dummy = new char[BUFSIZ];
	strcpy(dummy,str);

	for(tc=0;tc<MAXARG;tc++){
		tp=strtok((tc==0)?dummy:NULL," \t");
		if(tp!=NULL&&((strncmp(tp,"#",1))!=0)) strcpy(iargv[tc],tp);
		else break;
	}

	delete [] dummy;

	return tc;
}

#endif
