#include "statement.h"
#include "mpi.h"
#include "time.h"
#define pml 10


extern void input_ac2d_parameters(int *shotnum,int *tmax,float *dt,float *fd,float *fo,
                                  int *nx,int *nz,float *dx,float *dz,float *disx,
                                  char fname1[],char fname2[],char fname3[],
                                  char fname4[],char fname5[],
                                  char vpname[],char rhoname[],char Qpname[],
				  int *flagout,int *flagsource);

extern void input_vp_model(int,int,float *,char vpname[]);

extern void cal_bbx_bbz(int nxt,int nzt,float *rhot,float *bbxt,float *bbzt);


extern void rik_cal(int tmax,float dt,float fd,float *rik,int);

extern void get_pml_parameters(int,float fd,float *i_ax,float *i_bx,float *i_az,float *i_bz,
				float *h_ax,float *h_bx,float *h_az,float *h_bz,int nxx,int nzz,
				float dx,float dz,float dt);
				
extern void cal_shift_parameters(float *kxsmr,float *kxscr,float *kzsmr,float *kzscr,
                          float *kxsmi,float *kxsci,float *kzsmi,float *kzsci,
                          float dx,float dz,int nxx,int nzz);


extern void get_max_min(int NXZ,float *vp,float *vmax,float *vmin);

extern float check_stability_ps(float vmax,float dt,float dx,float dz);

extern void cal_k2d(float dx,float dz,int n0,int n1,float *k2d1,float *k2d2,float *k2d3,float *k2da);


extern void input_geophone_pos(int shotid,int nr,float *gxp,char fname[]);


extern int getmaxtracenumber(int N,int *a);

extern void transfermodel(
				int nx,int nz,float fd,float fo,
				float *vp,float *Qp,float *gam,float *lam,
				float *yinda,float *tao
			 );



						   
extern  void acforward(
				int snapnum,int GPU_N,int tmax,float dt,float *rik,int myid,int numprocs,
			  	int nshot,struct MultiGPU singpu[],float *gxp,float *gzp,
			  	int nrmx,int nx,int nz,int nxx,int nzz,
		          	int pmll,float dx,float dz,float *k2da,
				float *k2d1,float *k2d2,float *k2d3,
		          	float sroffmx,float disx,
		          	float *vpt,float *vwdt,float *taot,float *lamt,float *gamt,
		          	float *rhot,float *bbxt,float *bbzt,
			  	int nxt,int nzt,int *tracenum,int flagout,
			  	float *i_ax,float *i_bx,float *i_az,float *i_bz,
			 	float *h_ax,float *h_bx,float *h_az,float *h_bz,
			  	float *kxsmr,float *kxscr,float *kzsmr,float *kzscr,
			  	float *kxsmi,float *kxsci,float *kzsmi,float *kzsci,
                          	char fname4[40],char fname5[40]
                          );
                               
                               
                               

int main(int argc,char*argv[])
{
    int myid,numprocs;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Get_processor_name(processor_name,&namelen);


    int GPU_N;
    int devID;
    cudaGetDevice(&devID);
    cudaGetDeviceCount(&GPU_N);
    printf(" Process %d of %d on %s with %d GPUs running... \n\n ",myid,numprocs,processor_name,GPU_N);


    time_t begin_time, end_time;
    clock_t start, end;
    float runtime=0.0;
	
    

	
	// information of source
	int shotnum;
        int tmax;
        float dt;
        float fd;
        float fo;
        
        int nxt;
        int nzt;
        float dx;
        float dz;
        float disx;
       
        int *tracenum;
    
	float *shotxp, *shotzp;
        
	char vpname[80];
	char rhoname[80];
	//char recname[60];
	char Qpname[80];
	
        float *vpt;
        float *rhot;
        float *bbxt;
        float *bbzt;
        float *Qpt;
        float *lamt;
        float *gamt;
        float *vwdt;
        float *taot;
        
        
        char fname1[40];
        char fname2[40];
        char fname3[40];
        char fname4[40];
        char fname5[40];

	int snapnum = 20;
        
        
	int flagsource;
	int flagout;
	float vmax,vmin,delta;

	if(myid==0)
	{
	    /* input simulation parameters */
	    input_ac2d_parameters(&shotnum,&tmax,&dt,&fd,&fo,&nxt,&nzt,&dx,&dz,&disx,
                                  fname1,fname2,fname3,fname4,fname5,vpname,rhoname,Qpname,
                                  &flagout,&flagsource);
                                  
	    if(shotnum<GPU_N){printf("Warning! Shotnum is less than GPU_N !\n");}
	}

        MPI_Bcast(&shotnum,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&tmax,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&dt,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        MPI_Bcast(&fd,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Bcast(&fo,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        MPI_Bcast(&nxt,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&nzt,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&dx,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Bcast(&dz,1,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Bcast(&disx,1,MPI_FLOAT,0,MPI_COMM_WORLD);


        MPI_Bcast(&fname1,40,MPI_CHAR,0,MPI_COMM_WORLD);    // trace number filename in each shot
        MPI_Bcast(&fname2,40,MPI_CHAR,0,MPI_COMM_WORLD);    // shot x position filename
        MPI_Bcast(&fname3,40,MPI_CHAR,0,MPI_COMM_WORLD);    // shot z position filename
        MPI_Bcast(&fname4,40,MPI_CHAR,0,MPI_COMM_WORLD);    // receiver x position prefix filename
        MPI_Bcast(&fname5,40,MPI_CHAR,0,MPI_COMM_WORLD);    // receiver z position prefix filename

        MPI_Bcast(vpname,80,MPI_CHAR,0,MPI_COMM_WORLD);     // true vp filename
        MPI_Bcast(rhoname,80,MPI_CHAR,0,MPI_COMM_WORLD);    // true density filename
        MPI_Bcast(Qpname,80,MPI_CHAR,0,MPI_COMM_WORLD);     // true Qp filename

        //MPI_Bcast(recname,60,MPI_CHAR,0,MPI_COMM_WORLD);    // observed binary data filename
        MPI_Bcast(&flagsource,1,MPI_INT,0,MPI_COMM_WORLD);  // ricker source or not
	MPI_Bcast(&flagout,1,MPI_INT,0,MPI_COMM_WORLD);     // ricker source or not

	
	if(shotnum<GPU_N){GPU_N=shotnum;}

	//GPU_N=1;

        int nzxt;
        nzxt=nxt*nzt;
        dt=dt/1000.0;

        tracenum=(int*)malloc(sizeof(int)*shotnum);
        shotxp=(float*)malloc(sizeof(float)*shotnum);
        shotzp=(float*)malloc(sizeof(float)*shotnum);
        
        vpt=(float*)malloc(sizeof(float)*nxt*nzt);
        Qpt=(float*)malloc(sizeof(float)*nxt*nzt);

        gamt=(float*)malloc(sizeof(float)*nxt*nzt);
        lamt=(float*)malloc(sizeof(float)*nxt*nzt);
        taot=(float*)malloc(sizeof(float)*nxt*nzt);
	vwdt=(float*)malloc(sizeof(float)*nxt*nzt);
        
        
        rhot=(float*)malloc(sizeof(float)*nxt*nzt);
        bbxt=(float*)malloc(sizeof(float)*nxt*nzt);
        bbzt=(float*)malloc(sizeof(float)*nxt*nzt);
        

        FILE *fp1,*fp2,*fp3;
        if(myid==0)
        {
		fp1=fopen(fname1,"rb");
                if(fp1==0){printf("Cannot open the tracenumber file !\n");exit(0);}
                fread(tracenum,sizeof(int),shotnum,fp1);
                fclose(fp1);
                
                fp2=fopen(fname2,"rb");
                if(fp2==0){printf("Cannot open the shot x position file !\n");exit(0);}
                fread(shotxp,sizeof(float),shotnum,fp2);
                fclose(fp2);

                fp3=fopen(fname3,"rb");
                if(fp3==0){printf("Cannot open the shot z position file !\n");exit(0);}
                fread(shotzp,sizeof(float),shotnum,fp3);
                fclose(fp3);

                input_vp_model(nxt,nzt,vpt,vpname);
                input_vp_model(nxt,nzt,rhot,rhoname);
                input_vp_model(nxt,nzt,Qpt,Qpname);
                
                get_max_min(nzxt,vpt,&vmax,&vmin);
                delta=check_stability_ps(vmax,dt,dx,dz);
                printf(" CFL number= %f %f %f\n",delta,vmax,vmin);
                //if(delta-1.0>1e-4){printf(" Error: Unstable simulation !\n");exit(0);}
                
                transfermodel(nxt,nzt,fd,fo,vpt,Qpt,gamt,lamt,vwdt,taot);
         }

         	MPI_Bcast(tracenum,shotnum,MPI_INT,0,MPI_COMM_WORLD);
         	MPI_Bcast(shotxp,shotnum,MPI_FLOAT,0,MPI_COMM_WORLD);
         	MPI_Bcast(shotzp,shotnum,MPI_FLOAT,0,MPI_COMM_WORLD);

		MPI_Bcast(vpt,nzxt,MPI_FLOAT,0,MPI_COMM_WORLD); 
   		MPI_Bcast(rhot,nzxt,MPI_FLOAT,0,MPI_COMM_WORLD);   		
   		MPI_Bcast(vwdt,nzxt,MPI_FLOAT,0,MPI_COMM_WORLD);
         	MPI_Bcast(taot,nzxt,MPI_FLOAT,0,MPI_COMM_WORLD);
         	MPI_Bcast(gamt,nzxt,MPI_FLOAT,0,MPI_COMM_WORLD);
         	MPI_Bcast(lamt,nzxt,MPI_FLOAT,0,MPI_COMM_WORLD);
   		
   		
   		
   		cal_bbx_bbz(nxt,nzt,rhot,bbxt,bbzt);

		 // calculate source //
		 float *rik;
		 rik=(float*)malloc(sizeof(float)*tmax);
	         rik_cal(tmax,dt,fd,rik,flagsource);
         

		 float xmin=0.0;
		 float xmax=(nxt-1)*dx;
		 float zmin=0.0;
		 float zmax=(nzt-1)*dz;

		
       
		 // distribute shots to each GPU //
		 int ave,rid,nshot,ridg;
		 int offsp;
		 rid=shotnum%(numprocs*GPU_N);
		 ave=(shotnum-rid)/(numprocs*GPU_N);
		 ridg=rid/GPU_N;
		 if(myid<ridg)
		 	{
				nshot=ave+1;
				offsp=myid*(ave+1)*GPU_N;
			}
		 else
		 	{	
				nshot=ave;
				offsp=ridg*(ave+1)*GPU_N+(myid-ridg)*ave*GPU_N;
			}

		 //printf("myid=%d nshot=%d off=%d\n",myid,nshot,offsp);
		 
		 if(myid==numprocs-1){printf("Run total shot number is %d\n",offsp+nshot*GPU_N);}



		 int i,j,is;
		 struct MultiGPU singpu[GPU_N];
		 for(i=0;i<GPU_N;i++)
		 {
		 	singpu[i].shotid=(int*)malloc(sizeof(int)*nshot);   // shot id at each GPU
			singpu[i].sxp=(float*)malloc(sizeof(float)*nshot);  // shot x position
			singpu[i].szp=(float*)malloc(sizeof(float)*nshot);  // shot z position 
			
			for(j=0;j<nshot;j++)
			{
				singpu[i].shotid[j]=offsp+i+j*GPU_N;
				is=singpu[i].shotid[j];

				singpu[i].sxp[j]=shotxp[is];
				singpu[i].szp[j]=shotzp[is];

				if(shotxp[is]<xmin||shotxp[is]>xmax){printf("%d shot x exceeds model!\n",is);exit(0);}
			 	if(shotzp[is]<zmin||shotzp[is]>zmax){printf("%d shot z exceeds model!\n",is);exit(0);}
			}
		 }


		 int nrmx,nr;
		 float xlp,xrp;
		 float zlp,zrp;
		 float temp1,temp2;
		 float sroffset,sroffmx=-1.0;
                 float *gxp,*gzp;

		 // maximum trace number cross all shots on each cpu node // 
		 nrmx=getmaxtracenumber(shotnum,tracenum);
		 gxp=(float*)malloc(sizeof(float)*nrmx);
                 gzp=(float*)malloc(sizeof(float)*nrmx);

		 
		 
		 for(i=offsp;i<offsp+nshot*GPU_N;i++)
		 {
			nr=tracenum[i];
			input_geophone_pos(i,nr,gxp,fname4);
			input_geophone_pos(i,nr,gzp,fname5);

			get_max_min(nr,gxp,&xrp,&xlp);
			get_max_min(nr,gzp,&zrp,&zlp);
			if(xlp<xmin||xrp>xmax){printf("%d shot receiver x exceeds model!\n",i);exit(0);}
			if(zlp<zmin||zrp>zmax){printf("%d shot receiver z exceeds model!\n",i);exit(0);}
          
			temp1=fabs(xrp-shotxp[i]);
			temp2=fabs(shotxp[i]-xlp);
			sroffset=(temp1>temp2)?temp1:temp2;
			if(sroffset>sroffmx){sroffmx=sroffset;}
			
		 }
                 
		 // maximum grid number for all shots on each cpu node //
                 float simxlen=2.0*(sroffmx+disx);
                 int nx=(int)(simxlen/dx+1.5);
		 int nz=nzt;
		 int nxx=nx+2*pml;
		 int nzz=nz+2*pml;

		 // allocate host memory //
		 for(i=0;i<GPU_N;i++)
		 {
		 	singpu[i].h_rho=(float*)malloc(sizeof(float)*nxx*nzz);
		 	singpu[i].h_bbx=(float*)malloc(sizeof(float)*nxx*nzz);
		 	singpu[i].h_bbz=(float*)malloc(sizeof(float)*nxx*nzz);
		 	
			singpu[i].h_vp=(float*)malloc(sizeof(float)*nxx*nzz);
		 	singpu[i].h_vwd=(float*)malloc(sizeof(float)*nxx*nzz);
		 	singpu[i].h_tao=(float*)malloc(sizeof(float)*nxx*nzz);
		 	singpu[i].h_gam=(float*)malloc(sizeof(float)*nxx*nzz);
		 	singpu[i].h_lam=(float*)malloc(sizeof(float)*nxx*nzz);
			
		 	
			singpu[i].gxg=(int*)malloc(sizeof(int)*nrmx);
			singpu[i].gzg=(int*)malloc(sizeof(int)*nrmx);
			
			singpu[i].h_obs=(float*)malloc(sizeof(float)*nrmx*tmax);
			
			//singpu[i].h_img_pj=(float*)malloc(sizeof(float)*nxt*nzt);
			//singpu[i].h_ssg_pj=(float*)malloc(sizeof(float)*nxt*nzt);
			
			singpu[i].rindex=(int*)malloc(sizeof(int)*nxx*nzz);
			singpu[i].flagrec=(int*)malloc(sizeof(int)*nxx*nzz);
		 }

		 


		 // get pml parameters //
		 float *i_ax,*i_bx;
		 float *i_az,*i_bz;
		 float *h_ax,*h_bx;
		 float *h_az,*h_bz;
                 
                 i_ax=(float*)malloc(sizeof(float)*nxx);
                 i_bx=(float*)malloc(sizeof(float)*nxx);
                 i_az=(float*)malloc(sizeof(float)*nzz);
                 i_bz=(float*)malloc(sizeof(float)*nzz);
                 
                 h_ax=(float*)malloc(sizeof(float)*nxx);
                 h_bx=(float*)malloc(sizeof(float)*nxx);
                 h_az=(float*)malloc(sizeof(float)*nzz);
                 h_bz=(float*)malloc(sizeof(float)*nzz);
         
		 get_pml_parameters(pml,fd,i_ax,i_bx,i_az,i_bz,h_ax,h_bx,h_az,h_bz,nxx,nzz,dx,dz,dt);
		 
		 
		 
		 float *kxsmr,*kxscr,*kzsmr,*kzscr,*kxsmi,*kxsci,*kzsmi,*kzsci;
      		 kxsmr=(float*)malloc(sizeof(float)*nxx);
	         kxsmi=(float*)malloc(sizeof(float)*nxx);
	         kxscr=(float*)malloc(sizeof(float)*nxx);
	         kxsci=(float*)malloc(sizeof(float)*nxx);
	         kzsmr=(float*)malloc(sizeof(float)*nzz);
	         kzsmi=(float*)malloc(sizeof(float)*nzz);
	         kzscr=(float*)malloc(sizeof(float)*nzz);
	         kzsci=(float*)malloc(sizeof(float)*nzz);
	         cal_shift_parameters(kxsmr,kxscr,kzsmr,kzscr,kxsmi,kxsci,kzsmi,kzsci,dx,dz,nxx,nzz);
	         

		 // calculate PS operator //
		 float *k2d1;
		 float *k2d2;
		 float *k2d3;
		 float *k2da;
		 
		 //int nxh=(int)(nxx/2)+1;
		 //int nzh=(int)(nzz/2)+1;
		 
		 int Ns1=nxx*nzz;
		 //int Ns2=nxh*nzh;
		
		 k2d1=(float*)malloc(sizeof(float)*Ns1);
		 k2d2=(float*)malloc(sizeof(float)*Ns1);
		 k2d3=(float*)malloc(sizeof(float)*Ns1);
		 k2da=(float*)malloc(sizeof(float)*Ns1);

		 cal_k2d(dx,dz,nxx,nzz,k2d1,k2d2,k2d3,k2da);
		 
		 
		 
		 /*for(i=0;i<GPU_N;i++)
		 {
		 	singpu[i].host_U=(float*)malloc(sizeof(float)*Ns1*col);
		 	singpu[i].host_V=(float*)malloc(sizeof(float)*row*Ns2);
		 	singpu[i].host_A=(double*)malloc(sizeof(double)*col*row);
		 	singpu[i].host_Vt=(float*)malloc(sizeof(float)*row*Ns1);
		 	
		 	singpu[i].host_We=(float*)malloc(sizeof(float)*Ns2);
		 	singpu[i].host_Wa=(float*)malloc(sizeof(float)*Ns2);
		 }*/
		 
		 /*float *img, *ssg, *img_sp, *ssg_sp;
		 img=(float*)malloc(sizeof(float)*nxt*nzt);
		 ssg=(float*)malloc(sizeof(float)*nxt*nzt);
		 
		 img_sp=(float*)malloc(sizeof(float)*nxt*nzt);
		 ssg_sp=(float*)malloc(sizeof(float)*nxt*nzt);
		 

		for(i=0;i<nxt*nzt;i++){img_sp[i]=0.0;ssg_sp[i]=0.0;}*/

                if(myid==0) {start=clock();}

		 // each process runs its nshot forward modeling //
		 acforward(
			 	   snapnum,GPU_N,tmax,dt,rik,myid,numprocs,nshot,singpu,
				   gxp,gzp,nrmx,nx,nz,nxx,nzz,
				   pml,dx,dz,k2da,k2d1,k2d2,k2d3,sroffmx,disx,
				   vpt,vwdt,taot,lamt,gamt,rhot,bbxt,bbzt,
				   nxt,nzt,
		                   tracenum,flagout,
				   i_ax,i_bx,i_az,i_bz,h_ax,h_bx,h_az,h_bz,
				   kxsmr,kxscr,kzsmr,kzscr,kxsmi,kxsci,kzsmi,kzsci,
		                   fname4,fname5
                          );
                          
                //MPI_Reduce(img_sp,img,nxt*nzt, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
                //MPI_Reduce(ssg_sp,ssg,nxt*nzt, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if(myid==0) {
				end=clock();printf("The cost of the run time is %f seconds\n",
				(double)(end-start)/CLOCKS_PER_SEC);
				
				//fp1=fopen("vacimg.bin","wb");
				//fwrite(img,sizeof(float),nxt*nzt,fp1);
				//fwrite(ssg,sizeof(float),nxt*nzt,fp1);
				//fclose(fp1);
				
			    }

		 

		 for(i=0;i<GPU_N;i++)
		 {
			free(singpu[i].shotid);
			free(singpu[i].sxp);
			free(singpu[i].szp);

			free(singpu[i].h_vwd);
			free(singpu[i].h_tao);
			free(singpu[i].h_gam);
			free(singpu[i].h_lam);
			free(singpu[i].h_vp);
			
			free(singpu[i].h_rho);
			free(singpu[i].h_bbx);
			free(singpu[i].h_bbz);
			
			free(singpu[i].gxg);
			free(singpu[i].gzg);
			
			free(singpu[i].h_obs);
			free(singpu[i].rindex);
			free(singpu[i].flagrec);
			
			//free(singpu[i].h_img_pj);
			//free(singpu[i].h_ssg_pj);
			
			
			/*free(singpu[i].host_U);
			free(singpu[i].host_V);
			free(singpu[i].host_A);
			free(singpu[i].host_Vt);
			
			free(singpu[i].host_We);
			free(singpu[i].host_Wa);*/
		 }             
		
                 free(tracenum);free(shotxp);free(shotzp);
                 free(rhot);free(bbxt);free(bbzt);free(rik);
                 
                 free(vpt);free(Qpt);free(vwdt);free(taot);
                 free(lamt);free(gamt);
                 
                 //free(img);free(ssg);
                 //free(img_sp);free(ssg_sp);
                 
		 free(i_ax);free(i_bx);free(i_az);free(i_bz);
		 free(h_ax);free(h_bx);free(h_az);free(h_bz);
		 free(kxsmr);free(kxscr);free(kzsmr);free(kzscr);
		 free(kxsmi);free(kxsci);free(kzsmi);free(kzsci);
		 
		 free(k2da);free(k2d1);free(k2d2);free(k2d3);
		 free(gxp);free(gzp);

   
     MPI_Finalize();
    
    return 0;
}
