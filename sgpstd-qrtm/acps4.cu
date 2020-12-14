#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "statement.h"

#define BLOCK_WIDTH 16
#define BLOCK_HEIGHT 16

__global__ void storecheckpoints
		(
			int flag,int nc,int nxx,int nzz,
			float *d_p,float *d_vx,float *d_vz,
			float *d_fipx,float *d_fipz,
			float *d_fivxx,float *d_fivzz,
			float *store_p,float *store_vx,float *store_vz,
			float *store_fipx,float *store_fipz,
			float *store_fivxx,float *store_fivzz
		);
		
__global__ void addsource(int nx,int nz,int xs,int zs,float *check_d_ps,float *d_rik,float dt);

__global__ void initialize_kernel(
					int nxx,int nzz,int pml,
					float *d_p,float *d_vx,float *d_vz,
					float *d_fipx,float *d_fipz,
					float *d_fivxx,float *d_fivzz,
					float *d_img,float *d_ssg,
					float *d_pb,float *d_vxb,float *d_vzb,
					float *d_fipxb,float *d_fipzb,
					float *d_fivxxb,float *d_fivzzb
				);
				
__global__ void initialize_kernel_ins(
					int nxx,int nzz,int pml,
					float *d_p,float *d_vx,float *d_vz,
					float *d_fipx,float *d_fipz,
					float *d_fivxx,float *d_fivzz,float *d_aux
				);
				
                                  
__global__ void fuzhi_kernel(
					int nxx,int nzz,float *d_vx,float *d_vz,
					cufftComplex *d_dis0,cufftComplex *d_dis1
			    );
			    
			    
__global__ void fuzhi_kernel_p (int nxx,int nzz,float *d_p,cufftComplex *d_dis0);

__global__ void multt_kernel(float eta,float dt,int nt,float tao,
			     int NXZ,int nxx,int nzz,
			     float *k2da,cufftComplex *d_dis0
			     );

__global__ void fuzhi_kernel_inverse(int nxx,int nzz,float *d_aux,cufftComplex *d_dis0);

__global__ void shift_vxz_kernel(
					float *kxsmr,float *kxsmi,float *kzscr,float *kzsci,
					cufftComplex *d_outdatax,cufftComplex *d_outdataz,
					cufftComplex *d_indatax,cufftComplex *d_indataz,
					int nxx,int nzz
				);
				
__global__ void shift_pxz_kernel(
					float *kxscr,float *kxsci,float *kzsmr,float *kzsmi,
     					cufftComplex *d_dis0,
     					cufftComplex *d_dis1,cufftComplex *d_dis2,
     					int nxx,int nzz
     				);

__global__ void multv_kernel(int N1,int nxx,int nzz,int ig,float *d_V,cufftComplex *d_psk,cufftComplex *d_pskf);

__global__ void fuzhi_to_ifftrow(int ig,int nxx,int nzz,int N1,cufftComplex *d_pskf,cufftComplex *d_pifft);

__global__ void fuzhi_r_kernel(int nxx,int nzz,int jg,int N1,double *d_pst);

__global__ void sumlw_kernel(int nxx,int nzz,int jg,int N1,int row,int ig,double *d_pst,double *d_A,cufftComplex *d_pifft);

__global__ void sumlw2_kernel(int nxx,int nzz,int N1,int col,double *d_pst,float *d_U,float *d_psc);

__global__ void update_vxx_vzz(
       				int nxx,int nzz,float *d_fi_vx,float *d_fi_vz,
       				float *di_a_x,float *di_b_x,float *di_a_z,float *di_b_z,
       				cufftComplex *d_outdatax,cufftComplex *d_outdataz
       			       );

__global__ void addsa_kernel(
				int restart,int pml,int nxx,int nzz,
				float *d_vp,float *d_rho,float dt,
				float *d_p,cufftComplex *d_dis0,cufftComplex *d_dis3,
				int xs,int zs,int nt,float *d_Rik,float *d_ssg
			    );
			    
__global__ void addsa_b_kernel(
				int pml,int nxx,int nzz,
				float *d_vp,float *d_rho,float dt,
				float *d_p,cufftComplex *d_dis0,cufftComplex *d_dis3,
				int nt,int tmax,float *d_obs,
				int *rindex, int *flagrec,int nr
			    );
			    
__global__ void vx_cal_kernel(
     				int nxx,int nzz,float dt,float *d_bbx,
     				float *d_fi_px,float *dh_a_x,
     				float *dh_b_x,cufftComplex *d_outdatax,float *d_vx
     			     );

__global__ void vz_cal_kernel(
     				int nxx,int nzz,float dt,float *d_bbz,
     				float *d_fi_pz,float *dh_a_z,
     				float *dh_b_z,cufftComplex *d_outdataz,float *d_vz
     			     );
     			     
__global__ void rtm_image(
				int nx,int nz,int nxx,int nzz,
				float *d_p,float *d_pb,
				float *d_img,float *d_imgf
			 );
			 
__global__ void filter_image(
				int nx,int nz,float dx,float dz,
				float *d_img,float *d_imgf
			    );

__global__ void sumpxpz(
       				int nxx,int nzz,
       				cufftComplex *d_outdatax,cufftComplex *d_outdataz,
				cufftComplex *d_dis0
       			       );

__global__ void  sum_total(int nxx,int nzz,
			   cufftComplex *d_dis0,cufftComplex *d_dis2,
			   cufftComplex *d_dis3,cufftComplex *d_dis4,
			   float *d_total,float *d_vp,float *d_lam,
			   float *d_gam,float *d_vwd,float dt);

__global__ void multtt_kernel(int nxx,int nzz,float *d_k2d3,cufftComplex *d_dis0,cufftComplex *d_dis1);
     	
extern "C" float get_max_tao(int nxx,int nzz,float *tao);

extern "C" float get_min_tao(int nxx,int nzz,float *tao);		       
                         
extern "C" void  expand_pml_areas(int nx,int nz,int nxx,int nzz,int pml,float *vp);

extern "C" void  select_part_model(int,float mlx,float dx,int nxt,int nzt,int nx,int nz,int nxx,int nzz,float *vpt,float *vp);

extern "C" void input_geophone_pos(int shotid,int nr,float *gxp,char fname[]);

extern "C" void acsglowrank2d(int row,int col,int NXZ,int Nh,
			       float *vp,float *tao,float *lam,float *gam,
			       float *k2dh,float dt,float *U,float *V,
			       double *A);
			       
extern "C" void checka_error(int row,int col,int N1,int N2,float *U,double *A,float *V,
			     float vp,float gam,float tao,float dt,
                             float *k2dh,float *We,float *Wa,int is);
                             
extern "C" void consist_Vt(int row,int nxx,int nzz,int nxh,int nzh,int N1,int N2,float *V,float *Vt);


extern "C" void projectomodel(
				int pml,float mlx,float dx,int nxt,int nzt,
			      	int nx,int nz,int nxx,int nzz,float *vpt,float *vp
			     );
			     
extern "C" void removedirect(
				int nrmx,int tmax,int nr,float dt,float fd,float xs,float zs,
				float *gxp,float *gzp,float v0,float *recsim
			    );


void Alloc_host_device_mem(
				int num,int nrmx,int nx,int nz,
				int nxx,int nzz,int nxh,int nzh,int NXZ,int Nh,
			   	int tmax,int GPU_N,struct MultiGPU singpu[]
			  );
			   
extern "C" enum action revolve(int* check,int* capo,int* fine,int snaps,int* info);



extern "C" void acforward(
				int snapnum,int GPU_N,int tmax,float dt,float *rik,int myid,int numprocs,
			  	int nshot,struct MultiGPU singpu[],float *gxp,float *gzp,
			  	int nrmx,int nx,int nz,int nxx,int nzz,
		          	int pml,float dx,float dz,float *k2da,float *k2d1,float *k2d2,float *k2d3,
		          	float sroffmx,float disx,
		          	float *vpt,
		          	float *rhot,float *bbxt,float *bbzt,
			  	int nxt,int nzt,int *tracenum,
			  	float *i_ax,float *i_bx,float *i_az,float *i_bz,
			 	float *h_ax,float *h_bx,float *h_az,float *h_bz,
			  	float *kxsmr,float *kxscr,float *kzsmr,float *kzscr,
			  	float *kxsmi,float *kxsci,float *kzsmi,float *kzsci,
			  	float *img_sp,float *ssg_sp,float fobs,float v0,
                          	char fname4[40],char fname5[40],char recname[60]
                          )                                
{
    // GPU configuration//
    dim3 dimBlock(BLOCK_WIDTH,BLOCK_HEIGHT);
    dim3 dimGrid((nxx+dimBlock.x-1)/dimBlock.x,(nzz+dimBlock.y-1)/dimBlock.y);

    
    int i,j,t,is,ns,nt;
    int ip;
    int shotid;
    int NXZ=nxx*nzz;
    
    int nxh=(int)(nxx/2)+1;
    int nzh=(int)(nzz/2)+1;
    int Nh=nxh*nzh;

    //printf("%d\n",snapnum);

    // allocate device memory for variables //
    Alloc_host_device_mem(snapnum,nrmx,nx,nz,nxx,nzz,nxh,nzh,NXZ,Nh,tmax,GPU_N,singpu);
    
    // create multi-stream //
    Multistream plansm[GPU_N];
    for(i=0;i<GPU_N;i++)
    {
	cudaSetDevice(i);
	cudaStreamCreate(&plansm[i].stream);	
	cufftSetStream(singpu[i].plan2dforward,plansm[i].stream);
	cufftSetStream(singpu[i].plan2dinverse,plansm[i].stream);
    }

     
    // copy values from host to device //
    for(i=0;i<GPU_N;i++)
    {
	cudaSetDevice(i); 
	cudaMemcpyAsync(singpu[i].d_rik,rik,tmax*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);

	
	cudaMemcpyAsync(singpu[i].d_k2d1,k2d1,NXZ*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	

	cudaMemcpyAsync(singpu[i].d_i_az,i_az,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_i_bz,i_bz,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
        cudaMemcpyAsync(singpu[i].d_i_ax,i_ax,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_i_bx,i_bx,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	
	cudaMemcpyAsync(singpu[i].d_h_az,h_az,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_h_bz,h_bz,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
        cudaMemcpyAsync(singpu[i].d_h_ax,h_ax,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_h_bx,h_bx,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	
			  
	cudaMemcpyAsync(singpu[i].d_kxsmr,kxsmr,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_kxscr,kxscr,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_kxsmi,kxsmi,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_kxsci,kxsci,nxx*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	
	cudaMemcpyAsync(singpu[i].d_kzsmr,kzsmr,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_kzscr,kzscr,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_kzsmi,kzsmi,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
	cudaMemcpyAsync(singpu[i].d_kzsci,kzsci,nzz*sizeof(float),cudaMemcpyHostToDevice,plansm[i].stream);
    }

    //printf("%f\n",rik[0]);

    
    float xs,zs;


    //char snapname[60];
    char obsname[80];
    FILE *fr;
    //FILE *fs=fopen("snap0.bin","wb");
    //FILE *fss=fopen("snap00.bin","wb");
    //FILE *fb=fopen("snap1.bin","wb");
    
    //int ick;
    
   // char checkname[60];
   // FILE *fck;
    
    //int flags;
    
    // for wavefield reconstruction //
    int check, capo, fine, steps, snaps, info;
    enum action whatodo;
    capo = 0;
    snaps = snapnum;
    steps = tmax ;
    info  = 3 ;
    fine  = steps + capo;
    check = -1;
    
    int ts = 0;
    int te = tmax-1;
    int tea = 0 ;
    int restart = 0;
	
   //FILE *fv=fopen("data.bin","wb");
	


    for(ns=0;ns<nshot;ns++)   // shot iteration //
    {
    	//flags = 1;
    	restart = 0; ts =0; te = tmax-1; tea = 0 ;
	capo = 0; snaps = snapnum; steps = tmax; info  = 3; fine  = steps + capo; check = -1;
    	

	//---------------------------------------------------------------------------//
	// 1. get velocity for each gpu //
	// 2. get source x and z grid index for each gpu //
	// 3. get receiver x and z grid index for each gpu with the trace number nr //
        
	for(ip=0;ip<GPU_N;ip++)       // gpu iteration  //
	{
		shotid=singpu[ip].shotid[ns];
		xs=singpu[ip].sxp[ns];
		zs=singpu[ip].szp[ns];
		printf("myid=%d shotid=%d nx=%d nz=%d xs=%f zs=%f\n",myid,shotid,nx,nz,xs,zs);
		
		singpu[ip].lgest=10000;
		singpu[ip].rgest=-1;

		// input vp from vpt model //
		singpu[ip].mlx=xs-sroffmx-disx;

		select_part_model(pml,singpu[ip].mlx,dx,nxt,nzt,nx,nz,nxx,nzz,vpt,singpu[ip].h_vp);
               
                
                select_part_model(pml,singpu[ip].mlx,dx,nxt,nzt,nx,nz,nxx,nzz,rhot,singpu[ip].h_rho);
                select_part_model(pml,singpu[ip].mlx,dx,nxt,nzt,nx,nz,nxx,nzz,bbxt,singpu[ip].h_bbx);
                select_part_model(pml,singpu[ip].mlx,dx,nxt,nzt,nx,nz,nxx,nzz,bbzt,singpu[ip].h_bbz);
                
                expand_pml_areas(nx,nz,nxx,nzz,pml,singpu[ip].h_vp);
		
                
                expand_pml_areas(nx,nz,nxx,nzz,pml,singpu[ip].h_rho);
                expand_pml_areas(nx,nz,nxx,nzz,pml,singpu[ip].h_bbx);
                expand_pml_areas(nx,nz,nxx,nzz,pml,singpu[ip].h_bbz);
                
                

		//fwrite(singpu[ip].h_vp,sizeof(float),NXZ,fv);
                //fwrite(singpu[ip].h_vwd,sizeof(float),NXZ,fv);
                //fwrite(singpu[ip].h_tao,sizeof(float),NXZ,fv);
		//fwrite(singpu[ip].h_gam,sizeof(float),NXZ,fv);
		//fwrite(singpu[ip].h_lam,sizeof(float),NXZ,fv);
		//fwrite(singpu[ip].h_rho,sizeof(float),NXZ,fv);
		//fwrite(singpu[ip].h_bbx,sizeof(float),NXZ,fv);
		//fwrite(singpu[ip].h_bbz,sizeof(float),NXZ,fv);
		
		 // locate source grid index on the model //
                singpu[ip].zsg=(int)(zs/dz+0.5)+pml;
                singpu[ip].xsg=(int)((sroffmx+disx)/dx+0.5)+pml;
		if(singpu[ip].xsg<0||singpu[ip].xsg>=nxx||singpu[ip].zsg<0||singpu[ip].zsg>=nzz)
                {
			printf("%d shot grid index exceeds model!\n",shotid);
			exit(0);
		}
		
		// locate receiver's grid index on the model //
                singpu[ip].nr=tracenum[shotid];
                input_geophone_pos(shotid,singpu[ip].nr,gxp,fname4);
	        input_geophone_pos(shotid,singpu[ip].nr,gzp,fname5);
	        
	        for(i=0;i<NXZ;i++)
		{
			singpu[ip].flagrec[i]=0;
		}


                for(i=0;i<singpu[ip].nr;i++)
		{
		    singpu[ip].gzg[i]=(int)(gzp[i]/dz+0.5)+pml;
		    if(gxp[i]>=xs){singpu[ip].gxg[i]=singpu[ip].xsg+(int)((gxp[i]-xs)/dx+0.5);}
		    else{singpu[ip].gxg[i]=singpu[ip].xsg-(int)((xs-gxp[i])/dx+0.5);}

		    if(singpu[ip].gxg[i]<pml){printf("warning: receiver x grid index in left pml area!\n");singpu[ip].gxg[i]=pml;}
                    if(singpu[ip].gxg[i]>=nxx-pml){printf("warning: receiver x grid index in right pml area!\n");singpu[ip].gxg[i]=nxx-pml-1;}
                    
                    // mark receiver grid index //
                    is=singpu[ip].gxg[i]*nzz+singpu[ip].gzg[i];
                    singpu[ip].flagrec[is]=1;
                    singpu[ip].rindex[is]=i;
                    
                    if(singpu[ip].gxg[i]>singpu[ip].rgest){singpu[ip].rgest=singpu[ip].gxg[i];}
		    if(singpu[ip].gxg[i]<singpu[ip].lgest){singpu[ip].lgest=singpu[ip].gxg[i];}
                }
                
                for(i=0;i<nxt*nzt;i++)
		{	
			singpu[ip].h_img_pj[i]=0.0; singpu[ip].h_ssg_pj[i]=0.0;
		}
		
		
		// input observed data //
		strcpy(obsname,"");
                sprintf(obsname,"%s%d_%d.bin",recname,shotid,singpu[ip].nr); 
                fr=fopen(obsname,"rb");
		fread(&singpu[ip].h_obs[0],sizeof(float),tmax*singpu[ip].nr,fr);
		fclose(fr);
		
		
		
		removedirect(nrmx,tmax,singpu[ip].nr,dt,fobs,xs,zs,gxp,gzp,v0,singpu[ip].h_obs);
		
		

	} // end gpu iteration 
        //---------------------------------------------------------------------------//


	// give a zero initial value //
	for(ip=0;ip<GPU_N;ip++)
	{
		cudaSetDevice(ip);
		
		cudaMemcpyAsync(singpu[ip].d_vp,singpu[ip].h_vp,NXZ*sizeof(float),cudaMemcpyHostToDevice,plansm[ip].stream);
		
		cudaMemcpyAsync(singpu[ip].d_rho,singpu[ip].h_rho,NXZ*sizeof(float),cudaMemcpyHostToDevice,plansm[ip].stream);
		cudaMemcpyAsync(singpu[ip].d_bbx,singpu[ip].h_bbx,NXZ*sizeof(float),cudaMemcpyHostToDevice,plansm[ip].stream);
		cudaMemcpyAsync(singpu[ip].d_bbz,singpu[ip].h_bbz,NXZ*sizeof(float),cudaMemcpyHostToDevice,plansm[ip].stream);
		
		
		
		cudaMemcpyAsync(singpu[ip].d_flagrec,singpu[ip].flagrec,nxx*nzz*sizeof(int),cudaMemcpyHostToDevice,plansm[ip].stream);
		cudaMemcpyAsync(singpu[ip].d_rindex,singpu[ip].rindex,nxx*nzz*sizeof(int),cudaMemcpyHostToDevice,plansm[ip].stream);
		
		initialize_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				(
					nxx,nzz,pml,
					singpu[ip].d_p,singpu[ip].d_vx,singpu[ip].d_vz,
				 	singpu[ip].d_fipx,singpu[ip].d_fipz,
				 	singpu[ip].d_fivxx,singpu[ip].d_fivzz,
				 	singpu[ip].d_img,singpu[ip].d_ssg,
				 	singpu[ip].d_pb,singpu[ip].d_vxb,singpu[ip].d_vzb,
				 	singpu[ip].d_fipxb,singpu[ip].d_fipzb,
				 	singpu[ip].d_fivxxb,singpu[ip].d_fivzzb
				);
	}

	//for(ip=0;ip<GPU_N;ip++)
	//{
	//	cudaSetDevice(ip);
		cudaDeviceSynchronize();
	//}
	
	do
	{
		whatodo = revolve(&check, &capo, &fine, snaps, &info);
		
		// case 1: store checkpoint //
		if ( (whatodo == takeshot) && (info > 1) )
		{
			for(ip=0;ip<GPU_N;ip++)
			{
				cudaSetDevice(ip);
				storecheckpoints<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
						(
							1,check,nxx,nzz,
							singpu[ip].d_p,singpu[ip].d_vx,singpu[ip].d_vz,
							singpu[ip].d_fipx,singpu[ip].d_fipz,
							singpu[ip].d_fivxx,singpu[ip].d_fivzz,
							singpu[ip].store_p,singpu[ip].store_vx,singpu[ip].store_vz,
							singpu[ip].store_fipx,singpu[ip].store_fipz,
							singpu[ip].store_fivxx,singpu[ip].store_fivzz
						);
						
				if( capo==0 )
				{
					addsource<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
							(
								nxx,nzz,singpu[ip].xsg,singpu[ip].zsg,
								singpu[ip].store_p,singpu[ip].d_rik,dt
							);
				}
						
			}
			
		} // end case 1
		
		
		// advance source wavefield //
		if ( (whatodo == advance) && (info > 2) )
		{
			te = capo;
			for(t=ts;t<=te;t++)
			{	
				for(ip=0;ip<GPU_N;ip++)
				{
					cudaSetDevice(ip);
					
					fuzhi_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
								(
									nxx,nzz,singpu[ip].d_vx,singpu[ip].d_vz,singpu[ip].d_dis0,singpu[ip].d_dis1
								);
								
					
				  	cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis0,CUFFT_FORWARD);
				  	cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis1,singpu[ip].d_dis1,CUFFT_FORWARD);
				  	
				  	
				  	shift_vxz_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				  				(
				  					singpu[ip].d_kxsmr,singpu[ip].d_kxsmi,
				  					singpu[ip].d_kzscr,singpu[ip].d_kzsci,
				  					singpu[ip].d_dis0,singpu[ip].d_dis1,
				  				        singpu[ip].d_dis2,singpu[ip].d_dis3,
				  				        nxx,nzz
				  				);
				      
				        cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis2,singpu[ip].d_dis2,CUFFT_INVERSE); // ps to x
					cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis3,singpu[ip].d_dis3,CUFFT_INVERSE); // ps to z

					update_vxx_vzz<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				      				(
				      					nxx,nzz,singpu[ip].d_fivxx,singpu[ip].d_fivzz,
				                                     	singpu[ip].d_i_ax,singpu[ip].d_i_bx,
				                                     	singpu[ip].d_i_az,singpu[ip].d_i_bz,
				                                     	singpu[ip].d_dis2,singpu[ip].d_dis3
				                                ); 

					sumpxpz<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
								(
				      					nxx,nzz,
				                                     	singpu[ip].d_dis2,singpu[ip].d_dis3,singpu[ip].d_dis0
				                                ); 

					cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis1,CUFFT_FORWARD);

					
					
					multtt_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>(nxx,nzz,singpu[ip].d_k2d1,
											        singpu[ip].d_dis1,singpu[ip].d_dis3);


					cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis3,singpu[ip].d_dis3,CUFFT_INVERSE); 


				 	// add dispersion and attenuation terms: inside computational domain //
				  	addsa_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				  				(
				  					restart,pml,nxx,nzz,
				  					singpu[ip].d_vp,singpu[ip].d_rho,dt,
				  					singpu[ip].d_p,singpu[ip].d_dis0,singpu[ip].d_dis3,
				  					singpu[ip].xsg,singpu[ip].zsg,
				  					t,singpu[ip].d_rik,singpu[ip].d_ssg
				  				);

					//cudaMemcpy(singpu[ip].h_p,singpu[ip].d_p,NXZ*sizeof(float),cudaMemcpyDeviceToHost);

				  	//fwrite(singpu[ip].h_p,sizeof(float),NXZ,fs);

				  	
				  	fuzhi_kernel_p<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
								(
									nxx,nzz,singpu[ip].d_p,singpu[ip].d_dis0
								);
								
					cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis0,CUFFT_FORWARD);
					
					shift_pxz_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
								(
									singpu[ip].d_kxscr,singpu[ip].d_kxsci,
									singpu[ip].d_kzsmr,singpu[ip].d_kzsmi,
				                                        singpu[ip].d_dis0,
				                                        singpu[ip].d_dis1,singpu[ip].d_dis2,
				                                        nxx,nzz
				                                );

					cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis1,singpu[ip].d_dis1,CUFFT_INVERSE);
					cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis2,singpu[ip].d_dis2,CUFFT_INVERSE);
	
				      
				      
				      vx_cal_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				      				(
				      					nxx,nzz,dt,singpu[ip].d_bbx,
				      					singpu[ip].d_fipx,singpu[ip].d_h_ax,singpu[ip].d_h_bx,
				      					singpu[ip].d_dis1,singpu[ip].d_vx
				      				);
				      				
				      vz_cal_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				      				(
				      					nxx,nzz,dt,singpu[ip].d_bbz,
				      					singpu[ip].d_fipz,singpu[ip].d_h_az,singpu[ip].d_h_bz,
				      					singpu[ip].d_dis2,singpu[ip].d_vz
				      				);
				} // end gpu
				
			}// end time
			ts = te +1 ;
		}// end case 2
		
		// ----- indicating already finishing forward modeling ----- //
		if ( (whatodo == firsturn) && (info > 2) )
		{
			restart = 1;
			// compute receiver wavefield at t=tmax-1//
			for(ip=0;ip<GPU_N;ip++)
			{
				cudaSetDevice(ip);
				cudaMemcpyAsync(singpu[ip].d_obs,singpu[ip].h_obs,nrmx*tmax*sizeof(float),cudaMemcpyHostToDevice,plansm[ip].stream);
				
				        
				
				fuzhi_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
						(
							nxx,nzz,singpu[ip].d_vxb,singpu[ip].d_vzb,singpu[ip].d_dis0,singpu[ip].d_dis1
						);
						
			
		          	cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis0,CUFFT_FORWARD);
		          	cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis1,singpu[ip].d_dis1,CUFFT_FORWARD);
		          	
		          	
		          	shift_vxz_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		          				(
		          					singpu[ip].d_kxsmr,singpu[ip].d_kxsmi,
		          					singpu[ip].d_kzscr,singpu[ip].d_kzsci,
		          					singpu[ip].d_dis0,singpu[ip].d_dis1,
		          				        singpu[ip].d_dis2,singpu[ip].d_dis3,
		          				        nxx,nzz
		          				);
		          	
		          
		          	cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis2,singpu[ip].d_dis2,CUFFT_INVERSE); // ps to x
				cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis3,singpu[ip].d_dis3,CUFFT_INVERSE); // ps to z

				update_vxx_vzz<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				      				(
				      					nxx,nzz,singpu[ip].d_fivxxb,singpu[ip].d_fivzzb,
				                                     	singpu[ip].d_i_ax,singpu[ip].d_i_bx,
				                                     	singpu[ip].d_i_az,singpu[ip].d_i_bz,
				                                     	singpu[ip].d_dis2,singpu[ip].d_dis3
				                                ); 

					sumpxpz<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
								(
				      					nxx,nzz,
				                                     	singpu[ip].d_dis2,singpu[ip].d_dis3,singpu[ip].d_dis0
				                                ); 

					cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis1,CUFFT_FORWARD);

					
					multtt_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>(nxx,nzz,singpu[ip].d_k2d1,
											        singpu[ip].d_dis1,singpu[ip].d_dis3);

					

					
					cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis3,singpu[ip].d_dis3,CUFFT_INVERSE); 
					

					




		         	// add dispersion and attenuation terms: inside computational domain // 
		          	addsa_b_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		          				(
		          					pml,nxx,nzz,
		          					singpu[ip].d_vp,
		          					singpu[ip].d_rho,dt,
		          					singpu[ip].d_pb,singpu[ip].d_dis0,singpu[ip].d_dis3,
		          					0,tmax,singpu[ip].d_obs,
		          					singpu[ip].d_rindex,singpu[ip].d_flagrec,singpu[ip].nr
		          				);// please note here //


		          	fuzhi_kernel_p<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
							(
								nxx,nzz,singpu[ip].d_pb,singpu[ip].d_dis0
							);
							
				cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis0,CUFFT_FORWARD);
				
				shift_pxz_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
							(
								singpu[ip].d_kxscr,singpu[ip].d_kxsci,
								singpu[ip].d_kzsmr,singpu[ip].d_kzsmi,
		                                                singpu[ip].d_dis0,
		                                                singpu[ip].d_dis1,singpu[ip].d_dis2,
		                                                nxx,nzz
		                                        );

				cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis1,singpu[ip].d_dis1,CUFFT_INVERSE);
				cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis2,singpu[ip].d_dis2,CUFFT_INVERSE);
		                                        
		                                        
		                
		                                        
		                                        
		               
		              
		              
		              vx_cal_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		              				(
		              					nxx,nzz,dt,singpu[ip].d_bbx,
		              					singpu[ip].d_fipxb,singpu[ip].d_h_ax,singpu[ip].d_h_bx,
		              					singpu[ip].d_dis1,singpu[ip].d_vxb
		              				);
		              				
		              vz_cal_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		              				(
		              					nxx,nzz,dt,singpu[ip].d_bbz,
		              					singpu[ip].d_fipzb,singpu[ip].d_h_az,singpu[ip].d_h_bz,
		              					singpu[ip].d_dis2,singpu[ip].d_vzb
		              				);
		              				
		              				
		              rtm_image<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		              				(
		              					nx,nz,nxx,nzz,singpu[ip].d_p,singpu[ip].d_pb,
		              					singpu[ip].d_img,singpu[ip].d_imgf
		              				);
		              				
			}// end gpu 

		}// end case 3
		
		
		// the source wavefields are available at this time //
		if ( (whatodo == youturn) && (info > 2) )
		{
			tea = capo; nt = tea;
			// compute receiver wavefield at t=nt //
			for(ip=0;ip<GPU_N;ip++)
			{
				cudaSetDevice(ip);
				
				fuzhi_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
						(
							nxx,nzz,singpu[ip].d_vxb,singpu[ip].d_vzb,singpu[ip].d_dis0,singpu[ip].d_dis1
						);
						
			
		          	cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis0,CUFFT_FORWARD);
		          	cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis1,singpu[ip].d_dis1,CUFFT_FORWARD);
		          	
		          	
		          	shift_vxz_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		          				(
		          					singpu[ip].d_kxsmr,singpu[ip].d_kxsmi,
		          					singpu[ip].d_kzscr,singpu[ip].d_kzsci,
		          					singpu[ip].d_dis0,singpu[ip].d_dis1,
		          				        singpu[ip].d_dis2,singpu[ip].d_dis3,
		          				        nxx,nzz
		          				);
		          	
		          
		          	cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis2,singpu[ip].d_dis2,CUFFT_INVERSE); // ps to x
				cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis3,singpu[ip].d_dis3,CUFFT_INVERSE); // ps to z

				update_vxx_vzz<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
				      				(
				      					nxx,nzz,singpu[ip].d_fivxxb,singpu[ip].d_fivzzb,
				                                     	singpu[ip].d_i_ax,singpu[ip].d_i_bx,
				                                     	singpu[ip].d_i_az,singpu[ip].d_i_bz,
				                                     	singpu[ip].d_dis2,singpu[ip].d_dis3
				                                ); 

					sumpxpz<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
								(
				      					nxx,nzz,
				                                     	singpu[ip].d_dis2,singpu[ip].d_dis3,singpu[ip].d_dis0
				                                ); 

					cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis1,CUFFT_FORWARD);

					
					
					multtt_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>(nxx,nzz,singpu[ip].d_k2d1,
											        singpu[ip].d_dis1,singpu[ip].d_dis3);

					

					
					cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis3,singpu[ip].d_dis3,CUFFT_INVERSE); 
					

					




		         	// add dispersion and attenuation terms: inside computational domain // 
		          	addsa_b_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		          				(
		          					pml,nxx,nzz,
		          					singpu[ip].d_vp,
		          					singpu[ip].d_rho,dt,
		          					singpu[ip].d_pb,singpu[ip].d_dis0,singpu[ip].d_dis3,
		          					tmax-1-nt,tmax,singpu[ip].d_obs,
		          					singpu[ip].d_rindex,singpu[ip].d_flagrec,singpu[ip].nr
		          				);// please note here //


		         	

		          	fuzhi_kernel_p<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
							(
								nxx,nzz,singpu[ip].d_pb,singpu[ip].d_dis0
							);
							
				cufftExecC2C(singpu[ip].plan2dforward,singpu[ip].d_dis0,singpu[ip].d_dis0,CUFFT_FORWARD);
				
				shift_pxz_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
							(
								singpu[ip].d_kxscr,singpu[ip].d_kxsci,
								singpu[ip].d_kzsmr,singpu[ip].d_kzsmi,
		                                                singpu[ip].d_dis0,
		                                                singpu[ip].d_dis1,singpu[ip].d_dis2,
		                                                nxx,nzz
		                                        );

				cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis1,singpu[ip].d_dis1,CUFFT_INVERSE);
				cufftExecC2C(singpu[ip].plan2dinverse,singpu[ip].d_dis2,singpu[ip].d_dis2,CUFFT_INVERSE);
		                                        
		                                        
		                
		                                        
		                                        
		               
		              
		              
		              vx_cal_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		              				(
		              					nxx,nzz,dt,singpu[ip].d_bbx,
		              					singpu[ip].d_fipxb,singpu[ip].d_h_ax,singpu[ip].d_h_bx,
		              					singpu[ip].d_dis1,singpu[ip].d_vxb
		              				);
		              				
		              vz_cal_kernel<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		              				(
		              					nxx,nzz,dt,singpu[ip].d_bbz,
		              					singpu[ip].d_fipzb,singpu[ip].d_h_az,singpu[ip].d_h_bz,
		              					singpu[ip].d_dis2,singpu[ip].d_vzb
		              				);
		              				
		              				
		              rtm_image<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
		              				(
		              					nx,nz,nxx,nzz,singpu[ip].d_p,singpu[ip].d_pb,
		              					singpu[ip].d_img,singpu[ip].d_imgf
		              				);
		              				
				
			}// end gpu iteration
			
		} // end case 4
		
		
		// load stored snapshots //
		if ((whatodo == restore) && (info > 2)) 
		{
			for(ip=0;ip<GPU_N;ip++)
			{
				cudaSetDevice(ip);
				
				storecheckpoints<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
						(
							0,check,nxx,nzz,
							singpu[ip].d_p,singpu[ip].d_vx,singpu[ip].d_vz,
							singpu[ip].d_fipx,singpu[ip].d_fipz,
							singpu[ip].d_fivxx,singpu[ip].d_fivzz,
							singpu[ip].store_p,singpu[ip].store_vx,singpu[ip].store_vz,
							singpu[ip].store_fipx,singpu[ip].store_fipz,
							singpu[ip].store_fivxx,singpu[ip].store_fivzz
						);
				
			}
			ts = capo + 1;
		}// end case 5
		
		
		// return error information //
		if ( whatodo == error )
		{
			printf(" irregular termination of revolve \n");
		      	switch(info)
		       	{
				case 10: printf(" number of checkpoints stored exceeds checkup, \n");
				 	 printf(" increase constant 'checkup' and recompile \n");
				 	 break;
				case 11: printf(" number of checkpoints stored = %d exceeds snaps = %d, \n"
						 ,check+1,snaps);
					 printf(" ensure 'snaps' > 0 and increase initial 'fine' \n");
					 break;
				case 12: printf(" error occurs in numforw \n");
					 break;
				case 13: printf(" enhancement of 'fine', 'snaps' checkpoints stored, \n");
					 printf(" increase 'snaps'\n");
					 break;
				case 14: printf(" number of snaps exceeds snapsup, ");
					 printf(" increase constant 'snapsup' and recompile \n");
					 break;
				case 15: printf(" number of reps exceeds repsup, ");
					 printf(" increase constant 'repsup' and recompile \n");       
			}
		} // end case 6
		
		
	}while( (whatodo != terminate) && (whatodo != error) );

	
	
	for(ip=0;ip<GPU_N;ip++)
	{
		cudaSetDevice(ip);
		
		filter_image<<<dimGrid,dimBlock,0,plansm[ip].stream>>>
                      				(
                      					nx,nz,dx,dz,singpu[ip].d_img,singpu[ip].d_imgf
                      				);
                      				
                cudaMemcpyAsync(singpu[ip].h_img,singpu[ip].d_imgf,nx*nz*sizeof(float),cudaMemcpyDeviceToHost,plansm[ip].stream);
                
	        cudaMemcpyAsync(singpu[ip].h_ssg,singpu[ip].d_ssg, nx*nz*sizeof(float),cudaMemcpyDeviceToHost,plansm[ip].stream);
	        
	        cudaStreamSynchronize(plansm[ip].stream);
		
	}
	
	for(ip=0;ip<GPU_N;ip++)
	{
			singpu[ip].lgest=singpu[ip].lgest-pml;
			singpu[ip].rgest=singpu[ip].rgest-pml;

			for(j=0;j<nz;j++)
			{
				for(i=0;i<singpu[ip].lgest;i++){singpu[ip].h_img[i*nz+j]=0.0;}
				for(i=singpu[ip].rgest;i<nx;i++){singpu[ip].h_img[i*nz+j]=0.0;}
			}

			projectomodel(pml,singpu[ip].mlx,dx,nxt,nzt,nx,nz,nxx,nzz,singpu[ip].h_img_pj,singpu[ip].h_img);
			projectomodel(pml,singpu[ip].mlx,dx,nxt,nzt,nx,nz,nxx,nzz,singpu[ip].h_ssg_pj,singpu[ip].h_ssg);

			for(i=0;i<nxt*nzt;i++)
                 	{
				img_sp[i]=img_sp[i]+singpu[ip].h_img_pj[i];
				ssg_sp[i]=ssg_sp[i]+singpu[ip].h_ssg_pj[i];
		 	}
	}
	

    } // end shot iteration


    // free the device and host memory //
    for(i=0;i<GPU_N;i++)
    {
		cudaSetDevice(i);

		cudaFreeHost(singpu[i].h_p);
		cudaFreeHost(singpu[i].h_img);
		cudaFreeHost(singpu[i].h_ssg);

		cufftDestroy(singpu[i].plan2dforward);
		cufftDestroy(singpu[i].plan2dinverse);

		cudaFree(singpu[i].d_p);
		cudaFree(singpu[i].d_vx);cudaFree(singpu[i].d_vz);
		cudaFree(singpu[i].d_fivxx);cudaFree(singpu[i].d_fivzz);
		cudaFree(singpu[i].d_fipx);cudaFree(singpu[i].d_fipz);
		
		cudaFree(singpu[i].d_vp);

		cudaFree(singpu[i].d_rho);
		cudaFree(singpu[i].d_bbx);cudaFree(singpu[i].d_bbz);
		
		cudaFree(singpu[i].d_rik);

		cudaFree(singpu[i].d_obs);
		cudaFree(singpu[i].d_img);
		cudaFree(singpu[i].d_ssg);
		cudaFree(singpu[i].d_imgf);
		cudaFree(singpu[i].d_k2d1);

		
		
		cudaFree(singpu[i].d_i_ax);cudaFree(singpu[i].d_i_bx);
		cudaFree(singpu[i].d_i_az);cudaFree(singpu[i].d_i_bz);
		cudaFree(singpu[i].d_h_ax);cudaFree(singpu[i].d_h_bx);
		cudaFree(singpu[i].d_h_az);cudaFree(singpu[i].d_h_bz);
		
		cudaFree(singpu[i].d_dis0);cudaFree(singpu[i].d_dis1);
		cudaFree(singpu[i].d_dis2);cudaFree(singpu[i].d_dis3);cudaFree(singpu[i].d_dis4);
		
		
		cudaFree(singpu[i].d_flagrec);cudaFree(singpu[i].d_rindex);



		// free stored variables //
		cudaFree(singpu[i].d_pb);cudaFree(singpu[i].d_vxb);cudaFree(singpu[i].d_vzb);
		cudaFree(singpu[i].d_fivxxb);cudaFree(singpu[i].d_fivzzb);
		cudaFree(singpu[i].d_fipxb);cudaFree(singpu[i].d_fipzb);
		
		cudaFree(singpu[i].store_p);
		cudaFree(singpu[i].store_vx);cudaFree(singpu[i].store_vz);
		cudaFree(singpu[i].store_fipx);cudaFree(singpu[i].store_fipz);
		cudaFree(singpu[i].store_fivxx);cudaFree(singpu[i].store_fivzz);
     }
	

     // destroy the multi-streams //
     for(i=0;i<GPU_N;i++)
     {
	 cudaSetDevice(i);	
	 cudaStreamDestroy(plansm[i].stream);
     }

    return;
	                        	                        
}// end function




void Alloc_host_device_mem(
				int num,int nrmx,int nx,int nz,
				int nxx,int nzz,int nxh,int nzh,int NXZ,int Nh,
			   	int tmax,int GPU_N,struct MultiGPU singpu[]
			  )
{
	int i;
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);

		cudaMalloc((void**)&singpu[i].d_p,sizeof(float)*NXZ);
		cudaMalloc((void**)&singpu[i].d_vx,sizeof(float)*NXZ);
		cudaMalloc((void**)&singpu[i].d_vz,sizeof(float)*NXZ);
		
		//cudaMalloc((void**)&singpu[i].d_aux,sizeof(float)*NXZ);
		//cudaMalloc((void**)&singpu[i].d_total,sizeof(float)*NXZ);
    		
    		cudaMalloc((void**)&singpu[i].d_fipx,sizeof(float)*NXZ);
    		cudaMalloc((void**)&singpu[i].d_fipz,sizeof(float)*NXZ);
    		
    		cudaMalloc((void**)&singpu[i].d_fivxx,sizeof(float)*NXZ);
    		cudaMalloc((void**)&singpu[i].d_fivzz,sizeof(float)*NXZ);
    
		cudaMalloc((void**)&singpu[i].d_vp,sizeof(float)*NXZ);
    		
    		
    		cudaMalloc((void**)&singpu[i].d_rho,sizeof(float)*NXZ);
    		cudaMalloc((void**)&singpu[i].d_bbx,sizeof(float)*NXZ);
    		cudaMalloc((void**)&singpu[i].d_bbz,sizeof(float)*NXZ);
 
    		cudaMalloc((void**)&singpu[i].d_rik,sizeof(float)*tmax);
    		//cudaMalloc((void**)&singpu[i].d_k2da,sizeof(float)*NXZ);
		cudaMalloc((void**)&singpu[i].d_k2d1,sizeof(float)*NXZ);
		//cudaMalloc((void**)&singpu[i].d_k2d2,sizeof(float)*NXZ);
		//cudaMalloc((void**)&singpu[i].d_k2d3,sizeof(float)*NXZ);
		    
		cudaMalloc((void**)&singpu[i].d_i_ax,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_i_bx,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_i_az,sizeof(float)*nzz);
		cudaMalloc((void**)&singpu[i].d_i_bz,sizeof(float)*nzz);
		
		cudaMalloc((void**)&singpu[i].d_h_ax,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_h_bx,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_h_az,sizeof(float)*nzz);
		cudaMalloc((void**)&singpu[i].d_h_bz,sizeof(float)*nzz);
		
		cudaMalloc((void**)&singpu[i].d_kxsmr,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_kxscr,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_kxsmi,sizeof(float)*nxx);
		cudaMalloc((void**)&singpu[i].d_kxsci,sizeof(float)*nxx);
		
		cudaMalloc((void**)&singpu[i].d_kzsmr,sizeof(float)*nzz);
		cudaMalloc((void**)&singpu[i].d_kzscr,sizeof(float)*nzz);
		cudaMalloc((void**)&singpu[i].d_kzsmi,sizeof(float)*nzz);
		cudaMalloc((void**)&singpu[i].d_kzsci,sizeof(float)*nzz);

		cudaMalloc((void**)&singpu[i].d_dis0,sizeof(cufftComplex)*NXZ);
		cudaMalloc((void**)&singpu[i].d_dis1,sizeof(cufftComplex)*NXZ);
		cudaMalloc((void**)&singpu[i].d_dis2,sizeof(cufftComplex)*NXZ);
		cudaMalloc((void**)&singpu[i].d_dis3,sizeof(cufftComplex)*NXZ);
		cudaMalloc((void**)&singpu[i].d_dis4,sizeof(cufftComplex)*NXZ);


	       //host device memory //
               cudaMallocHost((void **)&singpu[i].h_p, sizeof(float)*NXZ);
               cudaMallocHost((void **)&singpu[i].h_img, sizeof(float)*nx*nz);
               cudaMallocHost((void **)&singpu[i].h_ssg, sizeof(float)*nx*nz);
               

	       cufftPlan2d(&singpu[i].plan2dforward,nxx,nzz,CUFFT_C2C);
	       cufftPlan2d(&singpu[i].plan2dinverse,nxx,nzz,CUFFT_C2C);
	       
	       
	       // for summing gradient //
	       cudaMalloc((void**)&singpu[i].d_img,sizeof(float)*nx*nz);
	       cudaMalloc((void**)&singpu[i].d_imgf,sizeof(float)*nx*nz);
	       
	       cudaMalloc((void**)&singpu[i].d_ssg,sizeof(float)*nx*nz);
	       cudaMalloc((void**)&singpu[i].d_obs,sizeof(float)*nrmx*tmax);
	       
	       cudaMalloc((void**)&singpu[i].d_flagrec,sizeof(int)*NXZ);
	       cudaMalloc((void**)&singpu[i].d_rindex, sizeof(int)*NXZ);
	       
	       
	        // for wavefield reconstruction //
	       cudaMalloc((void**)&singpu[i].d_pb,sizeof(cufftComplex)*NXZ);
	       cudaMalloc((void**)&singpu[i].d_vxb,sizeof(cufftComplex)*NXZ);
	       cudaMalloc((void**)&singpu[i].d_vzb,sizeof(cufftComplex)*NXZ);
	      // cudaMalloc((void**)&singpu[i].d_auxb,sizeof(cufftComplex)*NXZ);
	       
	       cudaMalloc((void**)&singpu[i].d_fipxb,sizeof(float)*NXZ);
    	       cudaMalloc((void**)&singpu[i].d_fipzb,sizeof(float)*NXZ);
    	       cudaMalloc((void**)&singpu[i].d_fivxxb,sizeof(float)*NXZ);
    	       cudaMalloc((void**)&singpu[i].d_fivzzb,sizeof(float)*NXZ);
    	       
    	       cudaMalloc((void**)&singpu[i].store_p,sizeof(float)*NXZ*num);
    	       cudaMalloc((void**)&singpu[i].store_vx,sizeof(float)*NXZ*num);
    	       cudaMalloc((void**)&singpu[i].store_vz,sizeof(float)*NXZ*num);
    	       
    	       cudaMalloc((void**)&singpu[i].store_fipx,sizeof(float)*NXZ*num);
    	       cudaMalloc((void**)&singpu[i].store_fipz,sizeof(float)*NXZ*num);
    	       cudaMalloc((void**)&singpu[i].store_fivxx,sizeof(float)*NXZ*num);
    	       cudaMalloc((void**)&singpu[i].store_fivzz,sizeof(float)*NXZ*num);
    	      // cudaMalloc((void**)&singpu[i].store_aux,sizeof(float)*NXZ*num);
	       

	}
	
	return;
	
}



__global__ void initialize_kernel(
					int nxx,int nzz,int pml,
					float *d_p,float *d_vx,float *d_vz,
					float *d_fipx,float *d_fipz,
					float *d_fivxx,float *d_fivzz,
					float *d_img,float *d_ssg,
					float *d_pb,float *d_vxb,float *d_vzb,
					float *d_fipxb,float *d_fipzb,
					float *d_fivxxb,float *d_fivzzb
				)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       int is=ix*nzz+iz;
       
       int nx=nxx-2*pml;
       int nz=nzz-2*pml;

       if(iz<nzz&&ix<nxx)
       {
           d_p[is]=0.0;
           d_vx[is]=0.0;
           d_vz[is]=0.0;
           d_fipx[is]=0.0;
           d_fipz[is]=0.0;
           d_fivxx[is]=0.0;
           d_fivzz[is]=0.0;
          
           
           d_pb[is]=0.0;
           d_vxb[is]=0.0;
           d_vzb[is]=0.0;
           d_fipxb[is]=0.0;
           d_fipzb[is]=0.0;
           d_fivxxb[is]=0.0;
           d_fivzzb[is]=0.0;
       }
       
       if( iz<nz&&ix<nx&&ix>=0&&iz>=0 )
       {
       		d_img[ix*nz+iz]=0.0;
       		d_ssg[ix*nz+iz]=0.0;
       }

       __syncthreads();
          
}


__global__ void initialize_kernel_ins(
					int nxx,int nzz,int pml,
					float *d_p,float *d_vx,float *d_vz,
					float *d_fipx,float *d_fipz,
					float *d_fivxx,float *d_fivzz,float *d_aux
				)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       int is=ix*nzz+iz;

       if(iz<nzz&&ix<nxx)
       {
           d_p[is]=0.0;
           d_vx[is]=0.0;
           d_vz[is]=0.0;
           d_fipx[is]=0.0;
           d_fipz[is]=0.0;
           d_fivxx[is]=0.0;
           d_fivzz[is]=0.0;
           
           d_aux[is]=0.0;
       }
       __syncthreads();  
}


__global__ void fuzhi_kernel(int nxx,int nzz,float *d_vx,float *d_vz, cufftComplex *d_dis0, cufftComplex *d_dis1)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       
       if(iz<nzz&&ix<nxx)
       {
       	  d_dis0[area1].x=d_vx[area1];
          d_dis0[area1].y=0.0;
       	  
          d_dis1[area1].x=d_vz[area1];
          d_dis1[area1].y=0.0;
       }
       __syncthreads();     
}


__global__ void fuzhi_kernel_p (int nxx,int nzz,float *d_p,cufftComplex *d_dis0)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       if(iz<nzz&&ix<nxx)
       {
       	  d_dis0[area1].x=d_p[area1];
          d_dis0[area1].y=0.0;
       }
       __syncthreads();     
}


__global__ void shift_vxz_kernel(
					float *kxsmr,float *kxsmi,float *kzscr,float *kzsci,
     					cufftComplex *d_outdatax,cufftComplex *d_outdataz,
     					cufftComplex *d_indatax,cufftComplex *d_indataz,
     					int nxx,int nzz
     				)
{    
     int bx=blockIdx.x;
     int by=blockIdx.y;
     
     int tx=threadIdx.x;
     int ty=threadIdx.y;
     
     int iz=by*BLOCK_HEIGHT+ty;
     int ix=bx*BLOCK_WIDTH+tx;
     
     //int area1=iz*nxx+ix;
     int area2=ix*nzz+iz;
     
     if(iz<nzz&&ix<nxx)
     {   
         if(ix<nxx/2+1)
         {
         
             d_indatax[area2].x=-kxsmr[ix]*d_outdatax[area2].y
                                    -kxsmi[ix]*d_outdatax[area2].x;
             d_indatax[area2].y=kxsmr[ix]*d_outdatax[area2].x
                                   -kxsmi[ix]*d_outdatax[area2].y;
         }
         if(ix>nxx/2)
         {   
             d_indatax[area2].x=kxsmr[ix]*d_outdatax[area2].y
                                   +kxsmi[ix]*d_outdatax[area2].x;
                                   
             d_indatax[area2].y=-kxsmr[ix]*d_outdatax[area2].x
                                    +kxsmi[ix]*d_outdatax[area2].y;
         }
         d_indatax[area2].x=d_indatax[area2].x/nxx/nzz;
         d_indatax[area2].y=d_indatax[area2].y/nxx/nzz;
         
         
         if(iz<nzz/2+1)
         {
             d_indataz[area2].x=-kzscr[iz]*d_outdataz[area2].y
                                    -kzsci[iz]*d_outdataz[area2].x;
             d_indataz[area2].y=kzscr[iz]*d_outdataz[area2].x
                                   -kzsci[iz]*d_outdataz[area2].y;
         }
         if(iz>nzz/2)
         {
             d_indataz[area2].x=kzscr[iz]*d_outdataz[area2].y
                                   +kzsci[iz]*d_outdataz[area2].x;
             d_indataz[area2].y=-kzscr[iz]*d_outdataz[area2].x
                                    +kzsci[iz]*d_outdataz[area2].y;
         }

         d_indataz[area2].x=d_indataz[area2].x/nxx/nzz;
         d_indataz[area2].y=d_indataz[area2].y/nxx/nzz;
     }
     __syncthreads();
}


__global__ void update_vxx_vzz(
       				int nxx,int nzz,float *d_fi_vx,float *d_fi_vz,
       				float *di_a_x,float *di_b_x,float *di_a_z,float *di_b_z,
       				cufftComplex *d_outdatax,cufftComplex *d_outdataz
       			       )
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
       float sumx,sumz;
       
       int area2=ix*nzz+iz;
       int area1=area2;
       
       if(iz<nzz&&ix<nxx)
       {
           // calculate fi_vx
           sumx=d_outdatax[area2].x;        
           d_fi_vx[area1]=di_b_x[ix]*d_fi_vx[area1]+di_a_x[ix]*sumx;
           d_outdatax[area2].x=d_outdatax[area2].x+d_fi_vx[area1];
           
           //calculate fi_vz
           sumz=d_outdataz[area1].x;        
           d_fi_vz[area1]=di_b_z[iz]*d_fi_vz[area1]+di_a_z[iz]*sumz;
           d_outdataz[area1].x=d_outdataz[area1].x+d_fi_vz[area1];
       }
       __syncthreads();
}


__global__ void shift_pxz_kernel(
					float *kxscr,float *kxsci,float *kzsmr,float *kzsmi,
     					cufftComplex *d_dis0,
     					cufftComplex *d_dis1,cufftComplex *d_dis2,
     					int nxx,int nzz
     				)
{   
     int bx=blockIdx.x;
     int by=blockIdx.y;
     
     int tx=threadIdx.x;
     int ty=threadIdx.y;
     
     int iz=by*BLOCK_HEIGHT+ty;
     int ix=bx*BLOCK_WIDTH+tx;
     
     int area2=ix*nzz+iz;
     int area1=area2;
     
     if(iz<nzz&&ix<nxx)
     {   
         if(ix<nxx/2+1)
         {
         
             d_dis1[area1].x=-kxscr[ix]*d_dis0[area1].y
                                    -kxsci[ix]*d_dis0[area1].x;
             d_dis1[area1].y=kxscr[ix]*d_dis0[area1].x
                                   -kxsci[ix]*d_dis0[area1].y;
         }
         if(ix>nxx/2)
         {   
             d_dis1[area1].x=kxscr[ix]*d_dis0[area1].y
                                   +kxsci[ix]*d_dis0[area1].x;
                                   
             d_dis1[area1].y=-kxscr[ix]*d_dis0[area1].x
                                    +kxsci[ix]*d_dis0[area1].y;
         }
         d_dis1[area1].x=d_dis1[area1].x/nxx/nzz;
         d_dis1[area1].y=d_dis1[area1].y/nxx/nzz;
         
         
         if(iz<nzz/2+1)
         {
             d_dis2[area2].x=-kzsmr[iz]*d_dis0[area2].y
                                    -kzsmi[iz]*d_dis0[area2].x;
             d_dis2[area2].y=kzsmr[iz]*d_dis0[area2].x
                                   -kzsmi[iz]*d_dis0[area2].y;
         }
         if(iz>nzz/2)
         {
             d_dis2[area2].x=kzsmr[iz]*d_dis0[area2].y
                                   +kzsmi[iz]*d_dis0[area2].x;
             d_dis2[area2].y=-kzsmr[iz]*d_dis0[area2].x
                                    +kzsmi[iz]*d_dis0[area2].y;
         }

         d_dis2[area2].x=d_dis2[area2].x/nxx/nzz;
         d_dis2[area2].y=d_dis2[area2].y/nxx/nzz;
     }
     __syncthreads();
} 


__global__ void vx_cal_kernel(
     				int nxx,int nzz,float dt,float *d_bbx,
     				float *d_fi_px,float *dh_a_x,
     				float *dh_b_x,cufftComplex *d_outdatax,float *d_vx)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
       float sumx;
       
       int area2=ix*nzz+iz;
       int area1=area2;
       
       if(iz<nzz&&ix<nxx)
       {    
           // calculate fi_px
           sumx=d_outdatax[area1].x;
           d_fi_px[area2]=dh_b_x[ix]*d_fi_px[area2]+dh_a_x[ix]*sumx;
           
           // calculate vx
           d_vx[area2]=d_vx[area2]+dt*d_bbx[area1]*(d_outdatax[area1].x+d_fi_px[area2]);
       }
       __syncthreads();   
}

__global__ void vz_cal_kernel(
				int nxx,int nzz,float dt,float *d_bbz,
				float *d_fi_pz,float *dh_a_z,
    				float *dh_b_z,cufftComplex *d_outdataz,float *d_vz)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
       float sumz;
       
       int area1=ix*nzz+iz;
       
       if(iz<nzz&&ix<nxx)
       {    
           // calculate fi_pz
           sumz=d_outdataz[area1].x;
           d_fi_pz[area1]=dh_b_z[iz]*d_fi_pz[area1]+dh_a_z[iz]*sumz;
           
           // calculate vz
           d_vz[area1]=d_vz[area1]+dt*d_bbz[area1]*(d_outdataz[area1].x+d_fi_pz[area1]);
       }
       __syncthreads();   
}



__global__ void multv_kernel(int N1,int nxx,int nzz,int ig,float *d_V,cufftComplex *d_psk,cufftComplex *d_pskf)
{        
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       int NXZ=N1;
       int area2=ig*NXZ+area1;

       if(iz<nzz&&ix<nxx)
       {
          d_pskf[area1].x=d_psk[area1].x*d_V[area2]/NXZ;
          d_pskf[area1].y=d_psk[area1].y*d_V[area2]/NXZ;

       }
       __syncthreads();     
}


__global__ void fuzhi_to_ifftrow(int ig,int nxx,int nzz,int N1,cufftComplex *d_pskf,cufftComplex *d_pifft)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       int area2=ig*N1+area1;

       if(iz<nzz&&ix<nxx)
       {
          d_pifft[area2].x=d_pskf[area1].x;
          d_pifft[area2].y=d_pskf[area1].y;
       }
       __syncthreads();     
}


__global__ void fuzhi_r_kernel(int nxx,int nzz,int jg,int N1,double *d_pst)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       int area2=jg*N1+area1;

       if(iz<nzz&&ix<nxx)
       {
          d_pst[area2]=0.0;
       }
       __syncthreads();     
}

__global__ void sumlw_kernel(int nxx,int nzz,int jg,int N1,int row,int ig,double *d_pst,double *d_A,cufftComplex *d_pifft)
{

       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       int area2=jg*N1+area1;
       int area3=ig*N1+area1;

       if(iz<nzz&&ix<nxx)
       {
          d_pst[area2]=d_pst[area2]+d_A[jg*row+ig]*d_pifft[area3].x;
       }
       __syncthreads();     

}

__global__ void sumlw2_kernel(int nxx,int nzz,int N1,int col,double *d_pst,float *d_U,float *d_psc)
{

       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       int j,is;
       double temp;

       if(iz<nzz&&ix<nxx)
       {
          temp=0.0;
          for(j=0;j<col;j++)
          {
              is=area1*col+j;
              temp=temp+d_U[is]*d_pst[j*N1+area1];
          }
          d_psc[area1]=temp;
          
       }
       __syncthreads();     

}




__global__ void addsa_kernel(
				int restart,int pml,int nxx,int nzz,
				float *d_vp,float *d_rho,float dt,
				float *d_p,cufftComplex *d_dis0,cufftComplex *d_dis3,
				int xs,int zs,int nt,float *d_Rik,float *d_ssg
			    )
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       
       int nx=nxx-2*pml;
       int nz=nzz-2*pml;
       
       int ixx=ix-pml;
       int izz=iz-pml;
       int area2=ixx*nz+izz;
       
       if( iz<nzz&&ix<nxx )
       {
	   d_p[area1]=d_p[area1]+d_rho[area1]*d_vp[area1]*d_vp[area1]*dt*(d_dis0[area1].x-1.0/12.0*dt*d_vp[area1]*dt*d_vp[area1]*d_dis3[area1].x);

          
          
          if(iz==zs&&ix==xs){d_p[area1]=d_p[area1]+d_Rik[nt]*dt*1000.0;}
          
          if( restart==0&& ixx>=0&&ixx<nx&&izz>=0&&izz<nz ){ d_ssg[area2] = d_ssg[area2] + d_p[area1]*d_p[area1]; }
       }
       __syncthreads();     
}


__global__ void addsa_b_kernel(
				int pml,int nxx,int nzz,
				float *d_vp,float *d_rho,float dt,
				float *d_p,cufftComplex *d_dis0,cufftComplex *d_dis3,
				int nt,int tmax,float *d_obs,
				int *rindex, int *flagrec,int nr
			    )
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       int area1=ix*nzz+iz;
       
       
       if( iz<nzz && ix<nxx && ix>=0 && iz>=0 )
       {

	  d_p[area1]=d_p[area1]+d_rho[area1]*d_vp[area1]*d_vp[area1]*dt*(d_dis0[area1].x-1.0/12.0*dt*d_vp[area1]*dt*d_vp[area1]*d_dis3[area1].x);
          
          if( flagrec[area1]==1 ){ d_p[area1]= d_obs[ (tmax-1-nt)*nr+rindex[area1] ]; }  
       }
       
       
       __syncthreads();     
	
}

__global__ void rtm_image(
				int nx,int nz,int nxx,int nzz,
				float *d_p,float *d_pb,
				float *d_img,float *d_imgf
			 )
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
       int pml=(nxx-nx)/2;
       
       int area1 = ix*nzz+iz;
       int area2 = (ix-pml)*nz+iz-pml;
       
       if ( ix>=pml&&ix<nxx-pml&&iz>=pml&&iz<nzz-pml )
       {
       		d_img[area2] = d_img[area2] + d_p[area1]*d_pb[area1];
       		
       		d_imgf[area2] = d_img[area2];
       }
	
}
                      				
__global__ void filter_image(int nx,int nz,float dx,float dz,float *d_img,float *d_imgf)
{
	int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
       int area1 = ix*nz+iz;
       
       if ( ix>=1&&ix<nx-1&&iz>=1&&iz<nz-1 )
       {
       		d_imgf[area1] = (d_img[area1+1] + d_img[area1-1] -2.0*d_img[area1])/dz/dz
       				+ (d_img[area1+nz] + d_img[area1-nz]-2.0*d_img[area1])/dx/dx;
       }
       
       if (ix==0||iz==0||ix==nx-1||iz==nz-1){d_imgf[area1]=0.0;}
}


__global__ void multt_kernel(float eta,float dt,int nt,float tao,
			     int NXZ,int nxx,int nzz,
			     float *k2da,cufftComplex *d_dis0
			     )
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       float tt=nt*dt;
       float tmp;

       if(iz<nzz&&ix<nxx)
       {
       	     
       	  tmp=k2da[area1]-2.0*eta*k2da[area1]*(1.0+tao*dt*k2da[area1])/10000.0/( expf(-tao*tt*k2da[area1]) + eta*expf(tao*dt*k2da[area1])/10000.0 );
       	  //tmp=k2da[area1];
       	  
          d_dis0[area1].x=d_dis0[area1].x*tmp/NXZ;
          d_dis0[area1].y=d_dis0[area1].y*tmp/NXZ;
       }
	
	__syncthreads();
}

__global__ void fuzhi_kernel_inverse(int nxx,int nzz,float *d_aux,cufftComplex *d_dis0)
{
	int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       if(iz<nzz&&ix<nxx)
       {
       	  d_aux[area1]=d_dis0[area1].x;
       }
       __syncthreads();     
	
}

__global__ void storecheckpoints
		(
			int flag,int nc,int nxx,int nzz,
			float *d_p,float *d_vx,float *d_vz,
			float *d_fipx,float *d_fipz,
			float *d_fivxx,float *d_fivzz,
			float *store_p,float *store_vx,float *store_vz,
			float *store_fipx,float *store_fipz,
			float *store_fivxx,float *store_fivzz
		)
{
	 int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
       int area=ix*nzz+iz;
       int area1=nc*nxx*nzz+area;
       
       if (ix>=0&&ix<nxx&&iz>=0&&iz<nzz)
       {
       		if(flag==1)
       		{
       			store_p[area1]  = d_p[area];
       			store_vx[area1] = d_vx[area];
       			store_vz[area1] = d_vz[area];
       			store_fipx[area1] = d_fipx[area];
       			store_fipz[area1] = d_fipz[area];
       			store_fivxx[area1] = d_fivxx[area];
       			store_fivzz[area1] = d_fivzz[area];
       		}
       		else
       		{
       			d_p[area]  = store_p[area1];
       			d_vx[area] = store_vx[area1];
       			d_vz[area] = store_vz[area1];
       			
       			d_fipx[area]  = store_fipx[area1];
       			d_fipz[area]  = store_fipz[area1];
       			d_fivxx[area] = store_fivxx[area1];
       			d_fivzz[area] = store_fivzz[area1];
       		}
       }
}

__global__ void addsource(int nx,int nz,int xs,int zs,float *check_d_ps,float *d_rik,float dt)
{
	int bx=blockIdx.x;
        int by=blockIdx.y;
     
	int tx=threadIdx.x;
        int ty=threadIdx.y;
     
        int iz=by*BLOCK_HEIGHT+ty;
        int ix=bx*BLOCK_WIDTH+tx;

	int is=ix*nz+iz;

	if(iz==zs&&ix==xs)
	{
		check_d_ps[is]=check_d_ps[is]+d_rik[0]*dt*1000.0;
	}
}

__global__ void sumpxpz(
       				int nxx,int nzz,
       				cufftComplex *d_outdatax,cufftComplex *d_outdataz,
				cufftComplex *d_dis0
       			       )
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
      
       
       int area2=ix*nzz+iz;
       //int area1=area2;
       
       if(iz<nzz&&ix<nxx)
       {
		d_dis0[area2].x=d_outdatax[area2].x+d_outdataz[area2].x;
		d_dis0[area2].y=d_outdatax[area2].y+d_outdataz[area2].y;
		//d_dis0[area2].y=0.0;
       }
       __syncthreads();
}

__global__ void multtt_kernel(int nxx,int nzz,float *d_k2d3,cufftComplex *d_dis0,cufftComplex *d_dis1)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;

       int area1=ix*nzz+iz;
       int NXZ=nxx*nzz;

       if(iz<nzz&&ix<nxx)
       {
          d_dis1[area1].x=d_dis0[area1].x*d_k2d3[area1]/NXZ;
          d_dis1[area1].y=d_dis0[area1].y*d_k2d3[area1]/NXZ;

       }
       __syncthreads();     
}

__global__ void  sum_total(int nxx,int nzz,
			   cufftComplex *d_dis0,cufftComplex *d_dis2,
			   cufftComplex *d_dis3,cufftComplex *d_dis4,
			   float *d_total,float *d_vp,float *d_lam,
			   float *d_gam,float *d_vwd,float dt)
{
       int bx=blockIdx.x;
       int by=blockIdx.y;
     
       int tx=threadIdx.x;
       int ty=threadIdx.y;
     
       int iz=by*BLOCK_HEIGHT+ty;
       int ix=bx*BLOCK_WIDTH+tx;
       
        int area2=ix*nzz+iz;
	float y1,y2,y3,y4;
	float dtf=dt*dt;
	float vpf;

	if(iz<nzz&&ix<nxx)
       {
		y1=1.0-32.0*d_gam[area2];
		y2=32.0*d_gam[area2]*d_vwd[area2];
		y3=1.0-64.0*d_gam[area2];
		y4=64.0*d_gam[area2]*d_vwd[area2];
		vpf=d_vp[area2]*d_vp[area2];
	
	d_total[area2]= y1*d_dis0[area2].x+y2*d_dis2[area2].x-1.0/12.0*dtf*vpf*d_lam[area2]*(y3*d_dis3[area2].x+y4*d_dis4[area2].x);
       }
       __syncthreads();
	
}
