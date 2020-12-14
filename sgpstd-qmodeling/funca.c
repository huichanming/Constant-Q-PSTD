/* This file inludes all auxilary c functions serveing for forward modeling */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define pi 3.1415926


extern void input_ac2d_parameters(int *shotnum,int *tmax,float *dt,float *fd,
                                  int *nx,int *nz,float *dx,float *dz,float *disx,
                                  char fname1[40],char fname2[40],char fname3[40],
                                  char fname4[40],char fname5[40],
                                  char vpname[80],char rhoname[80],
                                  int *flagout,int *flagsource
                                  )
{
	char strtmp[256];
	FILE *fp=fopen("../input/ac2drealmodeling.txt","r");
	if(fp==0){printf("Cannot open ../ac2drealmodeling.txt!\n");exit(0);}
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",shotnum);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",tmax);
	fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",dt);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",fd);
	fscanf(fp,"\n");
	

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",nx);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",nz);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",dx);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",dz);
	fscanf(fp,"\n");
  
        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",disx);
	fscanf(fp,"\n"); 


	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",fname1);
	fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",fname2);
	fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",fname3);
	fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",fname4);
	fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",fname5);
	fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",vpname);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%s",rhoname);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",flagout);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",flagsource);
	fclose(fp);
	
	return;
}



extern void input_2d_geometry(int shotnum,float *shotxlc,float *shotzlc,float *recdepth)
{
    float x0,xinter,shotdep;
	char strtmp[256];
	FILE *fp=fopen("../input/geometry2d.txt","r");
	if(fp==0){printf("Cannot open the geometry2d.txt!\n");exit(0);}

    int i;
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&x0);
	fscanf(fp,"\n");
    
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&xinter);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&shotdep);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",recdepth);
	fclose(fp);

    for(i=0;i<shotnum;i++)
	{
		shotxlc[i]=x0+i*xinter;
		shotzlc[i]=shotdep;
	}
	
	return;
}

extern void cal_shot_rec_grid_index(int shotnum,float recdepth,float dx,float dz,float *shotxlc,float *shotzlc,
                                    int *shotxg,int *shotzg,int *recg)
{
	int i;
	for(i=0;i<shotnum;i++)
	{
		shotxg[i]=(int)(shotxlc[i]/dx);
		shotzg[i]=(int)(shotzlc[i]/dz);
	}
	*recg=(int)(recdepth/dz);
	
	return;
}

extern void cal_bbx_bbz(int nxx,int nzz,float *dens,float *bbx,float *bbz)
{
	int i,j,is;
	
	for(j=0;j<nzz;j++)
    	{
		for(i=0;i<nxx-1;i++)
		{
		    is=i*nzz+j;
		    bbx[is]=2.0/(dens[is]+dens[is+nzz]);
		}
     	}
    	for(j=0;j<nzz;j++)
    	{
        	bbx[(nxx-1)*nzz+j]=1.0/dens[(nxx-1)*nzz+j];
    	}
    	for(i=0;i<nxx;i++)
    	{
		for(j=0;j<nzz-1;j++)
		{
			is=i*nzz+j;
		    bbz[is]=2.0/(dens[is]+dens[is+1]);
		}
    	}
    	for(i=0;i<nxx;i++)
    	{
        	bbz[i*nzz+nzz-1]=1.0/dens[i*nzz+nzz-1];
    	}
}



extern void rik_cal(int tmax,float dt,float fd,float *rik,int flag)
{
	int i;
	float t0=1.0/fd;
	float temp;
	if(flag==2){
    for(i=0;i<tmax;i++)
	{
		temp=pow(pi*fd*(i*dt-t0),2.0);
		rik[i]=(1.0-2.0*temp)*exp(-temp);
	}}
    else{
   // 1st-order derivative of Gaussian function as source, not ricker wavelet //
    for(i=0;i<tmax;i++)
    {
		temp=pow(2.0*pi*fd*(i*dt-t0),2.0);
	    rik[i]=(i*dt-t0)*exp(-temp);
     }}
	 return;
}


extern void get_pml_parameters(
					int pml,float fm,
					float *i_a_x,float *i_b_x,float *i_a_z,float *i_b_z,
					float *h_a_x,float *h_b_x,float *h_a_z,float *h_b_z,
                               		int nxx,int nzz,float dx,float dz,float dt
                               )
{
     	float alpha_max=pi*fm,dmaxx,dmaxz,widthx,widthz;
        float Re=1e-5,vmax=5000.0;
        int n1=2,n2=1,n3=1,i,j;
        float temp,temp1,temp2,temp3,temp4;
        
        float *int_dx,*int_dz,*int_alphax,*int_alphaz;
        float *half_dx,*half_dz,*half_alphax,*half_alphaz;
        
        int_dx=(float*)malloc(sizeof(float)*nxx);
        half_dx=(float*)malloc(sizeof(float)*nxx);
        int_alphax=(float*)malloc(sizeof(float)*nxx);
        half_alphax=(float*)malloc(sizeof(float)*nxx);
        
        int_dz=(float*)malloc(sizeof(float)*nzz);
        half_dz=(float*)malloc(sizeof(float)*nzz);
        int_alphaz=(float*)malloc(sizeof(float)*nzz);
        half_alphaz=(float*)malloc(sizeof(float)*nzz);
          
        widthx=pml*dx;
        widthz=pml*dz;
        dmaxx=(1+n1+n2)*vmax*log(1.0/Re)/(2.0*widthx);
        dmaxz=(1+n1+n2)*vmax*log(1.0/Re)/(2.0*widthz);
     
        // integer absorbing parameters
        for(i=0;i<pml;i++)
        {
            temp1=pow(1.0*(pml-1-i)/pml,n1+n2);
            int_dx[i]=dmaxx*temp1;
            int_dz[i]=dmaxz*temp1;
         
            temp3=pow(1.0*i/(pml-1),n3);
            int_alphax[i]=alpha_max*temp3;
            int_alphaz[i]=int_alphax[i];

            int_dx[nxx-1-i]=int_dx[i];
            int_dz[nzz-1-i]=int_dz[i];
            
            int_alphax[nxx-1-i]=int_alphax[i];
            int_alphaz[nzz-1-i]=int_alphaz[i];
        }
        for(i=pml;i<nxx-pml;i++)
        {
           int_dx[i]=0.0;
           int_alphax[i]=int_alphax[pml-1];
        }
        for(j=pml;j<nzz-pml;j++)
        {
            int_dz[j]=0.0;
            int_alphaz[j]=int_alphaz[pml-1];
        }
     
        // half absorbing parameters
        for(i=0;i<pml-1;i++)
        {
            temp2=pow(1.0*(pml-1.5-i)/pml,n1+n2);
            half_dx[i]=dmaxx*temp2;
            half_dz[i]=dmaxz*temp2;
            
            half_dx[nxx-2-i]=half_dx[i];
            half_dz[nzz-2-i]=half_dz[i];
         
            temp4=pow(1.0*(i+0.5)/(pml-1),n3);
            half_alphax[i]=alpha_max*temp4;
            half_alphaz[i]=half_alphax[i];
            
            half_alphax[nxx-2-i]=half_alphax[i];
            half_alphaz[nzz-2-i]=half_alphaz[i];
        }
        for(i=pml-1;i<nxx-pml;i++)
        {
             half_dx[i]=0.0;
             half_alphax[i]=half_alphax[pml-2];
        }
        half_dx[nxx-1]=0.0;
        half_alphax[nxx-1]=half_alphax[nxx-2];
     
        for(j=pml-1;j<nzz-pml;j++)
        {
            half_dz[j]=0.0;
            half_alphaz[j]=half_alphaz[pml-2];;
        }
        half_dz[nzz-1]=0.0;
        half_alphaz[nzz-1]=half_alphaz[nzz-2];
    
        for(i=0;i<nxx;i++)
        {
            temp=int_dx[i]+int_alphax[i];
            i_b_x[i]=exp(-dt*temp);
            i_a_x[i]=int_dx[i]/temp*(i_b_x[i]-1.0);
        
            temp=half_dx[i]+half_alphax[i];
            h_b_x[i]=exp(-dt*temp);
            h_a_x[i]=half_dx[i]/temp*(h_b_x[i]-1.0);
        }
     
        for(j=0;j<nzz;j++)
        {
            temp=int_dz[j]+int_alphaz[j];
            i_b_z[j]=exp(-dt*temp);
            i_a_z[j]=int_dz[j]/temp*(i_b_z[j]-1.0);
        
            temp=half_dz[j]+half_alphaz[j];
            h_b_z[j]=exp(-dt*temp);
            h_a_z[j]=half_dz[j]/temp*(h_b_z[j]-1.0);
        }
        
        free(int_dx);free(int_dz);free(int_alphax);free(int_alphaz);
        free(half_dx);free(half_dz);free(half_alphax);free(half_alphaz);
     
     return;
}

extern void cal_shift_parameters(float *kxsmr,float *kxscr,float *kzsmr,float *kzscr,
                          float *kxsmi,float *kxsci,float *kzsmi,float *kzsci,
                          float dx,float dz,int nxx,int nzz)
{
      
    float *kkx,*kkz,*shiftxr,*shiftxi,*shiftzr,*shiftzi;
    int i,j,ir;
    float temp,dkx,dkz;
    
    dkx=2*pi/(nxx*dx);dkz=2*pi/(nzz*dz);
    
    kkx=(float*)malloc(sizeof(float)*nxx);
    kkz=(float*)malloc(sizeof(float)*nzz);
    shiftxr=(float*)malloc(sizeof(float)*nxx);
    shiftxi=(float*)malloc(sizeof(float)*nxx);
    shiftzr=(float*)malloc(sizeof(float)*nzz);
    shiftzi=(float*)malloc(sizeof(float)*nzz);
    
    int nxh=(int)(nxx/2)+1;
    int nzh=(int)(nzz/2)+1;

    for(i=0;i<nxh;i++)
    {
        kkx[i]=i*dkx;
    }
    for(i=nxh;i<nxx;i++)
    {
        kkx[i]=kkx[nxx-i];
    }

    for(j=0;j<nzh;j++)
    {
        kkz[j]=j*dkz;
    }
    for(j=nzh;j<nzz;j++)
    {
        kkz[j]=kkz[nzz-j];
    }
    
    for(i=0;i<nxh;i++)
    {
        temp=i*pi/nxx;
        shiftxr[i]=cos(temp);
        shiftxi[i]=sin(temp);

        kxsmr[i]=kkx[i]*shiftxr[i];
        kxsmi[i]=-kkx[i]*shiftxi[i];
        kxscr[i]=kkx[i]*shiftxr[i];
        kxsci[i]=kkx[i]*shiftxi[i];
     }
     for(i=nxh;i<nxx;i++)
     {
         shiftxr[i]=shiftxr[nxx/2*2-i];
         shiftxi[i]=-shiftxi[nxx/2*2-i];

         kxsmr[i]=kkx[i]*shiftxr[i];
         kxsmi[i]=-kkx[i]*shiftxi[i];
         kxscr[i]=kkx[i]*shiftxr[i];
         kxsci[i]=kkx[i]*shiftxi[i];
      }
      for(j=0;j<nzh;j++)
      {
          temp=j*pi/nzz;
          shiftzr[j]=cos(temp);
          shiftzi[j]=sin(temp);

          kzsmr[j]=kkz[j]*shiftzr[j];
          kzsmi[j]=kkz[j]*shiftzi[j];
          kzscr[j]=kkz[j]*shiftzr[j];
          kzsci[j]=-kkz[j]*shiftzi[j];
       }
       for(j=nzh;j<nzz;j++)
       {
           shiftzr[j]=shiftzr[nzz/2*2-j];
           shiftzi[j]=-shiftzi[nzz/2*2-j];

           kzsmr[j]=kkz[j]*shiftzr[j];
           kzsmi[j]=kkz[j]*shiftzi[j];
           kzscr[j]=kkz[j]*shiftzr[j];
           kzsci[j]=-kkz[j]*shiftzi[j];
        }
       free(kkx);free(kkz);free(shiftxr);free(shiftxi);
       free(shiftzr);free(shiftzi);
       
       return;
}



extern void input_vp_model(int nx,int nz,float *vp,char vpname[])
{

    FILE *fp=fopen(vpname,"rb");
    if(fp==0){printf("Cannot open the velocity file!\n");exit(0);}
	fread(vp,sizeof(float),nx*nz,fp);
	fclose(fp);
	return;
}


extern void input_part_model(int nx,int nz,float *vp,int nxo,int pml,int xg1,int xg2,float *vpc)
{
	int i,j,is,is2;
	int nzz=nz+2*pml;
	int nxx=nxo+2*pml;

	for(i=xg1;i<=xg2;i++)
	{
		for(j=0;j<nz;j++)
		{
			is=i*nz+j;
			is2=(i-xg1+pml)*nzz+j+pml;
			vpc[is2]=vp[is];
		}

		/* fill top and bottom pml areas  */
		for(j=0;j<pml;j++)
		{
			is=(i-xg1+pml)*nzz+pml;
			is2=(i-xg1+pml)*nzz+j;
			vpc[is2]=vpc[is];

			is=(i-xg1+pml)*nzz+nzz-pml-1;
			is2=(i-xg1+pml)*nzz+nzz-pml+j;
			vpc[is2]=vpc[is];
		}
	}

	/* filling left and right PML areas */
	for(j=0;j<nzz;j++)
	{
		for(i=0;i<pml;i++)
		{
			is=i*nzz+j;
			is2=pml*nzz+j;
			vpc[is]=vpc[is2];

			is=(nxx-pml+i)*nzz+j;
			is2=(nxx-pml-1)*nzz+j;
			vpc[is]=vpc[is2];
		}
	}
	return;
}



//for obtaining shots filename//
extern void num2str(char asc[4],int num)
{
     char asc1,asc2,asc3,asc4;
     asc4='\0';
     if(num<10)
     {
        asc3='0';
        asc2='0';
        asc1=(char)(num+48);
      }
      if(num>=10&&num<=99)
      {
         asc2=(char)(num/10+48);
         asc1=(char)(num-num/10*10+48);
         asc3='0';
      }
      if(num>99)
      {
         asc3=(char)(num/100+48);
         asc2=(char)(num%100/10+48);
         asc1=(char)(num%10+48);
      }
      asc[0]=asc3;
      asc[1]=asc2;
      asc[2]=asc1;
      asc[3]=asc4;
      return;
}

extern void get_max_min(int NXZ,float *vp,float *vmax,float *vmin)
{
     int i;
     float mn=10000000.0,mx=-10000000.0;
     for(i=0;i<NXZ;i++)
     {
         if(vp[i]-mn<1.0e-4){mn=vp[i];}
         if(vp[i]-mx>1.0e-4){mx=vp[i];}
     }
     *vmax=mx;
     *vmin=mn;
     return;
}

extern float check_stability_ps(float vmax,float dt,float dx,float dz)
{
      float delta;
      delta=pi*vmax*dt/2.0*sqrt(1.0/dx/dx+1.0/dz/dz);
      return delta;
}

extern void cal_k2d(float dx,float dz,int n0,int n1,float *k2d1)
{
	int n0h=(int)(n0/2)+1;
	int n1h=(int)(n1/2)+1;
	
	float dkx=2*pi/(n0*dx);
	float dkz=2*pi/(n1*dz);

	float *kkx,*kkz;
	kkx=(float*)malloc(sizeof(float)*n0);
	kkz=(float*)malloc(sizeof(float)*n1);

	int i,j,is;
    
    for(i=0;i<n0h;i++){kkx[i]=i*dkx;}
    for(i=n0h;i<n0;i++){kkx[i]=kkx[n0-i];}
    for(j=0;j<n1h;j++){kkz[j]=j*dkz;}
    for(j=n1h;j<n1;j++){kkz[j]=kkz[n1-j];}

    for(i=0;i<n0;i++)
    {
        for(j=0;j<n1;j++)
        {
            is=i*n1+j;
           

	    k2d1[is]=kkx[i]*kkx[i]+kkz[j]*kkz[j];
	   
        }
    }

    free(kkx);free(kkz);
    return;
}


extern void input_geophone_pos(int shotid,int nr,float *gxp,char fname[40])
{
	FILE *fp;
	char nname[60];
	
        strcpy(nname,"");
        sprintf(nname,"%s%d.bin",fname,shotid);
        fp=fopen(nname,"rb");
        if(fp==0){printf("Cannot open receiver positions file!\n");exit(0);}
        fread(gxp,sizeof(float),nr,fp);
        fclose(fp);
	
	return;
}


extern void transform_pos_to_grid_index(float dx,float dz,float xlp2,int pml,int nr,
			                float xs,float zs,float *gxp,float *gzp,
			                int *shotxg,int *shotzg,
				        int *gxg,int *gzg)
{
	int i;
	*shotxg=(int)((xs-xlp2)/dx+0.5)+pml;
	*shotzg=(int)(zs/dz+0.5)+pml;

	for(i=0;i<nr;i++)
	{
		gzg[i]=(int)(gzp[i]/dz+0.5)+pml;
		gxg[i]=(int)((gxp[i]-xlp2)/dx+0.5)+pml;
	}
	
	return;
}


extern int getmaxtracenumber(int N,int *a)
{
	int i;
	int mx=-1;
	for(i=0;i<N;i++)
        {
		if(a[i]>mx){mx=a[i];}
	}
	return mx;
}


extern void select_part_model(int pml,float mlx,float dx,int nxt,int nzt,int nx,int nz,int nxx,int nzz,float *vpt,float *vp)
{
	int i,j,is,is2;
	int flagx0;
	int flagx1;

      
	// --------------------- input model -----------------//
	if(mlx>=0.0) // one case
	{
	    flagx0=(int)(mlx/dx);
	    flagx1=flagx0+nx-1;
	    if(flagx1<nxt)
	    {
		for(i=flagx0;i<=flagx1;i++)
		{
			for(j=0;j<nzt;j++)
			{
                            is=i*nzt+j;
			    is2=(i-flagx0+pml)*nzz+j+pml;
			    vp[is2]=vpt[is];
			}
		}
	    }
	    else
	    {
		for(i=flagx0;i<nxt;i++)
		{
			for(j=0;j<nzt;j++)
			{
			    is=i*nzt+j;
			    is2=(i-flagx0+pml)*nzz+j+pml;
			    vp[is2]=vpt[is];
			}
		}
		for(i=nxt;i<=flagx1;i++)
		{
		    	for(j=pml;j<nzz-pml;j++)
		        {
			    is2=(i-flagx0+pml)*nzz+j;
			    is=(nxt-1-flagx0+pml)*nzz+j;
			    vp[is2]=vp[is];
			}
		}
	    }				
	}


            
        if(mlx<0.0) // another case
	{
	    flagx0=(int)(fabs(mlx)/dx+0.5);
	    flagx1=nx-flagx0;
	    if(flagx1<=nxt)
	    {
		for(i=0;i<flagx1;i++)
		{
                     for(j=0;j<nzt;j++)
		     {
			  is=i*nzt+j;
			  is2=(i+flagx0+pml)*nzz+j+pml;
			  vp[is2]=vpt[is];
		     }	
	        }
		for(i=0;i<flagx0;i++)
		{
		     for(j=pml;j<nzz-pml;j++)
		     {
			  is2=(i+pml)*nzz+j;
			  is=(flagx0+pml)*nzz+j;
                          vp[is2]=vp[is];

		     }
		}
	    }
	    else
	    {
	        for(i=0;i<nxt;i++)
	        {
                    for(j=0;j<nzt;j++)
		    {
		        is=i*nzt+j;
                        is2=(i+flagx0+pml)*nzz+j+pml;
		        vp[is2]=vpt[is];
		    }
	        }
	        for(i=nxt;i<flagx1;i++)
	        {
                    for(j=pml;j<nzz-pml;j++)
		    {
		        is2=(i+flagx0+pml)*nzz+j;
		        is=(nxt-1+flagx0+pml)*nzz+j;
		        vp[is2]=vp[is];
		    }
	        }
	        for(i=0;i<flagx0;i++)
	        {
                    for(j=pml;j<nzz-pml;j++)
		    {
		        is2=(i+pml)*nzz+j;
		        is=(flagx0+pml)*nzz+j;
		        vp[is2]=vp[is];
		    }
	       }				
	   }
       }

         //-----------------------------------------------------//

	return;
}

extern void expand_pml_areas(int nx,int nz,int nxx,int nzz,int pml,float *vp)
{
	int i,j,is,is2;
	for(i=pml;i<nxx-pml;i++)
	{
		for(j=0;j<pml;j++)
		{
			is=i*nzz+j;
			is2=i*nzz+pml;
			vp[is]=vp[is2];
		}
		for(j=nzz-pml;j<nzz;j++)
		{
			is=i*nzz+j;
			is2=i*nzz+nzz-pml-1;
			vp[is]=vp[is2];
		}
	}
	for(j=0;j<nzz;j++)
	{
		for(i=0;i<pml;i++)
		{
			is=i*nzz+j;
			is2=pml*nzz+j;
            vp[is]=vp[is2];
		}
		for(i=nxx-pml;i<nxx;i++)
		{
			is=i*nzz+j;
            is2=(nxx-pml-1)*nzz+j;
			vp[is]=vp[is2];
		}
	}

	return;
}

extern void cal_k2dh(int nxh,int nzh,int n0,int n1,float dx,float dz,float *k2dh)
{
	int n0h=nxh;
	int n1h=nzh;
	
	float dkx=2*pi/(n0*dx);
	float dkz=2*pi/(n1*dz);

	float *kkx,*kkz;
	kkx=(float*)malloc(sizeof(float)*n0);
	kkz=(float*)malloc(sizeof(float)*n1);

	int i,j,is;
    
    	for(i=0;i<n0h;i++){kkx[i]=i*dkx;}
    	for(i=n0h;i<n0;i++){kkx[i]=kkx[n0-i];}
    	for(j=0;j<n1h;j++){kkz[j]=j*dkz;}
    	for(j=n1h;j<n1;j++){kkz[j]=kkz[n1-j];}
	
	for(i=0;i<nxh;i++)
	{
		for(j=0;j<nzh;j++)
		{
			is=i*nzh+j;
			k2dh[is]=sqrt( kkx[i]*kkx[i]+kkz[j]*kkz[j] );
		}
	}
	k2dh[0]=1.0e-6;
}


extern void projectomodel(int pml,float mlx,float dx,int nxt,int nzt,int nx,int nz,int nxx,int nzz,float *vpt,float *vp)
{
	int i,j,is,is2;
	int flagx0;
	int flagx1;

      
	// --------------------- input model -----------------//
	if(mlx>=0.0) // one case
	{
	    flagx0=(int)(mlx/dx);
	    flagx1=flagx0+nx-1;
	    if(flagx1<nxt)
	    {
		for(i=flagx0;i<=flagx1;i++)
		{
			for(j=0;j<nzt;j++)
			{
                            is=i*nzt+j;
			    is2=(i-flagx0)*nz+j;
                            vpt[is]=vp[is2];
			}
		}
	    }
	    else
	    {
		for(i=flagx0;i<nxt;i++)
		{
			for(j=0;j<nzt;j++)
			{
			    is=i*nzt+j;
			    is2=(i-flagx0)*nz+j;
                            vpt[is]=vp[is2];
			}
		}
		
	    }				
	}


            
        if(mlx<0.0) // another case
	{
	    flagx0=(int)(fabs(mlx)/dx+0.5);
	    flagx1=nx-flagx0;
	    if(flagx1<=nxt)
	    {
		for(i=0;i<flagx1;i++)
		{
                     for(j=0;j<nzt;j++)
		     {
			  is=i*nzt+j;
			  is2=(i+flagx0)*nz+j;
                          vpt[is]=vp[is2];
		     }	
	        }
	    }
	    else
	    {
	        for(i=0;i<nxt;i++)
	        {
                    for(j=0;j<nzt;j++)
		    {
		        is=i*nzt+j;
                        is2=(i+flagx0)*nz+j;
                        vpt[is]=vp[is2];
		    }
	        }		
	   }
       }
       //-----------------------------------------------------//

	return;
}

extern void removedirect(int nrmx,int tmax,int nr,float dt,float fd,float xs,float zs,
                         float *gxp,float *gzp,float v0,float *recobs)
{
	int i,j,is;
	int ts,te;
	float dis;
	for(i=0;i<nr;i++)
	{
		dis=sqrt(pow((xs-gxp[i]),2.0)+pow((zs-gzp[i]),2.0));
		ts=(int)(1.0*dis/v0/dt+0.5);
		te=ts+(int)(2.0/fd/dt+0.5);
                if(te>=tmax){te=tmax-1;}

		for(j=ts;j<=te;j++)
		{
			is=j*nr+i;
			recobs[is]=0.0;
		}
	}
	
	return;
}

extern void transfermodel(
				int nx,int nz,float fd,float fo,
				float *vp,float *Qp,float *gam,float *lam,
				float *vwd,float *tao
			 )
{
	int i;
	int N=nx*nz;
	float wd=2*pi*fd;
	
	for(i=0;i<N;i++)
	{
		gam[i]=1.0/pi/Qp[i];
		lam[i]=pow(fd/fo,2.0*gam[i]);
		
		vwd[i] = pow(vp[i]/wd,1.0/16.0);
		tao[i]   = lam[i]*vp[i]/Qp[i];
	}
}


extern float get_max_tao(int nxx,int nzz,float *tao)
{
	int i;
	float mn=0.0;
	for(i=0;i<nxx*nzz;i++)
	{
		if ( tao[i]>=mn ){mn=tao[i];}
	}
	
	return mn;
}

extern float get_min_tao(int nxx,int nzz,float *tao)
{
	int i;
	float mn=1000.0;
	for(i=0;i<nxx*nzz;i++)
	{
		if ( tao[i]<mn ){mn=tao[i];}
	}
	
	return mn;
}

