#define pi 3.141592653
#include "4GINV.C"
#define sinc(x) sin(x)/(x)

void random_select(int nx,int mm,int *r1)
{
     time_t t;
     srand((unsigned)time(&t));
     
     int n,i=0,flag=1;
     do
     {
        r1[i]=rand()%nx;
        flag=1;
        for(n=0;n<i;n++)
        {
            if(r1[i]==r1[n])
            {flag=0;break;}
        }
        if(flag){i=i+1;}
      }while(i<mm);
}

void MultMatrix(int n1,int n2,int n3,double *W1,double *W2,double *W3)
{
    int i,j,k;
    double sum;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n3;j++)
        {
            sum=0.0;
            for(k=0;k<n2;k++)
            {sum=sum+W1[i*n2+k]*W2[k*n3+j];}
            W3[i*n3+j]=sum;
        }
    }
}

float dotxy(float *tv1,float *tv2,int n)
{
     int i=0;
     float sumx=0.0;
     for(i=0;i<n;i++)
     {
         sumx=sumx+tv1[i]*tv2[i];
     }
     return sumx;
}

float dotxx(float *tv1,int n)
{
     int i=0;
     float sumx=0.0;
     for(i=0;i<n;i++)
     {
         sumx=sumx+tv1[i]*tv1[i];
     }
     return sumx;
}


extern void checka_error(int row,int col,int N1,int N2,float *U,double *A,float *V,
			float vp,float gam,float tao,float dt,
                         float *k2dh,float *We,float *Wa,int is)
{
    int i,j;
    float temp;
    
    for(i=0;i<N2;i++)
    {
    	temp=vp*dt*pow(k2dh[i],1.0+gam);
    	We[i]=sinc(0.5*temp)*pow(k2dh[i],gam)*sqrt( 1+0.5*dt*tao*pow(k2dh[i],1.01) );
    }

    float *tp;
    tp=(float*)malloc(sizeof(float)*row);
    double tm1;

    for(j=0;j<row;j++)
    {  
        tm1=0.0;
        for(i=0;i<col;i++)
        {
           tm1=tm1+U[is*col+i]*A[i*row+j];
        }
        tp[j]=tm1;
        
    }
    
    for(j=0;j<N2;j++)
    {
        Wa[j]=0.0;
        for(i=0;i<row;i++)
        {
            Wa[j]=Wa[j]+tp[i]*V[i*N2+j];
        }
    }
    free(tp);
    return;
}


extern void consist_Vt(int row,int nxx,int nzz,int nxh,int nzh,int N1,int N2,float *V,float *Vt)
{
     int n,i,j,is1,is2;
     for(n=0;n<row;n++)
     {
              
         for(i=0;i<nxh;i++)
         {
              for(j=0;j<nzh;j++)
              {
                  is1=n*N2+i*nzh+j;
                  is2=n*N1+i*nzz+j;
                  Vt[is2]=V[is1];  
              }
              // z symetry  
              for(j=nzh;j<nzz;j++)
              {
                      is1=n*N2+i*nzh+nzz-j;
                      is2=n*N1+i*nzz+j;
                      Vt[is2]=V[is1];
              }
         }
         // x symetry
         for(i=nxh;i<nxx;i++)
         {
             for(j=0;j<nzh;j++)
             {
                 is1=n*N2+(nxx-i)*nzh+j;
                 is2=n*N1+i*nzz+j;
                 Vt[is2]=V[is1];
             }
             for(j=nzh;j<nzz;j++)
             {
                 is1=n*N2+(nxx-i)*nzh+nzz-j;
                 is2=n*N1+i*nzz+j;
                 Vt[is2]=V[is1];
             }
         }
     
   }
   return;
}


extern void acsglowrank2d(int row,int col,int N1,int N2,
			  float *vp,float *tao,float *lam,float *gam,
			  float *k2dh,float dt,float *U,float *V,double *A)
{
     int i,j,n,ii,jj,is;
     float temp;

     int *jrr,*r1;
     int mm=(int)(col*log(N1));

     jrr=(int*)malloc(sizeof(int)*col);
     r1=(int*)malloc(sizeof(int)*mm);

     for(i=0;i<col;i++)
     {
         jrr[i]=0;
     }

     random_select(N1,mm,r1);   //0:N1-1 select mm


     float *Sr,*Sv;
     Sr=(float*)malloc(sizeof(float)*mm*N2);
     Sv=(float*)malloc(sizeof(float)*mm*col);     

     for(i=0;i<mm;i++)
     {
         for(j=0;j<N2;j++)
         {
             is=i*N2+j;
             temp = vp[r1[i]]*dt*pow(k2dh[j],1.0+gam[r1[i]]);
             Sr[is]=sinc(0.5*temp)*pow( k2dh[j], gam[r1[i]] ) * sqrt( 1.0+0.5*dt*tao[r1[i]]*pow(k2dh[j],1.01) );
         }
     }

     int *flag1;
     flag1=(int*)malloc(sizeof(int)*N2);
     for(i=0;i<N2;i++)
     {
         flag1[i]=1;
     }
     
     float mx=0.0,normc=0.0;

     float *tv1,*tv2;
     tv1=(float*)malloc(sizeof(float)*mm);
     tv2=(float*)malloc(sizeof(float)*mm);

     for(n=0;n<col;n++)
     {
         mx=0.0;
         for(j=0;j<N2;j++)
         {
              normc=0.0;
              for(i=0;i<mm;i++)
              {
                  is=i*N2+j;
                  normc=normc+Sr[is]*Sr[is];
              }
              
              if((normc-mx>0.0)&&flag1[j]!=0){mx=normc;jrr[n]=j;}
         }
         
         flag1[jrr[n]]=0;
         
         
         if(n<col-1)
         {
             for(i=0;i<mm;i++)
             {
                 Sv[i*col+n]=Sr[i*N2+jrr[n]];
             }

             for(j=0;j<N2;j++)
             {
                  if(flag1[j]!=0)
                  {
                      for(i=0;i<mm;i++) 
                      {
                         tv1[i]=Sr[i*N2+j];
                      }

                      for(jj=0;jj<=n;jj++)
                      {

                         for(i=0;i<mm;i++)
                         {
                             tv2[i]=Sv[i*col+jj];
                         }
                         temp=dotxy(tv2,tv1,mm)/dotxx(tv2,mm);

                         for(i=0;i<mm;i++)
                         {
                             Sr[i*N2+j]=Sr[i*N2+j]-temp*tv2[i];
                         } 
                      }

                  }
            }
         }// end if

     }// end n iteration

     
     for(i=0;i<N1;i++)
     {
         for(j=0;j<col;j++)
         {
             is=i*col+j;
             
             temp = vp[i]*dt*pow( k2dh[ jrr[j] ],1.0+gam[i] );
             
             U[is] = sinc(0.5*temp) *pow( k2dh[jrr[j]], gam[i] ) * sqrt( 1.0+0.5*dt*tao[i]*pow(k2dh[jrr[j]],1.01) );
         }
     }

     /*FILE *fpu=fopen("Us.dat","wb");
     fwrite(U,sizeof(float),N1*col,fpu);
     fclose(fpu);
     for(n=0;n<col;n++)
     {
          printf("Selected k: n=%d k=%f\n",jrr[n],k2dh[jrr[n]]);
     }*/


     int *irr,*r2;
     int nn=(int)(row*log(N2));

     irr=(int*)malloc(sizeof(int)*row);
     r2=(int*)malloc(sizeof(int)*nn);
     for(i=0;i<row;i++)
     {
         irr[i]=-1;
     }

     random_select(N2,nn,r2);   //0:N2-1 select nn

     

     float *Sc,*Su;
     Sc=(float*)malloc(sizeof(float)*N1*nn);
     Su=(float*)malloc(sizeof(float)*row*nn);    

     for(i=0;i<N1;i++)
     {
         for(j=0;j<nn;j++)
         {
             is=i*nn+j;
             
             temp = vp[i]*dt*pow( k2dh[r2[j]],1.0+gam[i] );
             
             Sc[is] = sinc(0.5*temp) *pow( k2dh[r2[j]], gam[i] ) * sqrt( 1.0+0.5*dt*tao[i]*pow(k2dh[r2[j]],1.01) );
         }
     }

     int *flag2;
     flag2=(int*)malloc(sizeof(int)*N1);
     for(i=0;i<N1;i++)
     {
         flag2[i]=1;
     }
     
     float *tv11,*tv22;
     tv11=(float*)malloc(sizeof(float)*nn);
     tv22=(float*)malloc(sizeof(float)*nn);

     for(n=0;n<row;n++)
     {
         mx=0.0;
         for(i=0;i<N1;i++)
         {
            normc=0.0;
            for(j=0;j<nn;j++)
            {
                is=i*nn+j;
                normc=normc+Sc[is]*Sc[is];
            }

            if((normc-mx>0.0)&&flag2[i]!=0){mx=normc;irr[n]=i;}
         }

         //--------------//
         //if(irr[n]==-1)
         //{
         //   irr[n]=n;
         //}
         //--------------//
     
         flag2[irr[n]]=0;
         

         if(n<row-1)
         {
            for(j=0;j<nn;j++)
            {
                is=n*nn+j;
                Su[is]=Sc[irr[n]*nn+j];
            }

            for(i=0;i<N1;i++)
            {
                if(flag2[i]!=0)
                {
                    for(j=0;j<nn;j++)
                    {
                        tv11[j]=Sc[i*nn+j];
                    }

                    for(ii=0;ii<=n;ii++)
                    {
                        for(j=0;j<nn;j++)
                        {
                            tv22[j]=Su[ii*nn+j];
                        }
                                     
                        temp=dotxy(tv22,tv11,nn)/dotxx(tv22,nn);

                        for(j=0;j<nn;j++)
                        {
                            Sc[i*nn+j]=Sc[i*nn+j]-temp*tv22[j];
                        }

                    }
                }
            }
                       
         }//end if

     }//end n iteration

     

     for(i=0;i<row;i++)
     {
         for(j=0;j<N2;j++)
         {
             	
             	temp = vp[irr[i]]*dt*pow( k2dh[j],1.0+gam[irr[i]] );
             
             V[i*N2+j] = sinc(0.5*temp) *pow( k2dh[j], gam[irr[i]] ) * sqrt( 1.0+0.5*dt*tao[irr[i]]*pow(k2dh[j],1.01) );
         }
     }
     
     
     /*FILE *fpv=fopen("Vs.dat","wb");
     fwrite(V,sizeof(float),row*N2,fpv);
     fclose(fpv);
     printf("\n");
     for(i=0;i<row;i++)
     {
         printf("Selected r: n=%d v=%f Q=%f\n",irr[i],vp[irr[i]],Qp[irr[i]]);
     }*/


     int np=(mm<nn)?mm:nn;
     int *r3,*r4;
     r3=(int*)malloc(sizeof(int)*np);
     r4=(int*)malloc(sizeof(int)*np);

     random_select(N1,np,r3);   //0:N1-1 select np
     random_select(N2,np,r4);   //0:N2-1 select np
     
 
     double *Us,*Vs,*Mes;
     Us=(double*)malloc(sizeof(double)*np*col);
     Vs=(double*)malloc(sizeof(double)*row*np);
     Mes=(double*)malloc(sizeof(double)*np*np);
     
     for(i=0;i<np;i++)
     {
         for(j=0;j<col;j++)
         {
             Us[i*col+j]=U[r3[i]*col+j];
         }
     }
     for(i=0;i<row;i++)
     {
         for(j=0;j<np;j++)
         {
             Vs[i*np+j]=V[i*N2+r4[j]];
         }
     }
     for(i=0;i<np;i++)
     {
         for(j=0;j<np;j++)
         {
              temp=vp[r3[i]]*dt*pow(k2dh[r4[j]],1.0+gam[r3[i]]);
              Mes[i*np+j]=sinc(0.5*temp)*pow( k2dh[r4[j]], gam[r3[i]] ) * sqrt( 1.0+0.5*dt*tao[r3[i]]*pow(k2dh[r4[j]],1.01) ) ;
         }
     }
     
     // pesuinverse of Us (np*col)
     double eps=1.0e-8; 
     double *u1,*v1,*Usv;
     u1=(double*)malloc(sizeof(double)*np*np);
     v1=(double*)malloc(sizeof(double)*col*col);
     Usv=(double*)malloc(sizeof(double)*col*np);
         
     int ka1=(np>col)?(np+1):(col+1);
     i=ginv(Us,np,col,Usv,eps,u1,v1,ka1);
     if(i<=0){printf("Fails to obtian pinv of Us!\n");}

     // pesuinverse of Vs (row*np)
     double *u11,*v11,*Vsv;
     u11=(double*)malloc(sizeof(double)*row*row);
     v11=(double*)malloc(sizeof(double)*np*np);
     Vsv=(double*)malloc(sizeof(double)*np*row);
     int ka2=(row>np)?(row+1):(np+1);
     i=ginv(Vs,row,np,Vsv,eps,u11,v11,ka2);
     if(i<=0){printf("Fails to obtian pinv of Vs!\n");}

     
     double *Wt1;
     Wt1=(double*)malloc(sizeof(double)*col*np);
     MultMatrix(col,np,np,Usv,Mes,Wt1);
     MultMatrix(col,np,row,Wt1,Vsv,A);

     /*FILE *fpa=fopen("Ms.dat","wb");
     fwrite(A,sizeof(double),col*row,fpa);
     fclose(fpa);*/

     free(Sr);free(Sv);free(Su);free(Sc);
     free(flag1);free(flag2);free(tv1);free(tv2);
     free(tv11);free(tv22);free(Wt1);free(Mes);
     free(u1);free(v1);free(u11);free(v11);
     free(Us);free(Vs);free(Usv);free(Vsv);
     
     free(jrr);free(irr);
     free(r1);free(r2);free(r3);free(r4);
     
     return;

}
