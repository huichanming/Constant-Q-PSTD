/* generate double-side geometry */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int *argc,char *argv)
{
	int		shotnum;
	int		nr;
	int		*tracenum;
	int		nx;
	int		nz;
	float	offset;
	float	dx;
	float	dz;
	float	*sxp;
	float	*szp;
	float	*rxp;
	float	*rzp;
	float	sx0;
	float	sxdx;
	float	sdep;
	float	rdep;
	float   trdx;

	char strtmp[256];
	FILE *fp=fopen("ac2dgeometry.txt","r");
	if(fp==0){printf("Cannot open ac2dgeometry.txt!\n");exit(0);}
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",&shotnum);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&sx0);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&sxdx);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&sdep);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&rdep);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&offset);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&trdx);
	fscanf(fp,"\n");

	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&dx);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%f",&dz);
        fscanf(fp,"\n");

        fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",&nx);
	fscanf(fp,"\n");
	
	fgets(strtmp,256,fp);
	fscanf(fp,"\n");
	fscanf(fp,"%d",&nz);
	fclose(fp);

	int i,j;
	tracenum=(int*)malloc(sizeof(int)*shotnum);
	sxp=(float*)malloc(sizeof(float)*shotnum);
	szp=(float*)malloc(sizeof(float)*shotnum);

	float widthx=offset;
	float xmin=0.0;
	float xmax=(nx-1)*dx;
	float zmin=0.0;
	float zmax=(nz-1)*dz;
	float xlp,xrp;

	char namerx[50];
	char namerz[50];
	char namesx[50];
	char namesz[50];
	char fname[50];

	FILE *ftd;

        for(i=0;i<shotnum;i++)
	{
		sxp[i]=sx0+i*sxdx;	
		szp[i]=sdep;
		if(sxp[i]<0.0||sxp[i]>xmax||szp[i]<0.0||szp[i]>zmax)
		{
			printf("Error! shot position exceeds model scope!!!\n");exit(0);
		}

		xlp=(sxp[i]-widthx>0.0)?sxp[i]-widthx:0.0;
		xrp=(sxp[i]+widthx<xmax)?sxp[i]+widthx:xmax;

		nr=(int)((xrp-xlp)/trdx)+1;
		rxp=(float*)malloc(sizeof(float)*nr);
		rzp=(float*)malloc(sizeof(float)*nr);
		for(j=0;j<nr;j++)
		{
			rxp[j]=xlp+j*trdx;
			rzp[j]=rdep;
			if(rxp[j]<0.0||rxp[j]>xmax||rzp[j]<0.0||rzp[j]>zmax)
			{
				printf("Error! Receiver position exceeds model scope!!!\n");exit(0);
			}
		}

		strcpy(namerx,"");sprintf(namerx,"gxp%d.bin",i);
		strcpy(namerz,"");sprintf(namerz,"gzp%d.bin",i);
		ftd=fopen(namerx,"wb");fwrite(rxp,sizeof(float),nr,ftd);fclose(ftd);
                ftd=fopen(namerz,"wb");fwrite(rzp,sizeof(float),nr,ftd);fclose(ftd);	
		tracenum[i]=nr;

		printf("id=%d nr=%d xlp=%f xrp=%f xdis=%f\n",i,nr,xlp,xrp,xrp-xlp);

		free(rxp);free(rzp);
	}

	strcpy(namesx,"");sprintf(namesx,"shotxp.bin");
	strcpy(namesz,"");sprintf(namesz,"shotzp.bin");
	ftd=fopen(namesx,"wb");fwrite(sxp,sizeof(float),shotnum,ftd);fclose(ftd);
	ftd=fopen(namesz,"wb");fwrite(szp,sizeof(float),shotnum,ftd);fclose(ftd);
	strcpy(fname,"");sprintf(fname,"tracenum.bin");
	ftd=fopen(fname,"wb");fwrite(tracenum,sizeof(int),shotnum,ftd);fclose(ftd);
	
	free(tracenum);free(sxp);free(szp);
	
	return 0;
}
