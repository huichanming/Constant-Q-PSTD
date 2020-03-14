#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "cufft.h"
#include "cuda.h"
#define pi 3.1415926

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "revolve.h"

//extern "C"
struct MultiGPU
{
	// host memory //
	int *shotid;

	float *sxp;
	float *szp;
	
	float *h_vp; 
	float *h_rho;
	float *h_bbx;
	float *h_bbz;
	
	float *h_yinda;
	float *h_tao;
	float *h_lam;
	float *h_gam;
	
	float *h_dvp;
	
	float *h_Qp;
	float *h_vpq;
	float *h_vwd;
	float *h_obs;
	
	int *rindex;
	int *flagrec;
	
	int lgest;
	int rgest;
	float mlx;
	float taomn;
	float taomx;
	float taoc;

	int nr;   // tracenumber
	int xsg;  // source x index
	int zsg;  // source z index

	int *gxg; // receiver x index
	int *gzg; // receiver z index


        // host_device memory //
	float *h_p;
	
	float *h_img;
	float *h_ssg;
	
	float *h_img_pj;
	float *h_ssg_pj;


	// device memory //	
	float *d_vp;
	float *d_rho;
	float *d_bbx;
	float *d_bbz;
	
	float *d_Qp;
	float *d_vpq;
	float *d_vwd;
	
	float *d_dvp;
	
	float *d_yinda;
	float *d_tao;
	float *d_lam;
	float *d_gam;


	float *d_p;
	float *d_ps;
	float *d_pf;
	
	
	float *d_px,*d_pxs,*d_pxf;
	float *d_pz,*d_pzs,*d_pzf;
	float *d_fixx,*d_fizz;
	
	float *d_vx;
	float *d_vz;
	float *d_fipx;
	float *d_fipz;
	float *d_fivxx;
	float *d_fivzz;
	
	
	// for born modeling //
	float *ds_p;
	float *ds_ps;
	float *ds_pf;
	
	float *ds_px,*ds_pxs,*ds_pxf;
	float *ds_pz,*ds_pzs,*ds_pzf;
	float *ds_fixx,*ds_fizz;
	
	
	float *d_rik;
	float *d_k2d;
	float *d_k2d2;
	float *d_k2d3;
	float *d_k2d1;
	
	float *d_k2da;
	
	float *d_kxsmr,*d_kxscr,*d_kzsmr,*d_kzscr;
	float *d_kxsmi,*d_kxsci,*d_kzsmi,*d_kzsci;
	
	// for lowrank decomposition //
	float *host_U;
	float *host_V;
	double *host_A;
	float *host_Vt;
	
	float *host_We;
	float *host_Wa;
	
	float *dev_U;
	float *dev_V;
	double *dev_A;
	
	
	cufftComplex *d_psk,*d_pskf;
        cufftComplex *d_pifft;

        float *d_psc;
        double *d_pst;
        
        float *d_psc_x;
        float *d_psc_z;
	

	float *d_atx, *d_btx, *d_ctx;
	float *d_atz, *d_btz, *d_ctz;
	float *d_ax,  *d_bx,  *d_az, *d_bz;
	
	float *d_i_ax, *d_i_bx;
	float *d_h_ax, *d_h_bx;
	float *d_i_az, *d_i_bz;
	float *d_h_az, *d_h_bz;
	
	float *d_aux;
	float *d_aua;
	
	float *d_img;
	float *d_imgf;
	float *d_ssg;
	float *d_obs;
	
	int *d_flagrec;
	int *d_rindex;

	cufftComplex *d_dis0;
	cufftComplex *d_dis1;
	cufftComplex *d_dis2;
	cufftComplex *d_dis3;
	cufftComplex *d_dis4;

	float *d_total;
   

	// cufft related handles //
        cufftHandle plan2dforward;
	cufftHandle plan2dinverse;
	
        cufftHandle plan2dforward2;
        cufftHandle plan2dinverse2;
        cufftHandle plan2dinverse3;
        
        // receiver wavefield //
        float *d_vxb, *d_vzb, *d_pb;
        float *d_fipxb, *d_fipzb, *d_fivxxb, *d_fivzzb;
        float *d_auxb;
        
        float *store_p,*store_vx,*store_vz;
        float *store_fipx,*store_fipz;
        float *store_fivxx,*store_fivzz;
        float *store_aux;
	
};

struct Multistream
{
	cudaStream_t stream;
};


