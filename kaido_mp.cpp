#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gdal_priv.h"
#include <omp.h>
#include "cpl_conv.h"

int main(int argc, char ** argv) 
{
    if(argc < 6){
     printf("###############################\n");
     printf("USAGE:\n");
     printf("kaido_mp.exe input.tif output.tif R type div\n");
     printf("type 1:chijyou 2:chika 3:onetani\n");
     printf("###############################\n");
     exit(1);
    }
    
    CPLSetConfigOption( "GDAL_CACHEMAX", "1024" );
    CPLSetConfigOption( "GDAL_DATA", "data" );

    GDALDataset *poDataset;     
    const float degrees_to_radians = 3.14159 / 180.0;
    const float radians_to_degrees = 180.0 / 3.14159;           
    double      adfGeoTransform[6];

    const char  *pszFilename = argv[1];
    const char  *pszKaidoFilename = argv[2];
    double R = atof(argv[3]);
    int type = atoi(argv[4]);
    int divmax = atoi(argv[5]);
    GDALAllRegister(); 

    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    

    GDALRasterBand  *poBand;       
    poBand = poDataset->GetRasterBand( 1 );
    poDataset->GetGeoTransform( adfGeoTransform );

    const double cellsizeY = adfGeoTransform[5];
    const double cellsizeX = adfGeoTransform[1];
    const float nullValue = -9999;
    const int   nXSize = poBand->GetXSize();
    const int   nYSize = poBand->GetYSize();
    const int celln = ceil(R*2/cellsizeX);
    const int pcell = (celln + 1)/2; 
    
    printf("sizeX=%f,sizeY=%f\n",cellsizeX,cellsizeY);
    printf("winsize=%d,center=%d\n",celln,pcell);
   

    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset      *poKaidoDS;    
    GDALRasterBand   *poKaidoBand; 
    char **papszOptions = NULL;
    
    papszOptions = CSLSetNameValue( papszOptions, "PROFILE", "GeoTIFF" );

    poKaidoDS = poDriver->Create(pszKaidoFilename,nXSize,nYSize,1,GDT_Float32, papszOptions );
    poKaidoDS->SetGeoTransform( adfGeoTransform );    
    poKaidoDS->SetProjection( poDataset->GetProjectionRef() );
    poKaidoBand = poKaidoDS->GetRasterBand(1);
    poKaidoBand->SetNoDataValue(-9999);   
  
    int percent[11];
    for(int i=0;i<=10;i++){
          percent[i] = nXSize*nYSize*0.1*i;
    }
    int start[100],end[100];
    for(int i=0;i<divmax;i++){
     start[i] = i*nYSize/divmax-celln;
     end[i] = (i+1)*nYSize/divmax-celln;
     if(i==0) start[i]=0;
//     printf("%d,%d,%d,%d\n",start[i],end[i],end[i]-start[i],nYSize);
    } 
    
 
    float win[10000];
    float kaidoBuf[50000];
       
    int p=0;
    int strp=0;
   
    for(int div=0;div < divmax;div++){
   
     float *bigwin  = (float *) VSIMalloc3((end[div]-start[div]+celln),nXSize,sizeof(float));
     if( bigwin == NULL ){
        printf("### Please gain div count. ###\n");
        exit(1); 
     }
     
     poBand->RasterIO( GF_Read, 0, start[div], nXSize,(end[div]-start[div]+celln) ,bigwin, nXSize, (end[div]-start[div]+celln), GDT_Float32, 0, 0 );
    
    #pragma omp parallel for firstprivate(win,kaidoBuf)

    for (int i = start[div]; i < end[div]; i++) 
    {
      //printf("%d\n",i);
        for (int j = 0; j < nXSize; j++) 
        {
            #pragma omp critical
            {
	     if(p >= percent[strp]){
              printf("%d%s",strp*10,"...");
              strp++;
	     }
             p++;
	    }
            if (i > nYSize-celln || j > nXSize-celln ) 
            {

                kaidoBuf[j] = nullValue;
                continue;
            }
            ////////////////////winに代入////////////////////////////////////
            for(int ii=0;ii<celln;ii++){
             for(int jj=0;jj<celln;jj++){
              win[ii*celln + jj] = bigwin[((i-start[div]+ii)*nXSize + j) + jj];
             }
            } 
            ////////////////////////////////////////////////////////

    
              int maxelev;
              int maxpos;
              int minelev;
              int minpos;
	      int start;
	      int last;
              int m;
              int center = celln*(pcell-1)+pcell-1;
              //left
              maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = pcell-1;
              last = 0;
	      for(int k=start;k>last;k--){
	       if(win[celln*(pcell-1)+k-1]>maxelev){
		 maxelev = win[celln*(pcell-1)+k-1];
                 maxpos = k;
	       }
               if(win[celln*(pcell-1)+k-1]<minelev){
		 minelev = win[celln*(pcell-1)+k-1];
                 minpos = k;
	       }
              }
              float k_left_max = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX))*radians_to_degrees;
	      float k_left_min = atan((minelev - win[center])/((pcell-minpos)*cellsizeX))*radians_to_degrees;
              float k_left = k_left_max > 0 ? k_left_max : k_left_min;
              //right
	      maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;

	      start = pcell+1;
              last = celln;
	      for(int k=start;k<last;k++){
	       if(win[celln*(pcell-1)+k-1]>maxelev){
		 maxelev = win[celln*(pcell-1)+k-1];
                 maxpos = k;
	       }
               if(win[celln*(pcell-1)+k-1]<minelev){
		 minelev = win[celln*(pcell-1)+k-1];
                 minpos = k;
	       }
              }
              float k_right_max = atan((maxelev - win[center])/((maxpos-pcell)*cellsizeX))*radians_to_degrees;
	      float k_right_min = atan((minelev - win[center])/((minpos-pcell)*cellsizeX))*radians_to_degrees;
	      float k_right = k_right_max > 0 ? k_right_max : k_right_min; 
             //top
             maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = 1;
              last = pcell;
	      for(int k=start;k<last;k++){
	       if(win[celln*(k-1)+pcell-1]>maxelev){
		 maxelev = win[celln*(k-1)+pcell-1];
                 maxpos = k;
	       }
               if(win[celln*(k-1)+pcell-1]<minelev){
		 minelev = win[celln*(k-1)+pcell-1];
                 minpos = k;
	       }
              }
              float k_top_max = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX))*radians_to_degrees;
              float k_top_min = atan((minelev - win[center])/((pcell-minpos)*cellsizeX))*radians_to_degrees;
              float k_top = k_top_max > 0 ? k_top_max : k_top_min; 
	      //bottom
              maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = pcell+1;
              last = celln+1;
	      for(int k=start;k<last;k++){
	       if(win[celln*(k-1)+pcell-1]>maxelev){
		 maxelev = win[celln*(k-1)+pcell-1];
                 maxpos = k;
	       }
               if(win[celln*(k-1)+pcell-1]<minelev){
		 minelev = win[celln*(k-1)+pcell-1];
                 minpos = k;
	       }
              }
              float k_bottom_max = atan((maxelev - win[center])/((maxpos-pcell)*cellsizeX))*radians_to_degrees;
              float k_bottom_min = atan((minelev - win[center])/((minpos-pcell)*cellsizeX))*radians_to_degrees;
              float k_bottom = k_bottom_max > 0 ? k_bottom_max : k_bottom_min;
	      //topleft
              maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = pcell-1;
              last = 0;
	      for(int k=start;k>last;k--){
	       if(win[celln*(k-1)+k-1]>maxelev){
		 maxelev = win[celln*(k-1)+k-1];
                 maxpos = k;
	       }
               if(win[celln*(k-1)+k-1]<minelev){
		 minelev = win[celln*(k-1)+k-1];
                 minpos = k;
	       }
               if((pcell-k) * cellsizeX * 1.4142 >= R) break;
              }
              float k_topleft_max = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX)*1.4142)*radians_to_degrees;
              float k_topleft_min = atan((minelev - win[center])/((pcell-minpos)*cellsizeX)*1.4142)*radians_to_degrees;
	      float k_topleft = k_topleft_max > 0 ? k_topleft_max : k_topleft_min;
              //bottomright
              maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = pcell+1;
              last = celln;
	      for(int k=start;k<last;k++){
	       if(win[celln*(k-1)+k-1]>maxelev){
		 maxelev = win[celln*(k-1)+k-1];
                 maxpos = k;
	       }
               if(win[celln*(k-1)+k-1]<minelev){
		 minelev = win[celln*(k-1)+k-1];
                 minpos = k;
	       }
               if ((k-pcell)*cellsizeX*1.4142 >= R) break;
              }
              float k_bottomright_max = atan((maxelev - win[center])/((maxpos-pcell)*cellsizeX)*1.4142)*radians_to_degrees;
              float k_bottomright_min = atan((minelev - win[center])/((minpos-pcell)*cellsizeX)*1.4142)*radians_to_degrees;
              float k_bottomright = k_bottomright_max > 0 ? k_bottomright_max : k_bottomright_min;
	      //topright
              maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = pcell-1;
              last = 0;
	      m = pcell + 1;
	      for(int k=start;k>last;k--){
	       if(win[celln*(k-1)+m-1]>maxelev){
		 maxelev = win[celln*(k-1)+m-1];
                 maxpos = k;
	       }
               if(win[celln*(k-1)+m-1]<minelev){
		 minelev = win[celln*(k-1)+m-1];
                 minpos = k;
	       }
               if ((pcell-k)*cellsizeX*1.4142 >= R) break;
               m++;
              }
              float k_topright_max = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX)*1.4142)*radians_to_degrees;
              float k_topright_min = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX)*1.4142)*radians_to_degrees;
	      float k_topright = k_topright_max > 0 ? k_topright_max : k_topright_min;
              //bottomleft
              maxelev=-9999;
              maxpos = 0;
              minelev=9999;
              minpos = 0;
	      start = pcell-1;
              last = 0;
	      m = pcell + 1;
	      for(int k=start;k>last;k--){
	       if(win[celln*(m-1)+k-1]>maxelev){
		 maxelev = win[celln*(m-1)+k-1];
                 maxpos = k;
	       }
               if(win[celln*(m-1)+k-1]<minelev){
		 minelev = win[celln*(m-1)+k-1];
                 minpos = k;
	       }
               if( (pcell-k)*cellsizeX*1.4142 >= R) break;
               m++;
              }
              float k_bottomleft_max = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX)*1.4142)*radians_to_degrees;
              float k_bottomleft_min = atan((maxelev - win[center])/((pcell-maxpos)*cellsizeX)*1.4142)*radians_to_degrees;
              float k_bottomleft = k_bottomleft_max > 0 ? k_bottomleft_max : k_bottomleft_min;
              //90-maxが地上開度,90+minが地下開度
	      if(type==1){
               ////地上開度
                kaidoBuf[j] = 90-(k_left_max + k_right_max + k_top_max + k_bottom_max + k_topleft_max + k_topright_max + k_bottomright_max + k_bottomleft_max)/8.0;
              }else if(type==2){ 
              ////地下開度
                kaidoBuf[j] = 90+(k_left_min + k_right_min + k_top_min + k_bottom_min + k_topleft_min + k_topright_min + k_bottomright_min + k_bottomleft_min)/8.0;
              }else if(type==3){
	      ////尾根谷度（(地上開度-地下階度)/2）
                kaidoBuf[j] = (-1.0*(k_left_min + k_right_min + k_top_min + k_bottom_min + k_topleft_min + k_topright_min + k_bottomright_min + k_bottomleft_min)-(k_left_max + k_right_max + k_top_max + k_bottom_max + k_topleft_max + k_topright_max + k_bottomright_max + k_bottomleft_max))/16.0;
              }
        }
         poKaidoBand->RasterIO( GF_Write, pcell-1, pcell-1+i, nXSize-(pcell-1), 1, kaidoBuf, nXSize-(pcell-1), 1, GDT_Float32, 0, 0 ); 
      }
          VSIFree(bigwin);
    }
    

    //delete poBand;
    delete poKaidoDS;
    delete poDataset;
    return 0;
}
