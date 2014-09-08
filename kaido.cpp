#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES //M_PI
#include <math.h>
#include "gdal_priv.h"
#include <omp.h>
#include "cpl_conv.h"

//メモリが固定なので、オーバーするとエラーになるかも
//その場合、逐次計算に変更が必要(TODO)

//type 1:地上開度 2:地下開度 3:尾根谷度
double getKaido(int x , int y,double cellsizeX,double R,int nXSize,int nYSize,float *data,int type,float nullValue){
   int iDifX[] = {0,1,1,1,0,-1,-1,-1};
   int iDifY[] = {1,1,0,-1,-1,-1,0,1};
   double iDifL[] = {1,1.41421356,1,1.41421356,1,1.41421356,1,1.41421356};
   double maxAngle[8],minAngle[8];
   double dDifHeight=0; 
   double dDist=0;
   double dKaido=0;

   for(int i=0;i<8;i++){
       int j = 1;
       maxAngle[i] = -100;
       minAngle[i] = 100;
       dDist = iDifL[i] * j * cellsizeX;
       int xx,yy;
       double dAngle;
       while(dDist < R){
           xx = x + iDifX[i]*j;
           yy = y + iDifY[i]*j;
           if( 0<=xx && xx<nXSize && 0<=yy && yy<nYSize && data[nXSize*yy+xx]!=nullValue){
               dDifHeight = data[nXSize*yy+xx] - data[nXSize*y+x];
               dAngle = atan(dDifHeight/dDist);
               if(dAngle > maxAngle[i]) maxAngle[i] = dAngle;
               if(dAngle < minAngle[i]) minAngle[i] = dAngle;
               j++;
               dDist = iDifL[i]*j*cellsizeX;
           }else{
               return -100;
           }    
       }
       if(type==1){
           dKaido = dKaido + M_PI_2 - maxAngle[i];
       }else if(type==2){
           dKaido = dKaido + M_PI_2 + minAngle[i];
       }else if(type==3){
           dKaido = dKaido - (maxAngle[i] + minAngle[i])/2.0;
       }
    }
    return dKaido/8.0;
}

int main(int argc, char ** argv){

    if(argc < 5){
     printf("###############################\n");
     printf("USAGE:\n");
     printf("kaido.exe input.tif output.tif R type\n");
     printf("type 1:chijyou 2:chika 3:onetani\n");
     printf("###############################\n");
     exit(1);
    }
    
    const char  *pszFilename = argv[1];
    const char  *pszKaidoFilename = argv[2];
    double R = atof(argv[3]);
    int type = atoi(argv[4]);
    
    const float radians_to_degrees = 180.0 / M_PI;

    CPLSetConfigOption( "GDAL_CACHEMAX", "1024" );
    CPLSetConfigOption( "GDAL_DATA", "data" );
    GDALDataset *poDataset,*poKaidoDS;
    GDALRasterBand  *poBand,*poKaidoBand;
    GDALDriver *poDriver;
    char **papszOptions = NULL;    
    double adfGeoTransform[6];

    GDALAllRegister(); 
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    poBand = poDataset->GetRasterBand( 1 );
    
    poDataset->GetGeoTransform( adfGeoTransform );

    const double cellsizeY = adfGeoTransform[5];
    const double cellsizeX = adfGeoTransform[1];
    const int   nXSize = poBand->GetXSize();
    const int   nYSize = poBand->GetYSize();
    const float nullValue = poBand->GetNoDataValue();

    poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    papszOptions = CSLSetNameValue( papszOptions, "PROFILE", "GeoTIFF" );
    poKaidoDS = poDriver->Create(pszKaidoFilename,nXSize,nYSize,1,GDT_Float32, papszOptions );
    poKaidoDS->SetGeoTransform( adfGeoTransform );    
    poKaidoDS->SetProjection( poDataset->GetProjectionRef() );
    poKaidoBand = poKaidoDS->GetRasterBand(1);
    poKaidoBand->SetNoDataValue(nullValue);   
      
    float *data  = (float *) VSIMalloc3(nYSize,nXSize,sizeof(float));
    poBand->RasterIO( GF_Read, 0, 0, nXSize,nYSize ,data, nXSize, nYSize, GDT_Float32, 0, 0 );

    double kaido;
    float kaidoBuf[50000];
    int p=0;
    int strp=0;
    int percent[11];
    
    for(int i=0;i<=10;i++){
          percent[i] = nXSize*nYSize*0.1*i;
    }
    
    
    #pragma omp parallel for private(kaido,kaidoBuf)
	for (int y = 0; y < nYSize; y++){
        #pragma omp critical
		{
	       if(p >= percent[strp]){
                printf("%d...",strp*10);
                fflush(stdout);
                strp++;
	       }
           p=p+nXSize;
		}
        for (int x = 0; x < nXSize; x++){
            kaido=getKaido(x,y,cellsizeX,R,nXSize,nYSize,data,type,nullValue);
            if (kaido==-100){
                kaidoBuf[x] = nullValue;
            }else{
                kaidoBuf[x] = kaido*radians_to_degrees;
            }
		}
        poKaidoBand->RasterIO( GF_Write, 0, y, nXSize, 1, kaidoBuf, nXSize, 1, GDT_Float32, 0, 0 ); 
    }
    printf("100\n");
    VSIFree(data);
    delete poKaidoDS;
    delete poDataset;
    return 0;
}
