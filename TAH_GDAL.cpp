#include <stdio.h>
#include <string.h>
#include <math.h>

//GDAL Libraries
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

float FlowConveyanceFromNLCD(int L)
{
	//These are adapted from McCuen (1998).  McCuen, R. H. 1998. Hydrologic analysis and design. New Jersey: Prentice-Hall.
	if(L==11){ return 25.0; } //Open Water
	if(L==21){ return 5.1; } //Developed, Open Space
	if(L==22){ return 10.2; } //Developed, Low Intensity
	if(L==23){ return 15.3; } //Developed, Med Intensity
	if(L==24){ return 20.4; } //Developed, High Intensity
	if(L==31){ return 2.0; } //Barren Land
	if(L==41){ return 1.6; } //Deciduous Forest
	if(L==42){ return 1.6; } //Evergreen Forest
	if(L==43){ return 1.6; } //Mixed Forest
	if(L==52){ return 7.0; } //Shrub/Scrub
	if(L==71){ return 7.0; } //Herbacueous
	if(L==81){ return 3.5; } //Hay/Pasture
	if(L==82){ return 2.0; } //Cultivated Crops
	if(L==90){ return 1.6; } //Woody Wetlands
	if(L==95){ return 1.6; } //Emergent Herbaceuous Wetlands
	return 2.0;
}

float RunoffCoefficientsFromNLCD(int L)
{
	//These are adapted  (Cr) values taken from the American Society of Civil Engineers (ASCE) and Water Pollution Control Federation (WPCF) (1982).
	//  American Society of Civil Engineers (ASCE) and Water Pollution Control Federation (WPCF) (Joint Committee). 1982.
	//  Design and construction of sanitary and storm sewers. ASCE (American Society of Civil Engineers) Manual on
	//  Engineering Practice No. 37 and WPCF (Water Pollution Control Federation) Manual of Practice No. 9., 5th ed.
	if(L==11){ return 1.0; } //Open Water
	if(L==21){ return 0.2; } //Developed, Open Space
	if(L==22){ return 0.3; } //Developed, Low Intensity
	if(L==23){ return 0.4; } //Developed, Med Intensity
	if(L==24){ return 0.5; } //Developed, High Intensity
	if(L==31){ return 0.15; } //Barren Land
	if(L==41){ return 0.15; } //Deciduous Forest
	if(L==42){ return 0.15; } //Evergreen Forest
	if(L==43){ return 0.15; } //Mixed Forest
	if(L==52){ return 0.15; } //Shrub/Scrub
	if(L==71){ return 0.15; } //Herbacueous
	if(L==81){ return 0.15; } //Hay/Pasture
	if(L==82){ return 0.15; } //Cultivated Crops
	if(L==90){ return 0.15; } //Woody Wetlands
	if(L==95){ return 0.15; } //Emergent Herbaceuous Wetlands
	return 0.15;
}

double Cellsize_Conversion( double cellsize, double lat_f, double &cellsize_X, double &cellsize_Y )
{
	//These are based on lat-long conversions from www.zodiacal.com/tools/lat_table.php
	double lat = fabs(lat_f);
	double lat_conv, lon_conv, Conv_Fact;
	double lat_up, lat_down, lon_up, lon_down, lat_base;

	if(cellsize > 0.5)	//This indicates that the DEM is projected, so no need to convert from geographic into projected.
	{
		cellsize_X = cellsize;
		cellsize_Y = cellsize;
		return cellsize;
	}

	if( lat<-1e-16 )
	{
		printf("\nPlease use lat and long values greater than or equal to 0.  Thanks!\n");
		printf("Now that you've learned your lesson, hit -99 and start over.   ");
		scanf("%lf", &Conv_Fact);
		return Conv_Fact;
	}
	if( lat>=0 && lat<=10 ){ lat_up=110.61; lat_down=110.57; lon_up=109.64; lon_down=111.32; lat_base=0.0; }
	else if( lat>10 && lat<=20 ){ lat_up=110.7; lat_down=110.61; lon_up=104.64; lon_down=109.64; lat_base=10.0; }
	else if( lat>20 && lat<=30 ){ lat_up=110.85; lat_down=110.7; lon_up=96.49; lon_down=104.65; lat_base=20.0; }
	else if( lat>30 && lat<=40 ){ lat_up=111.03; lat_down=110.85; lon_up=85.39; lon_down=96.49; lat_base=30.0; }
	else if( lat>40 && lat<=50 ){ lat_up=111.23; lat_down=111.03; lon_up=71.70; lon_down=85.39; lat_base=40.0; }
	else if( lat>50 && lat<=60 ){ lat_up=111.41; lat_down=111.23; lon_up=55.80; lon_down=71.70; lat_base=50.0; }
	else if( lat>60 && lat<=70 ){ lat_up=111.56; lat_down=111.41; lon_up=38.19; lon_down=55.80; lat_base=60.0; }
	else if( lat>70 && lat<=80 ){ lat_up=111.66; lat_down=111.56; lon_up=19.39; lon_down=38.19; lat_base=70.0; }
	else if( lat>80 && lat<=90 ){ lat_up=111.69; lat_down=111.66; lon_up=0.0; lon_down=19.39; lat_base=80.0; }
	else
	{
		printf("\nPlease use legit (0-90) lat and long values.  Thanks!\n");
		printf("Now that you've learned your lesson, hit -99 and start over.   ");
		scanf("%lf", &Conv_Fact);
		return Conv_Fact;
	}

	//Latitude Conversion
	lat_conv=lat_down + (lat_up-lat_down)*(lat-lat_base)/10;
	cellsize_Y = cellsize * lat_conv * 1000.0;		//Converts from degrees to m

	//Longitude Conversion
	lon_conv=lon_down + (lon_up-lon_down)*(lat-lat_base)/10;
	cellsize_X = cellsize * lon_conv * 1000.0;		//Converts from degrees to m

	if( (lat_conv<lat_down || lat_conv>lat_up) || (lon_conv<lon_up || lon_conv>lon_down) )
	{
		printf("\nProblem in conversion!! \n");
		printf("Hit -99 and Reprogram!!   ");
		scanf("%lf", &Conv_Fact);
		return Conv_Fact;
	}
	Conv_Fact = 1000.0*(lat_conv+lon_conv)/2;

	return Conv_Fact*cellsize;
}

float FindSlope(float **DEM, int **FAC, int ro, int co, int nrows, int ncols, double cellsize, double lat_f, double &cellsize_X, double &cellsize_Y)
{
	int r, c, rmin, rmax, cmin, cmax;
	double dist;
	float S=-9999.9;
	int r1, c1, r2, c2;
	int FAC1=-1, FAC2=-2;
	float DEM1=-1, DEM2=-2;

	rmin = ro-1;  if(rmin<0){rmin=0;} if(rmin>=(nrows-1)){rmin=nrows-1;}
	rmax = ro+1;  if(rmax<0){rmax=0;} if(rmax>=(nrows-1)){rmax=nrows-1;}
	cmin = co-1;  if(cmin<0){cmin=0;} if(cmin>=(ncols-1)){cmin=ncols-1;}
	cmax = co+1;  if(cmax<0){cmax=0;} if(cmax>=(ncols-1)){cmax=ncols-1;}

	for(r=rmin; r<=rmax; r++)
	{
		for(c=cmin; c<=cmax; c++)
		{
			if(r==ro && c==co){
				continue; }
			if(FAC[r][c]>FAC1)
			{
				r1=r;
				c1=c;
				FAC1=FAC[r][c];
				DEM1=DEM[r][c];
			}
			if(FAC[r][c]>FAC2 && FAC1!=FAC2)
			{
				r2=r;
				c2=c;
				FAC2=FAC[r][c];
				DEM2=DEM[r][c];
			}
		}
	}

	if(DEM1>0.0 && DEM2>0.0)
	{
		Cellsize_Conversion( cellsize, lat_f, cellsize_X, cellsize_Y );
		dist = sqrt( ((r1-r2)*cellsize_Y)*((r1-r2)*cellsize_Y) + ((c1-c2)*cellsize_X)*((c1-c2)*cellsize_X) );
		S = fabs(DEM2-DEM1)/dist;
	}
	if(S<0.01){
		S=0.03; }

	return S;
}

float Calculate_Distance(int ro1, int ro2, int co1, int co2, double cellsize_X, double cellsize_Y)
{
	float dist = 0.0;
	dist = sqrt( ((ro1-ro2)*cellsize_Y)*((ro1-ro2)*cellsize_Y) + ((co1-co2)*cellsize_X)*((co1-co2)*cellsize_X) );
	return dist;
}

int main(int argc, const char * argv[])
{
	int r, c, ro, co, ro1 , co1 , ro2, co2,  cf, ncols, nrows, max_fac;
	double h, distance, val, xll, yll, cellsize, ND, t=60.0;
	char JUNK[1024];
	bool quit=false;
	int WaterLC = 11;
	int RiverBuffer = 1;
	int RiverFAC_Threshold = 500;
	int numcells=0;
	int Max_FAC=-9999;
	double Slope_Avg=0.0;
	float MaxTT=0.0;
	float sd=0.0;

	float C_calibration = 0.3;  //Calibration Parameter
	C_calibration = 1.0;

	FILE *UH_File;

	// File Names
	char FAC_Name[1024], MSK_Name[1024], DEM_Name[1024], SLP_Name[1024], LC_Name[1024], TT_Name [1024], TPM_Name[1024], UH_Name[1024];
	strcpy (DEM_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\NED_30m\\DuchesneStrawberry_DEM.tif" );
	strcpy (FAC_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\DuchesneStrawberry_FAC_nc.tif" );
	strcpy (MSK_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\Watersheds\\DuchesneStrawberry_Watershed_09279000.tif" );                           //If you don't want to use the Mask, just put the Land Cover here again because it is an integer.
	strcpy (SLP_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\DuchesneStrawberry_FDIR_SLOPE.tif" );
	strcpy (LC_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\DuchesneStrawberry_LandCover.tif" );
	strcpy (TT_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\UnitHydrographs\\TOTALTIME_09279000.tif" );
	strcpy (TPM_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\UnitHydrographs\\TIME_PER_METER_09279000.tif" );
	strcpy (UH_Name, "C:\\Projects\\2022_UCRAF\\GeospatialDatasets\\HYDRO_30m\\UnitHydrographs\\UH_09279000.txt" );

	//FAC at outlet
	int FAC_Gage = -9999;  //If you know the FAC value of the gage, use this.  Otherwise leave as a negative value and program will use max fac value as the outlet (do this especially if you have a mask).
	//FAC_Gage = 687654;
	//FAC_Gage = 74084;


	/**************************************************/
	/*  Read in the first raster and allocate memory  */
	/**************************************************/

	printf("\n\nUSING GDAL TO OPEN RASTERS\n");
	GDALAllRegister();
	//Read In the DEM for Spatial Reference
	GDALDatasetH FACDataset, LCDataset, DEMDataset, SLPDataset, MSKDataset;
	GDALDriverH INDriver, MODDriver, MSKDriver;
	double INGeoTransform[6], FACGeoTransform[6], SLPGeoTransform[6], LCGeoTransform[6], MSKGeoTransform[6];
	char FACProjection[2048], DEMProjection[2048], SLPProjection[2048], LCProjection[2048], MSKProjection[2048];
	double xur, yur;
	int slashmark=0, strlength;

	//These are for the Output Raster
	const char *pszFormat = "GTiff";
    GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
    char **papszMetadata;
	double OutGeoTransform[6];
	papszMetadata = GDALGetMetadata( hDriver, NULL );
	if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) ){ printf( "Driver %s supports Create() method.\n", pszFormat ); }
	if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) ){ printf( "Driver %s supports CreateCopy() method.\n", pszFormat ); }

	//Read the FAC Raster Header data
	FACDataset = GDALOpen( FAC_Name, GA_ReadOnly );
	if( FACDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", FAC_Name ); printf( "\nHit a Key and Try Again\n", FAC_Name );  scanf( "%s", &FAC_Name ); return 1; }
	else{ printf("\nWorking on %s", FAC_Name );}
	GDALRasterBandH FACBand = GDALGetRasterBand( FACDataset, 1 );
	GDALGetGeoTransform( FACDataset, FACGeoTransform );
	strcpy( FACProjection, GDALGetProjectionRef( FACDataset));

	//Read the Land Cover Raster Header data
	LCDataset = GDALOpen( LC_Name, GA_ReadOnly );
	if( LCDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", LC_Name ); printf( "\nHit a Key and Try Again\n", LC_Name );  scanf( "%s", &LC_Name ); return 1; }
	else{ printf("\nWorking on %s", LC_Name );}
	GDALRasterBandH LCBand = GDALGetRasterBand( LCDataset, 1 );
	GDALGetGeoTransform( LCDataset, LCGeoTransform );
	strcpy( LCProjection, GDALGetProjectionRef( LCDataset));

	//Read the Mask Raster Header data
	MSKDataset = GDALOpen( MSK_Name, GA_ReadOnly );
	if( MSKDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", MSK_Name ); printf( "\nHit a Key and Try Again\n", MSK_Name );  scanf( "%s", &MSK_Name ); return 1; }
	else{ printf("\nWorking on %s", MSK_Name );}
	GDALRasterBandH MSKBand = GDALGetRasterBand( MSKDataset, 1 );
	GDALGetGeoTransform( MSKDataset, MSKGeoTransform );
	strcpy( MSKProjection, GDALGetProjectionRef( MSKDataset));

	//Read the SLP Raster Header data
	SLPDataset = GDALOpen( SLP_Name, GA_ReadOnly );
	if( SLPDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", SLP_Name ); printf( "\nHit a Key and Try Again\n", SLP_Name );  scanf( "%s", &SLP_Name ); return 1; }
	else{ printf("\nWorking on %s", SLP_Name );}
	GDALRasterBandH SLPBand = GDALGetRasterBand( SLPDataset, 1 );
	GDALGetGeoTransform( SLPDataset, SLPGeoTransform );
	strcpy( SLPProjection, GDALGetProjectionRef( SLPDataset));

	//Read the DEM Raster Header data
	DEMDataset = GDALOpen( DEM_Name, GA_ReadOnly );
	if( DEMDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", DEM_Name ); printf( "\nHit a Key and Try Again\n", DEM_Name );  scanf( "%s", &DEM_Name ); return 1; }
	else{ printf("\nWorking on %s", DEM_Name );}
	GDALRasterBandH DEMBand = GDALGetRasterBand( DEMDataset, 1 );
	GDALGetGeoTransform( DEMDataset, INGeoTransform );
	strcpy( DEMProjection, GDALGetProjectionRef( DEMDataset));

	ncols = GDALGetRasterBandXSize( FACBand );
	nrows = GDALGetRasterBandYSize( FACBand );
	cellsize = INGeoTransform[1];
	xll = INGeoTransform[0];
	yll = INGeoTransform[3] - nrows * fabs(INGeoTransform[5]);
	xur = xll + (ncols)*INGeoTransform[1];
	yur = INGeoTransform[3];
	printf("\n%lf, %lf, %lf, %lf\n", xll, yll, xur, yur);

	//Allocate menmory and read in the raster data
	printf("\nReading through the Raster Files\n");
	float **DEM = new float*[nrows];  //DEM (meters)
	float **SLP = new float*[nrows];  //Slope (decimal)
	//GInt16 **TT = new GInt16* [nrows]; //Total Time (min)
	//float **TT_Float = new float* [nrows]; //Total Time (min)
	int **TT = new int* [nrows]; //Total Time (min)
	float **Time_Per_Meter = new float*[nrows]; //Time to travel in minutes per meter.
	double **R = new double * [nrows]; //Rate of Tirr.e {rn.in/ m)
	int **D = new int* [nrows]; // Flow Accumulation Value
	int **LC = new int* [nrows]; // Land Cover Raster
	int **MSK = new int* [nrows]; // Mask Raster

	for (r=0; r<nrows; r++)
	{
		DEM[r] = new float[ncols];
		SLP[r] = new float[ncols];
		LC[r] = new int[ncols];
		MSK[r] = new int[ncols];
		D[r] = new int[ncols];
		//TT[r] = new GInt16[ncols];
		//TT_Float[r] = new float[ncols];
		TT[r] = new int[ncols];
		Time_Per_Meter[r] = new float[ncols];
		GDALRasterIO( FACBand, GF_Read, 0, (r), ncols, 1, D[r], ncols, 1, GDT_Int32, 0, 0 );
		GDALRasterIO( DEMBand, GF_Read, 0, (r), ncols, 1, DEM[r], ncols, 1, GDT_Float32, 0, 0 );
		GDALRasterIO( SLPBand, GF_Read, 0, (r), ncols, 1, SLP[r], ncols, 1, GDT_Float32, 0, 0 );
		GDALRasterIO( LCBand, GF_Read, 0, (r), ncols, 1, LC[r], ncols, 1, GDT_Int32, 0, 0 );
		GDALRasterIO( MSKBand, GF_Read, 0, (r), ncols, 1, MSK[r], ncols, 1, GDT_Int32, 0, 0 );

		for (c=0; c<ncols; c++)
		{
			if(D[r][c]>Max_FAC && MSK[r][c]>0){
				Max_FAC=D[r][c]; }
			TT[r][c]=0;
			if(D[r][c]>RiverFAC_Threshold){
				LC[r][c]=WaterLC; }
			/*
			if(DEM[r][c]>0 && D[r][c]>=0 && SLP[r][c]>=0.0 && LC[r][c]>0){
				Slope_Avg = Slope_Avg + SLP[r][c];
				numcells++; }
			*/
			if(DEM[r][c]>0.001 && MSK[r][c]>0)
			{
				numcells++;
				Slope_Avg = Slope_Avg + SLP[r][c];
			}
			else
			{
				D[r][c]=0;
				SLP[r][c]=0.0;
				LC[r][c]=0;
				TT[r][c]=-1;
			}
		}
	}
	Slope_Avg = Slope_Avg / numcells;
	printf("\nAverage Slope is %lf m/m \n", Slope_Avg);


	if(FAC_Gage<1){ FAC_Gage = Max_FAC; }
	printf("\nFAC at the Gage is %d\n", FAC_Gage);



	//Get the rate of time to cross each cell
	float K, S;
	double cellsize_X=0.0, cellsize_Y=0.0, lat_f;
	for (ro=0; ro<nrows; ro++)
	{
		lat_f = yur - ro*cellsize;
		for (co=0; co<ncols; co++)
		{
			if(MSK[ro][co]<=0){
				Time_Per_Meter[ro][co] = 0.0;
				continue; }

			if(D[ro][co]>=0 && DEM[ro][co]>=0.001)
			{
				K = RunoffCoefficientsFromNLCD(LC[ro][co]);
				//S = FindSlope(DEM, D, ro, co, nrows, ncols, cellsize, lat_f, cellsize_X, cellsize_Y);
				S = fabs(SLP[ro][co]);
				if(S<0.01)
				{
					S=Slope_Avg;
					//if(LC[ro][co]==WaterLC){
					//	S=0.01; }
					//else{
					//	S=0.03; }
				}
				Time_Per_Meter[ro][co] = (0.0547 / (C_calibration*K)) * pow(S,-0.5);
			}
			else
			{
				Time_Per_Meter[ro][co] = 0.0;
			}

			//If water is anywhere near the cell, just set the time
			for(r=ro-RiverBuffer; r<=(ro+RiverBuffer); r++)
			{
				for(c=co-RiverBuffer; c<=(co+RiverBuffer); c++)
				{
					if( c<0 || r<0 || c>=(ncols-1) || r>=(nrows-1) ){
						continue;}
					if(LC[r][c]==WaterLC && Time_Per_Meter[ro][co]>1.0)
					{
						Time_Per_Meter[ro][co]=(0.0547 / (C_calibration*25.0)) * pow(fabs(SLP[ro][co]),-0.5);
						if(Time_Per_Meter[ro][co]>1.0){
							Time_Per_Meter[ro][co] = 0.2; }
					}

				}
			}
		}
	}



	//Go througth the Each Cell
	printf ( " \nLooking at each cell and analyzing the Total Time\n") ;
	for ( ro= 0 ; ro<nrows; ro++)
	{
		lat_f = yur - ro*cellsize;
		Cellsize_Conversion( cellsize, lat_f, cellsize_X, cellsize_Y );
		for (co= 0; co<ncols; co ++ )
		{
			if(MSK[ro][co]<=0){
				continue; }
			if(D[ro][co]<0 || DEM[ro][co]<0.001){ TT[ro][co]=-1; continue; }
			quit= false; val=0.0; distance=0.0; ro1=ro; co1=co; 
			while(quit==false)
			{
				cf=0;
				// Look Diagonals First
				max_fac = D[ro1][co1];

				for(r=ro1-1; r<=(ro1+1); r++)
				{
					for(c=co1-1; c<=(co1+1); c++)
					{
						if((r==ro1 && c==co1) ||r<0 || c<0 || c>(ncols-1) || r>(nrows-1) ){
							continue; }

						if(D[r][c]>max_fac)
						{
							max_fac=D[r][c];
							ro2=r;
							co2=c;
							cf++;
							sd=Calculate_Distance(ro1, r, co1, c, cellsize_X, cellsize_Y);
						}
						else if(D[r][c]==max_fac)	//If by chance the two adjacent cells have the max_fac (very unlikely) then go with the closest one.
						{
							if( sd>Calculate_Distance(ro1, r, co1, c, cellsize_X, cellsize_Y) )
							{
								sd = Calculate_Distance(ro1, r, co1, c, cellsize_X, cellsize_Y);
								ro2=r;
								co2=c;
								cf++;
							}
						}
					}
				}


				/*
				r= (ro1)   ; c=(co1)   ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r= (ro1-1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r= (ro1-1) ; c=(co1+1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r= (ro1+1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r= (ro1+1) ; c=(co1+1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
		
				//Now Look at Cardinal Directions (preferred if there is a tie between this and a diaoonal direction
				r= (ro1-1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r =(ro1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r= (ro1) ; c=(co1+1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				r= (ro1+1) ; c=(co1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; cf++; } }
				*/
		
				//Determnine the distance across the cell
				//if ( ro1 == ro2 || co1 == co2 ) { h=1.0; } //This means the adjacent cell is in a cardirnal directiontion
				//else{ h =sqrt (2.0) ; } //This mean the adjacent cell is at a diagonal
				//h = h *cellsize; distance = distance + h ;

				//h = Calculate_Distance(ro1, ro2, co1, co2, cellsize_X, cellsize_Y);
				h = sd; //The short distance (sd) is now calculated in the loop above.
				distance = distance + h;


		
				if ( co1 >=ncols || ro1>=nrows || co1<0 || ro1<0 ) { val=-1; quit=true; }

				//This is the source cell, so this is the end of the line
				else if(D[ro1][co1]==FAC_Gage){ val=1; /*val=val + Time_Per_Meter[ro1][co1]*h;*/ quit=true;}
				else if(cf==0){ val=-1; quit=true; }
				else if(D[ro2][co2]==FAC_Gage){val=val + Time_Per_Meter[ro1][co1]*h; quit=true;}
				else if(D[ro2][co2]>FAC_Gage){val=-1; quit=true;}
				//else if(D[ro2] [co2]==D[ro1][co1]){val=-1; quit=true; }	//I changed this one to -1 because it doesn't drain to the outlet specified.
				else if(D[ro2] [co2]==D[ro1][co1]){val=val + Time_Per_Meter[ro1][co1]*h; quit=true;}

				// This cell has already been analyzed, so just use it's value
				else if(TT[ro2][co2]<-0.9){ val=-1; quit=true; } //This cell ultimately does not drain to the outlet speficied.
				else if(TT[ro2][co2]>0.0){ val = val + 0.5*h*(Time_Per_Meter[ro1][co1] + Time_Per_Meter[ro2][co2]) + TT[ro2][co2]; quit =true ; }
				else{ val = val+ 0.5*h* ( Time_Per_Meter[ro1][co1] + Time_Per_Meter[ro2][co2]); }

				//The adjacent cell is now the cell to evaluate in the next loop
				co1=co2;
				ro1=ro2;
			}
			if(val>=-1 && val<2000000000) //An integer has to be less than 2,147,483,647.  This is important for printing Total Time raster.
			{
				TT[ro][co] = int(val);
				if(TT[ro][co]>MaxTT){ MaxTT = TT[ro][co]; }
			}
		}
	}


	//Write the Output Rasters
	if( hDriver == NULL )
    exit( 1 );
	papszMetadata = GDALGetMetadata( hDriver, NULL );
	if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) ){ printf( "Driver %s supports Create() method.\n", pszFormat ); }
	if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) ){ printf( "Driver %s supports CreateCopy() method.\n", pszFormat ); }
	OutGeoTransform[0] = INGeoTransform[0];	// top-left x (CONTSANT)
	OutGeoTransform[1] = INGeoTransform[1];	// w-e pixel resolution (CONSTANT)
	OutGeoTransform[2] = INGeoTransform[2];	// 0 (CONSTANT)
	OutGeoTransform[3] = INGeoTransform[3];
	OutGeoTransform[4] = INGeoTransform[4];	// 0 (CONSTANT)
	OutGeoTransform[5] = INGeoTransform[5];	// n-s pixel resolution (negative value) (CONSTANT)

	//Write the Time per Meter Raster
	GDALDatasetH TPM_OUT;
	TPM_OUT = GDALCreate(hDriver, TPM_Name, ncols, nrows, 1, GDT_Float32, NULL);
	GDALSetGeoTransform( TPM_OUT, OutGeoTransform );
	GDALSetProjection(TPM_OUT, DEMProjection);
	GDALRasterBandH TPMBand = GDALGetRasterBand( TPM_OUT, 1 );
	for (r=0; r<nrows; r++)
	{
		GDALRasterIO(TPMBand,GF_Write,0,r,ncols, 1, Time_Per_Meter[r], ncols, 1, GDT_Float32, 0, 0);
	}
	GDALClose(TPM_OUT);

	//Write the Total Time Raster
	GDALDatasetH TT_OUT;
	TT_OUT = GDALCreate(hDriver, TT_Name, ncols, nrows, 1, GDT_UInt32, NULL);
	GDALSetGeoTransform( TT_OUT, OutGeoTransform );
	GDALSetProjection(TT_OUT, DEMProjection);
	GDALRasterBandH TTBand = GDALGetRasterBand( TT_OUT, 1 );
	for (r=0; r<nrows; r++)
	{
		GDALRasterIO(TTBand,GF_Write,0,r,ncols, 1, TT[r], ncols, 1, GDT_UInt32, 0, 0);
		//GDALRasterIO(TTBand,GF_Write,0,r,ncols, 1, TT[r], ncols, 1, GDT_Float32, 0, 0);
	}
	GDALClose(DEMDataset);
	GDALClose(FACDataset);
	GDALClose(SLPDataset);
	GDALClose(LCDataset);
	GDALClose(MSKDataset);
	GDALClose(TT_OUT);



	//Write the Unit Hydrograph
	int ordinate, num_ordinates;
	num_ordinates = int(MaxTT/60.0)+1;
	int *UH = new int [num_ordinates]; // Flow Accumulation Value

	for(r=0; r<num_ordinates; r++){UH[r]=0;}

	for(r=0; r<nrows; r++)
	{
		for(c=0; c<ncols; c++)
		{
			if(DEM[r][c]>0.001 && D[r][c]>=0 && SLP[r][c]>=0.0 && LC[r][c]>0 && TT[r][c]>=0.001)
			{
				ordinate = int(TT[r][c]/60);
				UH[ordinate]++;
			}
		}
	}
	
	UH_File = fopen(UH_Name, "w");
	for(r=0; r<num_ordinates; r++)
	{
		fprintf(UH_File, "%d %d\n", int(r*60), UH[r]);
	}
	fclose(UH_File);
	return 0;
}