#include <stdio.h>
#include <string.h>
#include <math.h>

//GDAL Libraries
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

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

float CurveNumbersFromNLCD(int L, char LMH[])
{
	//https://www.wheatfield.ny.us/DocumentCenter/View/1434/T1?bidId=
	//  1) Land Cover Classifications - National Land Cover Database 2011 (http://www.mrlc.gov/nlcd11_leg.php)
	//  2) Curve Numbers - Stormwater Conveyance Modeling and Design, Haestad and Durrans, Table 5.5 to 5.7

	if(strcmp(LMH,"Low")==0)
	{
		if(L==11){ return 100.0; } //Open Water
		if(L==21){ return 74.0; } //Developed, Open Space
		if(L==22){ return 77.0; } //Developed, Low Intensity
		if(L==23){ return 82.0; } //Developed, Med Intensity
		if(L==24){ return 88.0; } //Developed, High Intensity
		if(L==31){ return 88.0; } //Barren Land
		if(L==41){ return 73.0; } //Deciduous Forest
		if(L==42){ return 77.0; } //Evergreen Forest
		if(L==43){ return 75.0; } //Mixed Forest
		if(L==52){ return 77.0; } //Shrub/Scrub
		if(L==71){ return 74.0; } //Herbacueous
		if(L==81){ return 71.0; } //Hay/Pasture
		if(L==82){ return 82.0; } //Cultivated Crops
		if(L==90){ return 82.0; } //Woody Wetlands
		if(L==95){ return 77.0; } //Emergent Herbaceuous Wetlands
		return 78.357;
	}
	else if(strcmp(LMH,"High")==0)
	{
		if(L==11){ return 100.0; } //Open Water
		if(L==21){ return 79.0; } //Developed, Open Space
		if(L==22){ return 82.0; } //Developed, Low Intensity
		if(L==23){ return 86.0; } //Developed, Med Intensity
		if(L==24){ return 93.0; } //Developed, High Intensity
		if(L==31){ return 88.0; } //Barren Land
		if(L==41){ return 79.0; } //Deciduous Forest
		if(L==42){ return 85.0; } //Evergreen Forest
		if(L==43){ return 82.0; } //Mixed Forest
		if(L==52){ return 83.0; } //Shrub/Scrub
		if(L==71){ return 80.0; } //Herbacueous
		if(L==81){ return 78.0; } //Hay/Pasture
		if(L==82){ return 86.0; } //Cultivated Crops
		if(L==90){ return 86.0; } //Woody Wetlands
		if(L==95){ return 83.0; } //Emergent Herbaceuous Wetlands
		return 83.571;
	}
	else  //Should correlate to Average Values
	{
		if(L==11){ return 100.0; } //Open Water
		if(L==21){ return 75.75; } //Developed, Open Space
		if(L==22){ return 79.0; } //Developed, Low Intensity
		if(L==23){ return 83.5; } //Developed, Med Intensity
		if(L==24){ return 90.5; } //Developed, High Intensity
		if(L==31){ return 88.0; } //Barren Land
		if(L==41){ return 75.25; } //Deciduous Forest
		if(L==42){ return 80.0; } //Evergreen Forest
		if(L==43){ return 77.75; } //Mixed Forest
		if(L==52){ return 79.25; } //Shrub/Scrub
		if(L==71){ return 76.25; } //Herbacueous
		if(L==81){ return 73.5; } //Hay/Pasture
		if(L==82){ return 83.5; } //Cultivated Crops
		if(L==90){ return 83.5; } //Woody Wetlands
		if(L==95){ return 79.25; } //Emergent Herbaceuous Wetlands
		return 80.357;
	}
	return 80.357;
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

float Calculate_Distance(int ro1, int ro2, int co1, int co2, double cellsize_X, double cellsize_Y)
{
	float dist = 0.0;
	dist = sqrt( ((ro1-ro2)*cellsize_Y)*((ro1-ro2)*cellsize_Y) + ((co1-co2)*cellsize_X)*((co1-co2)*cellsize_X) );
	return dist;
}

int main(int argc, const char * argv[])
{
	int r, c, ncols, nrows;
	double xll, yll, cellsize, t=60.0;
	char JUNK[1024];
	int numcells=0;

	int TimeStep_min = 60*24;
	int m, n, num_rasters=0;
	float IncrementalFlow;
	double cellsize_X=0.0, cellsize_Y=0.0, lat_f;

	float RunoffCalibration = 1.0;

	//Curve Number Information    https://engineering.purdue.edu/mapserve/LTHIA7/documentation/scs.htm
	float Ia, P, S;

	int b, baseflow_mult;

	FILE *InRasterList_File, *Outflow_File, *UH_File;

	// File Names
	char TT_Name[1024], LC_Name[1024], UH_Name[1024], Rast_Name[1024], Outflow_Name[1024], InRasterList_File_Name[1024];
	strcpy (TT_Name, "C:\\Projects\\2022_GCUH\\BullLakeTestSite\\RasterData\\TOTALTIME.tif" );
	strcpy (LC_Name, "C:\\Projects\\2022_GCUH\\BullLakeTestSite\\RasterData\\BullLake_LANDCOVER.tif" );
	strcpy (InRasterList_File_Name, "C:\\Projects\\2022_GCUH\\BullLakeTestSite\\RunoffRasterList.txt" );
	strcpy (UH_Name, "C:\\Projects\\2022_GCUH\\BullLakeTestSite\\UH.txt" );
	strcpy (Outflow_Name, "C:\\Projects\\2022_GCUH\\BullLakeTestSite\\Runoff.txt" );


	//Open the Input Rast List to see how many rasters to process
	if( (InRasterList_File=fopen( InRasterList_File_Name, "r" )) != NULL ) { printf( "\n Opened %s with no problems \n", InRasterList_File_Name ); }
	else{ printf( "\n PROBLEMS OPENING %s !!!!! \n", InRasterList_File_Name ); printf( "\nHit a Key and Try Again\n", InRasterList_File_Name ); scanf( "%s", &InRasterList_File_Name ); return 1; }

	while((fscanf(InRasterList_File,"%s", &Rast_Name)) !=EOF )
	{
		num_rasters++;
	}
	printf("\nWe have %d runoff rasters to evaluate\n", num_rasters);
	rewind(InRasterList_File);

	//Open the Unit Hydrograph
	int i, num_uh_ordinates=0, val;
	float uh_sum=0.0;
	if( (UH_File=fopen( UH_Name, "r" )) != NULL ) { printf( "\n Opened %s with no problems \n", UH_Name ); }
	else{ printf( "\n PROBLEMS OPENING %s !!!!! \n", UH_Name ); printf( "\nHit a Key and Try Again\n", UH_Name ); scanf( "%s", &UH_Name ); return 1; }
	while((fscanf(UH_File,"%d %d", &val, &val)) !=EOF )
	{
		num_uh_ordinates++;
	}
	printf("\nWe have %d ordinates in the unit hydrograph rasters to evaluate\n", num_uh_ordinates);
	rewind(UH_File);
	int *uh_ts = new int[num_uh_ordinates];
	int *uh_val = new int[num_uh_ordinates];
	for(i=0; i<num_uh_ordinates; i++)
	{
		fscanf(UH_File,"%d %d", &uh_ts[i], &uh_val[i]);
		uh_sum += uh_val[i];
	}
	fclose(UH_File);


	//Open and Read the Raster Information
	printf("\n\nUSING GDAL TO OPEN RASTERS\n");
	GDALAllRegister();
	//Read In the DEM for Spatial Reference
	GDALDatasetH RUNOFFDataset, TTDataset, LCDataset;
	GDALDriverH INDriver;
	double INGeoTransform[6], RUNOFFGeoTransform[6];
	char RUNOFFProjection[2048], TTProjection[2048], LCProjection[2048];
	double xur, yur;
	int slashmark=0, strlength;

	//Read the Total Time (TT) Raster Header data
	TTDataset = GDALOpen( TT_Name, GA_ReadOnly );
	if( TTDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", TT_Name ); printf( "\nHit a Key and Try Again\n", TT_Name );  scanf( "%s", &TT_Name ); return 1; }
	else{ printf("\nWorking on %s", TT_Name );}
	GDALRasterBandH TTBand = GDALGetRasterBand( TTDataset, 1 );
	GDALGetGeoTransform( TTDataset, INGeoTransform );
	strcpy( TTProjection, GDALGetProjectionRef( TTDataset));

	ncols = GDALGetRasterBandXSize( TTBand );
	nrows = GDALGetRasterBandYSize( TTBand );
	cellsize = INGeoTransform[1];
	xll = INGeoTransform[0];
	yll = INGeoTransform[3] - nrows * fabs(INGeoTransform[5]);
	xur = xll + (ncols)*INGeoTransform[1];
	yur = INGeoTransform[3];
	printf("\n%lf, %lf, %lf, %lf\n", xll, yll, xur, yur);


	//Read the Land Cover (LC) Raster Header data
	LCDataset = GDALOpen( LC_Name, GA_ReadOnly );
	if( LCDataset == NULL ){ printf( "\n PROBLEMS OPENING %s !!!!! \n", LC_Name ); printf( "\nHit a Key and Try Again\n", LC_Name );  scanf( "%s", &LC_Name ); return 1; }
	else{ printf("\nWorking on %s", LC_Name );}
	GDALRasterBandH LCBand = GDALGetRasterBand( LCDataset, 1 );
	GDALGetGeoTransform( LCDataset, INGeoTransform );
	strcpy( LCProjection, GDALGetProjectionRef( LCDataset));


	/**************************************************/
	/*  Read in the first raster and allocate memory  */
	/**************************************************/

	//Allocate menmory and read in the raster data
	printf("\nReading through the Raster Files\n");
	int **TT = new int* [nrows]; //Total Time (min)
	int **LC = new int* [nrows]; //LandCover
	int **Runoff = new int * [nrows]; //Runoff data
	float **RunoffCoeff = new float * [nrows]; //Runoff Coefficient
	float **CurveNumber = new float * [nrows]; //Curve Number

	for (r=0; r<nrows; r++)
	{
		TT[r] = new int[ncols];
		LC[r] = new int[ncols];
		Runoff[r] = new int[ncols];
		RunoffCoeff[r] = new float[ncols];
		CurveNumber[r] = new float[ncols];
		GDALRasterIO( TTBand, GF_Read, 0, (r), ncols, 1, TT[r], ncols, 1, GDT_Int32, 0, 0 );
		GDALRasterIO( LCBand, GF_Read, 0, (r), ncols, 1, LC[r], ncols, 1, GDT_Int32, 0, 0 );
		for (c=0; c<ncols; c++)
		{
			Runoff[r][c]=0;
			RunoffCoeff[r][c] = RunoffCoefficientsFromNLCD(LC[r][c]) * RunoffCalibration;
			CurveNumber[r][c] = CurveNumbersFromNLCD(LC[r][c], "Med");
		}
	}
	GDALClose(TTDataset);
	GDALClose(LCDataset);

	//Allocate Memory for the outflow
	int ordinate, num_ordinates;
	num_ordinates = int(num_rasters*60*24 / TimeStep_min) + int(7*24*60 / TimeStep_min);  //The number of ordinates is based on the number of rasters.  Also go 7 days past the last raster input.
	float *Outflow = new float[num_ordinates];
	float *Outflow_from_UH = new float[num_ordinates];
	float *TotalRunoffPerDay = new float[num_ordinates];
	for(ordinate=0; ordinate<num_ordinates; ordinate++){
		Outflow[ordinate]=0.0; Outflow_from_UH[ordinate]=0.0; TotalRunoffPerDay[ordinate]=0.0; }


	//Get the Area by Latitude
	float *Area= new float[nrows];
	for(r=0; r<nrows; r++)
	{
		lat_f = yur - r*cellsize;
		Cellsize_Conversion( cellsize, lat_f, cellsize_X, cellsize_Y );
		Area[r] = cellsize_X * cellsize_Y;
	}


	//Now lets go through the raster list (assuming the raster list is daily
	for(n=0; n<num_rasters; n++)
	{
		fscanf(InRasterList_File,"%s",&Rast_Name);

		//Open the Runoff raster
		RUNOFFDataset = GDALOpen( Rast_Name, GA_ReadOnly );
		if( RUNOFFDataset == NULL ){
			printf( "\n PROBLEMS OPENING %s !!!!! \n", Rast_Name );
			continue;
		}
		else{ printf("\nWorking on %s", Rast_Name );}

		GDALRasterBandH RUNOFFBand = GDALGetRasterBand( RUNOFFDataset, 1 );
		GDALGetGeoTransform( RUNOFFDataset, RUNOFFGeoTransform );
		strcpy( RUNOFFProjection, GDALGetProjectionRef( RUNOFFDataset));

		for(r=0; r<nrows; r++)
		{
			GDALRasterIO( RUNOFFBand, GF_Read, 0, (r), ncols, 1, Runoff[r], ncols, 1, GDT_Int32, 0, 0 );
			for(c=0; c<ncols; c++)
			{
				if(Runoff[r][c]>0 && TT[r][c]>=0 && TT[r][c]<2000000)
				{
					//Rational Method
					//IncrementalFlow = float(Runoff[r][c]*Area[r]/100000) * RunoffCoeff[r][c];  //Have to divide by 100000 to get to meters (https://nsidc.org/data/g02158).  By mulitplying by Area you get cubic meters.

					//Curve Number Method (remember this method uses inches)
					/*
					S = (1000.0/CurveNumber[r][c])-10.0;
					Ia = 0.2*S;
					P = float(Runoff[r][c]*39.3701/100000.0);  //Have to divide by 100000 to get to meters (https://nsidc.org/data/g02158).  Multiply by 39.3701 to get to inches.
					IncrementalFlow = (P-Ia)*(P-Ia)/((P-Ia)+S);	//Runoff in inches.
					IncrementalFlow = (IncrementalFlow*Area[r]/39.3701);  //Convert back to meters and also multiply by Area to get cubic meters;
					*/


					//Straight Runoff (no infiltration)
					IncrementalFlow = float(Runoff[r][c]*Area[r]/100000);  //Have to divide by 100000 to get to meters (https://nsidc.org/data/g02158).  By mulitplying by Area you get cubic meters.


					TotalRunoffPerDay[n]+=IncrementalFlow;  //This is summing up the runoff.  Should be m3.
					IncrementalFlow = IncrementalFlow / (24*60);  //Incremental Flow is now cubic meters per minute

					for(m=0; m<(60*24/TimeStep_min); m++)   //Have to do this if your runoff data is daily and you want hourly flow rates
					{
						ordinate = int(n*60*24 / TimeStep_min) + int(TT[r][c] / TimeStep_min) + m;
						if(ordinate>=num_ordinates){ printf("\nError with Ordinate (%d) being higher than limit (%d)", ordinate, num_ordinates); continue;}
						IncrementalFlow = IncrementalFlow*TimeStep_min;  //Now we have volume of flow for that timestep.
						Outflow[ordinate]+=IncrementalFlow * RunoffCoeff[r][c];

						//Let's make-up baseflow and see what happens.  Try adding in the same runoff, but put it multiple times in the future
						baseflow_mult=4;
						ordinate = int(n*60*24 / TimeStep_min) + int(TT[r][c] / TimeStep_min)*baseflow_mult + m - 1;
						for(b=0; b<baseflow_mult; b++)
						{
							ordinate = ordinate+1;
							if(ordinate>=num_ordinates){ printf("\nError with Ordinate (%d) being higher than limit (%d) when using UH", ordinate, num_ordinates); continue;}
							Outflow[ordinate]+=IncrementalFlow * float(1.0/baseflow_mult) * (1.0-RunoffCoeff[r][c]);
						}

					}
				}
			}
		}
		GDALClose(RUNOFFDataset);

		//Now the unit hydrograph
		baseflow_mult=4;
		for(i=0; i<num_uh_ordinates; i++)
		{
			ordinate = int(n*60*24 / TimeStep_min) + int(uh_ts[i] / TimeStep_min);
			if(ordinate>=num_ordinates){ printf("\nError with Ordinate (%d) being higher than limit (%d) when using UH", ordinate, num_ordinates); continue;}
			Outflow_from_UH[ordinate]+=(TotalRunoffPerDay[n]*uh_val[i] * 0.2 /uh_sum)/(24*60*60);  //This gets to cubic meters per second

			//Let's make-up baseflow and see what happens.  Try adding in the same runoff, but put it multiple times in the future
			ordinate = int(n*60*24 / TimeStep_min) + int(uh_ts[i] / TimeStep_min)*baseflow_mult - 1;
			for(b=0; b<baseflow_mult; b++)
			{
				ordinate = ordinate+1;
				if(ordinate>=num_ordinates){ printf("\nError with Ordinate (%d) being higher than limit (%d) when using UH", ordinate, num_ordinates); continue;}
				Outflow_from_UH[ordinate]+=float(1.0/baseflow_mult)*(TotalRunoffPerDay[n]*uh_val[i] * 0.8 /uh_sum)/(24*60*60);  //This gets to cubic meters per second
			}
		}
	}
	fclose(InRasterList_File);


	//Now get the runoff data to cubic meters per second
	Outflow_File = fopen(Outflow_Name,"w");
	for(ordinate=0; ordinate<num_ordinates; ordinate++)
	{
		Outflow[ordinate] = Outflow[ordinate] / (TimeStep_min*60.0);
		fprintf(Outflow_File,"%f %f\n", Outflow[ordinate], Outflow_from_UH[ordinate]);
	}
	fclose(Outflow_File);
	return 0;
}