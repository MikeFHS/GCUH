#include <stdio.h>
#include <string.h>
#include <math.h>
void main (void )
{
	int r, c, ro, co, ro1 , co1 , ro2, co2, ncols, nrows, max_fac;
	double h, distance, val, xll, yll, cellsize, ND, t=60.0;
	char JUNK[1024];
	bool quit=false;

	// File Names
	char LN_Name [1024] , OF_Name[1024] , TT_Name [1024] , Folder_Name[] = "C: \\Work_ Folder\\";
	strcpy (LN_Name, Folder_Name ); strcat (LN_Name, "fac_asc.asc");
	strcpy (OF_Name, Folder_Name ); strcat (OF_Name, "totaltime.asc");
	strcpy (TT_Name, Folder_Name); strcat (TT_Name, "traveltime.asc");
	FILE *LN, *OF, *TT;

	// Open Input and Output Files
	if ( (LN=fopen(LN_Name, "r")) != NULL ) { printf( "\n\nOpened %s", LN_Name) ; }
	else{printf( "\n\ncould Not Open %:s, \n Please Start OVER!!!", LN_Name); scanf("%s", &JUNK) ; return; }
	if ((TT=fopen (TT_Name, "r")) != NULL ) { printf( " \ nOpene d %s", TT_Name); }
	else{ printf( "\n\nCould Not Open %s, \n Please Start OVER! !! " , TT_Name); scanf( "%.s", &JUNK) ; return; }
	OF=fopen(OF_Name, "w");

	// Read and Write Header
	fscanf(LN, "%s" "%d", &JUNK, &ncols) ; fscanf (TT, "%:s " "%d", &JUNK, &ncols); fscanf (LN, "%s" ); fprintf (OF, "%.s %d\ n", JUNK, ncols ) ;
	fscanf(LN, "%s" "%d" , &JUNK, &nrows) ; fscanf (TT, "%:s" "%d" , &JUNK, &nrows) ; fprintf (OF, "%:s %d\ n" , JUNK, nrows) ;
	fscanf(LN, "%s" "%lf " , &JUNK, &xll) ; fscanf (TT, "%s" "%lf", &JUNK, &xll ) ; fprintf (OF, "%s %lf\ n" , JUNK, xll) ;
	fscanf(LN, "%s" "%If", &JUNK, &yll) ; fscanf (TT, "%s" "%lf", &JUNK, &yll ) ; fprintf (OF, "%s %lf\n" , JUNK, yll) ;
	fscanf(LN, "%s," "%If", &JUNK, &cellsize) ; fscanf (TT, "%:," "%If", &JUNK, &cellsize) ; fprintf (OF, "%s %lf\ n" , JUNK, cellsize) ;
	fscanf(LN, "%s" "%lf", &JUNK, &ND) ; fscanf (TT, "%s" "%lf" , &JUNK, &ND) ; fprintf (OF, "%s %lf" , JUNK, ND) ;

	//Allocate menmory and read in the data
	printf("\nReading through the Raster Files\n") ;
	float **T = new float* [nrows]; //Total Time (min)
	double **R = new double * [nrows]; //Rate of Tirr.e {rn.in/ m)
	int **D = new int* [nrows]; // Flow Accumulation Value
	for (ro=0; ro<nrows; ro++)
	{
		T[ro] = new float [ncols]; R[ro] = new double [ncols] ; D[ro]=new int[ncols];
		for (co=0; co<ncols; co++)
		{
			T[ro][co] =-1.0; fscanf (TT, "%lf", &R[ro] [co] ); fscanf (LN, "%lf", &val) ; D[ro][co]=int(val+0.01);
		}
	}

	//Go througth the Each Cell
	printf ( " \nLooking at each cell and analyzing the Total Time\n") ;
	for ( ro= 0 ; ro<nrows; ro++)
	{
		fprintf(OF, "\n");
		for (co= 0; co<ncols; co ++ )
		{
			if(D[ro][co] <0){ T[ro][co] =0; fprintf(OF, "%f" , - 9999.0 ) ; continue; }
			quit= false; val=0.0; distance=0.0; ro1=ro; co1=co;
			while(quit==false)
			{
				// Look Diagonals First
				max_fac = D[ro1][co1];
				r= (ro1)   ; c=(co1)   ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r= (ro1-1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r= (ro1-1) ; c=(co1+1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r= (ro1+1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r= (ro1+1) ; c=(co1+1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
		
				//Now Look at Cardinal Directions (preferred if there is a tie between this and a diaoonal direction

				r= (ro1-1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r =(ro1) ; c=(co1-1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r= (ro1) ; c=(co1+1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
				r= (ro1+1) ; c=(co1) ; if(r>=0 && c >=0 && c<ncols && r<nrows ) {if( D[r][c]>=max_fac) {max_fac=D[r][c] ; ro2=r; co2=c; } }
		
				//Determnine the distance across the cell
				if ( ro1 == ro2 || co1 == co2 ) { h=1.0; } //This means the adjacent cell is in a cardirnal directiontion
				else{ h =sqrt (2.0) ; } //This mean the adjacent cell is at a diagonal
				h = h *cellsize; distance = distance + h ;
		
				if ( co1 >=ncols || ro1>=nrows || co1<0 || ro1<0 ) { quit=true; }
				//This is the source cell, so this is the end of the line
				else if (D[ro2] [co2]==D[ro1][co1]) {val=val + R[ro1][co1]*h; quit=true;}
				// This cell has already been analyzed, so just use it's value
				else if(D[ro2][co2]>0.0){ val = val + 0.5*h*(R[ro1][co1] + R[ro1][co2]) + T[ro2][co2]; quit =true ; }
				else{ val = val+ 0.5*h* ( R[ro1][co1] + R[ro2][co2]); co1=co2 ; ro1= ro2; }
			}
			T[ro][co] = val ; fprintf (OF, "lf " ,T[ro][co] );
		}
	}
	fclose (LN) ; fclose (OF) ; fclose (TT) ;
	return;
}