#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define SIZECOMPUESTOS 30/*maximum number of metabolites in the biochemical network*/
#define SIZEREACCIONES 30/*maximum number of reactions in the biochemical network*/


int lectura_compuestos(double a[][SIZECOMPUESTOS])
{
	
	int i=0;
	int j=0;
	
	FILE *g;
	
	g=fopen("input_compuestos.txt","rt");
	if(g==NULL)
	{
		printf("\n\nError al abrir el fichero con la informaci√≥n de compuestos.\n\n");
	}
	

	while(fscanf(g,"%lf %lf",&a[i][j],&a[i][j+1])!=EOF)
	{
		i++;
		a[i][2]=0;
				 
	}
	fclose(g);
	return(i);
}

int lectura_reacciones(double b[][SIZEREACCIONES])
{
	
	int i=0;
	int j=0;
	
	
	FILE *h;
	
	h=fopen("input_reacciones.txt","rt");
	if(h==NULL)
	{
		printf("\n\nError al abrir el fichero con la informaci√≥n de reacciones.\n\n");
	}
	
	while(fscanf(h,"%lf %lf %lf %lf %lf %lf",&b[i][j],&b[i][j+1],&b[i][j+2],&b[i][j+3],&b[i][j+4],&b[i][j+5])!=EOF)
	{
		
		i++;
		b[i][j+6]=0;
	}
	fclose(h);
	return (i);
	
	
}

int lectura_parametros(double c[][SIZEREACCIONES])
{
	FILE *m;
	m=fopen("input_parametros.txt","rt");
	if(m==NULL)
	{
		printf("\n\nError al abrir el fichero con la informaciÛn de par·metros.\n\n");
	}
	
	while(fscanf(m,"%lf %lf %lf",&c[0][0],&c[0][1],&c[0][2])!=EOF)
	{
		
	}
	fclose(m);
	
	
}	



double michaelis(double conc_sustrato,double km,double k2, double e0)
{
	double var;
	
	var=(k2*conc_sustrato*e0/(conc_sustrato+km));
	
	return(var);
}

int main()
{
    double tiempo=0;
    double almacen;
    int i;
    int j;
    double contador_impresion; 
    int bandera=0;
    double sustrato[2];
    double producto[2];
    int numcompuest;
    int numreacc;
    FILE *m;
    FILE *n;
	
    m=fopen("output1.txt","at");
	if(m==NULL)
	{
		printf("\n\nError al abrir el fichero de output 1.\n\n");
	}
    n=fopen("output2.txt","at");
	if(n==NULL)
	{
		printf("Error output2");
	}
	
	
    double compuestos[SIZECOMPUESTOS][SIZECOMPUESTOS]; 
	
	
    double reacciones[SIZEREACCIONES][SIZEREACCIONES]; 
	
    double PARAMETROS[SIZEREACCIONES][SIZEREACCIONES]; 
    contador_impresion=PARAMETROS[0][2]; 	
	
	numcompuest=lectura_compuestos(compuestos);
	numreacc=lectura_reacciones(reacciones);
	lectura_parametros(PARAMETROS);
	
	
	
    while(tiempo<=PARAMETROS[0][1])
    {
		    
	if(contador_impresion<PARAMETROS[0][2])
 	{
                contador_impresion=contador_impresion+1;
        }else
        {	
		fprintf(m,"%lf\t",tiempo);
		for(i=0;i<numcompuest;i++)
		{
			fprintf(m,"%lf ",compuestos[i][1]);
		}
		fprintf(m,"\n");
		
		fprintf(n,"%lf\t",tiempo);
		for(i=0;i<numcompuest;i++)
		{
			fprintf(n,"%lf ",compuestos[i][2]);
		}
		fprintf(n,"\n");
               contador_impresion=1;
        }
    
		
	
		for(i=0,j=0;i<numreacc;i++)
		{
			
			sustrato[0]=reacciones[i][1];
			producto[0]=reacciones[i][2];
			
			
			while(bandera==0)
			{
				if(sustrato[0]==compuestos[j][0])
				{
					sustrato[1]=compuestos[j][1];
					bandera=1;
					j=0;
					
				}else{
					j++;
				}
			}
			
			while(bandera==1)
			{
				if(producto[0]==compuestos[j][0])
				{
					producto[1]=compuestos[j][1];
					bandera=0;
					j=0;
				}else{
					j++;
				}
			}
			
			reacciones[i][6]=michaelis(sustrato[1],reacciones[i][3],reacciones[i][4],reacciones[i][5]);
			
		}
		
		for(i=0;i<numcompuest;i++)
		{
			
			compuestos[i][2]=0;
			
			almacen=compuestos[0][1];
			
			for(j=0;j<numreacc;j++)
			{	
				
				if((compuestos[i][0]==reacciones[j][1])) 
				{
					compuestos[i][1]=compuestos[i][1]-reacciones[j][6]*PARAMETROS[0][0];
					compuestos[i][2]=compuestos[i][2]-reacciones[j][6];
				}
				
				
				if((compuestos[i][0]==reacciones[j][2]))
				{
					compuestos[i][1]=compuestos[i][1]+reacciones[j][6]*PARAMETROS[0][0];
					compuestos[i][2]=compuestos[i][2]+reacciones[j][6];
				}
				
				
				
				
			}
			
			
			if(i==0)
			{
				compuestos[i][1]=almacen;
			}  
			
			
			if(compuestos[i][1]<0)
			{
				compuestos[i][1]=0;
			}
			
			
		}
		tiempo=tiempo+PARAMETROS[0][0];
	} 
    fclose(m);
    fclose(n);
    return 0;
}
