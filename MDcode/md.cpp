#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "para.h"
#include <math.h>
#include <time.h>
int main (int argc, char *argv[])
{
    srand(time(0));
    if(argc > 1) {
    ind = atoi(argv[1]); //index for differentiating the output files
	start=atoi(argv[2]); //starting configuration
	N = atoi(argv[3]); //Number of particles
	rho = atof(argv[4]); //density (reduced units)
	set_tem = atof(argv[5]); //temperature (reduced units)
    }


//  ind ;     //this value specifies the index for different runs of the program
              //for examole if I run the program for one set of condition then I 
              //would use a different ind value for that run. this is basically to 
              //differentiate between different runs.
//	start=0;  //this value starts with different initial conditions
              // see initialcondition() function for more details.




box = pow ((double(N)/rho), 1./3.);
                
cout<<"ind  "<<ind<<"start   "<<start<<"N   "<<N<<"rho   "<<rho<<"box   "<<box<<"\n";
initialcondition();
simulation();
rdfsample();
return 0;
}//end of main program
//----------------------------------------------------------------------
//-----------------------RDF Calculation--------------------------------
//----------------------------------------------------------------------
void rdfsample()
{
double delr = box/double(2*maxbin);//why 2*maxbin??
double pi = 3.1415926;
double boxinv = 1.0/box;//why are we calaculating this value?
double hbox = 0.5*box;//half box length
double rij2, r1;//rij2 is rij^square and r1 is the magnitude of the distance between the particles
int nsample = 0, bin;
Cordinates rij;
char IntStr[80];
ofstream ofrdf;
    sprintf( IntStr, "rdf%d.dat",ind);
    ofrdf.open (IntStr, ofstream::out | ofstream::app);
for(int i=0; i<maxbin; i++)
{
hist[i]=0;
}
cout<<"\n ------------------------------------------------ \n";
cout<<"\n --------------- RDF Calculation ----------------\n";
cout<<"\n ------------------------------------------------ \n";
for (int t=0; t<=5000; t++)
{
force();
integrate();
if(t%1 ==0) {
//cout<<t<<"\n";
nsample += 1;
for (int i = 0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
       rij.x = colloid[i].x - colloid[j].x ;
       rij.y = colloid[i].y - colloid[j].y ;
       rij.z = colloid[i].z - colloid[j].z ;
       
	   if (rij.x >  hbox) rij.x -= box;//MIC
       if (rij.x < -hbox) rij.x += box;//MIC
       if (rij.y >  hbox) rij.y -= box;//MIC
       if (rij.y < -hbox) rij.y += box;//MIC
       if (rij.z >  hbox) rij.z -= box;//MIC
       if (rij.z < -hbox) rij.z += box;//MIC
       rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
       r1 = sqrt(rij2);
       /*Instead of calculating square root of "r" outside the if condition
        *we can check whether r1^2 < hbox^2 and then if the condition satisfies
        *we calcuate the square root of "r". this will be computationally 
		*less intensive*/
       if(r1 <= hbox){
         bin = int(r1/delr);//Address of the particle i.e the bin number to which the particle j belongs wrt the ith partice. 
	hist[bin] += 2.0;//this is because there are pairs. jthe particle belongs to bin and vice versa fo jthe particle
} 
}
}
}
if(t%500 == 0)
{
    cout<<t<<"\n";
}
}

double idcon = 4./3.*pi*rho;
double rlow, rup, ideal, rbin;
if(ofrdf.is_open()) {
for(int i=0; i<maxbin; i++)
{
rlow = double(i)*delr;
rup = rlow + delr;
ideal = idcon*(rup*rup*rup - rlow*rlow*rlow);
hist[i] = hist[i]/double(N*nsample)/ideal;
rbin = rlow + delr/2.0;
            ofrdf << rbin <<"  "<<hist[i]<< "\n";
}
}
else {cerr << "unable to open file for rdf output \n";}
    ofrdf.close();
return;
}
//---------------------------------------------------------------------
//-----------------------Main Simulation--------------------------------
//----------------------------------------------------------------------
void simulation()
{
//int step = 10000;//my line of code
double E_msq = 0,E_rms,E,Eav;
for(int t=0; t<5000; t++) {
force();
integrate();
writeconf();
kinetic();
E = KE+U;
//cout<<E;

E_msq += RMS(E);//Average of E^2
Eav += EAV(E);//Average value of E
 
if(t%10 ==0) {
    printout(t);
    cout<<"time "<<t<<" "<<"\n";}

}

E_rms = sqrt(E_msq);//my line of code
Sigma_E = E_msq - (Eav*Eav);//standard deviation

cout<<"standard devaition for energy is "<<Sigma_E<<"\n";
cout<<"E_RMS value is  "<<E_rms;//my line of code



 


return;
}

//----------------------------------------------------------------------
//-----------------------Printing Values--------------------------------
//----------------------------------------------------------------------

void printout(double t)
{
char IntStr[80];
    ofstream ofout;
    sprintf( IntStr, "values%d.dat",ind);
    ofout.open (IntStr, ofstream::out | ofstream::app);
    if (ofout.is_open())
    {
            ofout << t<<"  "<<U<< "  " << KE << "  "<<U+KE<<"  "<<Temp<<"  "<<press<< "\n";
        }
        else {cerr << "unable to open file for config output \n";}
    ofout.close();
return;
}

//----------------------------------------------------------------------
//--------------------Kinetic Energy&Temp Calculation-------------------
//----------------------------------------------------------------------
void kinetic()
{
double vsq = 0.0;
for(int i=0; i<N; i++){
vsq += v[i].x*v[i].x +v[i].y*v[i].y+v[i].z*v[i].z;//sum of the squares of all the velocities
}
KE = vsq / 3.0;
Temp = vsq/ double(2*N);//This comes from the equipartition theory.

return;
}
//----------------------------------------------------------------------
//-----------------------Initialisation---------------------------------
//----------------------------------------------------------------------
void initialcondition()
{
    if (start ==0){   //simple cubic
	int index=0;
        int baseunit = ceil(pow(N, 1.0/3.0));/*ceil(x) : Returns the smallest integer that is greater
		                                     than or equal to x (i.e : rounds up the nearest integer).
											 But why use it here?*/
	double lengthcell = box/baseunit;/*the length of the cell or we can say the size of each lattice in the simple cubic structure is
	                                  basically the box length/(number of particles)^(1/3) that is what is done in the previous step 
									  when declaring baseunit = N^(1/3).Basically N^(1/3) ensures that along each length of the box 
									  the particles are arranged in simple cubic fashion*/
	for (int iz=0; iz<baseunit; iz++) {
	  for(int iy=0; iy<baseunit; iy++) {
	     for(int ix=0; ix<baseunit; ix++){
		colloid[index].x=double(ix)*lengthcell;//didnt undertand this?
		colloid[index].y=double(iy)*lengthcell;
		colloid[index].z=double(iz)*lengthcell;
		index+=1;
             }
	}
	}	 
	}
	else if(start == 1){                          //fcc crystal
        int index = 0;
        int baseunit;
	baseunit = ceil(pow(N/4, 1.0/3.0));//why N/4?
        double lengthcell = double(box)/double(baseunit);
	cout<<"length"<<lengthcell<<"\n";
	double facto = 1.0/sqrt(2.0); //what is this?//
        for (int iz=0; iz < baseunit ; iz++) {
            for (int iy=0; iy < baseunit ; iy++) {
                for (int ix = 0; ix < baseunit; ix++) {
                    colloid[index].x = double(ix)*lengthcell;
                    colloid[index].y = double(iy)*lengthcell;
                    colloid[index].z = double(iz)*lengthcell;
		    colloid[index+1].x = colloid[index].x+lengthcell/2.0;
                    colloid[index+1].y = colloid[index].y+lengthcell/2.0;
                    colloid[index+1].z = colloid[index].z;//why not add lenghtcell/2 ?
                    colloid[index+2].x = colloid[index].x+lengthcell/2.0;
                    colloid[index+2].y = colloid[index].y;
                    colloid[index+2].z = colloid[index].z+lengthcell/2.0;
	            colloid[index+3].x = colloid[index].x;
                    colloid[index+3].y = colloid[index].y+lengthcell/2.0;
                    colloid[index+3].z = colloid[index].z+lengthcell/2.0;
                    index+=4;
                }
            }
        }
    }               
        
    else if (start == 2) {                 //create random initial configuration// begin //
    
    colloid[0].x = (randnum()-0.5)*box;
    colloid[0].y = (randnum()-0.5)*box;            
    colloid[0].z = (randnum()-0.5)*box;
    double r1, r2;
    Cordinates rij, t; 
    double hbox = box/2.0;
    for (int i=1;i<N; i++) {
        do {
            t.x = (randnum()-0.5)*box;
            t.y = (randnum()-0.5)*box;
            t.z = (randnum()-0.5)*box;
            for (int j=0; j<i; j++) {
		rij.x = t.x - colloid[j].x;
		rij.y = t.y - colloid[j].y;
		rij.z = t.z - colloid[j].z;
                if (rij.x >  hbox) rij.x -= box;
        	if (rij.x < -hbox) rij.x += box;
        	if (rij.y >  hbox) rij.y -= box;
        	if (rij.y < -hbox) rij.y += box;
        	if (rij.z >  hbox) rij.z -= box;
        	if (rij.z < -hbox) rij.z += box;
                r2 = rij.x*rij.x+rij.y*rij.y+rij.z*rij.z;
                r1 = sqrt(r2);
                if(r1<1.0) break;
            }
        }while (r1<1.0);
        
        colloid[i].x = t.x;
        colloid[i].y = t.y;
        colloid[i].z = t.z;
    }
    }

   	else if (start == 3) 
    {               //read from file// begin //

        char dum[100];
        FILE * fin;
        char IntStr[80];
        sprintf( IntStr, "configin.%d.xyz", ind);
        fin = fopen (IntStr, "r");
        fscanf(fin, "%d", &N);
        fscanf(fin, " ");
        for (int i = 0; i<N; i++) {
            fscanf(fin, "%lf%lf%lf", &colloid[i].x, &colloid[i].y, &colloid[i].z);
        }
    }
    
    writeconf();

    //---------------------initialize velocities------------------------
    sumv.x = 0.0;
    sumv.y = 0.0;
    //sumv.z = 0.0;
    double vsq = 0.0;
    for(int i=0; i<N; i++){//why for loop initialized from i = 0 rather than i = 1 ?
    v[i].x = (randnum()-0.5);//why randnum() - 0.5??
    v[i].y = (randnum()-0.5);
    v[i].z = (randnum()-0.5);
    sumv.x += v[i].x;
    sumv.y += v[i].y;
    sumv.z += v[i].z;
    vsq += v[i].x*v[i].x +v[i].y*v[i].y +v[i].z*v[i].z;   
    }
    cout<<"x momentum"<<" "<<sumv.x<<" "<<"y momentum"<<" "<<sumv.y<<" "<<"z momentum"<<" "<<sumv.z<<"\n";
    sumv.x /= double(N);
    sumv.y /= double(N);
    sumv.z /= double(N);
    vsq /= double(N);

    double temp_init = vsq / double(2);
    cout<<"Temperature initialized"<<"  "<<temp_init<<"  "<<"\n";
    double scale = sqrt(2.0*set_tem/vsq);//From where does this scaling come from?
    sumvn.x = 0.0;
    sumvn.y = 0.0;
    sumvn.z = 0.0;
    for(int i=0; i<N; i++){
    	v[i].x = (v[i].x-sumv.x)*scale;
    	v[i].y = (v[i].y-sumv.y)*scale;
    	v[i].z = (v[i].z-sumv.z)*scale;
        sumvn.x += v[i].x;
        sumvn.y += v[i].y;
        sumvn.z += v[i].z;
        vsq += v[i].x*v[i].x +v[i].y*v[i].y+v[i].z*v[i].z;
	}
    vsq /= double(N);
    double temp_init_reset = vsq / double(2);
    cout<<"Temperature after reset"<<"  "<<temp_init_reset<<"  "<<"\n";
    cout<<"x momentum"<<" "<<sumvn.x<<" "<<"y momentum"<<" "<<sumvn.y<<" "<<"\n"<<"z momentum"<<" "<<sumvn.z<<"\n";

return;
}
//----------------------------------------------------------------------
double randnum()
{
    return double(rand())/double(RAND_MAX);
}
//----------------------------------------------------------------------
//-----------------------Force Calculation------------------------------
//----------------------------------------------------------------------
void force()
{
double boxinv = 1.0/box;
double boxinv3 = boxinv*boxinv*boxinv;
double hbox = 0.5*box;
double rcut = 2.5;
double rcutsq = rcut*rcut;
double pi = 3.1415926;
for (int i = 0; i<N; i++){//why for loop initialeized from i=0 rather than i=1?
f[i].x=0.0;
f[i].y=0.0;
f[i].z=0.0;
}
int ncut=0;
U = 0;
double W = 0;
Cordinates ri, rij, fij,fi0;
double r1, rij2,sr2,sr3,sr9, sr6, Uij,Wij,fr,UTC,WTC,ideal;
for(int i=0; i<N; i++){
	ri.x = colloid[i].x;
	ri.y = colloid[i].y;
	ri.z = colloid[i].z;
	fi0.x = f[i].x;
	fi0.y = f[i].y;
	fi0.z = f[i].z;
	for(int j = i+1; j<N; j++){
		rij.x = ri.x - colloid[j].x ;
		rij.y = ri.y - colloid[j].y ;
		rij.z = ri.z - colloid[j].z ;
		if (rij.x >  hbox) rij.x -= box; 
		if (rij.x < -hbox) rij.x += box;
		if (rij.y >  hbox) rij.y -= box;
		if (rij.y < -hbox) rij.y += box;
		if (rij.z >  hbox) rij.z -= box;
		if (rij.z < -hbox) rij.z += box;
		rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
		r1 = sqrt(rij2);
		if(r1 < rcut) {
			sr2 = 1.0/rij2;
			sr6 = sr2*sr2*sr2;
			Uij = sr6*(sr6 - 1.0);
			Wij = -Uij-sr6*sr6;
			U += Uij;
			W += Wij;
			fr = -Wij/rij2;
			fij.x = fr*rij.x;
			fij.y = fr*rij.y;
			fij.z = fr*rij.z;
			fi0.x += fij.x;
			fi0.y += fij.y;
			fi0.z += fij.z;
			f[j].x -= fij.x;
			f[j].y -= fij.y;
			f[j].z -= fij.z; 
			ncut += 1;
		}
	}
	f[i].x = fi0.x;
	f[i].y = fi0.y;
	f[i].z = fi0.z;
}
sr2 = 1.0/rcutsq;
sr6 = sr2*sr2*sr2;
sr3 = sr2/rcut;
sr9 = sr3*sr3*sr3;
Uij = sr6*(sr6 - 1.0);
U -= double(ncut)*Uij;

for(int i=0; i<N; i++){
f[i].x *= 24.0;
f[i].y *= 24.0;
f[i].z *= 24.0;
}

W = 24.0*W/3.0;
U *= 4.0;

UTC = 8./9.*pi*rho*double(N)*(sr9 - 3.0*sr3);
WTC = 16./9.*pi*rho*rho*double(N)*(2.*sr9 - 3.*sr3);
U += UTC;
W += WTC;
ideal = double(N)*boxinv3*set_tem;
press = ideal+W*boxinv3;
//cout<<"UTC     "<<UTC<<" "<<"\n";

return;
}
//----------------------------------------------------------------------
//---------------------------Integration--------------------------------
//----------------------------------------------------------------------
void integrate()
{
double dt2=dt/2.0;
double dtsq2 = dt*dt2;
double hbox = 0.5*box;



for(int i=0; i<N; i++) {
	dummy[i].x = colloid[i].x;
	dummy[i].y = colloid[i].y;
	dummy[i].z = colloid[i].z;
	colloid[i].x += (dt*v[i].x + dtsq2*f[i].x);
	colloid[i].y += (dt*v[i].y + dtsq2*f[i].y);
	colloid[i].z += (dt*v[i].z + dtsq2*f[i].z);
	//if(i>0)
	//{
	//	if((dummy[i].x>0&&colloid[i].x<0)||(dummy[i].x<0&&colloid[i].x>0))
		
	//}
	if(colloid[i].x >= box) colloid[i].x -= box;
	if(colloid[i].x <= 0.0) colloid[i].x += box;
 	if(colloid[i].y >= box) colloid[i].y -= box;
    if(colloid[i].y <= 0.0) colloid[i].y += box;
 	if(colloid[i].z >= box) colloid[i].z -= box;
    if(colloid[i].z <= 0.0) colloid[i].z += box;

	v[i].x += (dt2* f[i].x);
	v[i].y += (dt2* f[i].y);
	v[i].z += (dt2* f[i].z);
}

force();

for(int i=0; i<N; i++) {
	v[i].x += (dt2* f[i].x);
	v[i].y += (dt2* f[i].y);
	v[i].z += (dt2* f[i].z);
}

return;
}

//----------------------------------------------------------------------
//---------------------Printing coordinates-----------------------------
//----------------------------------------------------------------------
void writeconf()
{
    char IntStr[80];
    ofstream of;
    sprintf( IntStr, "config%d.xyz",ind);
    of.open (IntStr, ofstream::out | ofstream::app);
    if (of.is_open())
    {
        of << N<< "\n";
        of << "\n";
	if ((start==0) || (start ==2)){
        for (int i=0; i<N; i++)
        {
            of << 'h'<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< colloid[i].z << "\n";
        }    
    } 
	else if (start ==1) {
	 for (int i=0; i<N; i++)
        {
	   if (i%4 ==0) {
            of <<"he"<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< colloid[i].z << "\n";
	   }else {of << 'h'<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< colloid[i].z << "\n";}
	}
        }

        }
	else {cerr << "unable to open file for config output \n";}
    of.close();
    
    return;
}
//---------------------------------------------------------------------------------------//
//my lines of code---------------------
double RMS(double p)
{
	double E_msq;
	int step =5000;
	E_msq = (p*p)/step;
	
	
	return(E_msq);
}
double EAV(double q)
{
	double Eav;
	int step = 5000;
	Eav = q/step;
	
	return(Eav);
}

