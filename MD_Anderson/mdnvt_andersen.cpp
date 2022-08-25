
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "para.h"
#include <math.h>
#include <ctime>
int main (int argc, char *argv[])
{
    srand(time(0));///what is this?
    /*if(argc > 1)
    {
    ind = atoi(argv[1]);
	start=atoi(argv[2]);
	N = atoi(argv[3]);
	rho = atof(argv[4]);
	set_tem = atof(argv[5]);
    }*/
    ind=6; start=2; N=2; rho=0.05; set_tem=0.5;

box = pow ((double(N)/rho), 1./3.);
cout<<"ind  "<<ind<<"start   "<<start<<"N   "<<N<<"rho   "<<rho<<"box   "<<box<<"\n";
initialcondition();
simulation();
     //msdsample(10,0,5000); // or msdsample2()
rdfsample();
//anderson();
return 0;
}

double gauss()//what does this function do?
{
    double pr=2,v1,v2,l1;
    do{
        v1 = 2.*randnum()-1.0;
        v2 = 2.*randnum()-1.0;
        pr=v1*v1+v2*v2;
    }while(pr >=1.0);
        l1 = v1*sqrt(-2*log(pr)/(pr));
    return sqrt(set_tem)*l1;
}


void anderson(double nu)
{
    int i;
    double T_req=set_tem,ran[N];
    double u1_ran, u2_ran, u3_ran, a1, a2, a3;
    float n1_ran, n2_ran, n3_ran, pi=3.1415;
    for(i=0;i<N;i++)
    {
        if(randnum()<=nu*dt)
        {
            v[i].x=gauss();
            v[i].y=gauss();
            v[i].z=gauss();
            //cout<<v[i].x<<"\t"<<v[i].y<<"\t"<<v[i].z<<"\n";

        }
    }
}
/// ### Updated on 12th Mar 2015 Version 1.0
// function to calculate Mean Square Displacement V1.3
void msdsample(int delt, int stime, int etime)
{
	//Initialization
	char IntStr[80];
	int nbins = 0,timesteps = 1, dt = 0,time = etime - stime;
	double msd[1000] = {0}, msdave[1000] = {0};          // mean square displacement and average mean square displacement
	Cordinates rt1[N], rt0[N];

	//cout<<"Computing absolute positions \n";
	sprintf( IntStr, "absconfig%d.xyz",ind);
	absconfig(IntStr, delt, stime, etime);//what does this function do?

	cout<<"Computing <Mean Square Displacement> for "<<N<<" atoms \n";
	sprintf( IntStr, "absconfig%d.xyz",ind);
   	ifstream file(IntStr,std::ifstream::binary);
		for(dt = delt; dt <= time/2; (dt = dt + delt),timesteps++) 
		{		//cout<<"dt == "<<dt<<"\n";
		    nbins = 0;
			for(int t0 = 0; t0 < (time - dt); (t0 = t0 + dt),nbins++) 
			{	//cout<<"t == "<<t0<<"\t";
				readconf(file, true, t0/delt, &rt0[0]);
				readconf(file, false, (dt/delt)-1, &rt1[0]);
				for(int i = 0; i < N; i++) 
				{
					msd[nbins] +=  pow(( rt1[i].x - rt0[i].x ), 2) + pow(( rt1[i].y - rt0[i].y ), 2) + pow(( rt1[i].z - rt0[i].z ), 2) ;
				}
				
			msd[nbins] /= double(N);
			msdave[timesteps-1] += msd[nbins];   //cout<<"bin no "<<nbins+1<<"\t msd of bin  == "<<msd[nbins]<<"\n";
			}
		msdave[timesteps-1] /= nbins; 	//cout<<"Average msd [timesteps-1] "<<msdave[timesteps-1]<<"\n";
		}

	file.close();
	sprintf( IntStr, "msd%d.dat",ind);

	cout<<"Saving <Mean Square Displacement> to file "<<IntStr<<endl;
	ofstream ofmsd;
	ofmsd.open(IntStr, ofstream::out | ofstream::app);
	if(ofmsd.is_open()) {								//  Writing the MSDaverage to file
		for(int i = 1; i < timesteps; i++) {
			dt = delt * i;
			ofmsd<<dt<<"  "<<msdave[i-1]<<" "<<log(dt)<<" "<<log(msdave[i-1])<<" "<<dt<<" "<<sqrt(msdave[i-1])<<"\n";
		}
	}
	else {
	cerr << "unable to open file for msd average output \n";
	}
	ofmsd.close();

}
//Function converts config to absolute config
void absconfig(char* outputfilename,int delt, int stime, int etime)
{
	cout<<"Computing absolute positions \n";
	char IntStr[80];
	Cordinates r[N], rprev[N],rabs[N];
	sprintf( IntStr, "config%d.xyz",ind);
    ifstream file(IntStr,std::ifstream::binary);//reading config file as an input

	for(int t = stime; t <= etime; t++)	
	{
		//cout<<t<<endl;
		if(t > stime) 
		{
			readconf(file, false, 0, &r[0]);//what does this function do?
			for (int i = 0; i<N; i++) 
			{
			   rabs[i].x = r[i].x + nint( (rprev[i].x - r[i].x )/box) * box ;
			   rabs[i].y = r[i].y + nint( (rprev[i].y - r[i].y )/box) * box ;
			   rabs[i].z = r[i].z + nint( (rprev[i].z - r[i].z )/box) * box ;
			}
			copycoordinates(&rabs[0], &rprev[0]);//"size" argument missing for this function
		}else 
		{
			readconf(file, true, t, &rabs[0]);
			readconf(file, true, t, &rprev[0]);		
		}
		if(t%delt == 0) 
		{
			writecoord(&rabs[0], outputfilename); 				//Sampling absolute position every delt time
		}
	}
	file.close();
}
//Reads config file and return coordinates for time t
std::ifstream& readconf(std::ifstream& file,bool abs ,unsigned int time, Cordinates *r)//config.xyz file is passed as a parameter
{
	string str; double num;
	static streampos fpos = 0;
	if (abs) { file.seekg(std::ios::beg); }
	else{ file.seekg(fpos); }
	// skiping time
	for(int j = 0; j < time; j++) 
	{ 
		file>>str; 	
		for(int i=0; i < N; i++)
		{ 
		file>>str>>num>>num>>num; 
		} 	
	}
	//reading data for time "time"
	file>>str; 	
	for(int i = 0; i < N; i++) 
	{ 
	file >> str >> r[i].x >> r[i].y >> r[i].z;	
	}
	fpos = file.tellg();  //remembering the file read position
	return file;
}
//Writes coordinates to config file
void writecoord(Cordinates *input, char *fname) {
    ofstream of;
    of.open (fname, ofstream::out | ofstream::app);
    if (of.is_open())
    {
        of << N<< "\n";
        of << "\n";
		for (int i=0; i<N; i++)
		{
			of <<"h"<<" "<<input[i].x << "  " << input[i].y << "  "<< input[i].z << "\n";
		}
	}
	else {cerr << "unable to open file for config output \n";}
    of.close();
    return;
}
// Use an alternate for native "nint" function
//int nint(double d)
//{
//   if (d == 0) {    	return 0;  }
//   else if(d >0) {  	return ceil(d);    }
//   else if(d <0) {  	return floor(d);    }
//}
int nint(double d)
{
   return floor(d+0.5);
}
// Function to copy cordinates from one variable to another
void copycoordinates(Cordinates *a, Cordinates *b, int size ) 
{
	for(int i = 0; i< size; i++) 
	{
		b[i].x = a[i].x ;
		b[i].y = a[i].y ;
		b[i].z = a[i].z ;
	}
}
/// ###

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
	hist[i]=0;//what is hist?
	}

	for (int t=0; t<20000; t++)
	{
	force();
	integrate();
	anderson(1);
		if(t%1 ==0) 
		{
		cout<<t<<"\n";
		nsample += 1;
			for (int i = 0; i<N; i++) 
				{
    				for(int j=i+1; j<N; j++) 
    					{
       					rij.x = colloid[i].x - colloid[j].x ;
       					rij.y = colloid[i].y - colloid[j].y ;
     					rij.z = colloid[i].z - colloid[j].z ;
       
					if (rij.x >  hbox) rij.x -= box;//Periodic boundary condition
       					if (rij.x < -hbox) rij.x += box;//Periodic boundary condition
				       	if (rij.y >  hbox) rij.y -= box;//Periodic boundary condition
       					if (rij.y < -hbox) rij.y += box;//Periodic boundary condition
       					if (rij.z >  hbox) rij.z -= box;//Periodic boundary condition
       					if (rij.z < -hbox) rij.z += box;//Periodic boundary condition
       					rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
       					r1 = sqrt(rij2);
       /*Instead of calculating square root of "r" outside the if condition
        *we can check whether r1^2 < hbox^2 and then if the condition satisfies
        *we calcuate the square root of "r". this will be computationally 
		*less intensive*/
       						if(r1 <= hbox)
       						{
        					bin = int(r1/delr);//Address of the particle i.e the bin number to which the particle j belongs wrt the ith partice. 
						hist[bin] += 2.0;//this is because there are pairs. jthe particle belongs to bin and vice versa fo jthe particle
						} 
					}
				}
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

void simulation()
{

for(int t=0; t<200000; t++)
{
    if(t%100 ==0)
    {
    printout(t);
    cout<<"time "<<t<<" "<<"\n";
    }
force();
integrate();
anderson(1);
writeconf();
kinetic();

}

return;
}
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


void kinetic()
{
double vsq = 0.0;
for(int i=0; i<N; i++){
vsq += v[i].x*v[i].x +v[i].y*v[i].y+v[i].z*v[i].z;
}
KE = vsq / 2.0;
Temp = vsq/ double(3*N);

return;
}

void initialcondition()
{
    if (start ==0){   //simple cubic
	int index=0;
        int baseunit = ceil(pow(N, 1.0/3.0));
	double lengthcell = box/baseunit;
	for (int iz=0; iz<baseunit; iz++) {
	  for(int iy=0; iy<baseunit; iy++) {
	     for(int ix=0; ix<baseunit; ix++){
		colloid[index].x=double(ix)*lengthcell;
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
	baseunit = ceil(pow(N/4, 1.0/3.0));
        double lengthcell = double(box)/double(baseunit);
	cout<<"length"<<lengthcell<<"\n";
//	double facto = 1.0/sqrt(2.0);
        for (int iz=0; iz < baseunit ; iz++) {
            for (int iy=0; iy < baseunit ; iy++) {
                for (int ix = 0; ix < baseunit; ix++) {
                    colloid[index].x = double(ix)*lengthcell;
                    colloid[index].y = double(iy)*lengthcell;
                    colloid[index].z = double(iz)*lengthcell;
		    colloid[index+1].x = colloid[index].x+lengthcell/2.0;
                    colloid[index+1].y = colloid[index].y+lengthcell/2.0;
                    colloid[index+1].z = colloid[index].z;
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

   	else if (start == 3) {               //read from file// begin //

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
    // initialize velocities
    sumv.x = 0.0;
    sumv.y = 0.0;
    sumv.z = 0.0;
    double vsq = 0.0;
    for(int i=0; i<N; i++){
    v[i].x = (randnum()-0.5);
    v[i].y = (randnum()-0.5);
    v[i].z = (randnum()-0.5);
    sumv.x += v[i].x;
    sumv.y += v[i].y;
    sumv.z += v[i].z;
    vsq += v[i].x*v[i].x +v[i].y*v[i].y+v[i].z*v[i].z;
    }
    cout<<"x momentum"<<" "<<sumv.x<<" "<<"y momentum"<<" "<<sumv.y<<" "<<"z momentum"<<" "<<sumv.z<<"\n";
    sumv.x /= double(N);
    sumv.y /= double(N);
    sumv.z /= double(N);
    vsq /= double(N);

    double temp_init = vsq / double(3);
    cout<<"Temperature initialized"<<"  "<<temp_init<<"  "<<"\n";
    double scale = sqrt(3.0*set_tem/vsq);
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
    double temp_init_reset = vsq / double(3);
    cout<<"Temperature after reset"<<"  "<<temp_init_reset<<"  "<<"\n";
    cout<<"x momentum"<<" "<<sumvn.x<<" "<<"y momentum"<<" "<<sumvn.y<<" "<<"z momentum"<<" "<<sumvn.z<<"\n";

return;
}

double randnum()
{
    return double(rand())/double(RAND_MAX);
}

void force()
{
double boxinv = 1.0/box;
double boxinv3 = boxinv*boxinv*boxinv;
double hbox = 0.5*box;
double rcut = 2.5;
double rcutsq = rcut*rcut;
double pi = 3.1415926;
for (int i = 0; i<N; i++){
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

void integrate()
{
double dt2=dt/2.0;
double dtsq2 = dt*dt2;
double hbox = 0.5*box;

for(int i=0; i<N; i++) {
	colloid[i].x += (dt*v[i].x + dtsq2*f[i].x);
	colloid[i].y += (dt*v[i].y + dtsq2*f[i].y);
	colloid[i].z += (dt*v[i].z + dtsq2*f[i].z);
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
	else if (start ==1)
        {
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

