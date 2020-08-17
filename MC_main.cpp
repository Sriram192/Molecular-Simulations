#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include<cmath>
#include<vector>
#include<cstddef>
#include"param.h"
#include<math.h>
using namespace std;
int main(int argc, char *argv[])
{
    if(argc>1)
    {
    run_index   =   atoi(argv[1]); //index number for the output files
    N           =   atoi(argv[2]); // Number of particles
    Temp        =   atof(argv[3]);    //Temperature of the simulation
    rho         =   atof(argv[4]); //density of the system
    start_index =   atoi(argv[5]); //Initialization index for different starting configurations
    }
    

    box = pow(static_cast<double>(N/rho),1.0/3.0);
    std::cout<<"\n *************** sim parameters *************** \n";
    std::cout<<"index : "<<run_index<<" No. of Particles : "<<N<<" Temperature : "<<Temp<<" density : "<<rho<<" Box length : "<<box<<"\n";
    
    initial_config(start_index);
    
    /*for(auto pos : Particle)
    {
        std::cout<<pos.x<<"  "<<pos.y<<"  "<<pos.y<<"\n";
    }*/

}
void initial_config(int write_index)
{
	if(write_index==1)
	{
		
        sc(write_index); 
	}
	else if(write_index==2)
	{
        
		fcc(write_index);
	}
	else if(write_index == 3)
    {
        bcc(write_index); 
    }
    else
	{
		std::cerr<<" Enter a valid initial configuration index --- either 1 for sc or 2 for fcc ";
		exit(0);
	}
	
}
void sc(int write_index) // simple cubic
{
    cout<<"inside SC\n";
    int units = ceil(pow(N,1.0/3.0));
    double baselength = box/units;
    int particle_index = 0;
    cout<<units<<" "<<baselength<<" "<<box<<"\n";
    
    for(int ix = 0; ix < units; ix++)
    {
        for(int iy = 0; iy < units; iy++)
        {
            for(int iz = 0; iz < units; iz++)
            {
                Particle[particle_index].x = static_cast<double>(ix * baselength);
                Particle[particle_index].y = static_cast<double>(iy * baselength);
                Particle[particle_index].z = static_cast<double>(iz * baselength);
                particle_index += 1;
                
                //cout<<Particle[particle_index].x<<" "<<Particle[particle_index].y<<" "<<Particle[particle_index].z<<"\n";
            }
            
        }
        
    }
    Particle.resize(particle_index);cout<<"vecctor size : "<<Particle.size()<<"\n";
    write_config(write_index);
}
void fcc(int write_index) //face centered cubic
{
    int particle_index = 0;
    int units;
	units = ceil(pow(N/4, 1.0/3.0));
    double lengthcell = static_cast<double>(box/units);
//	
double facto = 1.0/sqrt(2.0);
        for (int iz=0; iz < units ; iz++) {
            for (int iy=0; iy < units ; iy++) {
                for (int ix = 0; ix < units; ix++) {
                    Particle[particle_index].x = double(ix)*lengthcell;
                    Particle[particle_index].y = double(iy)*lengthcell;
                    Particle[particle_index].z = double(iz)*lengthcell;
		            Particle[particle_index+1].x = Particle[particle_index].x+lengthcell/2.0;
                    Particle[particle_index+1].y = Particle[particle_index].y+lengthcell/2.0;
                    Particle[particle_index+1].z = Particle[particle_index].z;
                    Particle[particle_index+2].x = Particle[particle_index].x+lengthcell/2.0;
                    Particle[particle_index+2].y = Particle[particle_index].y;
                    Particle[particle_index+2].z = Particle[particle_index].z+lengthcell/2.0;
	                Particle[particle_index+3].x = Particle[particle_index].x;
                    Particle[particle_index+3].y = Particle[particle_index].y+lengthcell/2.0;
                    Particle[particle_index+3].z = Particle[particle_index].z+lengthcell/2.0;
                    particle_index+=4;
                    
                }
                
            }
            
        }
        Particle.resize(particle_index); cout<<"vecctor size : "<<Particle.size()<<"\n";
        write_config(write_index);
}
void bcc(int write_index) //body centered cubic
{
    int particle_index = 0;
    int units = ceil(pow(N/2,1.0/3.0));
    double baselength = static_cast<double>(box/units);
    for(int ix = 0; ix < units; ix++)
    {
        for(int iy = 0; iy < units; iy++)
        {
            for(int iz = 0; iz < units; iz++)
            {
                Particle[particle_index].x = static_cast<double>(ix * baselength);
                Particle[particle_index].y = static_cast<double>(iy * baselength);
                Particle[particle_index].z = static_cast<double>(iz * baselength);
                Particle[particle_index+1].x = Particle[particle_index].x + static_cast<double>(baselength/2.0);
                Particle[particle_index+1].y = Particle[particle_index].y + static_cast<double>(baselength/2.0);
                Particle[particle_index+1].z = Particle[particle_index].z + static_cast<double>(baselength/2.0);
                particle_index += 2;
                
            }            
        }        
    }
    Particle.resize(particle_index);cout<<"vecctor size : "<<Particle.size()<<"\n";
    write_config(write_index);
}
double potential(int particle_index, double part_x, double part_y, double part_z)
{
    double rij , xij , yij , zij , r_sq , r, r12, r6 , U = 0.0, U_tail = 0.0;
    double r_cut = 2.5;
    double r_cut_sq = r_cut*r_cut;
    for(int other_particle = 0; other_particle < N; other_particle++)
    {
        if(other_particle != particle_index)
        {
            xij = static_cast<double>(part_x - Particle[other_particle].x);
            yij = static_cast<double>(part_y - Particle[other_particle].y);
            zij = static_cast<double>(part_z - Particle[other_particle].z);
            xij = periodic(xij);
            yij = periodic(yij);
            zij = periodic(zij);
            r_sq = static_cast<double>(xij*xij + yij*yij + zij*zij);
            if(r_sq < r_cut_sq)
            {
                r = sqrt(r_sq);
                r12 = static_cast<double>(1 / pow(r,12));
                r6 = static_cast<double>(1 / pow(r,6));
                U += r12 - r6;
            }
        }
    }
double pi = 3.141592;
double r9 = static_cast<double>(1 / pow(r_cut,9));
double r3 = static_cast<double>(1 / pow(r_cut,3)); 
U_tail = 8./9.*pi*rho*double(N)*(r9 - 3.0*r3); //tail correction
return (4*U + U_tail);
}
double periodic(double pos)
{
    double hbox = static_cast<double>(box)/2.0;
	
	if(pos<hbox) pos += box_length;
	if(pos>hbox) pos -= box_length;
	
	
return pos;
}

double randnum()
{
    return double(rand())/double(RAND_MAX);
}
double presssure()
{
    double
}
void write_config(int write_index)
{
    if(write_index == 1 || write_index == 2 || write_index == 3)
    {
        cout<<"inside write_config\n";
        char IntStr[80];
        ofstream of;
        sprintf( IntStr, "config%d.xyz",run_index);
        of.open (IntStr, ofstream::out | ofstream::trunc);
        if (of.is_open())
        {
            of << N<< "\n";
            of << "\n";
	    
            for (int i=0; i<N; i++)
            {
                of << 'h'<<"  "<<Particle[i].x << "  " << Particle[i].y << "  "<< Particle[i].z << "\n";
            }    
        }
    of.close();
    }
    else
    {
        char OutStr[80];
        ofstream out;
        sprintf( OutStr, "Run_config%d.xyz",run_index);
        out.open (OutStr, ofstream::out | ofstream::trunc);
        if(out.is_open())
        {
            out << N<< "\n";
            out << "\n";
	    
            for (int i=0; i<N; i++)
            {
                out << 'h'<<"  "<<Particle[i].x << "  " << Particle[i].y << "  "<< Particle[i].z << "\n";
            }    
        }
    out.close();
    }
    
}     