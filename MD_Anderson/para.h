using namespace std;
const int maxcolloids = 2000;
const int maxbin = 500;
//const double rcut = 2.5;
class Cordinates
{
public:
    double x, y, z;
}vector;

Cordinates colloid[maxcolloids], v[maxcolloids],f[maxcolloids],sumv,sumvn;
void initialcondition();
double randnum(), gauss();
void printout(double), writeconf(),kinetic(), force(), simulation(), integrate(),rdfsample();
double box, rho, set_tem;
int N, start, ind, baseunit;
double hist[maxbin];
double U, press, dt=0.005, KE, Temp;

//Last Update on 12-Mar-2015
void copycoordinates(Cordinates *a,Cordinates *b, int size = N);
int nint(double d);
void writecoord(Cordinates *input, char *fname);
void msdsample(int delt, int stime ,int etime);

void anderson(double nu=2);

void absconfig(char* outputfilename,int delt = 10,int stime = 1,int etime = 5000);

std::ifstream& readconf(std::ifstream& file,bool abs ,unsigned int time, Cordinates *r);
//


