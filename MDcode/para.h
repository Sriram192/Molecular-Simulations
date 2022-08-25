using namespace std;
const int maxcolloids = 5000;//what is maxcolloids??
const int maxbin = 200;//what is maxbin?
const double rcut = 2.5;
class Cordinates{
public:
    double x, y, z;
}vector;

Cordinates colloid[maxcolloids], v[maxcolloids],f[maxcolloids],sumv,sumvn;
Cordinates dummy[maxcolloids];
void initialcondition();
double randnum(),RMS(double p),EAV(double q);
void printout(double), writeconf(),kinetic(), force(), simulation(), integrate(),rdfsample();
double box, rho, set_tem,t;
int N, start, ind, baseunit;
double hist[maxbin];//what is hist?
//const int step = 5000;
double U, press, dt=0.005, KE, Temp,E_msq,E_rms,Sigma_E;

