#include <iostream>
#include <string>
#include <cmath>
//#include <boost/serialization/array_wrapper.hpp>

#include <vector>
// Optimization
#define BOOST_UBLAS_NDEBUG
#include <fstream>
#include <exception>
#include <string>
#include <utility>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::math::tools;

// type definitions
typedef double value_type;// or typedef float value_type;
typedef boost::numeric::ublas::vector< value_type > state_type;
typedef boost::numeric::ublas::matrix< value_type > matrix_type;
typedef rosenbrock4< value_type > stepper_type;

// constants
const value_type pressure = 13.332237; //pascal soit 0.1 torr
const value_type Tg =0.02758 ;
const value_type L = 3e-2; //distance netre deux plaques en m
const value_type k_b = 1.38064852e-23; // constante de boltzman en J/K
const value_type K = 8.6173303e-5; //constqnte de boltzman en eV/k
const value_type D_Amet=0.005; //diffusion des Ar* en m2/s
const value_type pi = M_PI;
const value_type diff = pow((pi/L), 2);
const value_type DP = 1.58e25;//puissance totale du systeme
const value_type n_Ar =  (0.1/760)*2.69e25;
const value_type n_SiH4_ini = n_Ar/100.;
const value_type n_Arp_ini = 1.e16;
const int Nbr_espece=21;

const float C=1.35e21;

const int Nbr_K=47;
int jmax=Nbr_K;
int imax=9;
double **Tab;

value_type p1,p2,g1,g2,g3,g4,Tp,Tx,Tj;
state_type Kt(jmax, 0.0);


value_type k (int ind, value_type Tp) //K58 SiH3- + SiH3+ ->  Si2H4 + H2
{
    value_type K;
    K= Tab[6][ind]*pow(Tp,Tab[7][ind])*exp(-Tab[8][ind]/Tp);
    return K;
}
struct Condition
{
  value_type tol=1.e-9;
  bool operator() (value_type min, value_type max)  {
    return abs(min - max) <= tol;
  }
};



struct nsystem
{
  void operator()(const state_type &n, state_type &dndt, const value_type &t)
  {


    /*0=e, 1=Armet, 2=SiH3-, 3=SiH2-, 4=SiH3+, 5=SiH4, 6=SiH3,
    7=H, 8=SiH2, 9=H2, 10=H2+, 11=Si2H5, 12=Si2H2, 13=Si2H4-,
    14=Si2H6, 15=Si2H3-, 16=Si2H5-, 17=SiH-, 18=SiH, 19=Si, 20=Arp, 21=NP*/

for (int k=0;k<Nbr_espece;k++)
{
dndt[k]=0;
}


for (int j=0;j<jmax;j++)
{

 p1=Tab[0][j];
 p2=Tab[1][j];
 g1=Tab[2][j];
 g2=Tab[3][j];
 g3=Tab[4][j];
 g4=Tab[5][j]; 
 Tp=(p2==0 or g3==0)?Te:Tg;

 Kt[j]={Tab[6][j]*pow(Tp,Tab[7][j])*exp(-Tab[8][j]/Tp)};

if (p1==200)
{
Tx=n_Ar*n[p2]*Kt[j];
}
else if (p2==100)
{
Tx=n[p1]*Kt[j];
}
else
{
Tx=n[p1]*n[p2]*Kt[j];
}


if(p1!=200) {dndt[p1]=dndt[p1]-Tx;}
if(p2!=100) {dndt[p2]=dndt[p2]-Tx;}
if(g1!=200) {dndt[g1]=dndt[g1]+Tx;}
if(g2!=100) {dndt[g2]=dndt[g2]+Tx;}
if(g3!=100) {dndt[g3]=dndt[g3]+Tx;}
if(g4!=100) {dndt[g4]=dndt[g4]+Tx;}

dndt[5]=dndt[5]+C;

}
  }

  value_type Te;

};

struct jacobian
{
  void operator()(const state_type &n, matrix_type &jacobi,
                  const value_type &t, state_type &dfdt ) const
  {

for (int h=0;h<Nbr_espece;h++)

{
for (int p=0;p<Nbr_espece;p++)
{
jacobi(h,p)=0.0;
jacobi(h,p)=0.0;
jacobi(h,p)=0.0;
jacobi(h,p)=0.0;
jacobi(h,p)=0.0;
jacobi(h,p)=0.0;
}
}

for (int j=0;j<jmax;j++)
{

 p1=Tab[0][j];
 p2=Tab[1][j];
 g1=Tab[2][j];
 g2=Tab[3][j];
 g3=Tab[4][j];
 g4=Tab[5][j]; 
 Tp=(p2==0 or g3==0)?Te:Tg;

 Kt[j]={Tab[6][j]*pow(Tp,Tab[7][j])*exp(-Tab[8][j]/Tp)};
//cerr<<j<<endl;
for (int k=0;k<Nbr_espece;k++)
{

if (p1==200 and p2==k and p1!=p2)
	{Tj=n_Ar*Kt[j];
	jacobi(p2,k)=jacobi(p2,k)-Tj;
	if (g1!=200) {jacobi(g1,k)=jacobi(g1,k)+Tj;}
	if (g2!=100) {jacobi(g2,k)=jacobi(g2,k)+Tj;}
	if (g3!=100) {jacobi(g3,k)=jacobi(g3,k)+Tj;}
	if (g4!=100) {jacobi(g4,k)=jacobi(g4,k)+Tj;}
	}
if (p1!=200 and p2==k and p1!=p2) 
	{Tj=n[p1]*Kt[j];
	jacobi(p1,k)=jacobi(p1,k)-Tj;
	jacobi(p2,k)=jacobi(p2,k)-Tj;
	if (g1!=200) {jacobi(g1,k)=jacobi(g1,k)+Tj;}
	if (g2!=100) {jacobi(g2,k)=jacobi(g2,k)+Tj;}
	if (g3!=100) {jacobi(g3,k)=jacobi(g3,k)+Tj;}
	if (g4!=100) {jacobi(g4,k)=jacobi(g4,k)+Tj;}
	}
if (p2==100 and p1==k and p1!=p2) 
	{Tj=Kt[j];
	jacobi(p1,k)=jacobi(p1,k)-Tj;
	if (g1!=200) {jacobi(g1,k)=jacobi(g1,k)+Tj;}
	if (g2!=100) {jacobi(g2,k)=jacobi(g2,k)+Tj;}
	if (g3!=100) {jacobi(g3,k)=jacobi(g3,k)+Tj;}
	if (g4!=100) {jacobi(g4,k)=jacobi(g4,k)+Tj;}
	}
if (p1==k and p1!=p2 and p2!=100 ) 
	{Tj=n[p2]*Kt[j];
	jacobi(p1,k)=jacobi(p1,k)-Tj;
	jacobi(p2,k)=jacobi(p2,k)-Tj;
	if (g1!=200) {jacobi(g1,k)=jacobi(g1,k)+Tj;}
	if (g2!=100) {jacobi(g2,k)=jacobi(g2,k)+Tj;}
	if (g3!=100) {jacobi(g3,k)=jacobi(g3,k)+Tj;}
	if (g4!=100) {jacobi(g4,k)=jacobi(g4,k)+Tj;}
	}
if (p1==p2) 
	{Tj=2*n[p1]*Kt[j];
	if (p1!=200) {jacobi(p1,k)=jacobi(p1,k)-Tj;}
	if (p2!=100) {jacobi(p2,k)=jacobi(p2,k)-Tj;}
	if (g1!=200) {jacobi(g1,k)=jacobi(g1,k)+Tj;}
	if (g2!=100) {jacobi(g2,k)=jacobi(g2,k)+Tj;}
	if (g3!=100) {jacobi(g3,k)=jacobi(g3,k)+Tj;}
	if (g4!=100) {jacobi(g4,k)=jacobi(g4,k)+Tj;}

	}

}
}
    dfdt( 0 ) = 0.0;
    dfdt( 1 ) = 0.0;
    dfdt( 2 ) = 0.0;
    dfdt( 3 ) = 0.0;
    dfdt( 4 ) = 0.0;
    dfdt( 5 ) = 0.0;
    dfdt( 6 ) = 0.0;
    dfdt( 7 ) = 0.0;
    dfdt( 8 ) = 0.0;
    dfdt( 9 ) = 0.0;
    dfdt( 10 ) = 0.0;
    dfdt( 11 ) = 0.0;
    dfdt( 12 ) = 0.0;
    dfdt( 13 ) = 0.0;
    dfdt( 14 ) = 0.0;
    dfdt( 15 ) = 0.0;
    dfdt( 16 ) = 0.0;
    dfdt( 17 ) = 0.0;
    dfdt( 18 ) = 0.0;
    dfdt( 19 ) = 0.0;
    dfdt( 20 ) = 0.0;
  }
//stop modif
  value_type Te;
};


struct etemperature
{
  value_type operator()(value_type const& Te)
  
{

    return 
	-DP/n[0]
    +k(0,Te)*n_Ar*16.14 +k(1,Te)*n_Ar*12.31+ k(2,Te)*n[1]*5.39
    -k(3,Tg)*n[1]*n[1]*8.48
    -k(4,Te)*n[1]*12.31
    +k(5,Te)*n[5]*10.68 +k(6,Te)*n[5]*10.68 + k(7,Te)*n[5]*8.29 + k(8,Te)*n[5]*8.29
    +k(9,Te)*n[5]*24.1 + k(10,Te)*n[6]*1.94 + k(11,Te)*n[6]*1.30 + k(12,Te)*n[2]*1.16
    +k(13,Te)*n[3]*1.16 + k(14,Te)*n[8]*1.5*Te + k(15,Te)*n[9]*10.09
    +k(16,Te)*n[9]*16.05 +k(42,Te)*n[6]*1.5*Te +k(43,Te)*n[18]*1.5*Te
    +k(44,Te)*n[17]*1.25;
  }

  state_type n;
};

void write_density( const value_type t, const value_type Te, const state_type &n)
{
  cout << t  << '\t' <<Te <<'\t' << n[0] << '\t' << n[1] << '\t'
               << n[2] << '\t' << n[3] <<'\t'<< n[4] << '\t' << n[5] << '\t'
               << n[6] << '\t' << n[7] << '\t' << n[8] << '\t' << n[9] << '\t'
               << n[10] << '\t' << n[11] << '\t' << n[12] << '\t'
               << n[13] << '\t' << n[14] << '\t' << n[15] << '\t'
               << n[16] << '\t' << n[17] << '\t'<< n[18] << '\t'
               << n[19] << '\t' << n[20]  << '\t'<< n[13]+n[16]<<endl;
}

int main(int argc, char **argv)
{
  
ifstream fichier_k ("/home/cacot/Documents/CodeAutom/fichagarwal.dat");

Tab = new double*[imax];
Tab[0] = new double[imax*jmax];

	for(int i=1;i<imax;i++)
	{
	Tab[i]=Tab[i-1]+jmax;
	}

if(fichier_k)
{
    //Tout est prêt pour la lecture.
cerr<<"fichier ouvert"<<endl;


for(int j=0;j<jmax;j++)
{

           	fichier_k>>Tab[0][j]>>Tab[1][j]>>Tab[2][j]>>Tab[3][j]>>Tab[4][j]
		>>Tab[5][j]>>Tab[6][j]>>Tab[7][j]>>Tab[8][j];
       		
}
fichier_k.close();

}
else
{
    cerr << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
}

 //cerr << Tab[0][jmax-1]<<endl;



cout <<"t"<<'\t'<<"Te"<<'\t'<<"e"<<'\t'<<"Armet"<<'\t'<< "SiH3m"<<'\t'
               << "SiH2"<<'\t'<< "SiH3p"<<'\t'<< "SiH4"<<'\t'<< "SiH3"<<'\t'
               <<"H"<<'\t'<< "SiH2"<<'\t'<< "H2"<<'\t'<< "H2p"<<'\t'<< "Si2H5"
               <<'\t'<< "Si2H2"<<'\t'<<"Si2H4m"<<'\t'<<"Si2H6"<<'\t'<< "Si2H3m"
               <<'\t'<< "Si2H5m"<<'\t'<< "SiHm"<<'\t'<<"SiH"<<'\t'<< "Si"<<'\t'
               << "Arp"<<'\t'<<"NP"<<endl;
//clock_t t1,t2;

  // Time variables
  value_type t = 0.0;
  value_type dt = 1.0e-8;
  value_type Tmax = 20.e-3;
  value_type NT = Tmax/dt;

  // Root finding variables
  value_type min = 0.0001;
  value_type max = 20.0;
  boost::uintmax_t max_iter = 500;
  eps_tolerance<value_type> tol(30);

  // initial values
  value_type Te = 1.0;

  // Density vectors and initial condition
  state_type n_ini(Nbr_espece, 0.0); // initial conditions
  n_ini[0] = n_Arp_ini;
  n_ini[1] = n_Arp_ini;  // initial conditions
  n_ini[2] =  0.0;
  n_ini[3] = 0.0;
  n_ini[4] = 0.0;
  n_ini[5] = n_SiH4_ini;
  n_ini[6] = 0.0;
  n_ini[7] = 0.0;
  n_ini[8] = 0.0;
  n_ini[9] = 0.0;
  n_ini[10] = 0.0;
  n_ini[11] =0.0;
  n_ini[12] = 0.0;
  n_ini[13] = 0.0;
  n_ini[14] = 0.0;
  n_ini[15] = 0.0;
  n_ini[16] = 0.0;
  n_ini[17] = 0.0;
  n_ini[18] = 0.0;
  n_ini[19] = 0.0;
  n_ini[20] = n_Arp_ini;



  state_type n_new(Nbr_espece, 0.0);  // first step same as initial conditions
  n_new = n_ini;
  state_type n_err(Nbr_espece, 0.0); //error

  // declare the functor etemperature
  etemperature etemp;
  // assign initial values to functor etemp
  etemp.n = n_ini;

  //cerr << "\n[ii] Electrons  = " << etemp.n[0] << endl;
  //cerr << "\n[ii] Metastables  = " << etemp.n[1] << endl;
//   cout << "\n[ii] n  = " << etemp.n[2] << endl;

//t1=clock();
  // Find Te first calculation
  pair<value_type, value_type> pair_Te =\
                toms748_solve(etemp, min, max, tol, max_iter);

  Te = pair_Te.first;
  cerr << "\n[ii] Initial Temperature  = " << Te << endl;
//t2=clock()-t1;


//cout<<"timesec"<<(value_type )t2/CLOCKS_PER_SEC << endl;
/*cerr<<Tab[6][6]<<endl; 
cerr<<Tab[7][6]<<endl; 
cerr<<Tab[8][6]<<endl; 
value_type at={8.97e-9*pow(Te,-1)*exp(-10.68/Te)};
cerr<<at<<endl;
for (int j=0;j<jmax;j++)
{
 p2=Tab[1][j];
 g3=Tab[4][j];
 Tp=(p2==0 or g3==0)?Te:Tg;
 Kt[j]={Tab[6][j]*pow(Tp,Tab[7][j])*exp(-Tab[8][j]/Tp)};
cerr<<j<<'\t'<<Tp<<'\t'<<Kt[j]<<endl;
}

cerr<<k(44,Te)<<endl;
cerr<<k45(Te)<<endl;

for (int j=0;j<jmax;j++)
{p2=Tab[1][j];
 g3=Tab[4][j];
 Tp=(p2==0 or g3==0)?Te:Tg;
 Kt[j]={Tab[6][j]*pow(Tp,Tab[7][j])*exp(-Tab[8][j]/Tp)};
 cerr<<j<<'\t'<<Tp<<'\t'<<Kt[j]<<endl;
}*/

  // declare system and jacobian
  nsystem sys;
  jacobian jac;

  // declare stepper Rosenbrock
  stepper_type stepper;
cerr<<"patate1"<<endl;
  for (int i = 1; i <= NT; i++)
  {
    // update Te in system and jacobian
    sys.Te = Te;
    jac.Te = Te;
//cerr<<"patate"<<endl;
    // Integrate at least one step dt
    stepper.do_step( std::make_pair( sys, jac ), n_new, t, dt, n_err);
//cerr<<"patate2"<<endl;
    // assign values to functor etemp
    etemp.n = n_new;
    if (i%((int)(NT/50))==0)
    {
      write_density(t, Te, n_new);
    }
    // Find new Te
    pair<value_type, value_type> pair_Te =\
                  toms748_solve(etemp, min, max, tol, max_iter);
//cerr<<"patate3"<<endl;
    Te = pair_Te.first;

    t+= dt;
    n_ini = n_new;//update
  }

  value_type charge= (n_new[20]+n_new[4]+n_new[10]-n_new[0]-n_new[2]-n_new[3]-n_new[13]-n_new[15]-n_new[16]-n_new[17])/n_Arp_ini;

  cerr<<"charge/dArp="<<charge<<endl;

  value_type Si=(n_new[2]+n_new[3]+n_new[4]+n_new[13]*2+2*n_new[15]+n_new[16]*2
          +n_new[17]+n_new[5]+n_new[6]+n_new[8]+n_new[18]+2*n_new[11]+n_new[19]
          +n_new[12]*2+n_new[14]*2)/n_SiH4_ini;

  cerr<<"Si="<<Si<<endl;


  value_type H=(3*n_new[2]+2*n_new[3]+3*n_new[4]+2*n_new[10]+4*n_new[13]+3*n_new[15]
         +5*n_new[16]+n_new[17]+4*n_new[5]+3*n_new[6]+n_new[7]+2*n_new[8]
         +2*n_new[9]+n_new[18]+5* n_new[11]+2*n_new[12]+6*n_new[14])
            /(4*n_SiH4_ini);

  cerr<<"H="<<H<<endl;
   
 // Libération de la mémoire
    delete [] Tab[0];
    delete [] Tab;


  return 0;

}

