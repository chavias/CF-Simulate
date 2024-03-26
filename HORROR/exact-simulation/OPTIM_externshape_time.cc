/*

simulates phase transients
Can run up to 10 spins without the CSA tensors. 
Brute force integration of the MAS rotation (one cycle)
cosine modulated amplitude transient

*/

#include "gamma.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>      //added
#define NPROP 800000
#define MAXSPINS 10

//extern "C" int gethostname(char *,int);      // commeted out

using namespace std;

int main(int argc, char *argv[])

{
  spin_system ax;
  gen_op Ua,Ham, U, H[5], sigma1, sigmax,sigmay,sigmaz, detect[6];
  spin_T  Hdip[MAXSPINS][MAXSPINS];
  space_T Adip[MAXSPINS][MAXSPINS], Adip_R[MAXSPINS][MAXSPINS];
  space_T Acsa[MAXSPINS], Acsa_R[MAXSPINS];
  double D[MAXSPINS][MAXSPINS];
  double J[MAXSPINS][MAXSPINS];
  double iso_CSA[MAXSPINS];
  double eta_CSA[MAXSPINS];
  double delta_CSA[MAXSPINS];
  double phase[NPROP];
  double ampli[NPROP];
  int i,j,k ,k1,k2,count,qu,steps, n_sample;
  string name, names, nameB1;
  const double thetam=54.73561032;
  double mas_freq;
  double time, ltime, scale;
  double max_ampli;
  int nstart, nspins, n_max, n_file;
  double alpha,beta,gamma;
  double alpha_CSA[MAXSPINS],beta_CSA[MAXSPINS];
  double gamma_CSA[MAXSPINS];
  double alpha_D[MAXSPINS][MAXSPINS],beta_D[MAXSPINS][MAXSPINS];
  double gamma_D[MAXSPINS][MAXSPINS];
  struct rusage me;
  char hostname[200];
  fstream fpout;
  fstream fpin;

  int dooutput = 1;


  int value1[] = {1, 50, 100, 144, 200, 300, 538, 1154, 3000, 5000, 7000,10000};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107,  637, 1197, 1083, 1759};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271,  933, 1715, 1787, 3763};

  count = 1;
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
//setup for the spin system
  ax.read(names);
  nspins = ax.spins();
  if(nspins < 2 || nspins >= 10)
  { cerr << "This program is written for a two to nine spin system.\n";
    cerr << "Please change your spin system definition in \n";
    cerr << "the file " << names << ".\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { query_parameter(argc,argv,count++,"Scalar Coupling Constant    ? ", J[i][j]);
      query_parameter(argc,argv,count++,"Dipolar Coupling Constant   ? ", D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_D[i][j]);
    }
  }
  for(i=0;i<nspins;++i)
  { query_parameter(argc,argv,count++,"Isotropic chemical shift    ? ", iso_CSA[i]);
    query_parameter(argc,argv,count++,"anisotropy chemical shift   ? ", delta_CSA[i]);
    query_parameter(argc,argv,count++,"asymmetry chemical shift    ? ", eta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_CSA[i]);
  }
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"number of points            ? ", steps);
  query_parameter(argc,argv,count++,"Filename for shape (-1..1)  ? ", nameB1);
  query_parameter(argc,argv,count++,"number of total data points ? ", n_max);
  query_parameter(argc,argv,count++,"number of points per rotor  ? ", n_sample);
  query_parameter(argc,argv,count++,"MAS frequency		 ? ", mas_freq);
  query_parameter(argc,argv,count++,"Initial density operator?   ? ", nstart);
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);
  query_parameter(argc,argv,count++,"number of points in file    ? ", n_file);
  query_parameter(argc,argv,count++,"max amplitude                   ? ", max_ampli);


  if(nstart >= nspins)
  { cerr << "The starting density operator has a spin number which is\n";
    cerr << "larger than the number of spins.\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }

  // ampli is relative to max_ampli
  // the value of ampli has to be element of the interval [0,1]
  fpin.open(nameB1.c_str(),ios::in);
  if(fpin.is_open())
  { for(k=0;k<n_file;++k)
     fpin >> ampli[k] >> phase[k];
    fpin.close();
  }
  else
  { cerr << "File " << nameB1 << " does not exist\n";
    cerr << "Aborting ...\n\n";
  }

  if(dooutput==1)
  {
  string name3 = name+"_in.dat";
  fpout.open(name3.c_str(),ios::out);
  for(k=0;k<n_max;++k)
  { fpout << max_ampli*ampli[k%n_file] << "  " << phase[k%n_file] << "\n";
    phase[k]=phase[k]/180.0*PI;
  }
  fpout.flush();
  fpout.close();
  }

// #############################################################################################################



  gethostname(hostname,199);

  cout << "\n\nSimulation of arbitrary pulse Sequence with amplitude shaped pulses  \n";
  cout << "==========================================================\n\n";
  cout << "Program version: " << __FILE__ << " compiled at " << __DATE__ ", "
       << __TIME__ << "\n\n";
  cout << "running on machine " << hostname << "\n\n";
  cout << "Parameters:\n";
  cout << "rotation angle thetam          : " << thetam << " Degree\n";
  cout << "size of spin system            : " << nspins << " spins\n";
  cout << "initial density operator       : Iz(ax," << nstart << ")\n"; 
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { cout << "scalar  coupling constant (" << i << "," << j << ") : " << 
              J[i][j] << " Hz\n";
      cout << "dipolar coupling constant (" << i << "," << j << ") : " << 
              D[i][j] << " Hz\n";
      cout << "relativ orientation of D tensor: (" << alpha_D[i][j] << "," <<
               beta_D[i][j] << "," << gamma_D[i][j] << ")\n";
    }
  }
  for(i=0;i<nspins;++i)
  { cout << "isotropic chemical shift (" << i << ") : " << iso_CSA[i] << " Hz\n";
    cout << "anisotropy chemical shift (" << i << "): " << delta_CSA[i] << " Hz\n";
    cout << "asymmetry chemical shift (" << i << ") : " << eta_CSA[i] << "\n";
    cout << "relativ orientation of CSA tensor: (" << alpha_CSA[i] << "," <<
             beta_CSA[i] << "," << gamma_CSA[i] << ")\n";
  }
  cout << "Powder Quality Number:         " << qu << "  (" << value1[qu] <<
          " orientations)\n";
  cout << "MAS frequency:                 " << (mas_freq) << " Hz\n";
  cout << "# of data points:              " << n_max << "\n";
  cout << "# of data points/rotor:        " << steps << "\n";
  cout << "time resolution of calculation:" << 1.0/mas_freq/steps << " s\n";
  cout << "time resolution of output data:" << 1.0/mas_freq/steps*n_sample << " s\n";
  cout << "# of data points/B1 point:     " << n_sample << "\n";
  cout << "Output filename:               " << name << "\n";
  cout << "Max Amplitude:                 " << max_ampli << "\n";
  cout << "\n";
  cout.flush();

  time = 1.0/mas_freq/steps;
  block_2D data(18,n_max);
  for(k1=0;k1<18;++k1)
  for(k2=0;k2<n_max;++k2)
  { data(k1,k2) = 0;
  }

//setup for the hamiltonian
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { Hdip[i][j] = T_D(ax,i,j);
    }
  }


//setup for the space tensor
  matrix help(3,3,0);
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { help.put_h(-1.0/2.0,0,0);
      help.put_h(-1.0/2.0,1,1);
      help.put_h( 1.0,2,2);
      help   = - (complex) D[i][j] * help;
      Adip[i][j] = A2(help);
      Adip[i][j] = Adip[i][j].rotate(alpha_D[i][j],beta_D[i][j],gamma_D[i][j]);
    }
  }

  for(i=0;i<nspins;++i)
  { help.put(-1.0/2.0*(1.0+eta_CSA[i]),0,0);
    help.put(-1.0/2.0*(1.0-eta_CSA[i]),1,1);
    help.put( 1.0,2,2);
    help = (complex) delta_CSA[i] * help;
    Acsa[i] = A2(help);
    Acsa[i] = Acsa[i].rotate(alpha_CSA[i],beta_CSA[i],gamma_CSA[i]);
  }

  string name1 = name+".mat";
  string name2 = "result";

//here starts the powder loop
//reference JCP 59 (8), 3992 (1973).

  for(count=1; count<=value1[qu]; ++count)
  { beta  = 180.0 * count/value1[qu];
    alpha = 360.0 * ((value2[qu]*count) % value1[qu])/value1[qu];
    gamma = 360.0 * ((value3[qu]*count) % value1[qu])/value1[qu];
    if(count % 10 == 1)
    { getrusage(0, & me);
      cout << count << "\tbeta = " << beta << "\talpha = "
           << alpha << "\tgamma = " << gamma 
           << ",\ttime used: " << me.ru_utime.tv_sec << " seconds\n";
      cout.flush();
    }
    scale = sin(beta/180.0*PI);

    sigmax  = Ix(ax,nstart);
    sigmay  = Iy(ax,nstart);
    sigmaz  = Iz(ax,nstart);
    detect[0] = Ix(ax,0);
    detect[1] = Iy(ax,0);
    detect[2] = Iz(ax,0);
    detect[3] = Ix(ax,1);
    detect[4] = Iy(ax,1);
    detect[5] = Iz(ax,1);
  
//now we rotate the space tensor
    for(i=0;i<nspins;++i)
    { for(j=i+1;j<nspins;++j)
      { Adip_R[i][j] = Adip[i][j].rotate(alpha,beta,gamma);
      }
    }
    for(i=0;i<nspins;++i)
    { Acsa_R[i] = Acsa[i].rotate(alpha,beta,gamma);
    }

//zero all components
    for(i=0;i<5;++i)
      H[i] = gen_op();

//this is the dipolar part
    for(i=0;i<nspins;++i)
    { H[2] += iso_CSA[i]*Iz(ax,i);
    }


    for(k=-2;k<=2;++k)
    { for(i=0;i<nspins;++i)
      { H[k+2] += Acsa_R[i].component(2,k) * d2(k,0,thetam)*2.0/sqrt(6.0)*Iz(ax,i);
        for(j=i+1;j<nspins;++j)
        { if(ax.isotope(i) != ax.isotope(j))
          { H[k+2] += Adip_R[i][j].component(2,k) * d2(k,0,thetam) * 1.0/sqrt(6.0)*2*Iz(ax,i)*Iz(ax,j);
	    if(k==0)
	      H[k+2] += J[i][j] * Iz(ax,i)*Iz(ax,j);
          }
          else
          { H[k+2] += Adip_R[i][j].component(2,k) * d2(k,0,thetam) * Hdip[i][j].component(2,0);
	    if(k==0)
	      H[k+2] += J[i][j] * (Iz(ax,i)*Iz(ax,j)+Ix(ax,i)*Ix(ax,j)+Iy(ax,i)*Iy(ax,j));
          }
        }
      }
    }

     U = Ie(ax,0);

//now we calculate the propagator for the evolution

    detect[0].set_DBR();
    detect[1].set_DBR();

  for(k=0;k<n_max*n_sample;++k) // brute force over Fnp times R26 block
    { ltime=(k+0.5)*time;
      Ham = max_ampli * ampli[(k/(n_sample))%n_file]*( Fx(ax,"13C")*cos(phase[(k/(n_sample))%n_file])+Fy(ax,"13C")*sin(phase[(k/(n_sample))%n_file]) );
      for(i=-2;i<=2;++i)
        Ham += exp(complex(0,i*2.0*PI*ltime*(mas_freq))) * H[i+2];
//    cout << k << "  " << k/(n_sample) << "\n";
//    Ham.set_DBR();
//    cout << Ham << "\n";
      if(k%(n_sample)==0)
      { 
	i=k/(n_sample);
	sigma1=evolve(sigmax,U);	
        data(0,i) += proj(sigma1,detect[0])*scale;
        data(1,i) += proj(sigma1,detect[1])*scale;
        data(2,i) += proj(sigma1,detect[2])*scale;
        data(3,i) += proj(sigma1,detect[3])*scale; 
        data(4,i) += proj(sigma1,detect[4])*scale; 
        data(5,i) += proj(sigma1,detect[5])*scale;
	sigma1=evolve(sigmay,U);	
        data(6+0,i) += proj(sigma1,detect[0])*scale;
        data(6+1,i) += proj(sigma1,detect[1])*scale;
        data(6+2,i) += proj(sigma1,detect[2])*scale;
        data(6+3,i) += proj(sigma1,detect[3])*scale; 
        data(6+4,i) += proj(sigma1,detect[4])*scale; 
        data(6+5,i) += proj(sigma1,detect[5])*scale; 
	sigma1=evolve(sigmaz,U);	
        data(12+0,i) += proj(sigma1,detect[0])*scale;
        data(12+1,i) += proj(sigma1,detect[1])*scale;
        data(12+2,i) += proj(sigma1,detect[2])*scale;
        data(12+3,i) += proj(sigma1,detect[3])*scale; 
        data(12+4,i) += proj(sigma1,detect[4])*scale; 
        data(12+5,i) += proj(sigma1,detect[5])*scale; 
      }	
      U &= prop(Ham,time);
    }

  } // end of powder loop
  MATLAB(name1,name2,data,1);
  exit(0);
}
