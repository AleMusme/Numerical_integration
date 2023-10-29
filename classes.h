#ifndef __classes__
#define __classes__

#include "random.h"
#include "functions.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>

using namespace std;

class f{

    public:

        f(){}

        double Eval(double rtilde, double k){
            return boost::math::cyl_bessel_j(1, k*rtilde);
        }

};


class distro{

    public: 

        distro(){}

        distro(double r0, int t){

            m_r0 = r0;
            m_alpha = 1./r0;

            m_thickness = t;

            if(m_thickness==0){
                //m_zg = 0.196 * pow(r0, 0.633);
                m_zg = 1.2 * pow(r0, 0.633);
            }
            if(m_thickness==1){
                m_zd = 0.196 * pow(r0, 0.633);
            }
            if(m_thickness==2){
                m_zh = 0.196 * pow(r0, 0.633);
            }
        }

        double Eval(double k){

            if(m_thickness == 0){
                return (1. - exp(- m_zg * m_alpha * k)) * 1./pow( 1. + (k*k) , 3./2.);
            }
            if(m_thickness == 1){
                return k/(1 + m_zd*m_alpha*k) * 1./pow( 1. + (k*k)  , 3./2.);
            }
            if(m_thickness == 2){
                return (boost::math::beta(fabs(m_zh*m_alpha*k/2.), 0.0000000001) - 1.) * k/pow( 1. + (k*k) , 3./2.);
                //return (0.5 * (boost::math::digamma((m_zg*m_alpha*k/2. + 1.)*0.5) - boost::math::digamma((m_zg*m_alpha*k/2.) * 0.5)) - 1.) * k/pow( 1. + (k*k) , 3./2.);
            }
            else return 0;

        }

    private:
        int m_thickness;
        double m_r0;
        double m_alpha;
        double m_zd;
        double m_zg;
        double m_zh;

}; 

class Metropolis{

    public:

        Metropolis(){}

        Metropolis(distro& dis, Random& r, double x0 = 0., double d = 1.) :
            xn{x0}, xold{x0}, delta{d}, function{dis}, rnd{r} {

            accepted = 0;
            attempted = 0;
        }

        ~Metropolis(){}

        void SamplingStep(){

            xn = xold + delta*(rnd.Rannyu() - 0.5);

            double A = min(1. , function.Eval(xn)/function.Eval(xold));

            double r = rnd.Rannyu();

            if( r < A and xn > 0. ){

                xold = xn;
                accepted++;
            }
            attempted++;
        }

    public:
        int accepted, attempted;
        double xn, xold;
        double delta;

        distro function;
        Random rnd;

};

class BlockAverage{

    public:

        BlockAverage(f& fun, int n = 1, int b = 20):
            n_props{n}, nblk{b}, blk_av(n_props), glob_av(n_props), glob_av2(n_props), walker(n_props), function{fun} {
                blk_norm = 0.;
                stima = 0.;
                errore = 0.;
            }

        void Reset(int iblk){ //Reset block averages
    
            if(iblk == 1){

                for(int i=0; i<n_props; ++i)
                {
                    glob_av[i] = 0;
                    glob_av2[i] = 0;
                }
            }

            for(int i=0; i<n_props; ++i){

                blk_av[i] = 0;
            }

            blk_norm = 0;
        }

        void Measure(double r, double k){

            walker[iexp] = function.Eval(r, k);

        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk){ //Print results for current block
        
            const int wd=18;
            
            stima = blk_av[iexp]/blk_norm; 
            glob_av[iexp] += stima;
            glob_av2[iexp] += stima*stima;
            errore = error(glob_av[iexp],glob_av2[iexp],iblk);

        }


    public:
    
        int n_props, nblk;
        const int iexp=0;

        vector<double>  blk_av, glob_av, glob_av2;
        vector<double> walker;
        double stima;
        double errore;
        double blk_norm;

        f function;

};


#endif 