/*
 * File:   main.cpp
 * Author: Abuenameh
 *
 * Created on August 6, 2014, 11:21 PM
 */

#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <nlopt.h>
#include <complex>
#include <iostream>
#include <queue>
//#include <thread>
#include <nlopt.hpp>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time.hpp>
#include <boost/random.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/progress.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include <pagmo/src/pagmo.h>

//#include <bayesopt/bayesopt.hpp>

#include "mathematica.hpp"
#include "gutzwiller.hpp"

using namespace std;

//using boost::lexical_cast;
using namespace boost;
using namespace boost::random;
using namespace boost::filesystem;
using namespace boost::posix_time;

using namespace pagmo;
using namespace pagmo::problem;
using namespace pagmo::algorithm;

typedef boost::array<double, L> Parameter;

//template<typename T> void printMath(ostream& out, string name, T& t) {
//    out << name << "=" << ::math(t) << ";" << endl;
//}
//
//template<typename T> void printMath(ostream& out, string name, int i, T& t) {
//    out << name << "[" << i << "]" << "=" << ::math(t) << ";" << endl;
//}

double M = 1000;
double g13 = 2.5e9;
double g24 = 2.5e9;
double delta = 1.0e12;
double Delta = -2.0e10;
double alpha = 1.1e7;

double Ng = sqrt(M) * g13;

double JW(double W) {
    return alpha * (W * W) / (Ng * Ng + W * W);
}

double JWij(double Wi, double Wj) {
    return alpha * (Wi * Wj) / (sqrt(Ng * Ng + Wi * Wi) * sqrt(Ng * Ng + Wj * Wj));
}

Parameter JW(Parameter W) {
    Parameter v;
    for (int i = 0; i < L; i++) {
        v[i] = W[i] / sqrt(Ng * Ng + W[i] * W[i]);
    }
    Parameter J;
    for (int i = 0; i < L - 1; i++) {
        J[i] = alpha * v[i] * v[i + 1];
    }
    J[L - 1] = alpha * v[L - 1] * v[0];
    return J;
}

double UW(double W) {
    return -2 * (g24 * g24) / Delta * (Ng * Ng * W * W) / ((Ng * Ng + W * W) * (Ng * Ng + W * W));
}

Parameter UW(Parameter W) {
    Parameter U;
    for (int i = 0; i < L; i++) {
        U[i] = -2 * (g24 * g24) / Delta * (Ng * Ng * W[i] * W[i]) / ((Ng * Ng + W[i] * W[i]) * (Ng * Ng + W[i] * W[i]));
    }
    return U;
}

boost::mutex progress_mutex;
boost::mutex points_mutex;

struct Point {
    double x;
    double mu;
};

void norm2s(unsigned m, double *result, unsigned ndim, const double* x,
        double* grad, void* data);

double norm(const vector<double> x, vector<double>& norms) {
    const doublecomplex * f[L];
    for (int i = 0; i < L; i++) {
        f[i] = reinterpret_cast<const doublecomplex*> (&x[2 * i * idim]);
    }

    norms.resize(L);

    //    double norm = 1;
    for (int i = 0; i < L; i++) {
        double normi = 0;
        for (int n = 0; n <= nmax; n++) {
            normi += norm(f[i][n]);
        }
        //        norm *= normi;
        norms[i] = sqrt(normi);
    }
    //    return norm;
    return 0;
}

int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return a > b ? a : b;
}

struct PointResults {
    double W;
    double mu;
    double E0;
    double Eth;
    double E2th;
    double fs;
    vector<double> J;
    vector<double> U;
    double fmin;
    vector<double> fn0;
    vector<double> fmax;
    vector<double> f0;
    vector<double> fth;
    vector<double> f2th;
};

void phasepoints(Parameter& xi, phase_parameters pparms, queue<Point>& points, vector<PointResults>& pres, progress_display& progress) {

    mt19937 rng;
    rng.seed(time(NULL));
    uniform_real_distribution<> uni(-1, 1);

    int ndim = 2 * L * idim;

    vector<double> x(ndim);
    doublecomplex * f[L];
    for (int i = 0; i < L; i++) {
        f[i] = reinterpret_cast<doublecomplex*> (&x[2 * i * idim]);
    }

    vector<double> U(L), J(L), dU(L);

    vector<double> x0(ndim);
    doublecomplex * f0[L];
    for (int i = 0; i < L; i++) {
        f0[i] = reinterpret_cast<doublecomplex*> (&x0[2 * i * idim]);
    }

    vector<double> xabs(ndim / 2);
    double* fabs[L];
    for (int i = 0; i < L; i++) {
        fabs[i] = &xabs[i * idim];
    }

    vector<double> fn0(L);
    vector<double> fmax(L);

    vector<double> norms(L);

    funcdata data;
    data.canonical = pparms.canonical;
    //        parms.theta = theta;

    double theta = pparms.theta;

    funcdata2 fdata = {0, dU, U, J, 0, 0, vector<double>(6)};

    double scale = 1;

    for (;;) {
        Point point;
        {
            boost::mutex::scoped_lock lock(points_mutex);
            if (points.empty()) {
                break;
            }
            point = points.front();
            points.pop();
        }
        //        cout << "Got queued" << endl;

        PointResults pointRes;
        pointRes.W = point.x;
        pointRes.mu = point.mu;

        double W[L];
        for (int i = 0; i < L; i++) {
            W[i] = xi[i] * point.x;
        }
        double U0 = 1 / scale;
        fdata.U0 = U0;
        for (int i = 0; i < L; i++) {
            U[i] = UW(W[i]) / UW(point.x) / scale;
            dU[i] = U[i] - U0;
            J[i] = JWij(W[i], W[mod(i + 1)]) / UW(point.x) / scale;
        }
        pointRes.J = J;
        pointRes.U = U;

        nlopt::opt localopt(nlopt::LD_LBFGS, ndim);
        //        nlopt::opt localopt(nlopt::LD_MMA, ndim);
        //        nlopt::opt localopt(nlopt::LN_NELDERMEAD, ndim);
        localopt.set_lower_bounds(-1);
        localopt.set_upper_bounds(1);
        localopt.set_min_objective(energyfunc, &fdata);
        //        localopt.set_ftol_rel(1e-14);
        //        localopt.set_xtol_rel(1e-14);
        //        localopt.set_ftol_abs(1e-12);
        //        localopt.set_xtol_abs(1e-12);
        fdata.mu = point.mu / scale;

        vector<double> popx;

        de_1220 algo(2000);
        gsl_bfgs2 lalgo(100, 1e-8, 1e-8, 0.01, 1e-10);

        int npop = 20;

        ////        /*double*/ U0 = 0.7786907302065487;
        //        U0 = 0.8786907302065487;
        //        double Ua[] = {0.922213971320164, 0.856881508533313, 0.700324261083433, \
//0.571208141081931, 0.842825769013902};
        //        double dUa[] = {0.14352324111361536, 0.07819077832676435, -0.07836646912311562, \
//-0.2074825891246177, 0.06413503880735338};
        //////        double Ja[] = {0.07683663838672976, 0.07283001881815312, 0.07694395816589483, \
//////0.07347544296692873, 0.0761660621229901};
        //        double Ja[] = {0.07589911596407645, 0.07130662982267995, 0.07670884040551718, \
//0.07599776398921115, 0.0705576727237478};
        //////        double Ja[] = {0.07283001881815312, 0.07694395816589483, \
//////0.07347544296692873, 0.0761660621229901,0.07683663838672976};
        //        vector<double> dU(5);
        //        for(int i = 0; i < 5; i++) {
        //            U[i] = Ua[i];
        ////            dU[i] = U[i]-U0;//dUa[i];
        //            dU[i] = dUa[i];
        //            J[i] = Ja[i];
        //        }
        ////        
        //        energyprob prob(ndim, U0, dU, U, J, point.mu/scale, 0.1);
        //        fitness_vector fv(1);
        //        decision_vector dv(ndim);
        //        double dva[] = {0.12806249101263278, 0, 0.054426358315506206, 0, \
//0.5017203911948146, 0, 0.7032020484792338, 0, 0.4282632265309182, \
//0, 0.22585815635351836, 0, 0.41466066310273253, 0, \
//0.4858208770649655, 0, 0.3909971375850023, 0, 0.38811392766222463, \
//0, 0.5311079401305173, 0, 0.08029809068891339, 0, \
//0.5150540287549862, 0, 0.052190588891518695, 0, \
//0.4211423332437574, 0, 0.6842349703439841, 0, 0.23557373654220737, \
//0, 0.1759606349565957, 0, 0.5699253451548816, 0, \
//0.4367825914524495, 0, 0.33798191844475106, 0, \
//0.47171561640183807, 0, 0.38358047886589935, 0, \
//0.022906004326788187, 0, 0.5226675531791949, 0, \
//0.5421896797766332, 0, 0.5693170132769998, 0, 0.27577384468187405, \
//0, 0.012076767067972723, 0, 0.18036090673232205, 0};
        //        double dva2[2*L*dim];
        //        for (int i = 0; i < ndim; i++) {
        //            dv[i] = dva[i];
        //            dva2[i] = dva[i];
        //        }
        //        int nloops = 1;
        //        time_t start1 = time(NULL);
        //        for(int i = 0; i < nloops; i++) {
        //        prob.objfun_impl(fv, dv);
        //        }
        //        time_t end1 = time(NULL);
        //        cout << ::math(fv[0]) << endl;
        //        cout << (end1 - start1) << endl;
        //        double en, en2;
        //        time_t start2 = time(NULL);
        //        vector<double> grad(2*L*dim);
        //        double dx = 1e-6;
        //        int di = 12;
        //        dva2[di] += dx;
        //        for (int i = 0; i < nloops; i++) {
        //        en = energy(dva, Ua, U0, dUa, Ja, point.mu/scale, 0.1);
        //        en2 = energy(dva2, Ua, U0, dUa, Ja, point.mu/scale, 0.1);
        //        energygrad(dva, Ua, U0, dUa, Ja, point.mu/scale, 0.1, grad.data());
        //        }
        //        double den = (en2 - en)/dx;
        //        time_t end2 = time(NULL);
        //        cout << ::math(en) << endl;
        //        cout << "grad " << ::math(grad) << endl;
        //        cout << "den " << ::math(den) << endl;
        //        cout << (end2 - start2) << endl;
        //        exit(0);
        //        /*double*/ U0 = 2;//0.786907302065487;
        //                   fdata.U0 = U0;
        energyprob prob0(ndim, U0, dU, U, J, point.mu / scale, 0);
        population pop0(prob0, npop);
        //        algo.evolve(pop0);
        //        cout << "E0 PaGMO: " << pop0.champion().f << endl;

        double E0 = pop0.champion().f[0];
        popx = pop0.champion().x;
        fdata.theta = 0;
        try {
            localopt.optimize(popx, E0);
        } catch (std::exception& e) {
            printf("nlopt failed!: E0 refine: %f, %f\n", point.x, point.mu);
            cout << e.what() << endl;
            E0 = pop0.champion().f[0];
        }
        //        cout << ::math(fdata.Ei) << endl;
        //        cout << "E0 nlopt: " << ::math(E0) << endl;
        //        exit(0);

        norm(popx, norms);
        for (int i = 0; i < L; i++) {
            for (int n = 0; n <= nmax; n++) {
                x0[2 * (i * idim + n)] = popx[2 * (i * idim + n)] / norms[i];
                x0[2 * (i * idim + n) + 1] = popx[2 * (i * idim + n) + 1] / norms[i];
            }
            transform(f0[i], f0[i] + idim, fabs[i], std::ptr_fun<const doublecomplex&, double>(abs));
            fmax[i] = *max_element(fabs[i], fabs[i] + idim);
            fn0[i] = fabs[i][1];
        }

        pointRes.fmin = *min_element(fn0.begin(), fn0.end());
        pointRes.fn0 = fn0;
        pointRes.fmax = fmax;
        pointRes.f0 = x0;
        pointRes.E0 = E0;

        energyprob probth(ndim, U0, dU, U, J, point.mu / scale, theta);
        //        energy probth(ndim, 0, U, U, J, point.mu/scale, theta);
        population popth(probth, npop);
        //        algo.evolve(popth);
        //        cout << "Eth PaGMO: " << popth.champion().f << endl;
        //        cout << popth.champion().f << endl;

        double Eth = popth.champion().f[0];
        popx = popth.champion().x;
        fdata.theta = theta;
        try {
            localopt.optimize(popx, Eth);
        } catch (std::exception& e) {
            printf("nlopt failed!: Eth refine: %f, %f\n", point.x, point.mu);
            cout << e.what() << endl;
            Eth = popth.champion().f[0];
        }
        //        cout << ::math(fdata.Ei) << endl;
        //        cout << "Eth nlopt: " << ::math(Eth) << endl;

        //        lalgo.evolve(popth);
        //        cout << popth.champion().f << endl;
        //        cout << popth.champion().x << endl;

        //        double Eth = popth.champion().f[0];

        //        popx = popth.champion().x;
        norm(popx, norms);
        for (int i = 0; i < L; i++) {
            for (int n = 0; n <= nmax; n++) {
                x0[2 * (i * idim + n)] = popx[2 * (i * idim + n)] / norms[i];
                x0[2 * (i * idim + n) + 1] = popx[2 * (i * idim + n) + 1] / norms[i];
            }
        }

        pointRes.fth = x0;
        pointRes.Eth = Eth;

        energyprob prob2th(ndim, U0, dU, U, J, point.mu / scale, 2 * theta);
        population pop2th(prob2th, npop);
        //        algo.evolve(pop2th);
        //        cout << "E2th PaGMO: " << pop2th.champion().f << endl;
        //        cout << pop2th.champion().f << endl;

        double E2th = pop2th.champion().f[0];
        popx = pop2th.champion().x;
        fdata.theta = 2 * theta;
        try {
            localopt.optimize(popx, E2th);
        } catch (std::exception& e) {
            printf("nlopt failed!: E2th refine: %f, %f\n", point.x, point.mu);
            cout << e.what() << endl;
            E2th = pop2th.champion().f[0];
        }
        //        cout << ::math(fdata.Ei) << endl;
        //        cout << "E2th nlopt: " << ::math(E2th) << endl;

        //        lalgo.evolve(pop2th);
        //        cout << pop2th.champion().f << endl;
        //        cout << pop2th.champion().x << endl;

        //        double E2th = pop2th.champion().f[0];

        //        popx = pop2th.champion().x;
        norm(popx, norms);
        for (int i = 0; i < L; i++) {
            for (int n = 0; n <= nmax; n++) {
                x0[2 * (i * idim + n)] = popx[2 * (i * idim + n)] / norms[i];
                x0[2 * (i * idim + n) + 1] = popx[2 * (i * idim + n) + 1] / norms[i];
            }
        }

        pointRes.f2th = x0;
        pointRes.E2th = E2th;

        pointRes.fs = (E2th - 2 * Eth + E0) / (L * theta * theta);

        {
            boost::mutex::scoped_lock lock(points_mutex);
            pres.push_back(pointRes);
        }

        {
            boost::mutex::scoped_lock lock(progress_mutex);
            ++progress;
        }
    }

}

vector<double> nu;

/*
 *
 */
int main(int argc, char** argv) {

    mt19937 rng;
    uniform_real_distribution<> uni(-1, 1);

    int seed = lexical_cast<int>(argv[1]);
    int nseed = lexical_cast<int>(argv[2]);

    double xmin = lexical_cast<double>(argv[3]);
    double xmax = lexical_cast<double>(argv[4]);
    int nx = lexical_cast<int>(argv[5]);

    deque<double> x(nx);
    if (nx == 1) {
        x[0] = xmin;
    } else {
        double dx = (xmax - xmin) / (nx - 1);
        for (int ix = 0; ix < nx; ix++) {
            x[ix] = xmin + ix * dx;
        }
    }

    double mumin = lexical_cast<double>(argv[6]);
    double mumax = lexical_cast<double>(argv[7]);
    int nmu = lexical_cast<int>(argv[8]);

    deque<double> mu(nmu);
    if (nmu == 1) {
        mu[0] = mumin;
    } else {
        double dmu = (mumax - mumin) / (nmu - 1);
        for (int imu = 0; imu < nmu; imu++) {
            mu[imu] = mumin + imu * dmu;
        }
    }

    double D = lexical_cast<double>(argv[9]);
    double theta = lexical_cast<double>(argv[10]);

    int numthreads = lexical_cast<int>(argv[11]);

    int resi = lexical_cast<int>(argv[12]);

    double Wthresh = lexical_cast<double>(argv[13]);

    bool canonical = lexical_cast<bool>(argv[14]);

#ifdef AMAZON
    //    path resdir("/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/Gutzwiller Phase Diagram");
    path resdir("/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/Canonical Transformation Gutzwiller");
#else
    //    path resdir("/Users/Abuenameh/Dropbox/Amazon EC2/Simulation Results/Gutzwiller Phase Diagram");
    path resdir("/Users/Abuenameh/Documents/Simulation Results/Canonical Transformation Gutzwiller");
#endif
    if (!exists(resdir)) {
        cerr << "Results directory " << resdir << " does not exist!" << endl;
        exit(1);
    }
    for (int iseed = 0; iseed < nseed; iseed++, seed++) {
        ptime begin = microsec_clock::local_time();


        ostringstream oss;
        oss << "res." << resi << ".txt";
        path resfile = resdir / oss.str();
        while (exists(resfile)) {
            resi++;
            oss.str("");
            oss << "res." << resi << ".txt";
            resfile = resdir / oss.str();
        }
        if (seed < 0) {
            resi = seed;
            oss.str("");
            oss << "res." << resi << ".txt";
            resfile = resdir / oss.str();
        }

        Parameter xi;
        xi.fill(1);
        //        xi.assign(1);
        rng.seed(seed);

        int xiset = 0;
        double threshold = 0;

        while (true) {
            if (seed > -1) {
                for (int j = 0; j < L; j++) {
                    xi[j] = (1 + D * uni(rng));
                }
            }

            double W[L]; //, U[L], J[L];
            vector<double> U(L), J(L);
            for (int i = 0; i < L; i++) {
                W[i] = xi[i] * Wthresh;
            }
            for (int i = 0; i < L; i++) {
                U[i] = UW(W[i]) / UW(xmax);
                J[i] = JWij(W[i], W[mod(i + 1)]) / UW(xmax);
            }
            bool reject = false;
            //        for (int i = 0; i < L; i++) {
            //            double threshold = 1.2*(JWij(xmax,xmax)/UW(xmax));
            //                if (J[i]/U[i] > threshold || J[mod(i-1)]/U[i] > threshold) {
            //                    iseed--;
            //                    seed++;
            //                    reject = true;
            //                    break;
            //                }
            //        }
            for (int i = 0; i < L; i++) {
                U[i] = UW(W[i]) / UW(Wthresh);
            }
            for (int i = 0; i < L; i++) {
                int j1 = mod(i - 1);
                int j2 = mod(i + 1);
                for (int n = 0; n < nmax; n++) {
                    for (int m = 1; m <= nmax; m++) {
                        if (n != m - 1) {
                            if (fabs(eps(U, i, j1, n, m)) < threshold || fabs(eps(U, i, j2, n, m)) < threshold) {
                                reject = true;
                                break;
                            }
                        }
                    }
                    if (reject) {
                        break;
                    }
                }
                if (reject) {
                    break;
                }
            }
            //        if (reject) {
            //            iseed--;
            //            seed++;
            //            continue;
            //        }
            if (!reject) {
                break;
            }
            xiset++;
        }

        //        rng.seed(seed);
        nu = vector<double>(L, 0);
        //        for (int i = 0; i < L; i++) {
        //            nu[i] = 0.25 * uni(rng);
        //        }

        int Lres = L;
        int nmaxres = nmax;

        boost::filesystem::ofstream os(resfile);
        printMath(os, "canonical", resi, canonical);
        printMath(os, "Lres", resi, Lres);
        printMath(os, "nmaxres", resi, nmaxres);
        printMath(os, "seed", resi, seed);
        printMath(os, "theta", resi, theta);
        printMath(os, "Delta", resi, D);
        printMath(os, "xres", resi, x);
        printMath(os, "mures", resi, mu);
        printMath(os, "xires", resi, xi);
        printMath(os, "xiset", resi, xiset);
        printMath(os, "threshold", resi, threshold);
        os << flush;

        cout << "Res: " << resi << endl;

        //        multi_array<vector<double>, 2> Jres(extents[nx][nmu]);
        //        multi_array<vector<double>, 2> Ures(extents[nx][nmu]);
        //        //        multi_array<double, 2 > fcres(extents[nx][nmu]);
        //        multi_array<double, 2 > fsres(extents[nx][nmu]);
        //        //        multi_array<double, 2> dur(extents[nx][nmu]);
        //        //        multi_array<int, 2> iterres(extents[nx][nmu]);
        //        multi_array<double, 2 > fminres(extents[nx][nmu]);
        //        multi_array<vector<double>, 2> fn0res(extents[nx][nmu]);
        //        multi_array<vector<double>, 2> fmaxres(extents[nx][nmu]);
        //        multi_array<vector<double>, 2> f0res(extents[nx][nmu]);
        //        multi_array<vector<double>, 2> fthres(extents[nx][nmu]);
        //        multi_array<vector<double>, 2> f2thres(extents[nx][nmu]);
        //        multi_array<double, 2> E0res(extents[nx][nmu]);
        //        multi_array<double, 2> Ethres(extents[nx][nmu]);
        //        multi_array<double, 2> E2thres(extents[nx][nmu]);
        //        multi_array<double, 2> res0(extents[nx][nmu]);
        //        multi_array<double, 2> resth(extents[nx][nmu]);
        //        multi_array<double, 2> res2th(extents[nx][nmu]);
        //
        //        fresults fres = {fminres, fn0res, fmaxres, f0res, fthres, f2thres, f2thres};
        //        results results = {E0res, Ethres, E2thres, E2thres, fsres, res0, resth, res2th, res2th, Jres, Ures};

//        progress_display progress(nx * nmu);

        //        cout << "Queueing" << endl;
        double muwidth = 0.2;
        queue<Point> points;
        for (int ix = 0; ix < nx; ix++) {
//            double mu0 = x[ix] / 1e12 + 0.05;
            double mu0 = 7.142857142857143e-13*x[ix] + 0.08571428571428572;
            double mui = max(mumin, mu0 - muwidth);
            double muf = min(mumax, mu0 + muwidth);
            deque<double> mu(nmu);
            if (nmu == 1) {
                mu[0] = mui;
            } else {
                double dmu = (muf - mui) / (nmu - 1);
                for (int imu = 0; imu < nmu; imu++) {
                    mu[imu] = mui + imu * dmu;
                }
            }
            for (int imu = 0; imu < nmu; imu++) {
                Point point;
                point.x = x[ix];
                point.mu = mu[imu];
                points.push(point);
            }
        }
        for (int ix = 0; ix < nx; ix++) {
//            double mu0 = -3*x[ix] / 1e12 + 0.96;
            double mu0 = -2.142857142857143e-12*x[ix] + 0.942857142857143;
            double mui = max(mumin, mu0 - muwidth);
            double muf = min(mumax, mu0 + muwidth);
            deque<double> mu(nmu);
            if (nmu == 1) {
                mu[0] = mui;
            } else {
                double dmu = (muf - mui) / (nmu - 1);
                for (int imu = 0; imu < nmu; imu++) {
                    mu[imu] = mui + imu * dmu;
                }
            }
            for (int imu = 0; imu < nmu; imu++) {
                Point point;
                point.x = x[ix];
                point.mu = mu[imu];
                points.push(point);
            }
        }
        progress_display progress(points.size());
//        for (int imu = 0; imu < nmu; imu++) {
//            queue<Point> rowpoints;
//            for (int ix = 0; ix < nx; ix++) {
//                Point point;
//                point.x = x[ix];
//                point.mu = mu[imu];
//                points.push(point);
//            }
//        }

        phase_parameters parms;
        parms.theta = theta;
        parms.canonical = canonical;

        vector<PointResults> pointRes;

        //        cout << "Dispatching" << endl;
        thread_group threads;
        for (int i = 0; i < numthreads; i++) {
            //                        threads.emplace_back(phasepoints, std::ref(xi), theta, std::ref(points), std::ref(f0res), std::ref(E0res), std::ref(Ethres), std::ref(fsres), std::ref(progress));
            threads.create_thread(bind(&phasepoints, boost::ref(xi), parms, boost::ref(points), boost::ref(pointRes), boost::ref(progress)));
        }
        threads.join_all();

        vector<pair<double, double> > Wmu;
        vector<vector<double> > Js;
        vector<vector<double> > Us;
        vector<double> fs;
        vector<double> fmin;
        vector<vector<double> > fn0;
        vector<vector<double> > fmax;
        vector<vector<double> > f0;
        vector<vector<double> > fth;
        vector<vector<double> > f2th;
        vector<double> E0;
        vector<double> Eth;
        vector<double> E2th;

        for (vector<PointResults>::iterator iter = pointRes.begin(); iter != pointRes.end(); ++iter) {
            PointResults pres = *iter;
            Wmu.push_back(make_pair(pres.W, pres.mu));
            Js.push_back(pres.J);
            Us.push_back(pres.U);
            fs.push_back(pres.fs);
            fmin.push_back(pres.fmin);
            fn0.push_back(pres.fn0);
            fmax.push_back(pres.fmax);
            f0.push_back(pres.f0);
            fth.push_back(pres.fth);
            f2th.push_back(pres.f2th);
            E0.push_back(pres.E0);
            Eth.push_back(pres.Eth);
            E2th.push_back(pres.E2th);
        }

        printMath(os, "Wmu", resi, Wmu);
        printMath(os, "Js", resi, Js);
        printMath(os, "Us", resi, Us);
        printMath(os, "fs", resi, fs);
        printMath(os, "fn0", resi, fn0);
        printMath(os, "fmin", resi, fmin);
        printMath(os, "fmax", resi, fmax);
        printMath(os, "f0", resi, f0);
        printMath(os, "fth", resi, fth);
        printMath(os, "f2th", resi, f2th);
        printMath(os, "E0", resi, E0);
        printMath(os, "Eth", resi, Eth);
        printMath(os, "E2th", resi, E2th);

        ptime end = microsec_clock::local_time();
        time_period period(begin, end);
        cout << endl << period.length() << endl << endl;

        os << "runtime[" << resi << "]=\"" << period.length() << "\";" << endl;
    }

    //    time_t start = time(NULL);
    //


    //    time_t end = time(NULL);
    //
    //    printf("Runtime: %ld", end - start);

    return 0;
}

