// Uniform radial-velocity photon kernel. Angles are degrees, speeds in c.
// K(x) = < w_phi(mu) sqrt(1-beta^2)/(2 beta) >_beta, x=E/E0.
// Its integrated-photon weight is delta^2, as in RED-WINE/Nature Eq. (3).
#ifndef WINDKIN_UNIFORM_H
#define WINDKIN_UNIFORM_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <list>
#include <mutex>
#include <memory>
#include <stdexcept>
#include <vector>
#if __has_include("fftw3.h")
#include "fftw3.h"
#else
#include "fftw/fftw3.h"
#endif

namespace windkin_uniform {
constexpr double pi = 3.1415926535897932384626433832795;
constexpr double eps = std::numeric_limits<double>::epsilon();
inline double clip(double x, double a, double b) { return std::max(a, std::min(b, x)); }
struct Rule {
    std::vector<double> x, w, position, jacobian;
    explicit Rule(int n) : x(n), w(n), position(n), jacobian(n) {
        for (int i=0; i<(n+1)/2; ++i) {
            double z=std::cos(pi*(i+.75)/(n+.5)), derivative=0;
            for (int iteration=0; iteration<30; ++iteration) {
                double p=1, previous=0;
                for (int k=1; k<=n; ++k) {
                    const double old=p;
                    p=((2*k-1)*z*p-(k-1)*previous)/k; previous=old;
                }
                derivative=n*(z*p-previous)/(z*z-1);
                const double change=p/derivative; z-=change;
                if (std::abs(change)<2*eps) break;
            }
            x[i]=(1-z)/2; x[n-1-i]=(1+z)/2;
            w[i]=w[n-1-i]=1/((1-z*z)*derivative*derivative);
        }
        for (int i=0; i<n; ++i) {
            const double sine=std::sin(pi*x[i]/2);
            position[i]=sine*sine; jacobian[i]=pi/2*std::sin(pi*x[i]);
        }
    }
};
inline const Rule& rule4() { static const Rule r(4); return r; }
inline const Rule& rule8() { static const Rule r(8); return r; }
inline const Rule& rule16() { static const Rule r(16); return r; }

struct Cap {
    double ci, si, cin, cout, lower, upper, width;
    std::vector<double> tangencies;
    Cap(double incl, double inner, double outer, double center, double spread)
        : ci(incl==pi/2 ? 0 : std::cos(incl)), si(std::sin(incl)),
          cin(inner==pi/2 ? 0 : std::cos(inner)), cout(outer==pi/2 ? 0 : std::cos(outer)),
          lower(center-spread/2), upper(center+spread/2), width(spread) {
        const double limb=std::abs(pi/2-incl);
        if (si>4*eps) {
            if (std::abs(inner-limb)<=4*eps || std::abs(inner-(pi-limb))<=4*eps)
                cin=std::copysign(si,cin);
            if (std::abs(outer-limb)<=4*eps || std::abs(outer-(pi-limb))<=4*eps)
                cout=std::copysign(si,cout);
        }
        for (double angle : {incl-inner,incl+inner,incl-outer,incl+outer,0.,pi}) {
            double mu=std::cos(angle); if (std::abs(mu)<=4*eps) mu=0;
            tangencies.push_back(mu);
        }
        std::sort(tangencies.begin(),tangencies.end());
        tangencies.erase(std::unique(tangencies.begin(),tangencies.end()),tangencies.end());
    }
    double mu(double beta, double q) const {
        return q/beta+(1-q)*beta/(1+std::sqrt((1-beta)*(1+beta)));
    }
    double arcs(double mu) const {
        if (std::abs(mu)>1) return 0;
        const double center=ci*mu, sine=std::sqrt(std::max(0.,(1-mu)*(1+mu)));
        if (si*sine<=0) return (cout<=center && center<=cin) ? 1 : 0;
        const double radius=si*sine, correction=si*mu*mu/(1+sine);
        auto cap=[&](double cosine) {
            const double minus=std::abs(mu)<.25 ? (si-cosine)+center-correction : radius-(cosine-center);
            const double plus=std::abs(mu)<.25 ? (si+cosine)-center-correction : radius+(cosine-center);
            if (minus<=0) return 0.;
            if (plus<=0) return pi/2;
            return std::atan2(std::sqrt(minus),std::sqrt(plus));
        };
        return std::max(0.,2*(cap(cout)-cap(cin))/pi);
    }
    double density(double beta, double q) const {
        return arcs(mu(beta,q))*std::sqrt((1-beta)*(1+beta))/beta;
    }
    double quadrature(double a, double b, double q, double left, double right, const Rule& rule) const {
        double sum=0;
        for (size_t k=0; k<rule.x.size(); ++k) {
            const double t=left+(right-left)*rule.x[k];
            const double sine=left==0 && right==1 ? 0 : std::sin(pi*t/2);
            const double position=left==0 && right==1 ? rule.position[k] : sine*sine;
            const double jacobian=left==0 && right==1 ? rule.jacobian[k] : pi/2*std::sin(pi*t);
            sum+=rule.w[k]*density(a+(b-a)*position,q)*jacobian;
        }
        return sum*(right-left)*(b-a)/width;
    }
    double integrate(double a, double b, double q, double left=0, double right=1, int depth=0) const {
        const double low=quadrature(a,b,q,left,right,rule8());
        const double high=quadrature(a,b,q,left,right,rule16());
        if (std::abs(high-low)<=2e-9*((b-a)/width/((depth+1.)*(depth+1.))+std::abs(high))) return high;
        if (depth>=89 || left==(left+right)/2 || right==(left+right)/2)
            throw std::runtime_error("wuconv velocity quadrature did not converge");
        const double mid=(left+right)/2;
        return integrate(a,b,q,left,mid,depth+1)+integrate(a,b,q,mid,right,depth+1);
    }
    double evaluate(double logx, bool mono) const {
        const double q=-std::expm1(-logx);
        if (mono) return density((lower+upper)/2,q)/2;
        std::vector<double> cuts{lower,upper};
        // Unsquared Doppler roots in u=tanh(atanh(beta)/2).
        for (double m : tangencies) {
            const double discriminant=m*m-q*(2-q);
            if (discriminant<0) continue;
            const double large=m+std::copysign(std::sqrt(discriminant),m);
            const double roots[2]={large/(2-q),large!=0 ? q/large : 0};
            for (double u : roots) if (u>0 && u<1) cuts.push_back(clip(2*u/(1+u*u),lower,upper));
        }
        std::sort(cuts.begin(),cuts.end());
        double integral=0;
        for (size_t j=1; j<cuts.size(); ++j) {
            const double a=cuts[j-1], b=cuts[j];
            if (b<=a) continue;
            const double w=arcs(mu(a+(b-a)/2,q));
            if (w==0) continue;
            if (w==1 || si<=eps || (ci==0 && ((cin==1 && cout==0) || (cin==0 && cout==-1)))) {
                const double sa=std::sqrt((1-a)*(1+a)), sb=std::sqrt((1-b)*(1+b));
                const double ds=-(b-a)*(b+a)/(sb+sa);
                double primitive;
                if (a>=.8) {
                    const double product=sa*sb, x=-ds/(1-product);
                    const double remainder=x<.001 ? x*x*x*(1./3+x*x*(1./5+x*x/7)) : std::atanh(x)-x;
                    primitive=remainder+x*product;
                } else {
                    primitive=std::log1p((b-a)/a)-std::log1p(ds/(1+sa))+ds;
                }
                integral+=w*primitive/width;
            } else integral+=integrate(a,b,q);
        }
        return integral/2;
    }
};

using Pair=std::array<double,2>;
struct Panel {
    double left, right;
    std::array<Pair,4> coefficient; // Ascending powers of t=(log(x)-left)/(right-left).
    Pair cumulative;
};
struct Moments {
    std::vector<Cap> caps;
    std::vector<Panel> panels;
    bool mono;
    double lo, hi, fraction, tolerance;
    Pair total{{0,0}};
    Pair density(double logx) const {
        double value=0;
        for (const auto& cap : caps) value+=cap.evaluate(logx,mono);
        return {{value*std::exp(logx),value}};
    }
    void tabulate(double left, double right, int depth=0) {
        const double span=right-left;
        const auto& fit=rule4(); const auto& check=rule8();
        Panel p; p.left=left; p.right=right;
        for (auto& c : p.coefficient) c={{0,0}};
        // Four-node Lagrange interpolation, retaining normalized coordinates.
        for (int j=0; j<4; ++j) {
            std::array<double,4> polynomial{{1,0,0,0}};
            int degree=0;
            for (int k=0; k<4; ++k) if (k!=j) {
                const double divisor=fit.x[j]-fit.x[k];
                for (int d=degree+1; d>=0; --d)
                    polynomial[d]=((d ? polynomial[d-1] : 0)-fit.x[k]*(d<=degree ? polynomial[d] : 0))/divisor;
                ++degree;
            }
            const Pair value=density(left+span*fit.x[j]);
            for (int d=0; d<4; ++d) for (int m=0; m<2; ++m) p.coefficient[d][m]+=polynomial[d]*value[m];
        }
        Pair error{{0,0}}, scale{{0,0}};
        for (int k=0; k<8; ++k) {
            const Pair actual=density(left+span*check.x[k]);
            for (int m=0; m<2; ++m) {
                double fitted=0;
                for (int d=3; d>=0; --d) fitted=fitted*check.x[k]+p.coefficient[d][m];
                error[m]+=check.w[k]*std::abs(fitted-actual[m]);
                scale[m]+=check.w[k]*std::abs(actual[m]);
            }
        }
        bool accept=true;
        for (int m=0; m<2; ++m)
            accept=accept && span*error[m]<=tolerance*(span*scale[m]+fraction*(std::sqrt(span/(hi-lo))+1e-4));
        if (!accept) {
            const double mid=left+span/2;
            if (depth>=48 || mid==left || mid==right) throw std::runtime_error("wuconv kernel tabulation did not converge");
            tabulate(left,mid,depth+1); tabulate(mid,right,depth+1); return;
        }
        p.cumulative=total;
        for (int m=0; m<2; ++m) for (int d=0; d<4; ++d) total[m]+=span*p.coefficient[d][m]/(d+1);
        panels.push_back(p);
    }
    Moments(double center, double width, double inner, double outer, double incl, bool bicone,
            double lower, double upper, double tol=1e-7) : tolerance(tol) {
        const double speed=std::abs(center), spacing=std::nextafter(speed,INFINITY)-speed;
        mono=width<=4*spacing;
        fraction=std::sin((outer+inner)*pi/360)*std::sin((outer-inner)*pi/360)*(1+bicone);
        const double angle=(center<0 ? 180-incl : incl)*pi/180;
        caps.emplace_back(angle,inner*pi/180,outer*pi/180,speed,mono ? 0 : width);
        if (bicone) caps.emplace_back(pi-angle,inner*pi/180,outer*pi/180,speed,mono ? 0 : width);
        const double envelope=std::atanh(speed+(mono ? 0 : width/2));
        lo=std::max(lower,-envelope); hi=std::min(upper,envelope);
        if (lo>=hi || fraction==0) return;
        std::vector<double> knots{lo,hi,clip(0.,lo,hi)};
        for (const auto& cap : caps) for (double mu : cap.tangencies) {
            std::vector<double> speeds{cap.lower,cap.upper};
            if (cap.lower<mu && mu<cap.upper) speeds.push_back(mu);
            for (double beta : speeds) knots.push_back(clip(.5*(std::log1p(-beta)+std::log1p(beta))-std::log1p(-beta*mu),lo,hi));
        }
        std::sort(knots.begin(),knots.end());
        std::vector<double> unique;
        for (double knot : knots) if (unique.empty() || knot-unique.back()>8*eps*std::max(std::abs(unique.back()),std::min(1e-12,hi-lo))) unique.push_back(knot);
        const int intervals=static_cast<int>(std::ceil((hi-lo)/.12));
        for (int j=0; j<=intervals; ++j) unique.push_back(lo+(hi-lo)*j/intervals);
        std::sort(unique.begin(),unique.end());
        unique.erase(std::unique(unique.begin(),unique.end()),unique.end());
        for (size_t j=1; j<unique.size(); ++j) tabulate(unique[j-1],unique[j]);
    }
    Pair cumulative(double logx) const {
        if (panels.empty() || logx<=lo) return {{0,0}};
        if (logx>=hi) return total;
        auto it=std::upper_bound(panels.begin(),panels.end(),logx,[](double value,const Panel& p){return value<p.left;});
        const Panel& p=*(it-1);
        const double span=p.right-p.left, t=(logx-p.left)/span;
        Pair value=p.cumulative;
        for (int m=0; m<2; ++m) {
            double polynomial=0;
            for (int d=3; d>=0; --d) polynomial=polynomial*t+p.coefficient[d][m]/(d+1);
            value[m]+=span*t*polynomial;
        }
        return value;
    }
    double partial(double y, double a, double b) const {
        if (y<=0 || panels.empty()) return 0;
        const double tb=clip(std::log(y/b),lo,hi), ta=a>0 ? clip(std::log(y/a),lo,hi) : hi;
        if (a>0 && (b-a)/a<1e-5) {
            // Integrate the CDF over fractional source energy, split at its knots.
            const double ratio=(b-a)/b, origin=std::log(y/b);
            std::vector<double> cuts{0,1};
            for (const auto& p : panels) if (p.left>tb && p.left<ta) cuts.push_back(clip(-std::expm1(origin-p.left)/ratio,0.,1.));
            if (hi>tb && hi<ta) cuts.push_back(clip(-std::expm1(origin-hi)/ratio,0.,1.));
            std::sort(cuts.begin(),cuts.end());
            double result=0;
            for (size_t j=1; j<cuts.size(); ++j) for (int k=0; k<8; ++k) {
                const double u=cuts[j-1]+(cuts[j]-cuts[j-1])*rule8().x[k];
                result+=(cuts[j]-cuts[j-1])*rule8().w[k]*cumulative(origin-std::log1p(-ratio*u))[0];
            }
            return result;
        }
        const Pair lower=cumulative(tb), upper=cumulative(ta);
        return clip(lower[0]+(y*(upper[1]-lower[1])-a*(upper[0]-lower[0]))/(b-a),0.,total[0]);
    }
};

// FFTW execution is thread-safe on separate buffers; its planner is not.
inline std::mutex& plannerMutex() { static std::mutex mutex; return mutex; }

// FFTW plans and buffers persist per thread. No physics is cached across parameters.
struct FFT {
    int size;
    double *a, *b;
    fftw_complex *fa, *fb;
    fftw_plan forwarda, forwardb, backward;
    explicit FFT(int n) : size(n) {
        const std::lock_guard<std::mutex> lock(plannerMutex());
        a=fftw_alloc_real(n); b=fftw_alloc_real(n);
        fa=fftw_alloc_complex(n/2+1); fb=fftw_alloc_complex(n/2+1);
        if (!a || !b || !fa || !fb) throw std::bad_alloc();
        forwarda=fftw_plan_dft_r2c_1d(n,a,fa,FFTW_ESTIMATE);
        forwardb=fftw_plan_dft_r2c_1d(n,b,fb,FFTW_ESTIMATE);
        backward=fftw_plan_dft_c2r_1d(n,fa,a,FFTW_ESTIMATE);
        if (!forwarda || !forwardb || !backward) throw std::runtime_error("wuconv FFTW plan creation failed");
    }
    ~FFT() {
        const std::lock_guard<std::mutex> lock(plannerMutex());
        fftw_destroy_plan(forwarda); fftw_destroy_plan(forwardb); fftw_destroy_plan(backward);
        fftw_free(a); fftw_free(b); fftw_free(fa); fftw_free(fb);
    }
    FFT(const FFT&)=delete;
    FFT& operator=(const FFT&)=delete;
};

inline std::vector<double> convolve(const std::vector<double>& edges, const std::vector<double>& flux,
                                    const std::array<double,6>& p) {
    const size_t n=flux.size();
    if (n==0 || edges.size()!=n+1) throw std::invalid_argument("wuconv requires n+1 edges and n integrated fluxes");
    for (size_t j=0; j<=n; ++j) if (!std::isfinite(edges[j]) || edges[j]<0 || (j && edges[j]<=edges[j-1]))
        throw std::invalid_argument("wuconv requires finite nonnegative increasing energies");
    for (double v : p) if (!std::isfinite(v)) throw std::invalid_argument("wuconv parameters must be finite");
    for (double v : flux) if (!std::isfinite(v)) throw std::invalid_argument("wuconv flux must be finite");
    const double speed=std::abs(p[0]), width=p[1]; const bool bicone=p[5]!=0;
    if (speed>=1 || width<0 || speed-width/2<0 || speed+width/2>=1 || p[2]<0 || p[3]<p[2]
        || p[3]>(bicone ? 90 : 180) || p[4]<0 || p[4]>180)
        throw std::invalid_argument("wuconv: invalid speed interval or hollow-cone geometry");
    std::vector<double> result(n,0);
    const double fraction=std::sin((p[3]+p[2])*pi/360)*std::sin((p[3]-p[2])*pi/360)*(1+bicone);
    if (fraction==0) return result;
    if (speed==0 && width==0) {
        for (size_t j=0; j<n; ++j) result[j]=fraction*flux[j];
        return result;
    }
    size_t first=0,last=n;
    while (first<n && flux[first]==0) ++first;
    if (first==n) return result;
    while (flux[last-1]==0) --last;
    const double lower=edges[0]>0 ? std::max(-360.,std::log(edges[0]/edges[last])) : -360;
    const double upper=edges[first]>0 ? std::min(360.,std::log(edges[n]/edges[first])) : 360;
    // Cache only exact parameter/support matches, as in the Python reference.
    struct Entry { std::array<double,8> key; std::unique_ptr<Moments> moments; };
    static thread_local std::list<Entry> cache;
    const std::array<double,8> key{{p[0],p[1],p[2],p[3],p[4],p[5],lower,upper}};
    auto found=std::find_if(cache.begin(),cache.end(),[&](const Entry& entry){return entry.key==key;});
    if (found==cache.end()) {
        std::unique_ptr<Moments> calculated(new Moments(p[0],width,p[2],p[3],p[4],bicone,lower,upper));
        cache.push_front(Entry{key,std::move(calculated)});
        if (cache.size()>8) cache.pop_back();
    } else cache.splice(cache.begin(),cache,found);
    const Moments& moments=*cache.front().moments;
    if (moments.panels.empty()) return result;
    const double step=edges[0]>0 ? std::log(edges[n]/edges[0])/n : 0;
    bool geometric=edges[0]>0;
    for (size_t j=0; geometric && j<n; ++j)
        geometric=std::abs(std::log(edges[j+1]/edges[j])-step)<=1e-14+1e-10*std::abs(step);
    if (geometric) {
        int fftsize=1; while (fftsize<static_cast<int>(3*n-2)) fftsize*=2;
        static thread_local std::unique_ptr<FFT> fft;
        if (!fft || fft->size!=fftsize) fft.reset(new FFT(fftsize));
        std::fill(fft->a,fft->a+fftsize,0.); std::fill(fft->b,fft->b+fftsize,0.);
        std::copy(flux.begin(),flux.end(),fft->a);
        const double b=std::exp(step);
        double previous=moments.partial(std::exp(clip((-static_cast<double>(n)+1)*step,-360.,360.)),1,b);
        for (size_t j=0; j<2*n-1; ++j) {
            const int lag=static_cast<int>(j)-static_cast<int>(n)+1;
            const double next=moments.partial(std::exp(clip((lag+1)*step,-360.,360.)),1,b);
            fft->b[lag>=0 ? lag : fftsize+lag]=std::max(0.,next-previous); previous=next;
        }
        fftw_execute(fft->forwarda); fftw_execute(fft->forwardb);
        for (int j=0; j<=fftsize/2; ++j) {
            const double ar=fft->fa[j][0], ai=fft->fa[j][1], br=fft->fb[j][0], bi=fft->fb[j][1];
            fft->fa[j][0]=ar*br-ai*bi; fft->fa[j][1]=ar*bi+ai*br;
        }
        fftw_execute(fft->backward);
        const bool positive=std::all_of(flux.begin(),flux.end(),[](double v){return v>=0;});
        for (size_t j=0; j<n; ++j) result[j]=positive ? std::max(0.,fft->a[j]/fftsize) : fft->a[j]/fftsize;
    } else {
        // Exact finite-bin reference for irregular grids; use an explicit log grid for inference.
        const double minshift=std::exp(moments.lo), maxshift=std::exp(moments.hi);
        for (size_t j=first; j<last; ++j) if (flux[j]!=0) {
            const double a=edges[j],b=edges[j+1];
            const size_t begin=std::max<std::ptrdiff_t>(0,std::upper_bound(edges.begin(),edges.end(),a*minshift)-edges.begin()-1);
            const size_t end=std::min<size_t>(n,std::lower_bound(edges.begin(),edges.end(),b*maxshift)-edges.begin());
            double previous=moments.partial(edges[begin],a,b);
            for (size_t k=begin; k<end; ++k) {
                const double next=moments.partial(edges[k+1],a,b);
                result[k]+=flux[j]*std::max(0.,next-previous); previous=next;
            }
        }
    }
    return result;
}
} // namespace windkin_uniform
#endif
