//
//  runge_kutta4.hpp
//  opcltest
//
//  Created by Khan on 8/12/22.
//some from boost

#ifndef runge_kutta4_hpp
#define runge_kutta4_hpp

#include <stdio.h>
#include <vector>
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_TARGET_OPENCL_VERSION 120
#include <vexcl/vexcl.hpp>
#include <vexcl/vector.hpp>

struct container_algebra {
    
    template<class S1, class S2, class S3, class Op>
    static void for_each3(S1 &s1, S2 &s2, S3 &s3, Op op) {
        const size_t dim = s1.size();
        for(size_t n = 0; n < dim; ++n)
        op(s1[n], s2[n], s3[n]);
    }
    
    template<class S1, class S2, class S3,class S4, class S5, class S6, class Op>
    static void for_each6(S1 &s1, S2 &s2, S3 &s3,S4 &s4, S5 &s5, S6 &s6, Op op) {
        const size_t dim = s1.size();
        for(size_t n = 0; n < dim; ++n)
        op(s1[n], s2[n], s3[n],s4[n], s5[n], s6[n]);
    }
};



struct default_operations {
    template<class Fac1 = double, class Fac2 = Fac1>
    struct scale_sum2 {
        typedef void result_type;
        const Fac1 alpha1;
        const Fac2 alpha2;
        scale_sum2(Fac1 alpha1, Fac2 alpha2)
        : alpha1(alpha1), alpha2(alpha2) { }
        
        template<class T0, class T1, class T2>
        void operator()(T0 &t0, const T1 &t1, const T2 &t2) const {
            t0 = alpha1 * t1 + alpha2 * t2;
        }
    };
    //etc.

        
        
    template<class Fac1 = double, class Fac2 = Fac1, class Fac3 = Fac1, class Fac4 = Fac1,class Fac5 = Fac1>
    struct scale_sum5 {
        
        typedef void result_type;
        const Fac1 alpha1;
        const Fac2 alpha2;
        const Fac3 alpha3;
        const Fac4 alpha4;
        const Fac5 alpha5;
        
        scale_sum5(Fac1 alpha1, Fac2 alpha2, Fac3 alpha3, Fac4 alpha4, Fac5 alpha5)
        : alpha1(alpha1), alpha2(alpha2), alpha3(alpha3), alpha4(alpha4), alpha5(alpha5) { }
        
        template<class T0, class T1, class T2, class T3, class T4, class T5>
        void operator()(T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5) const {
            t0 = alpha1 * t1 + alpha2 * t2 + alpha3 * t3 + alpha4 * t4 + alpha5 * t5;
        }
    };
};

struct vector_space_algebra {
    template<class S1, class S2, class S3, class Op>
    static void for_each3(S1 &s1, S2 &s2, S3 &s3, Op op) {
        op(s1, s2, s3);
    }
    //etc.
    template<class S1, class S2, class S3, class S4, class S5, class S6, class Op>
    static void for_each6(S1 &s1, S2 &s2, S3 &s3, S4 &s4, S5 &s5, S6 &s6, Op op) {
        op(s1, s2, s3, s4, s5, s6);
    }
};

template<class T, size_t N>
void resize(const boost::array<T, N> &, boost::array<T,N>& ) {
}
template<class state_type>
void resize(const state_type &in, state_type &out) {
    out.resize(in.size());
}
template<class T>
void resize(const vex::vector<T> &in, vex::vector<T> &out) {
    out.resize(in.queue_list(), in.size());
}
template<class T, size_t N>
void resize(const vex::multivector<T,N> &in,
            vex::multivector<T,N> &out)
{
    out.resize(in.queue_list(), in.size());
}

template< class state_type, class value_type = double, class time_type = value_type,
class algebra = container_algebra, class operations = default_operations>
class runge_kutta4 {
public:
    runge_kutta4(size_t N)
    : N(N), x_tmp(N), k1(N), k2(N), k3(N), k4(N) { }
    template<typename System>
    void do_step(System system, state_type &x, double t, double dt)
    {
        adjust_size(x);
        const value_type one = 1;
        const time_type dt2 = dt/2, dt3 = dt/3, dt6 = dt/6;
        
        typedef typename operations::template scale_sum2<
        value_type, time_type> scale_sum2;
        typedef typename operations::template scale_sum5<
        value_type, time_type, time_type,
        time_type, time_type> scale_sum5;
        
        system(x, k1, t);
        algebra::for_each3(x_tmp, x, k1, scale_sum2(one, dt2));
        system(x_tmp, k2, t + dt2);
        algebra::for_each3(x_tmp, x, k2, scale_sum2(one, dt2));
        system(x_tmp, k3, t + dt2);
        algebra::for_each3(x_tmp, x, k3, scale_sum2(one, dt));
        system(x_tmp, k4, t + dt);
        algebra::for_each6(x, x, k1, k2, k3, k4,
                           scale_sum5(one, dt6, dt3, dt3, dt6));
    }
private:
    void adjust_size(const state_type &x) {
        resize(x, x_tmp);
        resize(x, k1);
        resize(x, k2);
        resize(x, k3);
        resize(x, k4);
    }
    const size_t N;
    state_type x_tmp, k1, k2, k3, k4;
};

typedef runge_kutta4< std::vector<double>, double, double,
container_algebra, default_operations > rk4_stepper_type;



typedef runge_kutta4< std::vector<double> > rk_stepper;

typedef runge_kutta4< std::vector<double>, double, double,
container_algebra, default_operations > rk4_stepper_type;

typedef runge_kutta4< vex::multivector<double, 3> , double, double,
vector_space_algebra, default_operations > rk4_ensemble_stepper_type;


typedef vex::vector<double> vex_type;

typedef runge_kutta4< vex::vector<double>, double, double,
vector_space_algebra, default_operations > gpu_stepper_type;


#endif /* runge_kutta4_hpp */

