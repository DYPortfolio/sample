//
//  Lattractor.cpp
//  opcltest
//
//  Created by Dmitry Yamolsky on 8/12/22.
//SDL visual demo for:
//Report on "Solving ordinary differential equations on GPUs"
//by Karsten Ahnert, Denis Demidov, and Mario Mulansky


#include <iostream>
#include <vector>
#include <chrono>

#define NO_SDL_GLEXT

#include "Lattractor.hpp"
#include "runge_kutta4.hpp"
#include "sawki_settings.h"


#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_TARGET_OPENCL_VERSION 120
#include <vexcl/vexcl.hpp>

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include <Opengl/glext.h>

using namespace std;

typedef std::vector<double> state_type;
struct lorenz {
    const double sigma, R, b;
    lorenz(const double sigma, const double R, const double b)
    : sigma(sigma), R(R), b(b) { }
    void operator()(const state_type &x,state_type &dxdt,double t)
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};



typedef vex::vector<double> vex3x_type;
struct vex_lorenz {
    const double sigma, R, b;
    vex_lorenz(const double sigma, const double R, const double b)
    : sigma(sigma), R(R), b(b) { }
    
    void operator()(const vex3x_type &x, vex3x_type &dxdt,
                    double t) const
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};


typedef vex::vector< double >    vector_type;
typedef vex::multivector< double, 3 > state_type_multi;
const double sigma = 10.0;
const double b = 8.0 / 3.0;

struct sys_func
{
    const vector_type &R;
    sys_func( const vector_type &_R ) : R( _R ) { }
    void operator()( const state_type_multi &x , state_type_multi &dxdt , double t ) const
    {
        dxdt(0) = -sigma * ( x(0) - x(1) );
        dxdt(1) = R * x(0) - x(1) - x(0) * x(2);
        dxdt(2) = - b * x(2) + x(0) * x(1);
    }
};

void lattractor() {
    
    
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    const int steps = 500000;
    const double dt = 0.003;
    const double t_max = .1;
    const size_t n = 2;
    double Rmin = 0.1 , Rmax = 50.0 , dR = ( Rmax - Rmin ) / double( n - 1 );
    
    std::deque<double> trackQueue;
    vex::Context ctx( vex::Filter::DoublePrecision );

    glPointSize(1.2f);
    glColor3f(1,1,1);
    glLineWidth(1.0f);
    std::vector<double> drawcoord(3);
    std::vector<double> drawcoord1(3);
    
    double pt_dist;

    rk4_stepper_type stepper(3);
    lorenz system(10.0, 28.0, 8.0/3.0);
    state_type x(3, 1.0);
    state_type x2(3, 1.0);
    x = {5.0, .1, .2};
    x2 = {5.0, .11, .2};

    gpu_stepper_type vex_stepper(3);
    vex_lorenz vex_system(10.0, 28.0, 8.0/3.0);//3 parameters
    vex_type vex_x(ctx,3);

    //initial state
    vex_x[0] = 5.0;
    vex_x[1] = .1;
    vex_x[2] = .2;
    
    for (unsigned int i = 0;i<2;i++)
    {
        drawcoord[i]  = (x[i]/20) ;
        drawcoord[i] = min(drawcoord[i],1.0);
        drawcoord[i] = max(drawcoord[i],-1.0);
     
        drawcoord1[i]  = (x2[i]/20) ;
        drawcoord1[i] = min(drawcoord1[i],1.0);
        drawcoord1[i] = max(drawcoord1[i],-1.0);
        
        
    }
    
    
    trackQueue.push_back(drawcoord[0] );
    trackQueue.push_back(drawcoord[1] );
    trackQueue.push_back(drawcoord1[0] );
    trackQueue.push_back(drawcoord1[1] );
     
    auto timer1 = high_resolution_clock::now();
    duration <double, std::milli> ms;
    
    for( size_t n=0 ; n<steps ; ++n ) {
        
        timer1 = high_resolution_clock::now();
        
        stepper.do_step(system, x, n*dt, dt);
        stepper.do_step(system, x2, n*dt, dt);
        //         vex_stepper.do_step(vex_system, vex_x, n*dt, dt);
        
        ms = ms + (high_resolution_clock::now()-timer1);
        
        
        
        for (unsigned int i = 0;i<2;i++)
        {
            drawcoord[i]  = (x[i]/20) ;
            drawcoord[i] = min(drawcoord[i],1.0);
            drawcoord[i] = max(drawcoord[i],-1.0);
            
            drawcoord1[i]  = (x2[i]/20) ;
            drawcoord1[i] = min(drawcoord1[i],1.0);
            drawcoord1[i] = max(drawcoord1[i],-1.0);
            
        }
        
        
        pt_dist = ((drawcoord1[0] - drawcoord[0]) * (drawcoord1[0] - drawcoord[0]) + (drawcoord1[1] - drawcoord[1]) * (drawcoord1[1] - drawcoord[1]));
        
        cout << "iteration #    " <<   "distance         " <<     "average time per do_step, in ms" <<endl;
        cout << n <<"          "<< pt_dist << "        " << ms.count()/n <<endl;
        
        SDL_Event event;
        SDL_PollEvent(&event);
        
        Uint8 *keys = SDL_GetKeyState(NULL) ;
        if( keys[SDLK_w] ) {
            break;
        }
        

        
        
        trackQueue.push_back(drawcoord[0] );
        trackQueue.push_back(drawcoord[1] );
        trackQueue.push_back(drawcoord1[0] );
        trackQueue.push_back(drawcoord1[1] );
        
        
        glAccum(GL_MULT, TAIL_LENGTH);
        glAccum(GL_RETURN, TAIL_BRIGHTNESS_RATIO);

        glPushMatrix();
        glBegin(GL_LINES);
        {
            for (unsigned int i = 0;i<trackQueue.size()/4;i++)
                glVertex2f(trackQueue[i*4],  trackQueue[i*4+1] );
        }
        glEnd();
        
        //

        glBegin(GL_LINES);
        {
            for (unsigned int i = 0;i<trackQueue.size()/4;i++)
                glVertex2f(trackQueue[i*4+2],  trackQueue[i*4+3] );
        }
        glEnd();
        
        
        glColor3f(1,0,0);
        glPointSize(3.0f);
        
glBegin(GL_POINTS);
{
        glVertex2f(drawcoord[0],  drawcoord[1]);
}
        
glBegin(GL_POINTS);
        {
                glVertex2f(drawcoord1[0],  drawcoord1[1]);
        }
        
        glEnd();
        glPointSize(1.2f);
        glColor3f(1,1,1);
        glPopMatrix();
        glAccum(GL_LOAD, TAIL_LENGTH);
        trackQueue.pop_front();
        trackQueue.pop_front();
        trackQueue.pop_front();
        trackQueue.pop_front();
        
        if (fmod(n,LATTRACTOR_REFRESH_PERIOD) == 0)
        SDL_GL_SwapBuffers();
    }
}
