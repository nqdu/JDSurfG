#pragma once 
#include "numerical.hpp"
#include <fstream>

// // useful function
// void smooth_grad(fvec &grad,int nx,int ny,int nz,float sigma_h,float sigma_v);
// void find_quadratic_min(const float *x,const float *y,int n,float &xmin,float &ymin);


/**
 * nonlinear operator
 * @tparam mod forward model operator, 3 functions should have \n 
 *   1. void mod.compute_grad(const fvec &x,fvec &dsyn,fvec &grad) \n 
 *   2. float mod.compute_misfit(const fvec &dsyn); \n
 *   3. void mod.forward(const fvec &x,fvec &dsyn,)
*/
template <class T> 
class NonlinOPT {
public:
    std::string result_dir;
    float MAX_STEP;
    int ITER_START;
    const int m_store = 5;
    const float Armijo_TOL = 0.01; 

    NonlinOPT(const std::string &outdir,float max_step,int iter_start)
    {
        this -> result_dir = outdir;;
        this -> MAX_STEP = max_step;
        this -> ITER_START = iter_start;
    }

private:
    void read_vector(fvec &v,const std::string &filename) const
    { 
        std::ifstream fp(filename,std::ios::binary);
        fp.read((char*)v.data(),sizeof(float)*v.size());
        fp.close();
    }

    void write_vector(const fvec &v,const std::string &filename) const 
    { 
        std::ofstream fp(filename,std::ios::binary);
        fp.write((const char*)v.data(),sizeof(float)*v.size());
        fp.close();
    }

    void get_cg_direction(int iter, const fvec &grad, fvec &direc) const 
    {
        std::string filename;

        // read previous gradient and search direction
        int n = grad.size();
        fvec direc_prev(n),grad_prev(n);
        if(iter > ITER_START) {
            filename = result_dir + "/grad" + std::to_string(iter-1) + ".bin";
            this -> read_vector(grad_prev,filename);
            filename = result_dir + "/direc"  + std::to_string(iter-1) + ".bin";
            this -> read_vector(direc_prev,filename);
        }

        if(iter == ITER_START) { 
            direc = -grad;
        }
        else{
            float beta = (grad * (grad - grad_prev)).sum() / grad_prev.square().sum();
            direc = -grad  + beta * direc_prev;
        }

        // save current search direction
        filename = result_dir + "/direc" + std::to_string(iter) + ".bin";
        this -> write_vector(direc,filename);
    }

    void get_lbfgs_direction(int iter, const fvec &x, 
                                fvec &direc) const 
    {
        std::string file;

        // save current log model
        int n = x.size();
        fvec x1(n),x0(n);
        file = result_dir + "/logmd" + std::to_string(iter) + ".bin";
        x1 = x.log();
        this -> write_vector(x1,file);
        x1.setZero();

        // check if iter_start = iter 
        if(ITER_START == iter) {
            file = result_dir + "/grad" + std::to_string(iter) + ".bin";
            this -> read_vector(direc,file);
            direc = - direc;
            return; 
        }

        // LBFGS arrays
        fvec a(1000),p(1000);
        fvec grad1(n),grad0(n);
        fvec grad_diff(n),x_diff(n), q_vec(n),r_vec(n);
        int iter_store = iter - m_store; 
        if(iter_store <= ITER_START) iter_store = ITER_START;

        // variables used
        float b,p_k_up,p_k_down,p_k;

        // set zero
        a.setZero(); p.setZero();
        r_vec.setZero(); q_vec.setZero();

        // start backward store
        file = result_dir + "/grad" + std::to_string(iter) + ".bin";
        this -> read_vector(q_vec,file);
        for(int istore = iter-1; istore >= iter_store; istore--) {
            file = result_dir + "/grad" + std::to_string(istore+1) + ".bin";
            this -> read_vector(grad1,file);
            file = result_dir + "/grad" + std::to_string(istore) + ".bin";
            this -> read_vector(grad0,file);
            file = result_dir + "/logmd" + std::to_string(istore+1) + ".bin";
            this -> read_vector(x1,file);
            file = result_dir + "/logmd" + std::to_string(istore) + ".bin";
            this -> read_vector(x0,file);

            // grad/logmd idff
            grad_diff = grad1 - grad0; 
            x_diff = x1 - x0;

            // coefs
            float p_tmp = (grad_diff * x_diff).sum();
            float a_tmp = (x_diff * q_vec).sum();
            p[istore] = 1. / p_tmp; 
            a[istore] = p[istore] * a_tmp;
            q_vec = q_vec - a[istore] * grad_diff;
        }
        int istore = iter - 1;
        file = result_dir + "/grad" + std::to_string(istore+1) + ".bin";
        this -> read_vector(grad1,file);
        file = result_dir + "/grad" + std::to_string(istore) + ".bin";
        this -> read_vector(grad0,file);
        file = result_dir + "/logmd" + std::to_string(istore+1) + ".bin";
        this -> read_vector(x1,file);
        file = result_dir + "/logmd" + std::to_string(istore) + ".bin";
        this -> read_vector(x0,file);
        grad_diff = grad1 - grad0; 
        x_diff = x1 - x0;

        // coefs
        p_k_up = (grad_diff * x_diff).sum();
        p_k_down = (grad_diff * grad_diff).sum();
        p_k = p_k_up / p_k_down;
        r_vec = p_k * q_vec;

        // forward store
        for(int istore = iter_store; istore <=iter-1; istore++) {
            file = result_dir + "/grad" + std::to_string(istore+1) + ".bin";
            this -> read_vector(grad1,file);
            file = result_dir + "/grad" + std::to_string(istore) + ".bin";
            this -> read_vector(grad0,file);
            file = result_dir + "/logmd" + std::to_string(istore+1) + ".bin";
            this -> read_vector(x1,file);
            file = result_dir + "/logmd" + std::to_string(istore) + ".bin";
            this -> read_vector(x0,file);

            // grad/logmd idff
            grad_diff = grad1 - grad0; 
            x_diff = x1 - x0;

            // coefs
            float b_tmp = (grad_diff * r_vec).sum();
            b = p[istore] * b_tmp;
            r_vec = r_vec + x_diff * (a[istore] - b);
        }

        //search direction
        direc = -r_vec;
    }

    float get_angle(const fvec &v1,const fvec &v2) {
        float l1 = std::sqrt(v1.square().sum());
        float l2 = std::sqrt(v2.square().sum());

        float a = (v1 / l1 * v2 / l2).sum();
        //printf("%f %f %f\n",l1,l2,a);
        float theta = std::acos(a) * 180. / M_PI;

        return theta;
    }

public:
    void update(int iter,T &mod,fvec &x,fvec &dsyn,
                const std::string &method = "CG")
    {
        int n = x.size(); // length 

        // compute current gradient, read previous one if required 
        fvec grad(n);
        mod.compute_grad(x,dsyn,grad);
        float chi = mod.compute_misfit(dsyn);
        int m = dsyn.size();

        // save current grad
        std::string file = result_dir + "/grad" + std::to_string(iter) + ".bin";
        this -> write_vector(grad,file);

        // get search direction
        fvec direc(n);
        if(method == "CG") {
            this -> get_cg_direction(iter,grad,direc);
        }
        else {
            this -> get_lbfgs_direction(iter,x,direc);
        }
        
        // check if the angle between 
        if(iter != ITER_START) {
            fvec invgrad = -grad;
            float theta = this -> get_angle(direc,invgrad);
            if(theta <= 90 && theta >= 0.) {
                printf("The search direction is accepted!\n");
            }
            else {
                printf("The search direction is not accepted! theta = %f\n",theta);
                printf("use negative grad instead\n");
                printf("ITER_START reset  = %d\n",iter);
                direc = invgrad;

                // clear information
                ITER_START = iter;
            }
        }

        // smoothing
        mod.smoothing(direc);

        // line search 
        printf("line searching ...\n");
        float dmax = direc.abs().maxCoeff();

        // estimate test line step
        float g1 = (grad * direc).sum();
        float step_fac; 
        if(method == "CG"){
            step_fac = -2. / g1 * chi;
        }
        else {
            step_fac = 1.;
        }
        
        // make sure the step fac is less than 3%
        if(std::abs(step_fac) * dmax > MAX_STEP || iter == ITER_START) {
            step_fac = MAX_STEP / dmax * step_fac / std::abs(step_fac);
        }
        printf("use trial step = %g, model relative variation in percent: %g\n",
                step_fac,step_fac * dmax);

        // compute misfit at the new point
        float chi1 = 0.;
        bool backtrack = true, interp_step = false;
        while(backtrack) {
            fvec x1 = x * (step_fac * direc).exp();
            fvec dsyn1;
            mod.forward(x1,dsyn1);
            chi1 = mod.compute_misfit(dsyn1);

            // check if Armijo condition is satisfied
            if(chi1 <= chi + Armijo_TOL * step_fac * g1) {
                backtrack = false;
                interp_step = false;
            }
            else if(chi1 < chi){
                backtrack = false;
                interp_step = true;
            }
            else { // chi1 > chi
                backtrack = true;
                interp_step = false;
            }

            if(!backtrack) {
                printf("line search success: two misfits = %g %g\n",chi,chi1);
            }
            else{
                step_fac *= 0.5;
                printf("line search failed: two misfits= %g %g\n",chi,chi1);
                printf("backtrack new step = %g relative variation = %g\n",step_fac,step_fac * dmax);
            }
        }

        // interp new step_fac if required 
        float alpha,chimax;
        if(interp_step) {
            float a = (chi1 - chi - g1 * step_fac) / (step_fac * step_fac);
            float b = g1, c = chi;
            alpha = -b / (2. * a);

    	    // check alpha
    	    if (alpha > step_fac || alpha < 0) {
    	    	alpha = step_fac;
    		      chimax = chi1;
    	    }
    	    else{
                chimax = a * alpha * alpha + b * alpha + c;
    	    }
        }
        else {
            alpha = step_fac;
            chimax = chi1;
        }
        
        printf("line search finished: step = %g relvar = %g misfit = %g\n",alpha,
                alpha*dmax,chimax);
        x = x * (alpha * direc).exp();
    }

};
