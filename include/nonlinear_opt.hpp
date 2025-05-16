#ifndef JDSURFG_NONLINEAR_OPT_H_
#define JDSURFG_NONLINEAR_OPT_H_

#include "numerical.hpp"
#include <fstream>


/**
 * nonlinear operator
 * @tparam mod forward model operator, 3 functions should have \n 
 *   1. void mod.compute_grad(const fvec &x,fvec &dsyn,fvec &grad) \n 
 *   2. float mod.compute_misfit(const fvec &dsyn); \n
 *   3. void mod.forward(const fvec &x,fvec &dsyn,)
*/
template <class T> 
class NonlinOPT {

private:
    std::string FLAG;
    int ITER_LS = 0;
    const float MSTORE = 5;

    // line search parameters
    const float WOLFE_M1 = 1.0E-4;
    const float WOLFE_BFGS_M2 = 0.9;
    const float WOLFE_CG_M2 = 0.1;
    float STEP_L, STEP_R;
    const float FACTOR = 10.;

public:
    std::string result_dir;
    float MAX_STEP;
    int ITER_START;
    int ITER;
    float STEP_FAC;

    NonlinOPT(const std::string &outdir,float max_step,int iter_start,int iter_current)
    {
        this -> result_dir = outdir;;
        this -> MAX_STEP = max_step;
        this -> ITER_START = iter_start;
        this -> ITER_LS = 0;
        this -> FLAG = "INIT";
        this -> STEP_FAC = -1.;
        this -> ITER = iter_current;
        STEP_L = 0.;
        STEP_R = 0.;
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

    void get_lbfgs_direction(int iter,fvec &direc) const 
    {
        std::string file;

        // save current log model
        int n = direc.size();
        fvec x1(n),x0(n);
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
        int iter_store = iter - MSTORE; 
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

    fvec get_next_model(T &mod,const fvec &grad,const fvec &x,float chi, 
                        const std::string &method, fvec &direc)
    {
        int iter = ITER;
        if(method == "CG") {
            this -> get_cg_direction(iter,grad,direc);
        }
        else {
            this -> get_lbfgs_direction(iter,direc);
        }

        // check if the angle between 
        if(iter != ITER_START) {
            fvec invgrad = -grad;
            float theta = this -> get_angle(direc,invgrad);
            if(theta > 92 ) {
                printf("The search direction is not accepted! theta = %f\n",theta);
                printf("use negative grad instead\n");
                printf("ITER_START reset  = %d\n",iter);
                direc = invgrad;

                // clear information
                ITER_START = iter;
            }
        }

        // smooth the search direction
        mod.smoothing(direc);

        // max of abs(direc)
        float dmax = direc.abs().maxCoeff();

        // initialize step_fac 
        float step_fac = STEP_FAC;
        if((iter == ITER_START) && (ITER_LS == 0)) {
            step_fac = -1.;
        }
        if((iter == ITER_START + 1) && (ITER_LS == 0) ) {
            if(method == "CG") {
                float g1 = (grad * direc).sum();
                step_fac = -2. / g1 * chi;
            }
            else {
                step_fac = 1.;
            }
        }
        
        // make sure the step fac is less than 3%
        if(step_fac < 0 || step_fac * dmax > MAX_STEP ) {
            step_fac = MAX_STEP / dmax;
        }
        printf("use trial step = %g, dmax = %g, relative change = %g\n",
                step_fac,dmax,step_fac * dmax);

        // reset STEP_FAC
        STEP_FAC = step_fac;

        // get next line search model
        fvec x1 = x * (STEP_FAC * direc).exp();

        // update flags
        if(FLAG != "LS") {
            FLAG = "LS";
            ITER_LS = 0;
        }

        return x1;
    }

public:
    void update(T &mod,fvec &x,fvec &dsyn,
                const std::string &method = "CG")
    {
        int n = x.size(); // length 
        fvec grad(n),direc(n); // gradient and search direction
        fvec x1;

        // run optimization
        if (FLAG == "INIT") {
            // compute gradient 
            mod.compute_grad(x,dsyn,grad);

            // save current grad,logmod and dsyn
            std::string file = result_dir + "/grad" + std::to_string(ITER) + ".bin";
            this -> write_vector(grad,file);

            // model
            file = result_dir + "/logmd" + std::to_string(ITER) + ".bin";
            x1 = x.log();
            this -> write_vector(x1,file);

            // read current dsyn
            file = result_dir + "/syn" + std::to_string(ITER) + ".bin";
            this -> read_vector(dsyn,file);
        }
        else if(FLAG == "GRAD") {
            // read current grad
            std::string file = result_dir + "/grad" + std::to_string(ITER) + ".bin";
            this -> read_vector(grad,file);

            // read current dsyn
            file = result_dir + "/syn" + std::to_string(ITER) + ".bin";
            this -> read_vector(dsyn,file);
        }

        // compute misift
        float fcost = mod.compute_misfit(dsyn);

        // get next model
        x1 = get_next_model(mod,grad,x,fcost,method,direc);

        // now do line search
        while (FLAG  == "LS") {
            printf("Line search begin, ITER_LS =  %d\n",ITER_LS);

            // compute gradient 
            fvec grad_next, dsyn_next;
            mod.compute_grad(x1,dsyn_next,grad_next);
            
            // reset LS parameters if required
            if (ITER_LS == 0) {
                STEP_L = 0.;
                STEP_R = 0.;
            }

            // get misifit for new model
            float fcost1 = mod.compute_misfit(dsyn_next);

            // compute inner prodcut
            float q = (grad * direc).sum();
            float q1 = (grad_next * direc).sum();

            // check WOLFE condition
            float m1 = WOLFE_M1;
            float m2 = WOLFE_BFGS_M2;
            if(method == "CG") m2 = WOLFE_CG_M2;
            bool cond1 = fcost1 <= (fcost + m1 * STEP_FAC * q);
            bool cond2 = q1 >= m2 * q;

            if (cond1 && cond2) {
                printf("Wolfe condition is satisfied!\n");
                printf("misfits current/next = %g %g\n",fcost,fcost1);
                FLAG = "GRAD";
                ITER += 1;
                ITER_LS = 0;

                // save gradient and data of next model 
                std::string file = result_dir + "/grad" + std::to_string(ITER) + ".bin";
                this -> write_vector(grad_next,file);
                file = result_dir + "/syn" + std::to_string(ITER) + ".bin";
                this -> write_vector(dsyn_next,file);

                // copy x1 to updated model
                x = x1;

                // save model
                file = result_dir + "/logmd" + std::to_string(ITER) + ".bin";
                x1 = x.log();
                this -> write_vector(x1,file);

                // end line search
                break;
            }
            else if(!cond1) {
                STEP_R = STEP_FAC;
                STEP_FAC = 0.5 * (STEP_R + STEP_L);
                ITER_LS += 1;

                printf("First Wolfe condition is not satisfied\n");
                printf("decrease step_fac to %g\n",STEP_FAC);
            }
            else if(!cond2) {
                STEP_L = STEP_FAC;
                if(STEP_R != 0.) {
                    STEP_FAC = 0.5 * (STEP_L + STEP_R);
                }
                else {
                    STEP_FAC *= FACTOR;
                }
                ITER_LS += 1;

                printf("Second Wolfe condition is not satisfied\n");
                printf("increase step_fac to %g\n",STEP_FAC);
            }

            // get next model if required
            x1 = get_next_model(mod,grad,x,fcost,method,direc);
        }
    }
};

#endif // end JDSURFG_NONLINEAR_OPT_H_
