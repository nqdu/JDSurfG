template<typename T>
class csr_matrix{
    private:
    int m,n;
    bool flag;

    public:
    restrict T *value;
    restrict int *indx;
    restrict int *indptr;
    int nonzeros;

    ~csr_matrix(){
        if(flag){
            delete[] value;
            delete[] indx;
            delete[] indptr;
        }
    }

    void mulvec(T *x,T *y,int m,int n){
        for(int i=0;i<m;i++){
            T s = 0.0;
            for(int j=indptr[i];j<indptr[i+1];j++){
                s += value[j];
            }
            y[i] = s;
        }
    }

};