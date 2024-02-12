void allocate_tasks(int ntasks,int nprocs,int myrank, int &start,int &end)
{
    // allocate jobs to each rank
    int sub_n = ntasks / nprocs;
    int num_larger_procs = ntasks - nprocs * sub_n;
    if (myrank < num_larger_procs){ sub_n = sub_n + 1;
        start = 0 + myrank * sub_n;
    }
    else if (sub_n > 0){ 
        start = 0 + num_larger_procs + myrank * sub_n;
    }
    else { // this process has only zero elements
        start = -1;
        sub_n = 0;
    }
    end = start + sub_n - 1;
}