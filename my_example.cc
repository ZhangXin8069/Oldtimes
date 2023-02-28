#include "my_vector.h"
#include <iostream>
#include <mpi.h>

// int main()
// {
//     MPI_Init(NULL, NULL);
//     int node_rank_, node_size_;
//     int lattice_x_(16), lattice_t_(16), lattice_spin_(2); // const lattice_spin=2
//     MPI_Comm_rank(MPI_COMM_WORLD, &node_rank_);
//     MPI_Comm_size(MPI_COMM_WORLD, &node_size_);
//     int backward_rank_ = (node_rank_ + node_size_ - 1) % node_size_;
//     int forward_rank_ = (node_rank_ + 1) % node_size_;
//     MPI_Request backward_send_req;
//     MPI_Request backward_recv_req;
//     MPI_Request forward_send_req;
//     MPI_Request forward_recv_req;
//     std::complex<double> *backward_send_vec_;
//     std::complex<double> *backward_recv_vec_;
//     std::complex<double> *forward_send_vec_;
//     std::complex<double> *forward_recv_vec_;
//     backward_send_vec_ = new std::complex<double>[lattice_t_];
//     backward_recv_vec_ = new std::complex<double>[lattice_t_];
//     forward_send_vec_ = new std::complex<double>[lattice_t_];
//     forward_recv_vec_ = new std::complex<double>[lattice_t_];

//     MPI_Isend(&backward_send_vec_, lattice_t_, MPI_DOUBLE_COMPLEX, backward_rank_, backward_rank_, MPI_COMM_WORLD, &backward_send_req);
//     MPI_Irecv(&forward_recv_vec_, lattice_t_, MPI_DOUBLE_COMPLEX,
//               forward_rank_, node_rank_, MPI_COMM_WORLD, &forward_recv_req);
//     std::cout << "backward_send:Rank# " << node_rank_ << "->Rank# " << backward_rank_ << std::endl;

//     MPI_Isend(&forward_send_vec_, lattice_t_, MPI_DOUBLE_COMPLEX,
//               forward_rank_, forward_rank_, MPI_COMM_WORLD, &forward_send_req);
//     MPI_Irecv(&backward_recv_vec_, lattice_t_, MPI_DOUBLE_COMPLEX, backward_rank_, node_rank_, MPI_COMM_WORLD, &backward_recv_req);
//     std::cout << "forward_send:Rank# " << node_rank_ << "->Rank# " << forward_rank_ << std::endl;

//     MPI_Wait(&forward_recv_req, MPI_STATUS_IGNORE);
//     std::cout << "forward_recv:Rank# " << node_rank_ << "<-Rank# " << forward_rank_ << std::endl;

//     MPI_Wait(&backward_recv_req, MPI_STATUS_IGNORE);
//     std::cout << "backward_recv:Rank# " << node_rank_ << "<-Rank# " << backward_rank_ << std::endl;

//     ComplexVector vector_0(4);
//     vector_0.clean_rand();
//     ComplexVector vector_1(vector_0);
//     std::cout << vector_0.norm2();
//     std::cout << vector_1.norm2X();
//     delete[] backward_send_vec_;
//     delete[] backward_recv_vec_;
//     delete[] forward_send_vec_;
//     delete[] forward_recv_vec_;
//     MPI_Finalize();
//     return 0;
// }

#include <iostream>
#include <mpi.h>
#include <thread>

void SendData(int rank, int size)
{
    int send_data = rank;
    MPI_Request request;
    MPI_Isend(&send_data, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
}

void RecvData(int rank, int size)
{
    int recv_data;
    MPI_Request request;
    MPI_Irecv(&recv_data, 1, MPI_INT, (rank + size - 1) % size, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    std::cout << "Process " << rank << " received " << recv_data << " from process " << (rank + size - 1) % size << std::endl;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int node_rank_, node_size_;
    int lattice_x_(16), lattice_t_(16), lattice_spin_(2); // const lattice_spin=2
    MPI_Comm_rank(MPI_COMM_WORLD, &node_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &node_size_);
    int backward_rank_ = (node_rank_ + node_size_ - 1) % node_size_;
    int forward_rank_ = (node_rank_ + 1) % node_size_;
    std::complex<double> *backward_send_vec_;
    std::complex<double> *backward_recv_vec_;
    std::complex<double> *forward_send_vec_;
    std::complex<double> *forward_recv_vec_;
    backward_send_vec_ = new std::complex<double>[lattice_t_];
    backward_recv_vec_ = new std::complex<double>[lattice_t_];
    forward_send_vec_ = new std::complex<double>[lattice_t_];
    forward_recv_vec_ = new std::complex<double>[lattice_t_];
    // std::thread send_thread(SendData, rank, size);
    // std::thread recv_thread(RecvData, rank, size);

    // send_thread.join();
    // recv_thread.join();
    for (int i(0); i < 3; i++)
    {
        SendData(node_rank_, node_size_);
        RecvData(node_rank_, node_size_);
    }
    delete[] backward_send_vec_;
    delete[] backward_recv_vec_;
    delete[] forward_send_vec_;
    delete[] forward_recv_vec_;
    MPI_Finalize();
    return 0;
}
