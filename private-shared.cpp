//CS - 770 Big Data Analytics
//HW3	OpenMP private
//Niko Grisel Todorov

#include <iostream>
#include <omp.h>
static int NUM = 4;
int main()
{
    int index = -1, tid = -1;
    std::cout << "begin index = " << index << "\n";
    std::cout << "begin tid   = " << tid << "\n";
    omp_set_num_threads(NUM);
#pragma omp parallel for private(tid) shared(index)
    for (int i=0; i<NUM; i++)
    {
        tid = omp_get_thread_num();
        index = tid;
        std::cout << "\nthread = " << tid << "| index = " << index << "\n";
    }
    std::cout << "\nserial index = " << index << " (private)\n";
    std::cout << "serial tid   = " << tid << " (public)\n";
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
