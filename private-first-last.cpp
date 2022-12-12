//CS - 770 Big Data Analytics
//HW3	OpenMP private-firstprivate-lastprivate
//Niko Grisel Todorov

#include <iostream>
#include <omp.h>
static int NUM = 4;
int main()
{
    int a = -1, b = -1, c = -1;
    std::cout << "begin a = " << a << "\n";
    std::cout << "begin b = " << b << "\n";
    std::cout << "begin c = " << c << "\n";
    omp_set_num_threads(NUM);
#pragma omp parallel for private(a) firstprivate(b) lastprivate(c)
    for (int i=0; i<NUM; i++)
    {
        a = omp_get_thread_num();
        b--;
        c++;
        std::cout << 
            "\na = " << a << 
            " | b = " << b << 
            " | c = " << c << "\n";
    }
    std::cout << "serial a = " << a << " (private)\n";
    std::cout << "serial b = " << b << " (firstprivate)\n";
    std::cout << "serial c = " << c << " (lastprivate)\n";
    return 0;
}
