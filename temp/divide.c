#include <stddef.h>

typedef func;

enum { SIMPLEX, HYPERRECTANGLE };

void divide(
    int dimens, int nrvert, int maxsub, bool uniform_subdiv, int numfun, double* verold, int* infold, double* rinfol, func integrand, int outsub, int num, double* vernew, int* infnew, double* rinfne, int ifail)
{
    if (maxsub == 1)
    {
        memcpy(vernew, verold, dimens*nrvert*sizeof(double));
        memcpy(infnew, infold, sizeof(int));
        memcpy(rinfne, rinfol, sizeof(double));
        outsub = 1;
        num = 0;
        ifail = 0;
        return;
    }

    int geometry = infold[0];
    switch (geometry)
    {
        case SIMPLEX:
        {
            if (dimens == 1)
            {
                memcpy(vernew, verold, dimens*nrvert*sizeof(double));
                memcpy(&vernew[dimens*nrvert], verold, dimens*nrvert*sizeof(double));
            }
            break;
        }
        case HYPERRECTANGLE: break;
        default: return;
    }

    double volume = rinfol[0]/outsub;

    for (size_t i = 0; i < outsub; ++i)
    {
        infnew[6*i] = geometry;
        infnew[6*i + 1] = infold[1] + 1;
        infnew[6*i + 2] = infold[2];
        infnew[6*i + 3] = 0;
        infnew[6*i + 4] = 0;
        rinfne[6*i] = volume;
        memcpy(&infnew[6*i + 5], &infold[5], sizeof(int));
        memcpy(&rinfne[6*i + 5], &rinfol[5], sizeof(double));
    }

    ifail = 0;
    return;
}

typedef vertex;

void divide_simplex(
    int dimens, int nf, int maxsub, double* verold, func integrand, int funcls, int outsub, double* vernew)
{
    double center = center_of(verold);
    double* frthdf;
    double* ewidth;
    for (size_t l = 0; l < dimens; ++l)
    {
        for (size_t k = l + 1; k < dimens; ++k)
        {
            vertex h = vertices[k] - vertices[l];
        }
    }
}