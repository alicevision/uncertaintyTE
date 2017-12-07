#include "compute.h"
#include "auxCmd.h"
#include "FactoryIO.h"
#include "IO.h"
#include "ColmapIO.h"
#include "JacobianIO.h"
#include "JacobianComposer.h"


void getCovariances(
    cov::Options &options,
    cov::Statistic &statistic,
    ceres::CRSMatrix &jacobian,
    double* points3D,
    cov::Uncertainty& uncertainties)
{
    JacobianComposer::findPts2Fix(options, options._numPoints, points3D);
    uncertainties.init(options);
    computeCovariances(options, statistic, jacobian, &uncertainties._camerasUnc[0], &uncertainties._pointsUnc[0]);
}

