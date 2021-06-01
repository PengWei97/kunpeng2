//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGrGrElasticEnergyLine.h"

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("PhaseFieldApp", ACGrGrElasticEnergyLine);

InputParameters
ACGrGrElasticEnergyLine::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Adds elastic energy contribution to the Allen-Cahn equation");
  params.addRequiredParam<MaterialPropertyName>(
      "D_tensor_name", "The elastic tensor derivative for the specific order parameter");
  params.addRequiredParam<MaterialPropertyName>(
      "D_elastic_energy_name", "The elastic tensor derivative for the specific order parameter");
  return params;
}

ACGrGrElasticEnergyLine::ACGrGrElasticEnergyLine(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    // _D_elastic_tensor(getMaterialProperty<RankFourTensor>("D_tensor_name")),
    _D_elastic_energy(getMaterialProperty<Real>("D_elastic_energy_name")),
    // _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
    // _lag_e_grgr(getMaterialPropertyByName<RankTwoTensor>("lage_grgr")),
     _D_elastic_tensor(getMaterialProperty<RankFourTensor>("D_tensor_name"))
{
}

Real
ACGrGrElasticEnergyLine::computeDFDOP(PFFunctionType type)
{
  // Access the heterogeneous strain calculated by the Solid Mechanics kernels
  // RankTwoTensor strain(_elastic_strain[_qp]);
  // RankTwoTensor strain(_lag_e_grgr[_qp]);

  // Compute the partial derivative of the stress wrt the order parameter
  // RankTwoTensor D_stress = _D_elastic_tensor[_qp] * strain;

  switch (type)
  {
    case Residual:
      // return 0.5 *
      //        D_stress.doubleContraction(strain); // Compute the deformation energy driving force
      return _D_elastic_energy[_qp];
      // return 0.5 * D_stress.doubleContraction(strain);


    case Jacobian:
      return 0.0;
  }

  mooseError("Invalid type passed in");
}

