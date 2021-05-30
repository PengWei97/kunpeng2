//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputerGrGrElasticEnergy.h"
// #include "RotationTensor.h"

registerMooseObject("kunpengApp", ComputerGrGrElasticEnergy);

InputParameters
ComputerGrGrElasticEnergy::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute an evolving elasticity energy coupled to a grain growth phase field model.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

ComputerGrGrElasticEnergy::ComputerGrGrElasticEnergy(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),  
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_energy_name(_base_name + "elasticity_energy"),
    _pk2_grgr(getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk2_grgr")),
    _lag_e_grgr(getMaterialPropertyByName<RankTwoTensor>(_base_name + "lage_grgr")),
    _piaolak2(declareProperty<RankTwoTensor>("piaolak2")),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _grain_tracker(getUserObject<GrainDataTrackerAdd<RankFourTensor,RealVectorValue>>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _elasticity_energy(declareProperty<Real>("elasticity_energy")),
    _D_elastic_energy(_op_num),
    _JtoeV(6.24150974e18)
{
  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    // declare elasticity tensor derivative properties
    _D_elastic_energy[op_index] = &declarePropertyDerivative<Real>(
        _elasticity_energy_name, getVar("v", op_index)->name());
  }
}

void
ComputerGrGrElasticEnergy::computerQpGrGrElasticityEnergy()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elasticity tensor
  _elasticity_energy[_qp] = 0;
  // _piaola2[_qp].zero();
  Real sum_h = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all elastic energy $ = $ 
    _elasticity_energy[_qp] += 0.5 * _pk2_grgr[_qp].doubleContraction(_lag_e_grgr[_qp])*h;
    sum_h += h;
  }
  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elasticity_energy[_qp] /=sum_h; // phi^e


  // Calculate elasticity tensor derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    (*_D_elastic_energy[op_index])[_qp] = 0;

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // 
    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;

    Real & EE_deriv = (*_D_elastic_energy[op_index])[_qp];

    EE_deriv = (0.5 * _pk2_grgr[_qp].doubleContraction(_lag_e_grgr[_qp])/sum_h - _elasticity_energy[_qp]) * dhdopi;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    // EE_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }
}
