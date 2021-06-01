//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeElasticityTensorCPAdd.h"
#include "RotationTensor.h"

registerMooseObject("TensorMechanicsApp", ComputeElasticityTensorCPAdd);

InputParameters
ComputeElasticityTensorCPAdd::validParams()
{
  InputParameters params = ComputeElasticityTensor::validParams();
  params.addClassDescription("Compute an elasticity tensor for crystal plasticity.");
  // params.addParam<UserObjectName>("read_prop_user_object",
  //                                 "The ElementReadPropertyFile "
  //                                 "GeneralUserObject to read element "
  //                                 "specific property values from file");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

ComputeElasticityTensorCPAdd::ComputeElasticityTensorCPAdd(const InputParameters & parameters)
  : ComputeElasticityTensor(parameters),
    // _read_prop_user_object(isParamValid("read_prop_user_object")
    //                            ? &getUserObject<ElementPropertyReadFile>("read_prop_user_object")
    //                            : nullptr),
    // _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles")),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _grain_tracker(getUserObject<GrainDataTrackerAdd<RankFourTensor,RealVectorValue>>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _D_elastic_tensor(_op_num),
    _crysrot(declareProperty<RankTwoTensor>("crysrot")),
    _JtoeV(6.24150974e18)
    // _R(_Euler_angles)
{
  // the base class guarantees constant in time, but in this derived class the
  // tensor will rotate over time once plastic deformation sets in
  revokeGuarantee(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME);

  // the base class performs a passive rotation, but the crystal plasticity
  // materials use active rotation: recover unrotated _Cijkl here
  // _Cijkl.rotate(_R.transpose());

  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    // declare elasticity tensor derivative properties
    _D_elastic_tensor[op_index] = &declarePropertyDerivative<RankFourTensor>(
        _elasticity_tensor_name, getVar("v", op_index)->name());
  }
}

// void
// ComputeElasticityTensorCPAdd::assignEulerAngles()
// {
//   if (_read_prop_user_object)
//   {
//     _Euler_angles_mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
//     _Euler_angles_mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
//     _Euler_angles_mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);
//   }
//   else
//     _Euler_angles_mat_prop[_qp] = _Euler_angles;
// }

// void
// ComputeElasticityTensorCPAdd::computeQpElasticityTensor()
// {
//   // Properties assigned at the beginning of every call to material calculation
//   assignEulerAngles();

//   _R.update(_Euler_angles_mat_prop[_qp]);

//   _crysrot[_qp] = _R.transpose();
//   _elasticity_tensor[_qp] = _Cijkl;
//   _elasticity_tensor[_qp].rotate(_crysrot[_qp]);
// }

void
ComputeElasticityTensorCPAdd::computeQpElasticityTensor()
{
    // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elasticity tensor
  _elasticity_tensor[_qp].zero();
  _crysrot[_qp].zero();
  Real sum_h = 0.0;
  
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all rotated elasticity tensors
    _elasticity_tensor[_qp] += _grain_tracker.getDataElasticity(grain_id) * h;
    // _elasticity_tensor[_qp] += _grain_tracker.getDataElasticity(grain_id);
    _crysrot[_qp] += RotationTensor(_grain_tracker.getDataRotation(grain_id)).transpose()* h;
    // _crysrot[_qp] += RotationTensor(_grain_tracker.getDataRotation(grain_id)).transpose();
    sum_h += h;
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  // sum_h = 1;
  _elasticity_tensor[_qp] /= sum_h;
  _crysrot[_qp] /= sum_h;


  // Calculate elasticity tensor derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    (*_D_elastic_tensor[op_index])[_qp].zero();

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    RankFourTensor & C_deriv = (*_D_elastic_tensor[op_index])[_qp];

    C_deriv = (_grain_tracker.getDataElasticity(grain_id) - _elasticity_tensor[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }
}
