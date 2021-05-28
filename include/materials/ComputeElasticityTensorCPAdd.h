//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeElasticityTensor.h"
// #include "ElementPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"
#include "GrainDataTrackerAdd.h"

/**
 * ComputeElasticityTensorCPAdd defines an elasticity tensor material object for crystal plasticity.
 */
class ComputeElasticityTensorCPAdd : public ComputeElasticityTensor
{
public:
  static InputParameters validParams();

  ComputeElasticityTensorCPAdd(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  // virtual void computeQpElasticityEnergy();

  // virtual void assignEulerAngles();

  std::string _elasticity_energy_name;
  
  Real _length_scale;
  Real _pressure_scale;

  const GrainDataTrackerAdd<RankFourTensor,RealVectorValue> & _grain_tracker;

  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;

  std::vector<MaterialProperty<RankFourTensor> *> _D_elastic_tensor;

  std::vector<MaterialProperty<Real> *> _D_elastic_energy;

  /**
   * Element property read user object
   * Presently used to read Euler angles -  see test
   */
  // const ElementPropertyReadFile * const _read_prop_user_object;

  // MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;

  /// Crystal Rotation Matrix

  MaterialProperty<RankTwoTensor> & _crysrot;
  MaterialProperty<Real> & _elasticity_energy;


  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<RankTwoTensor> & _lag_e;

  /// Rotation matrix
  // RotationTensor _R;
  const Real _JtoeV;
  
};



