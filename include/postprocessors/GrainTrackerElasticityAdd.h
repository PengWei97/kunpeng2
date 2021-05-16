//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainDataTrackerAdd.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

// #include "RealVectorValue.h"

class EulerAngleProvider;

/**
 * Manage a list of elasticity tensors for the grains
 */
class GrainTrackerElasticityAdd : public GrainDataTrackerAdd<RankFourTensor,RankTwoTensor>
{
public:
  static InputParameters validParams(); // 变量声明

  GrainTrackerElasticityAdd(const InputParameters & parameters); // 构造函数声明

protected:
  RankFourTensor newGrainElasticity(unsigned int new_grain_id);

  RankTwoTensor newGrainRotation(unsigned int new_grain_id);

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// unrotated elasticity tensor
  RankFourTensor _C_ijkl;

  /// object providing the Euler angles
  const EulerAngleProvider & _euler;


};

// /**
//  * Manage a list of elasticity tensors for the grains
//  */
// class GrainTrackerElasticityAdd : public GrainDataTrackerAdd<RankFourTensor,RealVectorValue>
// {
// public:
//   static InputParameters validParams(); // 变量声明

//   GrainTrackerElasticityAdd(const InputParameters & parameters); // 构造函数声明

// protected:
//   RankFourTensor newGrainElasticity(unsigned int new_grain_id);

//   RealVectorValue newGrainRotation(unsigned int new_grain_id);

//   /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
//   const bool _random_rotations;

//   /// unrotated elasticity tensor
//   RankFourTensor _C_ijkl;

//   /// object providing the Euler angles
//   const EulerAngleProvider & _euler;


// };