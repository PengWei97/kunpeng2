//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrainTrackerElasticityAdd.h"
#include "EulerAngleProvider.h"
#include "RotationTensor.h"

registerMooseObject("PhaseFieldApp", GrainTrackerElasticityAdd);

InputParameters
GrainTrackerElasticityAdd::validParams() // validParams的定义
{
  InputParameters params = GrainTracker::validParams(); // start with parent
  params.addParam<bool>("random_rotations",
                        true,
                        "Generate random rotations when the Euler Angle "
                        "provider runs out of data (otherwise error "
                        "out)");
  params.addRequiredParam<std::vector<Real>>("C_ijkl", "Unrotated stiffness tensor"); 
  params.addParam<MooseEnum>(
      "fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  return params;
}

GrainTrackerElasticityAdd::GrainTrackerElasticityAdd(const InputParameters & parameters) // 构造函数定义
  : GrainDataTrackerAdd<RankFourTensor,RankTwoTensor>(parameters), // 从基类中继承的变量
    _random_rotations(getParam<bool>("random_rotations")),
    _C_ijkl(getParam<std::vector<Real>>("C_ijkl"),
            getParam<MooseEnum>("fill_method").getEnum<RankFourTensor::FillMethod>()),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

RankFourTensor
GrainTrackerElasticityAdd::newGrainElasticity(unsigned int new_grain_id) // 对象定义
{
  EulerAngles angles;

  if (new_grain_id < _euler.getGrainNum())
    angles = _euler.getEulerAngles(new_grain_id);
  else
  {
    if (_random_rotations)
      angles.random();
    else
      mooseError("GrainTrackerElasticityAdd has run out of grain rotation data.");
  }

  RankFourTensor C_ijkl = _C_ijkl;
  C_ijkl.rotate(RotationTensor(RealVectorValue(angles)));

  return C_ijkl;
}

RankTwoTensor
GrainTrackerElasticityAdd::newGrainRotation(unsigned int new_grain_id) // 对象定义
{
  EulerAngles angles;

  if (new_grain_id < _euler.getGrainNum())
    angles = _euler.getEulerAngles(new_grain_id);
  else
  {
    if (_random_rotations)
      angles.random();
    else
      mooseError("GrainTrackerElasticityAdd has run out of grain rotation data.");
  }

  // // RankFourTensor C_ijkl = _C_ijkl;
  // // C_ijkl.rotate(RotationTensor(RealVectorValue(angles)));
  // RealVectorValue crysrot = RealVectorValue(angles);

  RankTwoTensor crysrot = RotationTensor(RealVectorValue(angles));

  return crysrot;
}

