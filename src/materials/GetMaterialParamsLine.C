// #include "GetMaterialParamsLineLine.h"
// #include "RankTwoTensor.h"
// #include "RankFourTensor.h"

// registerMooseObject("PhaseFieldApp", GetMaterialParamsLineLine);

// InputParameters
// GetMaterialParamsLineLine::validParams()
// {
//   InputParameters params = DerivativeFunctionMaterialBase::validParams();
//   params.addClassDescription("Free energy material for the elastic energy contributions.");
//   params.addParam<std::string>("base_name",
//                                "Optional parameter that allows the user to define "
//                                "multiple mechanics material systems on the same "
//                                "block, i.e. for multiple phases");
//   params.addRequiredParam<UserObjectName>(
//       "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
//   params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
//   params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
//   params.addRequiredCoupledVarWithAutoBuild(
//       "v", "var_name_base", "op_num", "Array of coupled variables");
//   return params;
// }

// GetMaterialParamsLineLine::GetMaterialParamsLineLine(const InputParameters & parameters)
//   : DerivativeFunctionMaterialBase(parameters),
//     _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
//     _stiff_tensor_name(_base_name + "stiff_tensor"),
//     _elasticity_energy_name(_base_name + "elasticity_energy"),
//     _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain"))
//     _length_scale(getParam<Real>("length_scale")),
//     _pressure_scale(getParam<Real>("pressure_scale")),
//     _grain_tracker(getUserObject<GrainDataTracker<RankFourTensor>>("grain_tracker")),
//     _op_num(coupledComponents("v")),
//     _vals(coupledValues("v")),
//     _elasticity_energy(declareProperty<Real>("elasticity_energy")),
//     _stiff_tensor(declareProperty<Real>("stiff_tensor")),
//     _D_elastic_energy(_op_num),
//     _JtoeV(6.24150974e18),
//     _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
//     _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain"))
// {
//     // Loop over variables (ops)
//   for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
//   {
//     // declare elasticity tensor derivative properties
//     _D_elastic_energy[op_index] = &declarePropertyDerivative<Real>(
//         _elasticity_energy_name, getVar("v", op_index)->name());
//   }

// }

// void
// GetMaterialParamsLineLine::initialSetup()
// {
//   validateCoupling<RankTwoTensor>(_base_name + "elastical_strain");
// }

// Real
// GetMaterialParamsLineLine::computeF()
// {
//     // Get list of active order parameters from grain tracker
//   const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

//   // Calculate elasticity tensor
//   _stiff_tensor[_qp].zero();
//   Real sum_h = 0.0;
//   RankTwoTensor strain(_elastic_strain[_qp]);

//   for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
//   {
//     auto grain_id = op_to_grains[op_index];
//     if (grain_id == FeatureFloodCount::invalid_id)
//       continue;

//     // Interpolation factor for elasticity tensors
//     Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

//     // Sum all rotated elasticity tensors
//     _stiff_tensor[_qp] += _grain_tracker.getData(grain_id) * h;
//     sum_h += h;
//   }

//   const Real tol = 1.0e-10;
//   sum_h = std::max(sum_h, tol);
//   _stiff_tensor[_qp] /= sum_h;

//   // Calculate elasticity tensor derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
//   for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
//     (*_D_elastic_energy[op_index])[_qp].zero();

//   for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
//   {
//     auto grain_id = op_to_grains[op_index];
//     if (grain_id == FeatureFloodCount::invalid_id)
//       continue;

//     Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
//     RankFourTensor & C_deriv = (*_D_elastic_energy[op_index])[_qp];

//     C_deriv = (_grain_tracker.getData(grain_id) - _stiff_tensor[_qp]) * dhdopi / sum_h * 0.5 * strain * doubleContraction(strain);

//     // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
//     C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;

//   return _elasticity_energy[_qp];
// }



#include "GetMaterialParamsLine.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

registerMooseObject("PhaseFieldApp", GetMaterialParamsLine);

InputParameters
GetMaterialParamsLine::validParams()
{
  InputParameters params = DerivativeFunctionMaterialBase::validParams();
  params.addClassDescription("Free energy material for the elastic energy contributions.");
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

GetMaterialParamsLine::GetMaterialParamsLine(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_energy_name(_base_name + "elasticity_energy"),
    _stiff_tensor_name(_base_name + "stiff_tensor"),
    _h_name(_base_name + "h"),
    _D_h_name(_base_name + "D_h"),
    _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _grain_tracker(getUserObject<GrainDataTracker<RankFourTensor>>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _elasticity_energy(declareProperty<Real>("elasticity_energy")),
    _stiff_tensor(declareProperty<RankFourTensor>("stiff_tensor")),
    _D_elastic_energy(_op_num),
    _h(_op_num),
    _D_h(_op_num),
    _JtoeV(6.24150974e18),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain"))
{
    // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    // declare elasticity tensor derivative properties
    _D_elastic_energy[op_index] = &declarePropertyDerivative<Real>(
        _elasticity_energy_name, getVar("v", op_index)->name());
    _h[op_index] = &declarePropertyDerivative<Real>(
        _h_name, getVar("v", op_index)->name());
    _D_h[op_index] = &declarePropertyDerivative<Real>(
        _D_h_name, getVar("v", op_index)->name());
  }

}

void
GetMaterialParamsLine::initialSetup()
{
  validateCoupling<RankTwoTensor>(_base_name + "elastic_strain");
}

Real
GetMaterialParamsLine::computeF()
{
// Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elasticity tensor
  _stiff_tensor[_qp].zero();
  _elasticity_energy[_qp] = 0;

  Real sum_h = 0.0;
  RankTwoTensor strain(_elastic_strain[_qp]);
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all rotated elasticity tensors
    _stiff_tensor[_qp] += _grain_tracker.getData(grain_id) * h;
    sum_h += h;
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _stiff_tensor[_qp] /= sum_h;

  // Calculate elasticity tensor derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    (*_D_elastic_energy[op_index])[_qp] = 0;

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    Real & C_deriv = (*_D_elastic_energy[op_index])[_qp];
    Real & HH = (*_h[op_index])[_qp];
    Real & HH_deriv = (*_D_h[op_index])[_qp];

    C_deriv = ((_grain_tracker.getData(grain_id) - _stiff_tensor[_qp]) * dhdopi / sum_h * 0.5 * _elastic_strain[_qp]).doubleContraction(_elastic_strain[_qp]); 

    HH = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;
    HH_deriv = dhdopi;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;

    _elasticity_energy[_qp] = 0.5 * (_stiff_tensor[_qp] * strain).doubleContraction(strain);
  }

  return _elasticity_energy[_qp];
}


