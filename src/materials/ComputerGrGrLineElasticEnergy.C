#include "ComputerGrGrLineElasticEnergy.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

registerMooseObject("PhaseFieldApp", ComputerGrGrLineElasticEnergy);

InputParameters
ComputerGrGrLineElasticEnergy::validParams()
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

ComputerGrGrLineElasticEnergy::ComputerGrGrLineElasticEnergy(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_energy_name(_base_name + "elasticity_energy"),
    _elastic_energy_name(_base_name + "elasticity_energy"),
    _h_name(_base_name + "h"),
    _D_h_name(_base_name + "D_h"),
    _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    // _pk2_grgr(getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk2_grgr")),
    // _lag_e_grgr(getMaterialPropertyByName<RankTwoTensor>(_base_name + "lage_grgr")),
    // _piaolak2(declareProperty<RankTwoTensor>("piaolak2")),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _grain_tracker(getUserObject<GrainDataTracker<RankFourTensor>>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _elasticity_energy(declareProperty<Real>("elasticity_energy")),
    _elastic_energy(declareProperty<Real>("elastic_energy")),
    _D_elastic_energy(_op_num),
    _h(_op_num),
    _D_h(_op_num),
    _JtoeV(6.24150974e18)
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
ComputerGrGrLineElasticEnergy::initialSetup()
{
  validateCoupling<RankTwoTensor>(_base_name + "elastic_strain");
}

Real
ComputerGrGrLineElasticEnergy::computeF()
{
    // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elasticity tensor
  _elasticity_energy[_qp] = 0;
  _elastic_energy[_qp] = 0;
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
    _elasticity_energy[_qp] += (0.5 * _elastic_strain[_qp].doubleContraction(_stress[_qp]))*h;
    sum_h += h;
  }
  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elasticity_energy[_qp] /= sum_h; // phi^e


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
    Real & HH = (*_h[op_index])[_qp];
    Real & HH_deriv = (*_D_h[op_index])[_qp];

    _elastic_energy[_qp] = (0.5 * _elastic_strain[_qp].doubleContraction(_stress[_qp]));

    EE_deriv = ((_elastic_energy[_qp] - _elasticity_energy[_qp])/sum_h) * dhdopi;

    HH = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;
    HH_deriv = dhdopi;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    // EE_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }

  return _elasticity_energy[_qp];
}