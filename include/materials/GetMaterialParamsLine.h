// #pragma once

// #include "DerivativeFunctionMaterialBase.h"
// #include "ADRankTwoTensorForward.h"
// #include "ADRankFourTensorForward.h"
// #include "GrainDataTracker.h"

// /**
//  * Material class to compute the elastic free energy and its derivatives
//  */
// class GetMaterialParamsLineLine : public DerivativeFunctionMaterialBase
// {
// public:
//   static InputParameters validParams();

//   GetMaterialParamsLineLine(const InputParameters & parameters);

//   virtual void initialSetup() override;

// protected:
//   // virtual void computerQpGrGrElasticityEnergy();
//   virtual Real computeF() override;
//   // virtual Real computeDF(unsigned int op_num) override;

//   const std::string _base_name;

//   std::string _stiff_tensor_name;

//   std::string _elasticity_energy_name;

//   const MaterialProperty<RankTwoTensor> & _elastic_strain;

//   Real _length_scale;
//   Real _pressure_scale;

//   const GrainDataTracker<RankFourTensor> & _grain_tracker;

//   /// Number of order parameters
//   const unsigned int _op_num;

//   /// Order parameters
//   const std::vector<const VariableValue *> _vals;

//   MaterialProperty<Real> & _elasticity_energy;

//   MaterialProperty<RankFourTensor> & _stiff_tensor;

//   // MaterialProperty<RankTwoTensor> & _piaola2;


//   std::vector<MaterialProperty<Real> *> _D_elastic_energy;
//   // std::vector<MaterialProperty<RankFourTensor> *> _D_elastic_tensor;


//   /// Conversion factor from J to eV
//   const Real _JtoeV;

//   /// Stress tensor
//   const MaterialProperty<RankTwoTensor> & _stress;

//   ///@{ Strain and derivatives
//   const MaterialProperty<RankTwoTensor> & _strain;
//   ///@}
// };


#pragma once

#include "DerivativeFunctionMaterialBase.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"
#include "GrainDataTracker.h"

/**
 * Material class to compute the elastic free energy and its derivatives
 */
class GetMaterialParamsLine : public DerivativeFunctionMaterialBase
{
public:
  static InputParameters validParams();

  GetMaterialParamsLine(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  // virtual void computerQpGrGrElasticityEnergy();
  virtual Real computeF() override;
  // virtual Real computeDF(unsigned int op_num) override;

  const std::string _base_name;

  std::string _elasticity_energy_name;

  std::string _stiff_tensor_name;

  std::string _h_name;

  std::string _D_h_name;

  const MaterialProperty<RankTwoTensor> & _elastic_strain;

  Real _length_scale;
  Real _pressure_scale;

  const GrainDataTracker<RankFourTensor> & _grain_tracker;

  /// Number of order parameters
  const unsigned int _op_num;

  /// Order parameters
  const std::vector<const VariableValue *> _vals;

  MaterialProperty<Real> & _elasticity_energy;

  MaterialProperty<RankFourTensor> & _stiff_tensor;

  // MaterialProperty<RankTwoTensor> & _piaola2;


  std::vector<MaterialProperty<Real> *> _D_elastic_energy;


  std::vector<MaterialProperty<Real> *> _h;

  std::vector<MaterialProperty<Real> *> _D_h;

  /// Conversion factor from J to eV
  const Real _JtoeV;

  /// Stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;

  ///@{ Strain and derivatives
  const MaterialProperty<RankTwoTensor> & _strain;
  ///@}
};

