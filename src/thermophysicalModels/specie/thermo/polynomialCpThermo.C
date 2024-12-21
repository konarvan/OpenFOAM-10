#include "polynomialCpThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace species
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addToRunTimeSelectionTable
(
    thermo,
    polynomialCpThermo<hePsiThermo<pureMixture<constIsoThermo<specie>>, sensibleEnthalpy>>,
    hePsiThermo
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// ************************************************************************* //
