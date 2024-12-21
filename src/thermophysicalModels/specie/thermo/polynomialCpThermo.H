#ifndef polynomialCpThermo_H
#define polynomialCpThermo_H

#include "thermo.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace species
{

template<class Thermo, template<class> class Type>
class polynomialCpThermo
:
    public Thermo,
    public Type<polynomialCpThermo<Thermo, Type>>
{
    // Private Data

        Polynomial<8> CpCoeffs_;


public:

    // Constructors

        polynomialCpThermo(const dictionary& dict)
        :
            Thermo(dict),
            CpCoeffs_(dict.subDict("CpCoeffs"))
        {}


    // Member Functions

        scalar Cp(const scalar p, const scalar T) const
        {
            return CpCoeffs_.value(T);
        }

        scalar Hs(const scalar p, const scalar T) const
        {
            return CpCoeffs_.integral(T) - CpCoeffs_.integral(T0_);
        }

        scalar Ha(const scalar p, const scalar T) const
        {
            return Hs(p, T) + Hf();
        }

        scalar Es(const scalar p, const scalar T) const
        {
            return Hs(p, T) - p/(rho(p, T)*Cv(p, T));
        }

        scalar Ea(const scalar p, const scalar T) const
        {
            return Es(p, T) + Hf();
        }

        scalar S(const scalar p, const scalar T) const
        {
            return CpCoeffs_.integral(T)/T - CpCoeffs_.integral(T0_)/T0_;
        }

        scalar HE(const scalar p, const scalar T) const
        {
            return Ha(p, T);
        }

        scalar G(const scalar p, const scalar T) const
        {
            return Ha(p, T) - T*S(p, T);
        }

        scalar A(const scalar p, const scalar T) const
        {
            return Ea(p, T) - T*S(p, T);
        }

        scalar cp(const scalar p, const scalar T) const
        {
            return Cp(p, T)*this->W();
        }

        scalar ha(const scalar p, const scalar T) const
        {
            return Ha(p, T)*this->W();
        }

        scalar hs(const scalar p, const scalar T) const
        {
            return Hs(p, T)*this->W();
        }

        scalar hc() const
        {
            return Hf()*this->W();
        }

        scalar s(const scalar p, const scalar T) const
        {
            return S(p, T)*this->W();
        }

        scalar he(const scalar p, const scalar T) const
        {
            return HE(p, T)*this->W();
        }

        scalar cv(const scalar p, const scalar T) const
        {
            return Cp(p, T) - this->R();
        }

        scalar es(const scalar p, const scalar T) const
        {
            return Es(p, T)*this->W();
        }

        scalar ea(const scalar p, const scalar T) const
        {
            return Ea(p, T)*this->W();
        }

        scalar g(const scalar p, const scalar T) const
        {
            return G(p, T)*this->W();
        }

        scalar a(const scalar p, const scalar T) const
        {
            return A(p, T)*this->W();
        }

        scalar THE(const scalar H, const scalar p, const scalar T0) const
        {
            return CpCoeffs_.T(H - Hf() + CpCoeffs_.integral(T0), T0);
        }

        scalar THs(const scalar Hs, const scalar p, const scalar T0) const
        {
            return CpCoeffs_.T(Hs + CpCoeffs_.integral(T0), T0);
        }

        scalar THa(const scalar H, const scalar p, const scalar T0) const
        {
            return THE(H, p, T0);
        }

        scalar TEs(const scalar E, const scalar p, const scalar T0) const
        {
            return CpCoeffs_.T(E + p/(rho(p, T0)*Cv(p, T0)) + CpCoeffs_.integral(T0), T0);
        }

        scalar TEa(const scalar E, const scalar p, const scalar T0) const
        {
            return TEs(E, p, T0);
        }

        scalar dKcdTbyKc(const scalar p, const scalar T) const
        {
            return 0;
        }

        scalar dcpdT(const scalar p, const scalar T) const
        {
            return CpCoeffs_.derivative(T);
        }

        void write(Ostream& os) const
        {
            Thermo::write(os);
            CpCoeffs_.writeEntry("CpCoeffs", os);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
