/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::species::phaseChangeCpThermo

Description
    Thermodynamics package extending polynomialCpThermo to include
    latent heat effects during phase changes.

SourceFiles
    phaseChangeCpThermo.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeCpThermo_H
#define phaseChangeCpThermo_H

#include "thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace species
{

/*---------------------------------------------------------------------------*\
                      Class phaseChangeCpThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class phaseChangeCpThermo
:
    public Thermo
{
    // Private member typedefs
    typedef phaseChangeCpThermo<Thermo> phaseChangeCpThermoType;

    // Private data

        //- Phase change temperatures
        List<scalar> phaseTemp_;

        //- Latent heats for each phase change
        List<scalar> latentHeat_;

        //- Temperature ranges for phase changes
        List<scalar> deltaT_;


public:

    //- Runtime type information
    TypeName("phaseChangeCpThermo");

    // Constructors

        //- Construct from components
        inline phaseChangeCpThermo(const Thermo& t)
        :
            Thermo(t),
            phaseTemp_(),
            latentHeat_(),
            deltaT_()
        {}

        //- Construct as copy
        inline phaseChangeCpThermo(const phaseChangeCpThermo& pt)
        :
            Thermo(pt),
            phaseTemp_(pt.phaseTemp_),
            latentHeat_(pt.latentHeat_),
            deltaT_(pt.deltaT_)
        {}

        //- Construct from dictionary
        phaseChangeCpThermo(const dictionary& dict)
        :
            Thermo(dict),
            phaseTemp_(dict.lookup("phaseTemp")),
            latentHeat_(dict.lookup("latentHeat")),
            deltaT_(dict.lookupOrDefault<List<scalar>>("deltaT", List<scalar>(phaseTemp_.size(), 1.0)))
        {
            // Check sizes match
            if (latentHeat_.size() != phaseTemp_.size())
            {
                FatalIOErrorInFunction(dict)
                    << "Size of latentHeat " << latentHeat_.size()
                    << " does not match size of phaseTemp " << phaseTemp_.size()
                    << exit(FatalIOError);
            }

            if (deltaT_.size() != phaseTemp_.size())
            {
                FatalIOErrorInFunction(dict)
                    << "Size of deltaT " << deltaT_.size()
                    << " does not match size of phaseTemp " << phaseTemp_.size()
                    << exit(FatalIOError);
            }

            // Check temperatures are in ascending order
            for (label i = 1; i < phaseTemp_.size(); i++)
            {
                if (phaseTemp_[i] <= phaseTemp_[i-1])
                {
                    FatalIOErrorInFunction(dict)
                        << "Phase temperatures must be in ascending order. "
                        << "Temperature " << phaseTemp_[i]
                        << " at position " << i
                        << " is not greater than previous temperature "
                        << phaseTemp_[i-1]
                        << exit(FatalIOError);
                }
            }

            // Check all latent heats are positive
            forAll(latentHeat_, i)
            {
                if (latentHeat_[i] <= 0)
                {
                    FatalIOErrorInFunction(dict)
                        << "Latent heat must be positive. "
                        << "Value " << latentHeat_[i]
                        << " at position " << i
                        << " is not positive"
                        << exit(FatalIOError);
                }
            }

            // Check all deltaT are positive
            forAll(deltaT_, i)
            {
                if (deltaT_[i] <= 0)
                {
                    FatalIOErrorInFunction(dict)
                        << "Temperature range must be positive. "
                        << "Value " << deltaT_[i]
                        << " at position " << i
                        << " is not positive"
                        << exit(FatalIOError);
                }
            }
        }

        //- Clone
        inline autoPtr<phaseChangeCpThermo> clone() const
        {
            return autoPtr<phaseChangeCpThermo>
            (
                new phaseChangeCpThermo(*this)
            );
        }


    //- Destructor
    virtual ~phaseChangeCpThermo()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "phaseChangeCpThermo<" + Thermo::typeName() + '>';
        }

        //- Return the name of the thermo physics
        word thermoName() const
        {
            return typeName();
        }

        //- Return the molecular weight [kg/kmol]
        scalar W() const
        {
            return Thermo::W();
        }

        //- Return the formation enthalpy [J/kmol]
        scalar Hf() const
        {
            return Thermo::Hf();
        }

        //- Return entropy [J/(kmol K)]
        scalar S(const scalar p, const scalar T) const
        {
            return Thermo::S(p, T);
        }

        //- Return compressibility rho/p [s^2/m^2]
        scalar psi(const scalar p, const scalar T) const
        {
            return Thermo::psi(p, T);
        }

        //- Return compression factor []
        scalar Z(const scalar p, const scalar T) const
        {
            return Thermo::Z(p, T);
        }

        //- Return (Cp - Cv)/Cp = (gamma - 1)/gamma
        scalar CpByCv(const scalar p, const scalar T) const
        {
            return Thermo::CpByCv(p, T);
        }

        //- Heat capacity at constant pressure [J/(kmol K)]
        virtual scalar Cp(const scalar p, const scalar T) const
        {
            const scalar baseCp = Thermo::Cp(p, T);

            // Add latent heat contributions using Gaussian distributions
            // around each phase change point to smooth the transitions
            scalar additionalCp = 0.0;
            
            forAll(phaseTemp_, i)
            {
                if (mag(T - phaseTemp_[i]) < deltaT_[i])
                {
                    scalar alpha = 1.0/(deltaT_[i]*sqrt(2.0*constant::mathematical::pi));
                    scalar beta = -0.5*sqr((T - phaseTemp_[i])/deltaT_[i]);
                    
                    additionalCp += latentHeat_[i]*alpha*exp(beta);
                }
            }
            
            return baseCp + additionalCp;
        }

        //- Sensible enthalpy [J/kmol]
        scalar Hs(const scalar p, const scalar T) const
        {
            const scalar baseHs = Thermo::Hs(p, T);

            // Add latent heat contributions
            scalar additionalH = 0.0;
            
            forAll(phaseTemp_, i)
            {
                if (T > phaseTemp_[i] + deltaT_[i])
                {
                    additionalH += latentHeat_[i];
                }
                else if (mag(T - phaseTemp_[i]) <= deltaT_[i])
                {
                    // Smooth transition using error function
                    scalar x = (T - phaseTemp_[i])/(deltaT_[i]*sqrt(2.0));
                    additionalH += 0.5*latentHeat_[i]*(1.0 + erf(x));
                }
            }
            
            return baseHs + additionalH;
        }

        //- Temperature from enthalpy
        scalar THE
        (
            const scalar H,
            const scalar p,
            const scalar T0
        ) const
        {
            // Use Newton-Raphson iteration to find temperature
            scalar T = T0;
            scalar Hn = Ha(p, T);
            scalar dT = T;
            label iter = 0;
            const label maxIter = 100;
            const scalar tolerance = 1e-6;

            while (mag(dT) > tolerance && iter < maxIter)
            {
                const scalar Cp = this->Cp(p, T);
                dT = (H - Hn)/Cp;
                T += dT;
                Hn = Ha(p, T);
                iter++;
            }

            if (iter == maxIter)
            {
                FatalErrorInFunction
                    << "Failed to converge temperature from enthalpy value " 
                    << H << " in " << maxIter << " iterations"
                    << abort(FatalError);
            }

            return T;
        }

        //- Override other temperature calculations to use THE
        scalar THs
        (
            const scalar Hs,
            const scalar p,
            const scalar T0
        ) const
        {
            return THE(Hs + this->Hf(), p, T0);
        }

        scalar THa
        (
            const scalar Ha,
            const scalar p,
            const scalar T0
        ) const
        {
            return THE(Ha, p, T0);
        }

        //- Write to Ostream
        void write(Ostream& os) const
        {
            Thermo::write(os);
            
            os.writeKeyword("phaseTemp") << phaseTemp_ << token::END_STATEMENT << nl;
            os.writeKeyword("latentHeat") << latentHeat_ << token::END_STATEMENT << nl;
            os.writeKeyword("deltaT") << deltaT_ << token::END_STATEMENT << nl;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
