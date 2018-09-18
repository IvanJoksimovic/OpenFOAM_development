/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
Class
    readProfileFromFile

Description

SourceFiles
    readProfileFromFile.C
\*---------------------------------------------------------------------------*/
#ifndef readProfileFromFileFvPatchField_H
#define readProfileFromFileFvPatchField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include <iostream>
#include <fstream>
#include <vector>

#include <stdio.h>  // Required to acces current working directory
#include <unistd.h> // Required to acces current working directory

namespace Foam
{
/*---------------------------------------------------------------------------*\
    Declaration of class: readProfileFromFileFvPatchField
\*---------------------------------------------------------------------------*/
class readProfileFromFileFvPatchField
:
public fixedValueFvPatchVectorField
{
    private :

        // Direction of the flow :
        vector flowNormal_;

        // Direction of the wall normal coordinate :
        vector wallNormal_;

        // Read from specific file :
        string readFromFile_;

        // Check if the profile is fully developed
        bool isFullyDeveloped_;

        // Tutbulent boundary layer thickness if the flow is not fully developed:
        scalar delta_;

    public:

    //- Runtime type information
        TypeName("readProfileFromFile");


        //- Construct from patch and internal field
        readProfileFromFileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        readProfileFromFileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a developedTurbulentVelocityProfileFvPatchVectorField onto a new patch
        readProfileFromFileFvPatchField
        (
            const readProfileFromFileFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new readProfileFromFileFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        readProfileFromFileFvPatchField
        (
            const readProfileFromFileFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new readProfileFromFileFvPatchField(*this, iF)
            );
        }

    // Member functions

        //- Return flow direction
        vector& flowNormal()
        {
            return flowNormal_;
        }
        //- Return wall normal coordinate
        vector& wallNormal()
        {
            return wallNormal_;
        }

        //-Return name of the file to read from
        string& readFromFile()
        {
            return readFromFile_;
        }
        //-Return the true value if the flow is fully developed
        bool& isFullyDeveloped()
        {
          return isFullyDeveloped_;
        }

        //- Return boundary layer thickness
        scalar& delta()
        {
            return delta_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

template <class T>
void linearyInterpolate(std::vector<T>& x,std::vector<T>& y, T& x_value, T& y_value);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


#endif
