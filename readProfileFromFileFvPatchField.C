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
\*---------------------------------------------------------------------------*/
#include "readProfileFromFileFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include <iostream>
#include <fstream>
#include <vector>

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

readProfileFromFileFvPatchField::
readProfileFromFileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowNormal_(pTraits<vector>::zero),
    wallNormal_(pTraits<vector>::zero),
    readFromFile_("foo.dat"),
    isFullyDeveloped_(true),
    delta_(0)
{}

readProfileFromFileFvPatchField::
readProfileFromFileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowNormal_(dict.lookup("flowNormal")),
    wallNormal_(dict.lookup("wallNormal")),
    readFromFile_(dict.lookup("readFromFile")),
    isFullyDeveloped_(dict.readIfPresent("delta",delta_)) // read boundary layer thickness only if "delta" keyword is found in dictionary
{
   evaluate();
}

readProfileFromFileFvPatchField::
readProfileFromFileFvPatchField
(
    const readProfileFromFileFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flowNormal_(ptf.flowNormal_),
    wallNormal_(ptf.wallNormal_),
    readFromFile_(ptf.readFromFile_),
    isFullyDeveloped_(ptf.isFullyDeveloped_),
    delta_(ptf.delta_)
{}

readProfileFromFileFvPatchField::
readProfileFromFileFvPatchField
(
    const readProfileFromFileFvPatchField& blpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(blpvf, iF),
    flowNormal_(blpvf.flowNormal_),
    wallNormal_(blpvf.wallNormal_),
    readFromFile_(blpvf.readFromFile_),
    isFullyDeveloped_(blpvf.isFullyDeveloped_),
    delta_(blpvf.delta_)
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void readProfileFromFileFvPatchField::updateCoeffs()
{
    const fvMesh& patchMesh = patch().boundaryMesh().mesh(); // fvMesh needed to retreive nodal values

    const vectorField& faceCenters = patch().Cf();                      // Retreive coordinates of face centers of the patch
    const scalarField& wallNormalCoordinate(faceCenters & wallNormal_); // Obtain normal distance to the wall by the dot profuct with wallNormal vector

    scalar Cmin=min(patchMesh.points() & wallNormal_);  // Minimal nodal value
    scalar Cmax=max(patchMesh.points() & wallNormal_);  // Maximal nodal value

    scalar channelH = abs(Cmax-Cmin);              // Calculate Heigth of the channel
    scalarField Unormal(wallNormalCoordinate.size());          // Scalar field of wall normal welocity

    std::fstream inputFile;                        //

    inputFile.open(readFromFile_.typeName);       // use .typeName member of readFromFile_ class to open file

    inputFile <<"Uspesno upisao!" <<endl;
    Info << "Uspeno upisao!" << endl;

    inputFile.close();


    forAll(wallNormalCoordinate, i) // Loop around all the face-center coordinates and interpolate velocity
    {
      if (wallNormalCoordinate[i]<delta_)
      {
        scalar yCenter=wallNormalCoordinate[i];
        Unormal[i]=0.5;

      }
      else if ((channelH-wallNormalCoordinate[i])<delta_)
      {
        scalar yCenter=channelH-wallNormalCoordinate[i];
        Unormal[i]=0.5;

      }
      else
      {
        Unormal[i]=0.5;
      }
    }

    vectorField::operator=(flowNormal_ *Unormal);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void readProfileFromFileFvPatchField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("flowNormal")       << flowNormal_       << token::END_STATEMENT << nl;
    os.writeKeyword("wallNormal")       << wallNormal_       << token::END_STATEMENT << nl;
    os.writeKeyword("readFromFile")     << readFromFile_     << token::END_STATEMENT << nl;
    os.writeKeyword("isFullyDeveloped") << isFullyDeveloped_ << token::END_STATEMENT << nl;
    os.writeKeyword("delta")            << delta_            << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    readProfileFromFileFvPatchField
);

} // End namespace Foam*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
