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
    // Input data from file ::
    std::vector<scalar> psiFromFile; //std::vector to be read from input file
    std::vector<scalar> yFromFile;  // std::vector to be read from input file
    scalar SIZE;

    std::ifstream inputFile(readFromFile_.typeName); // use .typeName member of readFromFile_ class to open file
    if(!inputFile.is_open())
    {
      //code error
    }
    inputFile >> SIZE;
    yFromFile.assign(SIZE,0);
    psiFromFile.assign(SIZE,0);

    for(int i=0;i<SIZE;i++)
    {
        inputFile >> yFromFile.at(i) >> psiFromFile.at(i);
    }
    Info << "Data input successfully finished!" << endl;
    inputFile.close();

    // Declare additional fields ::
    const fvMesh& patchMesh = patch().boundaryMesh().mesh();   // patchMesh needed to retreive nodal values
    const scalar Cmin=min(patchMesh.points() & wallNormal_);   // Find wall-normal coordinates of all nodes, then find the lovest one
    const scalar Cmax=max(patchMesh.points() & wallNormal_);   // Find wall-normal coordinates of all nodes, then find the biggest one
    const scalar channelH = abs(Cmax-Cmin);                    // Calculate Heigth of the channel

    const vectorField& faceCenters = patch().Cf();             // Retreive coordinates of face centers of the patch
    scalarField wallNormalDistance(faceCenters & wallNormal_);  // Obtain normal coordinate to the wall by the dot profuct with wallNormal vector
    forAll(wallNormalDistance, i)
    {
      wallNormalDistance[i]=wallNormalDistance[i]-Cmin;        // Now wall normal coordinate becomes distance, perpendicular to the wall
    }
    scalarField psiNormal(wallNormalDistance.size());        // Initialize a scalar field to be written to the field
    scalar yCenter;
    if(isFullyDeveloped_)
    {
      forAll(wallNormalDistance, i)
      {
        if(wallNormalDistance[i]<0.5*channelH)
        {
          yCenter=wallNormalDistance[i];
        }
        else
        {
          yCenter=channelH-wallNormalDistance[i];
        }
        linearyInterpolate<scalar>(yFromFile,psiFromFile,yCenter,psiNormal[i]);
      }

    }
    vectorField::operator=(flowNormal_ *psiNormal);

    fixedValueFvPatchVectorField::updateCoeffs();
}

template <class T>
void linearyInterpolate(std::vector<T>& x,std::vector<T>& y, T& x_value,T& y_value)
{
  for(unsigned int i=1;i<x.size();i++)
  {
    if(x_value>=x.at(i-1) && x_value <=x.at(i))
    {
      y_value= y.at(i-1)+(y.at(i)-y.at(i-1)*(x_value-x.at(i-1))/(x.at(i)-x.at(i-1)));
      return;
    }

  }
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
