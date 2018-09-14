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

#include <stdio.h>  // Required to acces current working directory
#include <unistd.h> // Required to acces current working directory

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
  /* ----------------------------------------------------------------------- *\
                         Input data from file ::
  \* ----------------------------------------------------------------------- */
	  Info << "Data will be read from the file ::" << readFromFile_<< endl;

    std::vector<scalar> psiFromFile; //std::vector, to be read from input file
    std::vector<scalar> yFromFile;  // std::vector, to be read from input file
    scalar SIZE;                    // Foam::scalar, size to be read from input file

    std::ifstream inputFile(readFromFile_);

    if(!inputFile.is_open())    // If file not opened, print error message and exit
    {
      char buff[FILENAME_MAX];
      getcwd( buff, FILENAME_MAX );  // Ignore compiler warning if any
      std::string current_working_dir(buff);

      FatalErrorInFunction
        << "File " << readFromFile_ << " not found!" << "\n"
        << "Check if : " << readFromFile_ << " is present in the current working directory ("
        << current_working_dir<<")"
        <<exit(FatalError);
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
  /* ----------------------------------------------------------------------- *\
            Access geometrical information of the current patch ::
  \* ----------------------------------------------------------------------- */

    label patchID=patch().boundaryMesh().findPatchID(patch().name()); // Retreive the current patch ID ;

    const fvMesh& mesh = patch().boundaryMesh().mesh();                // Access mesh

    const List<vector> & currentPatchPoints = mesh.boundaryMesh()[patchID].localPoints(); // Retreive list of all points that belong to the boundary patch

    scalar minPoint=min(currentPatchPoints & wallNormal_);  // Calculate minimal nodal value with scalar prod. with wall normal
    scalar maxPoint=max(currentPatchPoints & wallNormal_);  // Calculate maximal nodal value with scalar prod. with wall normal
    scalar channelH = maxPoint-minPoint;                            // Calculate heigth of the channel

    Info << "minPoin= " << minPoint << endl;
    Info << "maxPoin= " << maxPoint << endl;
    Info << "channelH= " << channelH << endl;

    const vectorField& patchFaceCentres = mesh.Cf().boundaryField()[patchID]; // Retreive centres of the faces belongigng to the patch;

    scalarField wallNormalCoordinate(patchFaceCentres & wallNormal_); // Coordinate normal to the wall

    scalar yCenter;

  /* ----------------------------------------------------------------------- *\
            Interpolate psi values  ::
  \* ----------------------------------------------------------------------- */
  Info <<"Number of faces :: "<< wallNormalCoordinate.size() << endl; // Ovde nastaje problem, nece da stampa pravilo ako je const scalarField&
  if(isFullyDeveloped_) // Ima problem da prihvati da je fully developed
  {
    forAll(wallNormalCoordinate, i) // Loop around all the face-center coordinates and interpolate velocity
    {
      yCenter = wallNormalCoordinate[i]-minPoint;
      Info << i <<" | " << yCenter << endl;

    }
  }
  else
  {

  }






    FatalError << "Namerno sam prekinuo ovde! "<< exit(FatalError);

   // vectorField::operator=(patch.nf & Unormal);

   // fixedValueFvPatchVectorField::updateCoeffs();
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

/*template <class T>
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
}*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    readProfileFromFileFvPatchField
);

} // End namespace Foam*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
