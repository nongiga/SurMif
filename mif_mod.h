// IsoMIF is a program to identify molecular interaction field similarities between proteins
// Copyright (C) 2015 - Matthieu Chartier (under the supervision or Rafael Najmanovich)

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <math.h> 
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iterator>

using namespace std;

#define NEW(p,type)     if ((p=(type *) malloc (sizeof(type))) == NULL) { \
        printf ("Could not allocate ptr!\n");                                    \
        exit(0);                                                        \
}

#define PI 3.14159265



//Atom struc
struct atom{
  float x,y,z,xr,yr,zr;
  string resn;
  string atomn;
  string chain;
  string alt;
  string pseudo;
  string at;
  int resnb; //residue ID in PDB
  int atomnb; //atomnb in PDB
  int dir;
  int rDir;
  int h; //if its a hydrogen 0 or 1
  int bs; //if its a binding site atom
  int mif; //if it should be used to calculated MIFS

  atom(float x, float y, float z, string resn, string atomn, string chain, string alt, int resnb, int atomnb, int h, int bs, int mif):  x(x), y(y), z(z), resn(resn), atomn(atomn), chain(chain), alt(alt), resnb(resnb), atomnb(atomnb), h(x), bs(bs), mif(mif)
  {
    if (alt.empty()) alt="-";
    if (chain.empty()) chain=" ";
  }
  void setRadial(float x, float y, float z, float d) { xr=x; yr=y;zr=z;dir=d;}
};

//Grid point struc
struct vertex{
  float x,y,z;
  int p;
  int bu;
  vector<int> ints;
  vector<float> nrgs;
  vector<float> angles;
  int grid[4];
  map<string,float> env;
  int id;
  int modulo;
  vertex(float x, float y, float z): x(x), y(y), z(z) {}
  vertex(float x, float y, float z, float p, float b, float m): x(x), y(y), z(z), p(p), bu(b), modulo(m) {}
  vertex() {}
};

struct angRef{
    string r1;
    string r2;
    int rDir;
    int ring;
};

struct coord{
  float x;
  float y;
  float z;
};

struct pseudovrtx{
    float x,y,z;
    float dist;
    string type;
    pseudovrtx(float x, float y, float z, float d, string t): x(x), y(y), z(z), dist(d), type(t) {}
};

struct pbVrtx{
  float dist;
  float nrg;
  pbVrtx(float d, float n): dist(d), nrg(n) {}
  pbVrtx(){} 
};

struct pwRun{
  string pdbF;
  string cleftF;
  string rnc;
  string ligF;
  pwRun(string p, string c, string r, string l): pdbF(p), cleftF(c), rnc(r), ligF(l) {}
  pwRun() {}
  pwRun(vector<string> vec) {
    for(int i=0; i<vec.size(); i+=2){
      if(vec[i].compare("-p")==0) pdbF=vec[i+1];
      if(vec[i].compare("-g")==0) cleftF=vec[i+1];
      if(vec[i].compare("-l")==0) rnc=vec[i+1];
      if(vec[i].compare("-lf")==0) ligF=vec[i+1];
    }
  }
};
vector<pwRun> pw;

string pairwiseF="";
string cmdLine="";
string cleftFile="";
string gridFile="";
string proteinFile="";
string outBase="";
string outGridBase="";
string tag="";
string chain="";
string basePath="";
string resnumc="";
string resnumcShort="";
string matrixF="";
string probesF="";
string ligFile="";
string vcFile="";
string ff="original";
string statsF="original";
string resnumcList[20];
string fields[13];
float distpbV=2.0;
float gridStep=0.5;
float stepsize=0.5;
float maxGridDist=4.0;
float minGridDist=2.5;
float atmPbMaxDist=8.0;
float gridLigDist=3.0;
float caT=5.0;
int uID=0;
int smoothDist=0;
int printDetails=0;
int printAtoms=1;
int* probesList;
float* pbDistTmax;
float* pbDistTmin;
int nbOfAts=0;
int nbOfAtoms=0;
int nbOfProbes=0;
int ss[4];
int ssm[4];
int zip=-1;
int bul=14;
int buD=40;
int nbResnumc =1;
int nbArgumentNecessaire =5;
float min_x, min_y, min_z, max_x, max_y, max_z;

float* epsilons;
float angThresh=60.0;
map<string,string> atomTypes;
map<string,angRef> atomRef;
map<int,string> idAt;
map<string,string> pseudoC;
map<string,string> ligAt;
map<string,int> eps;
map<string,int> hyd;
map<string,int> arm;
map<string,int> don;
map<string,int> acc;
map<string,int> chr;
map<string,float> minD;
map<string,float> maxD;
map<string,float> nrgT;
vector<string> probes;
vector<string> aa;
vector<pseudovrtx> pseudoList;
coord gmin, gmax;
int width, height, depth;

class  Protein{
  public:

    Protein(string);
    ~Protein(void);
    void readPDB(string);
    void getAtomDir();
    int getRefAtom(float&, float&, float&, string, int, string, string, int, float, float, float, string);
    void clearAll();
    vector<atom> PROTEIN;
    vector<float> LIGAND;
    vector<atom> LIGATOMS;

  private:

};

class  Grid{
  public:

    Grid(string, Protein&);
    ~Grid(void);
    int readGetCleft(string, vector<atom>&, vector<float>&);
    int generateID(int, int, int, int, int);
    int buildGrid(Protein&);
    void getBuriedness();
    int getDiag(int, int, int, int&);
    int inGridRes(vertex&, float);
    void smooth();
    void writeMif(vector<atom>&,vector<atom>&);
    void clearAll();
    map<int,vertex> GRID;
    vector<int> vrtxIdList;

  private:
    int vrtx050,vrtx100,vrtx150,vrtx200;
};

int ProcessCmdLine(int, char**);
float dist_3d(float, float, float, float, float, float);
void getMif(map<int,vertex>&, vector<atom>&,vector<int>&);
void getEnv(map<int,vertex>&, vector<atom>&,vector<int>&);
void getPseudo(map<int,vertex>&, vector<atom>&,vector<int>&);
int is_coord_in_cube(float,float,float,float,float,float,float);
void stripSpace(string &);
void getStats(map<int,vertex>&, vector<atom>&, vector<int>&);
float roundCoord(float, int);
float roundCoordUp(float);
float roundCoordDown(float);
double calcNrg(vertex&, atom&, int, int&, float&);
void getAtomRef();
void getPseudoC();
void getEpsilons();
void getAtomTypes();
void readLigFile();
void getProbes();
void getaa();
bool compByNrg(pbVrtx,pbVrtx);
void getPairwise();


//miscellaneous methods
//added by Noga
void readReferenceFiles();
void printMenu();
void printToCmdline();
int readCmdLine(int, char**);
void initializeInfo(int);
void initializeBounds();
void squareDistances();
void fetchFields(string line);