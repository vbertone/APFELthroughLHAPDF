/*
  Author: Valerio Bertone
 */

// Standard libs
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <glob.h>
#include <iomanip>
#include <locale>
#include <stdio.h>
#include <stdlib.h>
#include <functional>

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"

// APFEL++ libs
#include "apfel/dglapbuilder.h"
#include "apfel/grid.h"
#include "apfel/timer.h"
#include "apfel/alphaqcd.h"
#include "apfel/tabulateobject.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"
#include "apfel/dglap.h"
#include "apfel/rotations.h"

using namespace std;
using namespace apfel;

int main() {
  //const string set = "CT14nnlo";
  const string set = "MSTW2008nnlo68cl";

  // Open LHAPDF set
  //LHAPDF::PDF* dist = LHAPDF::mkPDF(set);
  LHAPDF::GridPDF* dist = new LHAPDF::GridPDF{set, 0};

  // Final scale
  double mu = 100;

  // Retrieve evolution parameters
  const int    pto   = dist->orderQCD();
  const double Qref  = 91.1876;
  const double asref = dist->alphasQ(Qref);
  const double mc    = dist->quarkThreshold(4);
  const double mb    = dist->quarkThreshold(5);
  const double mt    = dist->quarkThreshold(6);
  const double Qin   = dist->qMin();

  // Initialize APFEL++
  const Grid g{{SubGrid{100,1e-5,3}, SubGrid{60,1e-1,3}, SubGrid{50,6e-1,3}, SubGrid{50,8e-1,5}}};
  const vector<double> Masses = {0, 0, 0, mc, mb, mt};
  AlphaQCD a{asref, Qref, Masses, pto};
  const TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
  const auto InPDFs   = [&] (double const& x, double const& Q) -> map<int,double>{ return PhysToQCDEv(dist->xfxQ(x,Q)); };
  const auto DglapObj = InitializeDglapObjectsQCD(g, Masses);
  auto EvolvedPDFs    = BuildDglap(DglapObj, InPDFs, Qin, pto, as);
  const TabulateObject<Set<Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 0.95, 1000, 3};

  // Print results
  cout << scientific;

  cout << "\nmu = " << mu << " GeV\n" << endl;
  cout << "APFEL++: AlphaQCD(Q) = " << Alphas.Evaluate(mu) << endl;
  cout << "LHAPDF:  AlphaQCD(Q) = " << dist->alphasQ(mu) << endl;

  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  cout << endl;
  for (auto i = 2; i < 11; i++)
    {
      cout.precision(1);
      cout << xlha[i];
      cout.precision(4);
      cout << " APFEL++: " << TabulatedPDFs.EvaluatexQ(0,xlha[i],mu)
	   << " LHAPDF:  " << dist->xfxQ(0,xlha[i],mu)
	   << endl;
    }
  cout << "      " << endl;

  // Tests
/*
  const vector<double> xg = dist->xKnots();
  for (auto const& x : xg)
    cout << x << endl;

  cout << "\n";

  const vector<double> q2g = dist->q2Knots();
  for (auto const& q2 : q2g)
    cout << q2 << endl;
*/

  std::map<double, LHAPDF::KnotArrayNF> ka1 = dist->knotarrays();

  std::map<double, LHAPDF::KnotArrayNF>& ka = dist->knotarrays();
  for (auto const& m : ka)
    {
      cout << m.second.size() << endl;
      cout << sqrt(m.first) << endl;
      const vector<double> q2sg = m.second.get_first().q2s();
      for (auto const& q2 : q2sg)
	      cout << "- " << sqrt(q2) << endl;
    }
  cout << "\n";

  for (auto const& m : ka)
    {
      const vector<double> xsg = m.second.get_first().xs();
      for (auto const& x : xsg)
	      cout << "- " << x << endl;
    }

  return 0;
}
