//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <cmath>
#include <iostream>

#include <apfel/timer.h>
#include <apfel/evolutionsetup.h>
#include <apfel/initialiseevolution.h>

#include <LHAPDF/LHAPDF.h>
#include "LHAPDF/GridPDF.h"

using namespace apfel;
using namespace LHAPDF;
using namespace std;

int main()
{
  // Default EvolutionSetup object
  EvolutionSetup es{};

  // Feed it to the initialisation class
  InitialiseEvolution ev{es};

  // Get knot array from APFEL
  const map<double, map<int, LHKnotArray>> ka = ev.KnotArray();

  // Knot array to be fed to LHAPDF
  map<double, KnotArrayNF> knotarrays;

  // Vector of flavour indices to be filled
  vector<int> flavors;

  // Construct AlphaS object.
  vector<double> q2as;
  vector<double> als;

  // Loop over Q subgrids
  for (auto const& sg : ka)
    {
      // NF-dependent knot array for a single subgrid in Q
      KnotArrayNF kaNF;

      // Loop over flavour IDs and fill in the know array
      for (auto const& id : sg.second)
	kaNF.set_pid((id.first == 0 ? 21 : id.first), KnotArray1F{id.second.xs, id.second.q2s, id.second.xfs});

      // Fill in Q2 grid for alphas using the first element of the
      // map. Remove last element but do not do it fot the last
      // subgrid (this is why pop_back appears before).
      q2as.pop_back();
      q2as.insert(q2as.end(), sg.second.begin()->second.q2s.begin(), sg.second.begin()->second.q2s.end());

      // Fill in the flavour ID vector
      if (flavors.empty())
	for (auto const& id : sg.second)
	  flavors.push_back((id.first == 0 ? 21 : id.first));

      knotarrays.insert({sg.first, kaNF});
    }

  // Fill in alpha_s vector
  for (auto const& q2 : q2as)
    als.push_back(ev.Alphas(sqrt(q2)));

  // Define an LHAPDF GridPDF object
  GridPDF dist;

  // Pass knot arrays to LHAPDF
  dist.knotarrays() = knotarrays;

  // Set interpolator
  dist.setInterpolator((string) "logcubic");

  // Set extrapolator
  dist.setExtrapolator((string) "continuation");

  // Set flavours
  dist.setFlavors(flavors);

  // Set alpha_s
  AlphaS_Ipol* as = new AlphaS_Ipol{};
  as->setQ2Values(q2as);
  as->setAlphaSValues(als);
  dist.setAlphaS(as);

  cout << dist.xfxQ(0, 0.1, 100) << endl;
  cout << dist.alphasQ(100) << endl;
}
