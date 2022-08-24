//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/GridPDF.h>

// Function that construct an LHAPDF::GridPDF object using an
// apfel::InitialiseEvolution object as an input in the place of a
// standard LHAPDF grid.
LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const& ev)
{
  // Get knot array from APFEL
  const std::map<double, std::map<int, apfel::LHKnotArray>> ka = ev.KnotArray();

  // Knot array to be fed to LHAPDF
  std::map<double, LHAPDF::KnotArrayNF> knotarrays;

  // Vector of flavour indices to be filled
  std::vector<int> flavors;

  // Vectors needed to construct the LHAPDF::AlphaS object
  std::vector<double> q2as;
  std::vector<double> als;

  // Loop over Q subgrids
  for (auto const& sg : ka)
    {
      // NF-dependent knot array for a single subgrid in Q
      LHAPDF::KnotArrayNF kaNF;

      // Loop over flavour IDs and fill in the know array (use 21 for
      // the gluon).
      for (auto const& id : sg.second)
	kaNF.set_pid((id.first == 0 ? 21 : id.first), LHAPDF::KnotArray1F{id.second.xs, id.second.q2s, id.second.xfs});

      // Fill in the Q2 grid for alphas using the first element of the
      // map. Remove the last element but do not do it for the last
      // subgrid (this is why pop_back appears before).
      if (!q2as.empty())
	q2as.pop_back();
      q2as.insert(q2as.end(), sg.second.begin()->second.q2s.begin(), sg.second.begin()->second.q2s.end());

      // Fill in the flavour ID vector
      if (flavors.empty())
	for (auto const& id : sg.second)
	  flavors.push_back((id.first == 0 ? 21 : id.first));

      // Fill in the LHAPDF:KnotArrayNF
      knotarrays.insert({sg.first, kaNF});
    }

  // Fill in alpha_s vector
  for (auto const& q2 : q2as)
    als.push_back(ev.Alphas(sqrt(q2)));

  // Define an LHAPDF GridPDF object
  LHAPDF::GridPDF* dist = new LHAPDF::GridPDF{};

  // Pass knot arrays to LHAPDF
  dist->knotarrays() = knotarrays;

  // Set interpolator
  dist->setInterpolator((std::string) "logcubic");

  // Set extrapolator
  dist->setExtrapolator((std::string) "continuation");

  // Set flavours
  dist->setFlavors(flavors);

  // Set alpha_s
  LHAPDF::AlphaS_Ipol* as = new LHAPDF::AlphaS_Ipol{};
  as->setQ2Values(q2as);
  as->setAlphaSValues(als);
  dist->setAlphaS(as);
  
  return dist;
}

int main()
{
  // Open LHAPDF set
  LHAPDF::PDF* distLH = LHAPDF::mkPDF("CT18NNLO");

  // APFEL++ default EvolutionSetup object
  apfel::EvolutionSetup es{};

  // Construct vectors of thresholds and masses
  std::vector<double> ths, mss;
  for (auto const& f : distLH->flavors())
    if (f > 0 && f < 7)
      {
	ths.push_back(distLH->quarkThreshold(f));
	mss.push_back(distLH->quarkMass(f));
      }

  // Adjust evolution parameters to match those of the input set
  es.Q0                = distLH->qMin();
  es.Qmin              = distLH->qMin();
  es.Qmax              = distLH->qMax();
  es.PerturbativeOrder = distLH->orderQCD();
  es.QQCDRef           = 91.1876;
  es.AlphaQCDRef       = distLH->alphasQ(es.QQCDRef);
  es.Thresholds        = ths;
  es.Masses            = mss;
  es.InSet             = {[=] (double const& x, double const& Q) -> std::map<int,double> { return apfel::PhysToQCDEv(distLH->xfxQ(x, Q)); }};

  // Feed it to the initialisation class of APFEL++ and construct
  // pointer to LHAPDF::PDF object
  LHAPDF::PDF* distAP = mkPDF(apfel::InitialiseEvolution{es});

    // Initialize APFEL++
  const apfel::Grid g{{{100, 1e-5, 3}, {100, 1e-1, 3}, {100, 6e-1, 3}, {80, 8.5e-1, 5}}};
  apfel::AlphaQCD a{es.AlphaQCDRef, es.QQCDRef, es.Thresholds, es.PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 1, 1000, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
  const auto EvolvedPDFs = apfel::BuildDglap(apfel::InitializeDglapObjectsQCD(g, es.Thresholds),
					     [&] (double const& x, double const& Q) -> std::map<int, double>{ return apfel::PhysToQCDEv(distLH->xfxQ(x, Q)); },
					     es.Q0, es.PerturbativeOrder, as);
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Print test results
  std::cout << std::scientific;
  std::cout.precision(4);

  const double mu = 100;
  const std::vector<double> xlha{1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << "\nmu = " << mu << " GeV" << std::endl;

  std::cout << "\nLHAPDF (tabulated) evolution:" << std::endl;
  std::cout << "AlphaQCD(Q) = " << distLH->alphasQ(mu) << std::endl;
  std::cout << "   x    "
	    << "   u-ubar   "
	    << "   d-dbar   "
	    << " 2(ubr+dbr) "
	    << "   c+cbar   "
	    << "    gluon   "
	    << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> PDFmap = distLH->xfxQ(x, mu);
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << PDFmap.at(2) - PDFmap.at(-2)
		<< "  " << PDFmap.at(1) - PDFmap.at(-1)
		<< "  " << 2 * ( PDFmap.at(-2) + PDFmap.at(-1) )
		<< "  " << PDFmap.at(4) + PDFmap.at(-4)
		<< "  " << PDFmap.at(21)
		<< std::endl;
    }
  std::cout << "\nAPFEL++ evolution through LHAPDF:" << std::endl;
  std::cout << "AlphaQCD(Q) = " << distAP->alphasQ(mu) << std::endl;
  std::cout << "   x    "
	    << "   u-ubar   "
	    << "   d-dbar   "
	    << " 2(ubr+dbr) "
	    << "   c+cbar   "
	    << "    gluon   "
	    << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> PDFmap = distAP->xfxQ(x, mu);
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << PDFmap.at(2) - PDFmap.at(-2)
		<< "  " << PDFmap.at(1) - PDFmap.at(-1)
		<< "  " << 2 * ( PDFmap.at(-2) + PDFmap.at(-1) )
		<< "  " << PDFmap.at(4) + PDFmap.at(-4)
		<< "  " << PDFmap.at(21)
		<< std::endl;
    }

  std::cout << "\nAPFEL++ evolution stand alone:" << std::endl;
  std::cout << "AlphaQCD(Q) = " << as(mu) << std::endl;
  std::cout << "   x    "
	    << "   u-ubar   "
	    << "   d-dbar   "
	    << " 2(ubr+dbr) "
	    << "   c+cbar   "
	    << "    gluon   "
	    << std::endl;
  const std::map<int, apfel::Distribution> tpdfs = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(mu).GetObjects());
  for (auto const& x : xlha)
    {
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << tpdfs.at(2).Evaluate(x) - tpdfs.at(-2).Evaluate(x)
		<< "  " << tpdfs.at(1).Evaluate(x) - tpdfs.at(-1).Evaluate(x)
		<< "  " << 2 * ( tpdfs.at(-2).Evaluate(x) + tpdfs.at(-1).Evaluate(x) )
		<< "  " << tpdfs.at(4).Evaluate(x) + tpdfs.at(-4).Evaluate(x)
		<< "  " << tpdfs.at(0).Evaluate(x)
		<< std::endl;
    }
  std::cout << "\n";
}
