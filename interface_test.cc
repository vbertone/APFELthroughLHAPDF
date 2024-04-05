//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/mkpdf.h>
#include <LHAPDF/LHAPDF.h>
/*
#include <LHAPDF/GridPDF.h>

// Function that constructs a pointer to an LHAPDF::PDF object using
// an apfel::InitialiseEvolution object as an input in the place of a
// standard LHAPDF grid.
LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const& ev)
{
  // Get knot array from APFEL
  const std::map<double, std::map<int, apfel::LHKnotArray>> ka = ev.KnotArray();

  // Knot array to be fed to LHAPDF
  LHAPDF::KnotArray data;

  // Loop over Q subgrids to accumulate x-grid, q2-grid, and flavour
  // IDs.
  for (auto const& sg : ka)
    {
      // Get x-grid from first flavour of the first subgrid
      if (data.xs().empty())
	data.setxknots() = sg.second.begin()->second.xs;

      // Accumulate q2 grid using the first flavour
      data.setq2knots().insert(data.setq2knots().end(), sg.second.begin()->second.q2s.begin(), sg.second.begin()->second.q2s.end());

      // Fill in the flavour ID vector (use 21 for the gluon)
      if (data.setPids().empty())
	for (auto const& id : sg.second)
	  data.setPids().push_back((id.first == 0 ? 21 : id.first));
    }

  // Set up the knots of the Knotarray
  data.setShape() = std::vector<size_t>{data.xs().size(), data.q2s().size(), data.setPids().size()};
  data.fillLogKnots();
  data.initPidLookup();

  // Set size of data vector
  data.setGrid().resize(data.shape(0) * data.shape(1) * data.shape(2), 0.);

  // Loop over Q subgrids to accumulate data
  int nqprev = 0;
  for (auto const& sg : ka)
    {
      // Initialise flavour index
      int ifl = 0;

      // Get q2-subgrid size
      const int q2size = (int) sg.second.begin()->second.xfs.size() / data.shape(0);

      // Loop over flavours
      for (auto const& id : sg.second)
	{
	  // Get vector of distributions
	  const std::vector<double> f = id.second.xfs;

	  // Loop over x-grid
	  for (int ix = 0; ix < (int) data.shape(0); ix++)
	    {
	      int iq = nqprev;
	      // Loop over q2-sugrid
	      for (int iqs = 0; iqs < q2size; iqs++)
		data.setGrid()[ix * data.shape(2) * data.shape(1) + iq++ * data.shape(2) + ifl] = f[ix * q2size + iqs];
	    }
	  ifl++;
	}
      nqprev += q2size;
    }

  // Fill in alpha_s vector
  std::vector<double> als;
  for (auto const& q2 : data.q2s())
    als.push_back(ev.Alphas(sqrt(q2)));

  // Construct alpha_s object
  LHAPDF::AlphaS_Ipol* as = new LHAPDF::AlphaS_Ipol{};
  as->setQ2Values(data.q2s());
  as->setAlphaSValues(als);

  // Define an LHAPDF GridPDF object
  LHAPDF::GridPDF* dist = new LHAPDF::GridPDF{};

  // Pass knot arrays to LHAPDF
  dist->Data() = data;

  // Set interpolator
  dist->setInterpolator((std::string) "logcubic");

  // Set extrapolator
  dist->setExtrapolator((std::string) "continuation");

  // Set flavours
  dist->setFlavors(data.setPids());

  // Set alpha_s
  dist->setAlphaS(as);

  // Return object
  return dist;
}
*/
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
  es.QQCDRef           = apfel::ZMass;
  es.AlphaQCDRef       = distLH->alphasQ(es.QQCDRef);
  es.Thresholds        = ths;
  es.Masses            = mss;
  es.InSet             = {[=] (double const& x, double const& Q) -> std::map<int,double> { return apfel::PhysToQCDEv(distLH->xfxQ(x, Q)); }};

  // Feed it to the initialisation class of APFEL++ and construct
  // pointer to LHAPDF::PDF object
  apfel::InitialiseEvolution ev{es};
  LHAPDF::PDF* distAP = mkPDF(ev);

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

  std::cout << "\nAPFEL++ evolution through InitializeEvolution:" << std::endl;
  std::cout << "AlphaQCD(Q) = " << distAP->alphasQ(mu) << std::endl;
  std::cout << "   x    "
	    << "   u-ubar   "
	    << "   d-dbar   "
	    << " 2(ubr+dbr) "
	    << "   c+cbar   "
	    << "    gluon   "
	    << std::endl;
  const std::map<int, apfel::Distribution> td = apfel::QCDEvToPhys(ev.TabulatedDistributions().Evaluate(mu).GetObjects());
  for (auto const& x : xlha)
    {

      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << td.at(2).Evaluate(x) - td.at(-2).Evaluate(x)
		<< "  " << td.at(1).Evaluate(x) - td.at(-1).Evaluate(x)
		<< "  " << 2 * ( td.at(-2).Evaluate(x) + td.at(-1).Evaluate(x) )
		<< "  " << td.at(4).Evaluate(x) + td.at(-4).Evaluate(x)
		<< "  " << td.at(0).Evaluate(x)
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
