#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Functional/CompiledModel.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Statistics/MatSample.hpp>
#include <LatAnalyze/Core/Math.hpp>
#include <LatAnalyze/Numerical/MinuitMinimizer.hpp>
#include <LatAnalyze/Numerical/NloptMinimizer.hpp>
#include <LatAnalyze/Core/Plot.hpp>
#include <LatAnalyze/Statistics/XYSampleData.hpp>

using namespace std;
using namespace Latan;

int main(int argc, char *argv[])
{
    // parse arguments /////////////////////////////////////////////////////////
    OptParser            opt;
    bool                 parsed, imag, negate;
    string               plotFileName, outFileName, xName, yName, title, save;
    vector<string>       inFileName;
    double               xLow, xHigh, spacing;

    opt.addOption("i" , "imag"   , OptParser::OptType::trigger, true,
                "plot imaginary");
    opt.addOption("n" , "neg"   , OptParser::OptType::trigger, true,
                "negate input");
    opt.addOption("o", "output", OptParser::OptType::value  , true,
                  "output file", "");
    opt.addOption("x", "xAxis", OptParser::OptType::value  , true,
                  "x-axis name", "");
    opt.addOption("l", "xLow", OptParser::OptType::value  , true,
                  "x-axis lower bound", "0");
    opt.addOption("y", "yAxis", OptParser::OptType::value  , true,
                  "y-axis name", "");
    opt.addOption("", "spacing", OptParser::OptType::value  , true,
                  "spacing between points", "1");
    opt.addOption("s", "save", OptParser::OptType::value, true,
                    "saves the source and .pdf", "");
    opt.addOption("t", "title", OptParser::OptType::value  , true,
                  "plot title", "");
    opt.addOption("", "help"      , OptParser::OptType::trigger, true,
                   "show this help message and exit");
    parsed = opt.parse(argc, argv);
    if (!parsed or opt.gotOption("help") or opt.getArgs().size() != 1)
    {
        cerr << "usage: " << argv[0] << " <.h5/manifest file> <options> " << endl;
        cerr << endl << "Possible options:" << endl << opt << endl;
        
        return EXIT_FAILURE;
    }
    
    plotFileName = opt.getArgs().front();
    imag         = opt.gotOption("i");
    negate       = opt.gotOption("n");
    xName        = opt.optionValue("x");
    xLow         = opt.optionValue<double>("l");
    yName        = opt.optionValue("y");
    spacing      = opt.optionValue<double>("spacing");
    save         = opt.optionValue("s");
    title        = opt.optionValue("t");
    outFileName  = opt.optionValue<string>("o");


    if(plotFileName.find(".h5") == string::npos)
    {
        inFileName = readManifest(plotFileName);
    }
    
    // load and plot file(s) /////////////////////////////////////////////////////////
    DMatSample tmp;
    Index      nt,nSample;
    Plot       p;
    DVec       tAxis;

    bool coshModel=0;
    bool sinhModel=0;
    bool linearModel=0;
    bool constModel=0;



    if(inFileName.size() == 0)
    {
        tmp     = Io::load<DMatSample>(plotFileName);
        nt      = tmp[central].rows();
        nSample = tmp.size();
        if (negate)
            tmp *= -1;
        if(imag)
        {
            tmp     = tmp.block(0, 1, nt, 1);
        }
        else
        {
            tmp     = tmp.block(0, 0, nt, 1);
        }
        xHigh= xLow+spacing*(nt-1);
        tAxis.setLinSpaced(nt, xLow, xHigh);
        p  << PlotData(tAxis, tmp); 


        {
                Plot       p;
                DMatSample effMass(nSample);
                DVec       effMassT;
                Index      maxT = nt - 1;
                
                // // Log meff
                effMass.resizeMat(maxT, 1);
                effMassT.setLinSpaced(maxT, 0, maxT-1);
                FOR_STAT_ARRAY(effMass, s)
                {
                    for (Index t = 1; t < nt; ++t)
                    {
                        effMass[s](t - 1) = log(tmp[s](t-1)/tmp[s](t));
                    }
                }

                // SemiLept meff
                // effMass.resizeMat(maxT-1, 1);
                // effMassT.setLinSpaced(maxT-1, 1, maxT-1);
                // FOR_STAT_ARRAY(effMass, s)
                // {
                //     effMass[s](0 ) = 0;
                //     for (Index t = 0; t < nt-2; ++t)
                //     {
                //         effMass[s](t ) = 0.5*(tmp[s](t+2)-tmp[s](t));
                //     }
                // }

                p.reset();
                p << PlotRange(Axis::x, 1, maxT);
                p << Color("rgb 'red'") << PlotData(effMassT, effMass);
                p << Caption("Effective Mass");
                p.display();
                // if(savePlot != "")
                // {
                //     p.save(savePlot + "_effMass");
                // }
            }
    
    }
    else
    {
        tmp = Io::load<DMatSample>(inFileName[0]);
        nt      = tmp[central].rows();
        xHigh= xLow+spacing*(nt-1);
        tAxis.setLinSpaced(nt, xLow, xHigh);

        for(unsigned long i = 0; i < inFileName.size(); i++)
        {
            plotFileName = inFileName[i];
            tmp = Io::load<DMatSample>(plotFileName);
        if (negate)
            tmp *= -1;
        if(imag)
        {
            tmp     = tmp.block(0, 1, nt, 1);
        }
        else
        {
            tmp     = tmp.block(0, 0, nt, 1);
        }   
            p  << PlotData(tAxis, tmp);
            
        }
        
    }

    p << Label(xName, Axis::x);
    p << Label(yName, Axis::y);
    p << PlotRange(Axis::x, xLow, xHigh);
    p << Caption(title);

    if(save != "")
    {
        cout << "Saving plot and source code to " << save << endl;
        p.save(save + "/" + title);
    }
    
    else
    {
        cout << "Displaying plot..." << endl;
        p.display();    
    }
    
    

    // output //////////////////////////////////////////////////////////////////
    if (!outFileName.empty())
    {
        Io::save(tmp, outFileName);
        cout << "File saved as: " << outFileName << endl;
    }

    return EXIT_SUCCESS;
}
