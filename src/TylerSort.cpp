/* #region Includes */

// C++ Includes
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// ROOT Includes
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

// CASort Includes
#include <CASort/CAAddBack.hpp>
#include <CASort/CACalibration.hpp>
#include <CASort/CAConfiguration.hpp>
#include <CASort/CACrosstalkCorrection.hpp>
#include <CASort/CAGainCorrection.hpp>
#include <CASort/CAUtilities.hpp>

// Project Includes
#include "Histograms.hpp"

/* #endregion Includes */

// Main function
int main(int argc, char* argv[])
{

    /* #region Configuration Setup */

    auto args = CAUtilities::ParseArguments(argc, argv);

    std::cout << "=============== Welcome to TylerSort! ==================" << std::endl;
    CAUtilities::PrintConfiguration(args);

    printf("[INFO] ROOT multithreading enabled with %d threads\n", ROOT::GetThreadPoolSize());

    /* #endregion Configuration Setup */

    /* #region Calibration Setup */

    // Gain match functions
    std::vector<std::function<double(double)>> ccGainMatch, cbGainMatch, psGainMatch, ceGainMatch;
    auto funcsGainMatch = CAGainCorrection::MakeCorrections(args.gainShiftFile);

    // Cross-talk Correction Functions
    std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>> xTalkCorrection;
    try
    {
        xTalkCorrection = CACrosstalkCorrection::MakeCorrections("/mnt/Data_0/70Ge_Kowalewski/CAXTalkCorr/70Ge_xtk.caxt");
        printf("[INFO] Crosstalk correction functions retrieved\n");
    }
    catch (const std::exception& e)
    {
        printf("[WARN] Crosstalk correction functions not found, proceeding without crosstalk correction\n");
        xTalkCorrection = std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>>(4, [](std::array<double, 4> x)
                                                                                                   { return x; });
    }

    // Calibraration functions
    std::vector<std::function<double(double)>> ccECalibrate, cbECalibrate, psECalibrate, ceECalibrate;

    for (int det : {1, 3, 5, 7})
    {
        for (int xtal = 1; xtal <= 4; xtal++)
        {
            std::string calFileName = Form("%s/C%iE%i.cal_params.txt",
                                           args.calibrationDir.c_str(), det, xtal);
            ccECalibrate.push_back(CACalibration::MakeCalibration(calFileName));
        }
    }
    try
    {
        ccGainMatch = funcsGainMatch.at(0);
        printf("[INFO] Clover Cross Gain Match and Calibration functions retrieved\n");
    }
    catch (const std::exception& e)
    {
        printf("[WARN] Clover Cross Gain Match functions not found, proceeding without gain matching\n");
        ccGainMatch = std::vector<std::function<double(double)>>(16, [](double x)
                                                                 { return x; });
    }

    for (int det : {1, 2, 3, 5})
    {
        for (int xtal = 1; xtal <= 4; xtal++)
        {
            std::string calFileName = Form("%s/B%iE%i.cal_params.txt",
                                           args.calibrationDir.c_str(), det, xtal);
            cbECalibrate.push_back(CACalibration::MakeCalibration(calFileName));
        }
    }
    try
    {
        cbGainMatch = funcsGainMatch.at(1);
        printf("[INFO] Clover Back Gain Match and Calibration functions retrieved\n");
    }
    catch (const std::exception& e)
    {
        printf("[WARN] Clover Back Gain Match functions not found, proceeding without gain matching\n");
        cbGainMatch = std::vector<std::function<double(double)>>(16, [](double x)
                                                                 { return x; });
    }

    /* #endregion Calibration Setup */

    /* #region Event Loop Setup*/

    auto inputFile = TFile::Open(args.runFileName.c_str());
    if (!inputFile || inputFile->IsZombie())
    {
        throw std::runtime_error(Form("[ERROR] Error opening input file: %s", args.runFileName.c_str()));
    }
    printf("[INFO] Opened file %s\n", args.runFileName.c_str());

    // Peak at the TTree to get the number of events
    TTree* tree;
    inputFile->GetObject("clover", tree);
    if (!tree)
    {
        throw std::runtime_error("[ERROR] Error opening TTree");
    }

    ULong64_t nEntries = tree->GetEntries();

    printf("[INFO] Opened TTree \"clover\" and counted %llu events\n", nEntries);

    delete tree;
    inputFile->Close();

    // Atomic counter for processed entries
    std::atomic<uint64_t> processedEntries(0);

    // Start the progress bar in a separate thread
    std::thread progressBarThread(CAUtilities::DisplayProgressBar, std::ref(processedEntries), nEntries);

    printf("[INFO] Processing events...\n");

    // Limit max tasks to reduce TBB work-stealing related race conditions
    // This ensures better alignment between TBB tasks and ROOT's TThreadedObject slots
    ROOT::TTreeProcessorMT::SetTasksPerWorkerHint(1);

    // Create a TTreeReader to read the TTree
    ROOT::TTreeProcessorMT EventProcessor(args.runFileName.c_str(), "clover");

    /* #endregion Event Loop Setup */

    // Fill Function
    auto fillHistograms = [&](TTreeReader& eventReader)
    {
        /* #region Set the branch addresses for the TTree */

        // Clover Cross
        TTreeReaderArray<double> cc_amp_val(eventReader, "clover_cross.amplitude");
        TTreeReaderArray<double> cc_cht_val(eventReader, "clover_cross.channel_time");
        TTreeReaderArray<double> cc_mdt_val(eventReader, "clover_cross.module_timestamp");
        TTreeReaderArray<double> cc_plu_val(eventReader, "clover_cross.pileup");
        TTreeReaderArray<double> cc_trt_val(eventReader, "clover_cross.trigger_time");

        // Clover Back
        TTreeReaderArray<double> cb_amp_val(eventReader, "clover_back.amplitude");
        TTreeReaderArray<double> cb_cht_val(eventReader, "clover_back.channel_time");
        TTreeReaderArray<double> cb_mdt_val(eventReader, "clover_back.module_timestamp");
        TTreeReaderArray<double> cb_plu_val(eventReader, "clover_back.pileup");
        TTreeReaderArray<double> cb_trt_val(eventReader, "clover_back.trigger_time");

#if PROCESS_POS_SIG
        // Pos Sig
        TTreeReaderArray<double> ps_amp_val(eventReader, "pos_sig.amplitude");
        TTreeReaderArray<double> ps_cht_val(eventReader, "pos_sig.channel_time");
        TTreeReaderArray<double> ps_mdt_val(eventReader, "pos_sig.module_timestamp");
        TTreeReaderArray<double> ps_plu_val(eventReader, "pos_sig.pileup");
        TTreeReaderArray<double> ps_trt_val(eventReader, "pos_sig.trigger_time");
#endif // PROCESS_POS_SIG

#if PROCESS_CEBR_ALL
        // CEBR All
        TTreeReaderArray<double> ce_inL_val(eventReader, "cebr_all.integration_long");
        TTreeReaderArray<double> ce_cht_val(eventReader, "cebr_all.channel_time");
        TTreeReaderArray<double> ce_mdt_val(eventReader, "cebr_all.module_timestamp");
        TTreeReaderArray<double> ce_ins_val(eventReader, "cebr_all.integration_short");
        TTreeReaderArray<double> ce_trt_val(eventReader, "cebr_all.trigger_time");
#endif // PROCESS_CEBR_ALL

        /* #endregion */

        /* #region Get Histogram pointers*/

        // Use histograms defined in Histograms.hpp

        // Clover Cross Hists
#if FILL_RAW_HISTS
        auto cc_amp = Histograms::cc_amp.GetThreadLocalPtr();
        auto cc_cht = Histograms::cc_cht.GetThreadLocalPtr();
        auto cc_plu = Histograms::cc_plu.GetThreadLocalPtr();
        auto cc_trt = Histograms::cc_trt.GetThreadLocalPtr();
        auto cc_mdt = Histograms::cc_mdt.GetThreadLocalPtr();
#endif // FILL_RAW_HISTS

#if FILL_CAL_HISTS
        auto cc_chE = Histograms::cc_chE.GetThreadLocalPtr();
        auto cc_sum = Histograms::cc_sum.GetThreadLocalPtr();
        auto cc_abE = Histograms::cc_abE.GetThreadLocalPtr();
#endif // FILL_CAL_HISTS

#if FILL_XTALK_CORR_HISTS
        auto cc_abM = Histograms::cc_abM.GetThreadLocalPtr();
        auto c1_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::c1_xtk[0].GetThreadLocalPtr(),
            Histograms::c1_xtk[1].GetThreadLocalPtr(),
            Histograms::c1_xtk[2].GetThreadLocalPtr(),
            Histograms::c1_xtk[3].GetThreadLocalPtr(),
            Histograms::c1_xtk[4].GetThreadLocalPtr(),
            Histograms::c1_xtk[5].GetThreadLocalPtr(),
        };
        auto c3_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::c3_xtk[0].GetThreadLocalPtr(),
            Histograms::c3_xtk[1].GetThreadLocalPtr(),
            Histograms::c3_xtk[2].GetThreadLocalPtr(),
            Histograms::c3_xtk[3].GetThreadLocalPtr(),
            Histograms::c3_xtk[4].GetThreadLocalPtr(),
            Histograms::c3_xtk[5].GetThreadLocalPtr(),
        };
        auto c5_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::c5_xtk[0].GetThreadLocalPtr(),
            Histograms::c5_xtk[1].GetThreadLocalPtr(),
            Histograms::c5_xtk[2].GetThreadLocalPtr(),
            Histograms::c5_xtk[3].GetThreadLocalPtr(),
            Histograms::c5_xtk[4].GetThreadLocalPtr(),
            Histograms::c5_xtk[5].GetThreadLocalPtr(),
        };
        auto c7_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::c7_xtk[0].GetThreadLocalPtr(),
            Histograms::c7_xtk[1].GetThreadLocalPtr(),
            Histograms::c7_xtk[2].GetThreadLocalPtr(),
            Histograms::c7_xtk[3].GetThreadLocalPtr(),
            Histograms::c7_xtk[4].GetThreadLocalPtr(),
            Histograms::c7_xtk[5].GetThreadLocalPtr(),
        };
        std::array<decltype(c1_xtk), 4> cc_xtk = {c1_xtk, c3_xtk, c5_xtk, c7_xtk};
#endif // FILL_XTALK_CORR_HISTS

// Clover Back Hists
#if FILL_RAW_HISTS
        auto cb_amp = Histograms::cb_amp.GetThreadLocalPtr();
        auto cb_cht = Histograms::cb_cht.GetThreadLocalPtr();
        auto cb_plu = Histograms::cb_plu.GetThreadLocalPtr();
        auto cb_trt = Histograms::cb_trt.GetThreadLocalPtr();
        auto cb_mdt = Histograms::cb_mdt.GetThreadLocalPtr();
#endif // FILL_RAW_HISTS

#if FILL_CAL_HISTS
        auto cb_xtE = Histograms::cb_xtE.GetThreadLocalPtr();
        auto cb_sum = Histograms::cb_sum.GetThreadLocalPtr();
        auto cb_abE = Histograms::cb_abE.GetThreadLocalPtr();
#endif // FILL_CAL_HISTS

#if FILL_XTALK_CORR_HISTS
        auto cb_abM = Histograms::cb_abM.GetThreadLocalPtr();
        auto b1_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::b1_xtk[0].GetThreadLocalPtr(),
            Histograms::b1_xtk[1].GetThreadLocalPtr(),
            Histograms::b1_xtk[2].GetThreadLocalPtr(),
            Histograms::b1_xtk[3].GetThreadLocalPtr(),
            Histograms::b1_xtk[4].GetThreadLocalPtr(),
            Histograms::b1_xtk[5].GetThreadLocalPtr(),
        };
        auto b2_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::b2_xtk[0].GetThreadLocalPtr(),
            Histograms::b2_xtk[1].GetThreadLocalPtr(),
            Histograms::b2_xtk[2].GetThreadLocalPtr(),
            Histograms::b2_xtk[3].GetThreadLocalPtr(),
            Histograms::b2_xtk[4].GetThreadLocalPtr(),
            Histograms::b2_xtk[5].GetThreadLocalPtr(),
        };
        auto b3_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::b3_xtk[0].GetThreadLocalPtr(),
            Histograms::b3_xtk[1].GetThreadLocalPtr(),
            Histograms::b3_xtk[2].GetThreadLocalPtr(),
            Histograms::b3_xtk[3].GetThreadLocalPtr(),
            Histograms::b3_xtk[4].GetThreadLocalPtr(),
            Histograms::b3_xtk[5].GetThreadLocalPtr(),
        };
        auto b5_xtk = std::array<std::shared_ptr<TH2D>, 6>{
            Histograms::b5_xtk[0].GetThreadLocalPtr(),
            Histograms::b5_xtk[1].GetThreadLocalPtr(),
            Histograms::b5_xtk[2].GetThreadLocalPtr(),
            Histograms::b5_xtk[3].GetThreadLocalPtr(),
            Histograms::b5_xtk[4].GetThreadLocalPtr(),
            Histograms::b5_xtk[5].GetThreadLocalPtr(),
        };
        std::array<decltype(b1_xtk), 4> cb_xtk = {b1_xtk, b2_xtk, b3_xtk, b5_xtk};
#endif // FILL_XTALK_CORR_HISTS

#if PROCESS_POS_SIG

#endif // PROCESS_POS_SIG

#if PROCESS_CEBR_ALL

#endif // PROCESS_CEBR_ALL

        /* #endregion */

        // Loop over the entries in the tree
        while (eventReader.Next())
        {
#if FILL_RAW_HISTS
            // Module Time
            cc_mdt->Fill(cc_mdt_val[0] * Histograms::kNsPerBin);

            cb_mdt->Fill(cb_mdt_val[0] * Histograms::kNsPerBin);

            // Trigger Times
            cc_trt->Fill(cc_trt_val[0] * Histograms::kNsPerBin, 0);
            cc_trt->Fill(cc_trt_val[1] * Histograms::kNsPerBin, 1);

            cb_trt->Fill(cb_trt_val[0] * Histograms::kNsPerBin, 0);
            cb_trt->Fill(cb_trt_val[1] * Histograms::kNsPerBin, 1);
#endif // FILL_RAW_HISTS

            // Detector Loop
            for (size_t det = 0; det < 4; det++)
            {
                std::array<double, 4> cc_xtal_E = {NAN, NAN, NAN, NAN}, cb_xtal_E = {NAN, NAN, NAN, NAN};
                std::array<double, 4> cc_xtal_T = {NAN, NAN, NAN, NAN}, cb_xtal_T = {NAN, NAN, NAN, NAN};

                // Crystal Loop
                for (size_t xtal = 0; xtal < 4; xtal++)
                {
                    auto ch = det * 4 + xtal; // Channel number 0-15

// Raw Histograms
#if FILL_RAW_HISTS
                    cc_amp->Fill(cc_amp_val[ch], ch);
                    cc_cht->Fill(cc_cht_val[ch], ch);
                    cc_plu->Fill(cc_plu_val[ch], ch);

                    cb_amp->Fill(cb_amp_val[ch], ch);
                    cb_cht->Fill(cb_cht_val[ch], ch);
                    cb_plu->Fill(cb_plu_val[ch], ch);
#endif // FILL_RAW_HISTS

                    // Calibrated Histograms
                    if (!std::isnan(cc_amp_val[ch]) &&
                        !std::isnan(cc_cht_val[ch]))
                    {
                        // std::cout << "Channel: " << ch << ", ";
                        double energy = ccECalibrate[ch](ccGainMatch[ch](cc_amp_val[ch])); // Gain-match, then calibrate
                        double cht = cc_cht_val[ch] * Histograms::kNsPerBin;
                        cc_xtal_E[xtal] = energy;
                        cc_xtal_T[xtal] = cht;
#if FILL_CAL_HISTS
                        cc_chE->Fill(energy, ch);
                        cc_sum->Fill(energy, det); // ch / 4 is the detector number
#endif                                             // FILL_CAL_HISTS
                    }

                    if (!std::isnan(cb_amp_val[ch]) &&
                        !std::isnan(cb_cht_val[ch]))
                    {
                        // std::cout << "Channel: " << ch << ", ";
                        double energy = cbECalibrate[ch](cbGainMatch[ch](cb_amp_val[ch])); // Gain-match, then calibrate
                        double cht = cb_cht_val[ch] * Histograms::kNsPerBin;
                        cb_xtal_E[xtal] = energy;
                        cb_xtal_T[xtal] = cht;
#if FILL_CAL_HISTS
                        cb_xtE->Fill(energy, ch);
                        cb_sum->Fill(energy, det); // ch / 4 is the detector number
#endif                                             // FILL_CAL_HISTS
                    }
                } // End Crystal Loop

                // Add-Back Histograms
                if (std::any_of(cc_xtal_E.begin(), cc_xtal_E.end(), [](double x)
                                { return x > CAAddBack::kAddBackThreshold; }))
                {
#if FILL_XTALK_CORR_HISTS
                    // printf("[DEBUG] cc_xtal_E: [%.2f, %.2f, %.2f, %.2f]\n", cc_xtal_E[0], cc_xtal_E[1], cc_xtal_E[2], cc_xtal_E[3]);
                    unsigned int cc_mult = std::count_if(cc_xtal_E.begin(), cc_xtal_E.end(), [](double x)
                                                         { return !std::isnan(x); });
                    cc_abM->Fill(cc_mult, det);
                    if (cc_mult == 2)
                        CACrosstalkCorrection::FillXTalkHistograms(cc_xtk[det], cc_xtal_E, cc_xtal_T);
#endif // FILL_XTALK_CORR_HISTS
#if FILL_CAL_HISTS
                    auto energies_corr = xTalkCorrection[4 + det](cc_xtal_E);

                    // printf("[DEBUG] cc_xtal_C: [%.2f, %.2f, %.2f, %.2f]\n", energies_corr[0], energies_corr[1], energies_corr[2], energies_corr[3]);

                    auto energy_ab_corr = CAAddBack::GetAddBackEnergy(energies_corr, cc_xtal_T);
                    cc_abE->Fill(energy_ab_corr, det);
// std::cout << "[DEBUG] Done with add-back event" << std::endl;
#endif // FILL_CAL_HISTS
                }
                if (std::any_of(cb_xtal_E.begin(), cb_xtal_E.end(), [](double x)
                                { return x > CAAddBack::kAddBackThreshold; }))
                {
#if FILL_XTALK_CORR_HISTS
                    unsigned int cb_mult = std::count_if(cb_xtal_E.begin(), cb_xtal_E.end(), [](double x)
                                                         { return !std::isnan(x); });
                    cb_abM->Fill(cb_mult, det);
                    if (cb_mult == 2)
                        CACrosstalkCorrection::FillXTalkHistograms(cb_xtk[det], cb_xtal_E, cb_xtal_T);
                    // printf("[DEBUG] cb_xtal_E: [%.2f, %.2f, %.2f, %.2f]\n", cb_xtal_E[0], cb_xtal_E[1], cb_xtal_E[2], cb_xtal_E[3]);
#endif // FILL_XTALK_CORR_HISTS
#if FILL_CAL_HISTS

                    auto energies_corr = xTalkCorrection[det](cb_xtal_E);

                    // printf("[DEBUG] cc_xtal_C: [%.2f, %.2f, %.2f, %.2f]\n", energies_corr[0], energies_corr[1], energies_corr[2], energies_corr[3]);

                    auto energy_ab_corr = CAAddBack::GetAddBackEnergy(energies_corr, cb_xtal_T);
                    cb_abE->Fill(energy_ab_corr, det);
// std::cout << "[DEBUG] Done with add-back event" << std::endl;
#endif // FILL_CAL_HISTS
                }

            } // End Detector Loop

#if PROCESS_POS_SIG

#endif // PROCESS_POS_SIG

#if PROCESS_CEBR_ALL

#endif // PROCESS_CEBR_ALL

            processedEntries++;
        } // End Event Loop
    };

    // Loop over the entries in the TTree and fill the histograms appropriately
    TStopwatch timer;
    timer.Start();
    EventProcessor.Process(fillHistograms);
    timer.Stop();

    progressBarThread.join();

    printf("[INFO] Processed events in %.2f seconds (%.2f events/second)\n",
           timer.RealTime(),
           static_cast<double>(processedEntries) / timer.RealTime());

    // Save the histograms to a new ROOT file
    TFile* outfile = new TFile(args.outputFileName.c_str(), "RECREATE");
    if (!outfile || outfile->IsZombie())
    {
        throw std::runtime_error(Form("[ERROR] Error creating output file: %s", args.outputFileName.c_str()));
    }

    /* #region Write Histograms */

    // Use histograms defined in Histograms.hpp

    // Clover Cross Histograms
    auto cc_dir = outfile->mkdir("clover_cross");
    cc_dir->cd();
#if FILL_RAW_HISTS
    Histograms::cc_amp.Write();
    Histograms::cc_cht.Write();
    Histograms::cc_plu.Write();
    Histograms::cc_trt.Write();
    Histograms::cc_mdt.Write();
#endif // FILL_RAW_HISTS
#if FILL_CAL_HISTS
    Histograms::cc_chE.Write();
    Histograms::cc_sum.Write();
    Histograms::cc_abE.Write();
#endif // FILL_CAL_HISTS

#if FILL_XTALK_CORR_HISTS
    Histograms::cc_abM.Write();
    for (int i = 0; i < 6; i++)
    {
        Histograms::c1_xtk[i].Write();
        Histograms::c3_xtk[i].Write();
        Histograms::c5_xtk[i].Write();
        Histograms::c7_xtk[i].Write();
    }
#endif // FILL_XTALK_CORR_HISTS
    outfile->cd();

    // Clover Back Histograms
    auto cb_dir = outfile->mkdir("clover_back");
    cb_dir->cd();
#if FILL_RAW_HISTS
    Histograms::cb_amp.Write();
    Histograms::cb_cht.Write();
    Histograms::cb_plu.Write();
    Histograms::cb_trt.Write();
    Histograms::cb_mdt.Write();
#endif // FILL_RAW_HISTS
#if FILL_CAL_HISTS
    Histograms::cb_xtE.Write();
    Histograms::cb_sum.Write();
    Histograms::cb_abE.Write();
#endif // FILL_CAL_HISTS

#if FILL_XTALK_CORR_HISTS
    Histograms::cb_abM.Write();
    for (int i = 0; i < 6; i++)
    {
        Histograms::b1_xtk[i].Write();
        Histograms::b2_xtk[i].Write();
        Histograms::b3_xtk[i].Write();
        Histograms::b5_xtk[i].Write();
    }
#endif // FILL_XTALK_CORR_HISTS
    outfile->cd();

// Positive Signal Histograms
#if PROCESS_POS_SIG

#endif // PROCESS_POS_SIG

    // CeBr All Histograms
#if PROCESS_CEBR_ALL

#endif // PROCESS_CEBR_ALL

    /* #endregion */

    printf("[INFO] Saved histograms to file: %s\n", args.outputFileName.c_str());

    outfile->Close();
    delete outfile;

    printf("[INFO] Done!\n");

    return 0;
}