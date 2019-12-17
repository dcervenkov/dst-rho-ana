int extract_result() {
    TFile file("../plots/test/fit_results.root");
    TMacro* macro = (TMacro*)file.Get("all;1");
    /* cout << macro->GetLineWith("n_signal_plus_cross_model")->GetString(); */
    macro->Print();
    macro->SaveSource("current_result");
    return 0;
}
