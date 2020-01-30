int extract_result(const char* dir) {
    TString path(dir);
    path += "/fit_results.root";
    TFile file(path);
    TMacro* macro = (TMacro*)file.Get("all;1");
    /* cout << macro->GetLineWith("n_signal_plus_cross_model")->GetString(); */
    macro->Print();
    macro->SaveSource("current_result");
    return 0;
}
