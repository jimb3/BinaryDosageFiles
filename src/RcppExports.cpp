// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// VCF_to_BinaryDosage
int VCF_to_BinaryDosage(std::string vcfFilename, std::string outBaseFilename, unsigned int initSub);
RcppExport SEXP BinaryDosageFiles_VCF_to_BinaryDosage(SEXP vcfFilenameSEXP, SEXP outBaseFilenameSEXP, SEXP initSubSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type vcfFilename(vcfFilenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type outBaseFilename(outBaseFilenameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type initSub(initSubSEXP);
    rcpp_result_gen = Rcpp::wrap(VCF_to_BinaryDosage(vcfFilename, outBaseFilename, initSub));
    return rcpp_result_gen;
END_RCPP
}
// ExtractDosages
Rcpp::List ExtractDosages(std::string bdosageFilename, std::string mapFilename, unsigned int numSub, unsigned int numSNPs);
RcppExport SEXP BinaryDosageFiles_ExtractDosages(SEXP bdosageFilenameSEXP, SEXP mapFilenameSEXP, SEXP numSubSEXP, SEXP numSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bdosageFilename(bdosageFilenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type mapFilename(mapFilenameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numSNPs(numSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(ExtractDosages(bdosageFilename, mapFilename, numSub, numSNPs));
    return rcpp_result_gen;
END_RCPP
}
// ExtractMoreDosages
Rcpp::List ExtractMoreDosages(Rcpp::List inputs);
RcppExport SEXP BinaryDosageFiles_ExtractMoreDosages(SEXP inputsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type inputs(inputsSEXP);
    rcpp_result_gen = Rcpp::wrap(ExtractMoreDosages(inputs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"BinaryDosageFiles_VCF_to_BinaryDosage", (DL_FUNC) &BinaryDosageFiles_VCF_to_BinaryDosage, 3},
    {"BinaryDosageFiles_ExtractDosages", (DL_FUNC) &BinaryDosageFiles_ExtractDosages, 4},
    {"BinaryDosageFiles_ExtractMoreDosages", (DL_FUNC) &BinaryDosageFiles_ExtractMoreDosages, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BinaryDosageFiles(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
