#include <RcppArmadillo.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include "BinaryDosage.h"

using namespace Rcpp;

void SaveSNP(CBinaryDosage &bd, unsigned int n, arma::mat &d, arma::mat &p0, arma::mat &p1, arma::mat &p2,
             std::vector<std::string> &chromosome, std::vector<std::string> &snpName, std::vector<int> &basepairs,
             std::vector<std::string> &allele1, std::vector<std::string> &allele2) {
  d.col(n) = bd.Dosage();
  if (bd.Probabilities()) {
    p0.col(n) = bd.Probs().col(0);
    p1.col(n) = bd.Probs().col(1);
    p2.col(n) = bd.Probs().col(2);
  }
  chromosome.push_back(bd.Chromosome());
  snpName.push_back(bd.SNPName());
  basepairs.push_back(bd.Location());
  allele1.push_back(bd.FirstAllele());
  allele2.push_back(bd.SecondAllele());
}
//' Function to convert a VCF file to a binary dosage file for GxEScan
//' 
//' Function to convert a VCF file to a binary dosage file for GxEScan. 
//' The binary formatted file is smaller that the VCF file and reads in much quicker.
//' 
//' @param vcfFilename
//' Name of VCF file
//' @param outBaseFilename
//' Base filename of output files. Three files are output.
//' A family file with extension .fam.
//' A map file with extension .bim.
//' A binary dosage file with extenstion .bdosage.
//' @param initSub
//' Amount of memory to allocate for subjects. Should be the number of subjects.
//' Additional memory will be allocated if not enough is initially allocated.
//' Having to do this can slow down the speed of the routine because of memory
//' allocations/deallocations.
//' @return
//' 0 failure
//' otherwise number of subjects read 
//' @importFrom Rcpp evalCpp
//' @useDynLib GxEScanR
//' @export
// [[Rcpp::export]]
int VCF_to_BinaryDosage(std::string vcfFilename, std::string outBaseFilename, unsigned int initSub) {
  std::ifstream vcfFile;
  std::ofstream bdosefile;
  std::ofstream bimfile;
  std::ofstream famfile;
  std::string readLine;
  std::string junk;
  std::string colName;
  std::istringstream iss;
  std::vector<std::string> iid;
  std::string riid;
  std::string genotype;
  std::string chromosome, snpName, refAllele, altAllele;
  unsigned int location;
  double dosage, dosagec;
  double p0, p1, p2;
  double ptest;
  double psum;
  unsigned int numSub;
  unsigned int numSNPs;
  unsigned int ui;
  unsigned int num1;
  short *sdosage = NULL;
  short *sp1 = NULL;
  unsigned int oversum;
  
  vcfFile.open(vcfFilename.c_str());
  if (!vcfFile.good()) {
    std::cerr << "Unable to open VCF file" << std::endl;
    return 0;
  }

  getline(vcfFile, readLine);
  ui = 1;
  // Should check if first line indicates this is a vcf file
  // and in a supported format
  while (readLine[1] == '#' && ui < 25) {
    ++ui;
    getline(vcfFile, readLine);
  }
  iid.reserve(initSub);
  iss.str(readLine.substr(1));
  for (ui = 0; ui < 9; ++ui) {
    iss >> colName;
  }
  
  famfile.open((outBaseFilename + ".fam").c_str());
  numSub = 0;
  iss >> riid;
  while (!iss.fail()) {
    ++numSub;
    famfile << numSub << '\t' << riid << "\t0\t0\t9\t9" << std::endl;
    iid.push_back(riid);
    iss >> riid;
  }
  famfile.close();

  sdosage = new short[numSub];
  sp1 = new short[numSub + numSub + numSub];
  
  numSNPs = 0;
  bdosefile.open((outBaseFilename + ".bdosage").c_str(), std::ios::out | std::ios::binary);
  char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x3, 0x0, 0x1 };
  bdosefile.write(header, 8);
  bdosefile.write((const char *)&numSub, 4);
  
  bimfile.open((outBaseFilename + ".bim").c_str());
  vcfFile >> chromosome >> location >> snpName >> refAllele >> altAllele >> junk >> junk >> junk >> junk;
  
  while (!vcfFile.fail()) {
    ++numSNPs;
    bimfile << chromosome << '\t' << snpName << "\t0\t" << location << '\t' << refAllele << '\t' << altAllele << std::endl;
    std::memset(sdosage, 0, numSub * sizeof(short));
    std::memset(sp1, 0, 3 * numSub * sizeof(short));
    num1 = 0;
    oversum = 0;
    for (ui = 0; ui < numSub; ++ui) {
      vcfFile >> readLine;
      std::replace(readLine.begin(), readLine.end(), ':', ' ');
      std::replace(readLine.begin(), readLine.end(), ',', ' ');
      iss.str(readLine);
      iss.clear();
      iss >> genotype;
      iss >> dosage >> p0 >> p1 >> p2;
      psum = p0 + p1 + p2;
      sdosage[ui] = (short)((dosage + 0.00001) * 10000);
      ptest = sdosage[ui] / 10000.;
      if (ptest != dosage) {
        std::cerr << "Multiplication failure, dosage\t" << dosage << '\t' << sdosage[ui] << '\t' << p0 << '\t' << p1 << '\t' << p2 << '\t' << ptest << '\t' << sp1[num1] << std::endl;
        return 0;
      }
      dosagec = p1 + p2 + p2;
      if ((p2 != 0 && p0 != 0 && p1 != 0) || psum != 1 || dosagec != dosage) {
        sdosage[ui] |= 0x8000;
        sp1[num1] = (short)((p1 + 0.00001) * 10000);
        ptest = sp1[num1] / 10000.;
        if (ptest != p1) {
          std::cerr << "Multiplication failure, p1\t" << dosage << '\t' << sdosage[ui] << '\t' << p0 << '\t' << p1 << '\t' << p2 << '\t' << ptest << '\t' << sp1[num1] << std::endl;
          return 0;
        }
        if (psum != 1 || dosagec != dosage) {
          ++oversum;
          sp1[num1] |= 0x8000;
          ++num1;
          sp1[num1] = short((p0 + 0.00001) * 10000);
          ptest = sp1[num1] / 10000.;
          if (ptest != p0) {
            std::cerr << "Multiplication failure, p0\t" << dosage << '\t' << sdosage[ui] << '\t' << p0 << '\t' << p1 << '\t' << p2 << '\t' << ptest << '\t' << sp1[num1] << std::endl;
            return 0;
          }
          ++num1;
          sp1[num1] = short((p2 + 0.00001) * 10000);
          ptest = sp1[num1] / 10000.;
          if (ptest != p2) {
            std::cerr << "Multiplication failure, p2\t" << dosage << '\t' << sdosage[ui] << '\t' << p0 << '\t' << p1 << '\t' << p2 << '\t' << ptest << '\t' << sp1[num1] << std::endl;
            return 0;
          }
        }
        ++num1;
      }
      if (iss.fail()) {
        std::cerr << "Read failure" << std::endl;
        return 0;
      }
    }
    if ((numSNPs % 1000) == 0)
      std::cout << numSNPs << std::endl;
    bdosefile.write((const char *)sdosage, numSub + numSub);
    bdosefile.write((const char *)sp1, num1 + num1);
    vcfFile >> chromosome >> location >> snpName >> refAllele >> altAllele >> junk >> junk >> junk >> junk;
  }

  if (sdosage)
    delete[] sdosage;
  if (sp1)
    delete[] sp1;
  
  bdosefile.close();
  bimfile.close();
  vcfFile.close();
  return numSub;
}

//' Function to extract SNPs from a binary dosage file
//' 
//' Function to extract SNPs from a binary dosage file
//' 
//' @param bdosageFilename
//' Name of binary dosage file
//' @param mapFilename
//' Name of map file associated with dosage file
//' @param numSub
//' Number of subjects with data in dosage file
//' @param numSNPs
//' Number of SNPs to read in
//' @return
//' List with a vector of dosages and a matrix of probabilities
//' and a list of the input values
//' @importFrom Rcpp evalCpp
//' @useDynLib GxEScanR
//' @export
// [[Rcpp::export]]
Rcpp::List ExtractDosages(std::string bdosageFilename, std::string mapFilename, unsigned int numSub, unsigned int numSNPs) {
  CBinaryDosage bd;
  unsigned int ui;
  int ret;
  Rcpp::List inputs;
  Rcpp::List mapData;
  unsigned int numRead;
  std::streampos x;
  int *y;

  arma::uvec snpIDs;
  arma::mat dosages;
  arma::mat p0;
  arma::mat p1;
  arma::mat p2;
  arma::ivec stp;

  if (bd.ReadFile(bdosageFilename, numSub, mapFilename))
    return inputs;
  
  snpIDs.zeros(numSNPs);
  dosages.zeros(numSub, numSNPs);
  if (bd.Probabilities() == true) {
    p0.zeros(numSub, numSNPs);
    p1.zeros(numSub, numSNPs);
    p2.zeros(numSub, numSNPs);
  }
  
  for (ui = 0; ui < numSNPs;) {
    if (ui == 0)
      ret = bd.GetFirst();
    else
      ret = bd.GetNext();
    if (ret != 0)
      break;
    snpIDs(ui) = bd.CurrentSNP() + 1;
    dosages.col(ui) = bd.Dosage();
    if (bd.Probabilities() == true) {
      p0.col(ui) = bd.Probs().col(0);
      p1.col(ui) = bd.Probs().col(1);
      p2.col(ui) = bd.Probs().col(2);
    }
    ++ui;
  }
  numRead = ui;
  stp.set_size(sizeof(x) / sizeof(int));
  x = bd.FilePosition();
  y = (int *)&x;
  for (ui = 0; ui < sizeof(x) / sizeof(int); ++ui)
    stp(ui) = y[ui];

  mapData = Rcpp::List::create(
    Rcpp::Named("Chromosome") = bd.MapFile().Chromosome(),
    Rcpp::Named("SNP") = bd.MapFile().SNP(),
    Rcpp::Named("BasePairs") = bd.MapFile().BasePairs(),
    Rcpp::Named("Allele1") = bd.MapFile().FirstAllele(),
    Rcpp::Named("Allele2") = bd.MapFile().SecondAllele(),
    Rcpp::Named("Skipped") = bd.MapFile().Skipped() );
  inputs = Rcpp::List::create(
    Rcpp::Named("filename") = bdosageFilename,
    Rcpp::Named("StreamPos") = stp,
    Rcpp::Named("NumSub") = numSub,
    Rcpp::Named("NumSNPs") = numSNPs,
    Rcpp::Named("CurrentSNP") = bd.CurrentSNP(),
    Rcpp::Named("Version") = bd.Version(),
    Rcpp::Named("MapData") = mapData);
  
  if (bd.Probabilities() == true) {
    return Rcpp::List::create(
      Rcpp::Named("SNPID") = snpIDs,
      Rcpp::Named("Dosages") = dosages,
      Rcpp::Named("NumRead") = numRead,
      Rcpp::Named("P0") = p0,
      Rcpp::Named("P1") = p1,
      Rcpp::Named("P2") = p2,
      Rcpp::Named("Inputs") = inputs);
  }
  return Rcpp::List::create(
    Rcpp::Named("SNPID") = snpIDs,
    Rcpp::Named("Dosages") = bd.Dosage(),
    Rcpp::Named("NumRead") = numRead,
    Rcpp::Named("Inputs") = inputs);
}

//' Function to extract more SNPs from a binary dosage file
//' 
//' Function to extract more SNPs from a binary dosage file
//' 
//' @param inputs
//' List of inputs returned from ExtractDosages
//' @return
//' List with a vector of dosages and a matrix of probabilities
//' and a list of the input values
//' @importFrom Rcpp evalCpp
//' @useDynLib GxEScanR
//' @export
// [[Rcpp::export]]
Rcpp::List ExtractMoreDosages(Rcpp::List inputs) {
  CBinaryDosage bd;
  std::string filename;
  arma::ivec iLocation;
  arma::ivec iSkipped;
  std::vector<bool> skipped;
  Rcpp::List mapData;
  std::streampos loc;
  int *lp;
  unsigned int numSub;
  unsigned int numSNPs;
  unsigned int currentSNP;
  unsigned int numRead;
  double version;
  unsigned int ui;
  std::streampos x;
  int *y;
  arma::uvec snpIDs;
  arma::mat dosages;
  arma::mat p0;
  arma::mat p1;
  arma::mat p2;
  arma::ivec stp;
  
  filename = as<std::string>(inputs["filename"]);
  iLocation = as<IntegerVector>(inputs["StreamPos"]);
  numSub = as<unsigned int>(inputs["NumSub"]);
  numSNPs = as<unsigned int>(inputs["NumSNPs"]);
  currentSNP = as<unsigned int>(inputs["CurrentSNP"]);
  version = as<double>(inputs["Version"]);
  mapData = as<Rcpp::List>(inputs["MapData"]);
  iSkipped = as<IntegerVector>(mapData["Skipped"]);
  skipped.reserve(iSkipped.size());
  for (ui = 0; ui < iSkipped.size(); ++ui) {
    if (iSkipped(ui) == 1)
      skipped.push_back(true);
    else
      skipped.push_back(false);
  }
  lp = (int *)&loc;
  for (ui = 0; ui < sizeof(std::streampos) / sizeof(int); ++ui)
    lp[ui] = iLocation(ui);
  
  if (bd.ReOpen(filename, numSub, skipped, loc, currentSNP, version) != 0) {
    std::cout << "Failed to reopen" << std::endl;
    return inputs;
  }
  snpIDs.zeros(numSNPs);
  dosages.zeros(numSub, numSNPs);
  if (bd.Probabilities() == true) {
    p0.zeros(numSub, numSNPs);
    p1.zeros(numSub, numSNPs);
    p2.zeros(numSub, numSNPs);
  }
  
  for (ui = 0; ui < numSNPs;) {
    if (bd.GetNext() != 0)
      break;
    snpIDs(ui) = bd.CurrentSNP() + 1;
    dosages.col(ui) = bd.Dosage();
    if (bd.Probabilities() == true) {
      p0.col(ui) = bd.Probs().col(0);
      p1.col(ui) = bd.Probs().col(1);
      p2.col(ui) = bd.Probs().col(2);
    }
    ++ui;
  }
  numRead = ui;
  stp.set_size(sizeof(x) / sizeof(int));
  x = bd.FilePosition();
  y = (int *)&x;
  for (ui = 0; ui < sizeof(x) / sizeof(int); ++ui)
    stp(ui) = y[ui];
  
  inputs = Rcpp::List::create(
    Rcpp::Named("filename") = filename,
    Rcpp::Named("StreamPos") = stp,
    Rcpp::Named("NumSub") = numSub,
    Rcpp::Named("NumSNPs") = numSNPs,
    Rcpp::Named("CurrentSNP") = bd.CurrentSNP(),
    Rcpp::Named("Version") = bd.Version(),
    Rcpp::Named("MapData") = mapData);
  
  if (bd.Probabilities() == true) {
    return Rcpp::List::create(
      Rcpp::Named("SNPID") = snpIDs,
      Rcpp::Named("Dosages") = dosages,
      Rcpp::Named("NumRead") = numRead,
      Rcpp::Named("P0") = p0,
      Rcpp::Named("P1") = p1,
      Rcpp::Named("P2") = p2,
      Rcpp::Named("Inputs") = inputs);
  }
  return Rcpp::List::create(
    Rcpp::Named("SNPID") = snpIDs,
    Rcpp::Named("Dosages") = bd.Dosage(),
    Rcpp::Named("NumRead") = numRead,
    Rcpp::Named("Inputs") = inputs);
}

//' Function to extract a SNP from a binary dosage file
//' 
//' Function to extract a SNP from a binary dosage file
//' 
//' @param bdosageFilename
//' Name of binary dosage file
//' @param mapFilename
//' Name of map file associated with dosage file
//' @param numSub
//' Number of subjects with data in dosage file
//' @param snpName
//' Name of SNP to extract
//' @param flanking
//' Number of flanking SNPs on either side to include
//' @return
//' List with a vector of dosages and a matrix of probabilities
//' and a list of the input values
//' @importFrom Rcpp evalCpp
//' @useDynLib GxEScanR
//' @export
// [[Rcpp::export]]
Rcpp::List ExtractSNPDosages(std::string bdosageFilename, std::string mapFilename, unsigned int numSub, std::string snpName, unsigned int flanking) {
  CBinaryDosage bd;
  Rcpp::List dosages;
  unsigned int ui, uj, uk;
  unsigned int firstSNP;
  unsigned int lastSNP;
  unsigned int snpNum;
  unsigned numOut;
  Rcpp::DataFrame df;
  std::vector<std::string> chromosome;
  std::vector<std::string> snp;
  std::vector<int> location;
  std::vector<std::string> allele1;
  std::vector<std::string> allele2;
  
  arma::mat d;
  arma::mat p0;
  arma::mat p1;
  arma::mat p2;
  
  if (bd.ReadFile(bdosageFilename, numSub, mapFilename))
    return dosages;
  
  uj = 0;
  for (ui = 0; ui < bd.MapFile().NumSNPs(); ++ui) {
    if (bd.MapFile().Skipped()[ui] == false)
      ++uj;
    if (bd.MapFile().SNP()[ui] == snpName)
      break;
  }
  
  if (ui == bd.MapFile().NumSNPs())
    return dosages;
  if (bd.MapFile().Skipped()[ui] == true)
    return dosages;
  
  snpNum = uj;

  firstSNP = snpNum < flanking ? 0 : snpNum - flanking;
  lastSNP = snpNum + flanking  + 1 > bd.MapFile().NumUsed() ? bd.MapFile().NumUsed() : snpNum + flanking + 1;

  d.zeros(numSub, lastSNP - firstSNP);
  if (bd.Probabilities()) {
    p0.zeros(numSub, lastSNP - firstSNP);
    p1.zeros(numSub, lastSNP - firstSNP);
    p2.zeros(numSub, lastSNP - firstSNP);
  }
  
  bd.GetFirst();
  for (ui = 0; ui < firstSNP; ++ui)
    bd.GetNext();
  
  uj = 0;
  for (uj = 0; ui < lastSNP; ++ui, ++uj, bd.GetNext())
    SaveSNP(bd, uj, d, p0, p1, p2, chromosome, snp, location, allele1, allele2);

  df = DataFrame::create(_["Chromosome"] = chromosome,
                         _["SNPName"] = snp,
                         _["BasePairs"] = location,
                         _["Allele1"] = allele1,
                         _["Allele2"] = allele2);
  if (bd.Probabilities())
    dosages = Rcpp::List::create(Named("Dosage") = d,
                                 Named("P0") = p0,
                                 Named("P1") = p1,
                                 Named("P2") = p2,
                                 Named("SNPInfo") = df);
  else
  dosages = Rcpp::List::create(Named("Dosage") = d,
                               Named("SNPInfo") = df);
  return dosages;
}

// ********************************************************
//                CBinaryDosage
// Class to read in data form binary dosage files
// ********************************************************
// Constructor
CBinaryDosage::CBinaryDosage() : CBinaryGeneticData() {
	m_bMeasured = false;
	m_headerSize = 8;
	m_numSNPsUsed = 0;
	m_versionNumber = 0;
}
// Destructor
CBinaryDosage::~CBinaryDosage() {
}

// Read the file head and determine if file is valid
int CBinaryDosage::ReadFile(const std::string &_filename) {
//	unsigned long long actualSize;
//	unsigned long long expectedSize;
	const char bose[4] = { 'b', 'o', 's', 'e' };
	const char version_1_1[4] = { 0x0, 0x1, 0x0, 0x1 };
	const char version_1_2[4] = { 0x0, 0x1, 0x0, 0x2 };
	const char version_3_1[4] = { 0x0, 0x3, 0x0, 0x1 };
	char header[8];
	unsigned int expNumSub;

	if (OpenFile(_filename, true))
		return 1;

	m_infile.read(header, 4);
	if (memcmp(header, bose, 4)) {
		std::cerr << "File is not a binary dosage file" << std::endl;
		return 1;
	}
	m_infile.read(header, 4);
	if (memcmp(header, version_1_1, 4) == 0) {
		m_bProbabilities = false;
		m_versionNumber = 1.1;
	}
	else if (memcmp(header, version_1_2, 4) == 0) {
		m_bProbabilities = true;
		m_versionNumber = 1.2;
	}
	else if (memcmp(header, version_3_1, 4) == 0) {
	  m_bProbabilities = true;
	  m_versionNumber = 3.1;
	  m_infile.read((char *)&expNumSub, 4);
	  if (expNumSub != m_numSubjects) {
	    std::cerr << "Number of subjects in file doesn't agree with number of subjects indicated" << std::endl;
	    return 1;
	  }
	  m_headerSize = 12;
	}
	else {
		std::cerr << "File is not a binary dosage file" << std::endl;
		return 1;
	}
// Commented out because the code is OS specific
/*
	if (m_versionNumber == 1.1)
		expectedSize = (unsigned long long)(m_numSubjects)* m_numSNPs * 2 + 8;
	else
		expectedSize = (unsigned long long)(m_numSubjects)* m_numSNPs * 4 + 8;

	m_infile.seekg(0, ios_base::end);
	actualSize = FileSize(_filename.c_str());
//	cout << expectedSize << '\t' << actualSize << endl;
	if (actualSize != expectedSize) {
		cout << "Input file not of expected size" << endl;
		return 1;
	}
*/
  m_Dosage.set_size(m_numSubjects);
	
	if (m_versionNumber == 1.1)
		m_Probs.set_size(0, 0);
	else
	  m_Probs.set_size(m_numSubjects, 3);

	m_arraySize = m_numSubjects * sizeof(unsigned short);
	if (m_versionNumber == 1.2)
		m_arraySize += m_numSubjects * sizeof(unsigned short);
	
	m_SNPArray.resize((int)m_arraySize / 2);
	
	return 0;
}

// Read file and Map file
int CBinaryDosage::ReadFile(const std::string &_filename, unsigned int _numSub, const std::string &_mapFileName) {
	if (m_mapFile.ReadFile(_mapFileName)) {
		m_errorString = m_mapFile.ErrorString();
		return 1;
	}
  m_Skipped = m_mapFile.Skipped();
	m_numSubjects = _numSub;
	m_numSNPs = m_mapFile.NumSNPs();
	m_numSNPsUsed = m_mapFile.NumUsed();

	if (m_numSNPsUsed == 0) {
		m_errorString = "No SNPs to read";
		return 1;
	}

	return ReadFile(_filename);
}

int CBinaryDosage::ReOpen(const std::string &_filename, unsigned int _numSub, std::vector<bool> &_skipped,
           std::streampos _stp, unsigned int _currentSNP, double _versionNumber) {
  m_numSubjects = _numSub;
  m_Skipped = _skipped;
  m_numSNPs = m_Skipped.size();
  m_currentSNP = _currentSNP;
  m_versionNumber = _versionNumber;
  if (OpenFile(_filename, true))
    return 1;
  
  m_Dosage.set_size(m_numSubjects);
  
  if (m_versionNumber == 1.1) {
    m_Probs.set_size(0, 0);
    m_bProbabilities = false;
  } else {
    m_Probs.set_size(m_numSubjects, 3);
    m_bProbabilities = true;
  }

  m_arraySize = m_numSubjects * sizeof(unsigned short);
  if (m_versionNumber == 1.2)
    m_arraySize += m_numSubjects * sizeof(unsigned short);
  
  m_SNPArray.resize((int)m_arraySize / 2);

  m_infile.seekg(_stp);
  return 0;
}

// Process SNP if Version 1.1
void CBinaryDosage::ProcessSNP11() {
  m_Dosage = (arma::conv_to<arma::colvec>::from(m_SNPArray));
  m_Dosage.replace(65535, NA_REAL);
  m_Dosage /= 32767;
}

// Process SNP if Version 1.2
void CBinaryDosage::ProcessSNP12() {
  arma::Mat<unsigned short> rd(&m_SNPArray[0], 2, m_numSubjects, false, true);
	
	m_Probs.zeros();
	m_Probs.col(1) = arma::conv_to<arma::colvec>::from(rd.row(0));
	m_Probs.col(2) = arma::conv_to<arma::colvec>::from(rd.row(1));
	m_Probs.replace(65535, NA_REAL);
	m_Probs.col(0) = 65535 - (m_Probs.col(1) + m_Probs.col(2));
	m_Probs /= 65535.;
	m_Dosage = m_Probs.col(1) + m_Probs.col(2) + m_Probs.col(2);
}
// Process SNP if Version 3.1
void CBinaryDosage::ProcessSNP31() {
  unsigned int ui;
  unsigned short sp1;

  m_Probs.zeros();
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if ((m_SNPArray(ui) & 0x8000) != 0) {
      m_SNPArray(ui) &= 0x7FFF;
      m_Dosage(ui) = m_SNPArray(ui) / 10000.;
      m_infile.read((char *)&sp1, 2);
      if ((sp1 & 0x8000) != 0) {
        sp1 &= 0x7FFF;
        m_Probs(ui, 1) = sp1 / 10000.;
        m_infile.read((char *)&sp1, 2);
        m_Probs(ui, 0) = sp1 / 10000.;
        m_infile.read((char *)&sp1, 2);
        m_Probs(ui, 2) = sp1 / 10000.;
      }
      else {
        m_Probs(ui, 1) = sp1 / 10000.;
        m_Probs(ui, 2) = (m_Dosage[ui] - m_Probs(ui, 1)) / 2;
        m_Probs(ui, 0) = 1 - m_Probs(ui, 1) - m_Probs(ui, 2);
      }
    }
    else {
      m_Dosage(ui) = m_SNPArray(ui) / 10000.;
      if (m_Dosage(ui) < 1) {
        m_Probs(ui, 1) = m_Dosage(ui);
        m_Probs(ui, 0) = 1 -  m_Probs(ui, 1);
        m_Probs(ui, 2) = 0;
      }
      else {
        m_Probs(ui, 2) = m_Dosage(ui) - 1;
        m_Probs(ui, 1) = 1 -  m_Probs(ui, 2);
        m_Probs(ui, 0) = 0;
      }
    }
  }
}
// Process the data read in for the last SNP
void CBinaryDosage::ProcessSNP() {
	if (m_versionNumber == 1.1)
		ProcessSNP11();
	else if (m_versionNumber == 1.2)
		ProcessSNP12();
	else
	  ProcessSNP31();
}