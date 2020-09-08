//     Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte


//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.



// #include "OncoSimul.h"
// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
#include "debug_common.h"
#include "common_classes.h"
#include "bnb_common.h"
#include "new_restrict.h"
#include <cfloat>
#include <limits>
#include <Rcpp.h>
#include <iostream>
#include <random>
#include <set>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <ctime>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <time.h>
using namespace std;

using namespace Rcpp;
using std::vector;

// To track if mutation is really much smaller than birth/death
#define MIN_RATIO_MUTS_NR
#ifdef MIN_RATIO_MUTS_NR
// There is really no need for these to be globals?
// Unless I wanted to use them inside some function. So leave as globals.
double g_min_birth_mut_ratio_nr = DBL_MAX;
double g_min_death_mut_ratio_nr = DBL_MAX;
double g_tmp1_nr = DBL_MAX;
#endif

void sleep(unsigned int mseconds){
  clock_t goal = mseconds + clock();
  while (goal > clock())
    ;
}

void nr_fitness(spParamsP& tmpP,
		const spParamsP& parentP,
		const Genotype& ge,
		const fitnessEffectsAll& F,
		const TypeModel typeModel,
		std::vector<Genotype>& Genotypes,
		std::vector<spParamsP>& popParams) {
		       // const double& genTime,
		       // const double& adjust_fitness_B,
		       // const double& adjust_fitness_MF) {

  // We want a way to signal immediate non-viability of a clone. For
  // "birth-based" models that happens when any s = -1, as the fitness is
  // 0. By setting birth = 0.0 we ensure this clone does not get added and
  // we never reach into algo2, etc, leading to numerical problems.

  // With Bozic models, which are "death-based", it is different. For
  // bozic2, birth is bounded, so any death > 2 would lead to birth <
  // 0. For bozic1, deaths of around 50 lead to numerical issues.  The
  // general rule is: set those mutations to -inf, so prodDeathFitness
  // returns an inf for death, and that is recognized as "no
  // viability" (anything with death > 99)

  // The ones often used are bozic1, exp, mcfarlandlog

	if(F.frequencyDependentFitness){
		popParams.push_back(tmpP);
		Genotypes.push_back(ge);
	}

  if(typeModel == TypeModel::bozic1) {
    tmpP.death = prodDeathFitness(evalGenotypeFitness(ge, F, Genotypes, popParams));
    if( tmpP.death > 99) {
      tmpP.birth = 0.0;
    } else {
      tmpP.birth = 1.0;
    }
  // } else if (typeModel == TypeModel::bozic2) {
  //   double pp = prodDeathFitness(evalGenotypeFitness(ge, F));
  //   tmpP.birth = std::max(0.0, (1.0/genTime) * (1.0 - 0.5 * pp ));
  //   tmpP.death = (0.5/genTime) * pp;
  } else {
    double fitness = prodFitness(evalGenotypeFitness(ge, F, Genotypes, popParams));
    if( fitness <= 0.0) {
      tmpP.absfitness = 0.0;
      tmpP.death = 1.0;
      tmpP.birth = 0.0;
    } else{
      // Set appropriate defaults and change only as needed
      tmpP.death = parentP.death;
      tmpP.absfitness = parentP.absfitness;
      tmpP.birth = fitness;
      // exp, mcfarland, and mcfarlandlog as above. Next are the two exceptions.
      // if(typeModel == TypeModel::beerenwinkel) {
      // 	tmpP.absfitness = fitness;
      // 	tmpP.birth = adjust_fitness_B * tmpP.absfitness;
      // } else if(typeModel == TypeModel::mcfarland0) {
      // 	tmpP.absfitness = fitness;
      // 	tmpP.birth = adjust_fitness_MF * tmpP.absfitness;
      // }
    }
  }
	if(F.frequencyDependentFitness){
		popParams.pop_back();
		Genotypes.pop_back();
	}
  // Exp and McFarland and McFarlandlog are also like Datta et al., 2013
  // An additional driver gene mutation increases a cell’s fitness by a
  // factor of (1+sd), whereas an additional housekeeper gene mutation
  // decreases fitness by a factor of (1-sh) and the effect of multiple
  // mutations is multiplicative
}


// is this really any faster than the one below?
inline void new_sp_v(unsigned int& sp,
		     const Genotype& newGenotype,
		     const std::vector<Genotype> Genotypes) {
  sp = 0;
  for(sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] )
      break;
  }
}

inline unsigned int new_sp(const Genotype& newGenotype,
		    const std::vector<Genotype> Genotypes) {
  for(unsigned int sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] ) {
      return sp;
    }
  }
  return Genotypes.size();
}

void remove_zero_sp_nr(std::vector<int>& sp_to_remove,
			      std::vector<Genotype>& Genotypes,
			      std::vector<spParamsP>& popParams,
			      std::multimap<double, int>& mapTimes) {
  std::vector<spParamsP>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype>::iterator Genotypes_begin = Genotypes.begin();
  std::vector<int>::reverse_iterator r = sp_to_remove.rbegin();
  int remove_this;
  while(r != sp_to_remove.rend() ) {
    remove_this = *r;
    mapTimes.erase(popParams[remove_this].pv);
    popParams.erase(popParams_begin + remove_this);
    Genotypes.erase(Genotypes_begin + remove_this);
    ++r;
  }
}


inline void driverCounts(int& maxNumDrivers,
			 int& totalPresentDrivers,
			 std::vector<int>& countByDriver,
			 std::vector<int>& presentDrivers,
			 Rcpp::IntegerMatrix& returnGenotypes,
			 const vector<int>& drv){
  // Fill up the "countByDriver" table, how many genotypes each driver is
  // present.  Return the maximum number of mutated drivers in any
  // genotype, the vector with just the present drivers, and the total
  // number of present drivers.

  // We used to do count_NumDrivers and then whichDrivers
  maxNumDrivers = 0;
  int tmpdr = 0;
  int driver_indx = 0; // the index in the driver table
  for(int j = 0; j < returnGenotypes.ncol(); ++j) {
    tmpdr = 0;
    driver_indx = 0;
    for(int i : drv) {
      tmpdr += returnGenotypes(i - 1, j);
      countByDriver[driver_indx] += returnGenotypes(i - 1, j);
      ++driver_indx;
    }
    if(tmpdr > maxNumDrivers) maxNumDrivers = tmpdr;
  }
  if(returnGenotypes.ncol() > 0) {
    STOPASSERT(driver_indx == static_cast<int>( countByDriver.size()));
  } else {
    STOPASSERT(driver_indx <= static_cast<int>( countByDriver.size()));
  }
  for(size_t i = 0; i < countByDriver.size(); ++i) {
    if(countByDriver[i] > 0) {
      presentDrivers.push_back(i + 1);
      ++totalPresentDrivers;
    }
  }
}



// FIXME: why not keep the number of present drivers in the genotype? We
// call often the getGenotypeDrivers(ge, drv).size()


void nr_totPopSize_and_fill_out_crude_P(int& outNS_i,
					double& totPopSize,
					double& lastStoredSample,
					std::vector<Genotype>& genot_out,
					//std::vector<unsigned long long>& sp_id_out,
					std::vector<double>& popSizes_out,
					std::vector<int>& index_out,
					std::vector<double>& time_out,
					std::vector<double>& sampleTotPopSize,
					std::vector<double>& sampleLargestPopSize,
					std::vector<int>& sampleMaxNDr,
					std::vector<int>& sampleNDrLargestPop,
					bool& simulsDone,
					bool& reachDetection,
					int& lastMaxDr,
					double& done_at,
					const std::vector<Genotype>& Genotypes,
					const std::vector<spParamsP>& popParams,
					const double& currentTime,
					const double& keepEvery,
					const double& detectionSize,
					const double& finalTime,
					// const double& endTimeEvery,
					const int& detectionDrivers,
					const int& verbosity,
					const double& minDetectDrvCloneSz,
					const double& extraTime,
					const vector<int>& drv,
					const double& cPDetect,
					const double& PDBaseline,
					const double& checkSizePEvery,
					double& nextCheckSizeP,
					std::mt19937& ran_gen,
					const bool& AND_DrvProbExit,
					const std::vector<std::vector<int> >& fixation_l,
					const double& fixation_tolerance,
					const int& min_successive_fixation,
					const double& fixation_min_size,
					int& num_successive_fixation,
					POM& pom,
					const std::map<int, std::string>& intName,
					const fitness_as_genes& genesInFitness,
					const double& fatalPopSize = 1e15
					) {
  // Fill out, but also compute totPopSize
  // and return sample summaries for popsize, drivers.
  
  // This determines if we are done or not by checking popSize, number of
  // drivers, etc
  
  // static int lastMaxDr = 0; // preserves value across calls to Algo5 from R.
  // so can not use it.
  bool storeThis = false;
  bool checkSizePNow = false;
  totPopSize = 0.0;
  
  // this could all be part of popSize_over_m_dr, with a better name
  int tmp_ndr = 0;
  int max_ndr = 0;
  double popSizeOverDDr = 0.0;

  for(size_t i = 0; i < popParams.size(); ++i) {
    totPopSize += popParams[i].popSize;
    tmp_ndr = getGenotypeDrivers(Genotypes[i], drv).size();
    if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
    if(tmp_ndr >= detectionDrivers) popSizeOverDDr += popParams[i].popSize;
  }
  lastMaxDr = max_ndr;

  // Until fixation done here. Recall we use an OR operation for exiting
  // below.  Could be added to loop above.
  // And we call allGenesinGenotype also above, inside getGenotypeDrivers.
  // So room for speed ups?
  
  // Since we convert each genotype to a sorted allGenesinGenotype, iterate
  // over that first. Add that pop size if the combination is present in genotype.
  bool fixated = false;
  if(totPopSize > 0) { // Avoid silly things
    if( fixation_l.size() ) {
      std::vector<double> popSize_fixation(fixation_l.size());
      for(size_t i = 0; i < popParams.size(); ++i) {
	std::vector<int> thisg = allGenesinGenotype(Genotypes[i]);
	for(size_t fc = 0; fc != popSize_fixation.size(); ++fc) {
	  // Yes, fixation_l is sorted in R.
	  // if fixation_l[fc] starts with a -9, we are asking
	  // for exact genotype equality
	  if(fixation_l[fc][0] == -9) {
	    // // exact genotype identity?
	    std::vector<int> this_fix(fixation_l[fc].begin() + 1,
				      fixation_l[fc].end());
	    if(thisg == this_fix) {
	      popSize_fixation[fc] = popParams[i].popSize;
	    }
	  } else {
	  // gene combination in fixation element present in genotype?
	    if(std::includes(thisg.begin(), thisg.end(),
			     fixation_l[fc].begin(), fixation_l[fc].end()) ) {
	      popSize_fixation[fc] += popParams[i].popSize;
	    }
	  }
	}
      }
      // Any fixated? But avoid trivial of totPopSize of 0!
      // Now check of > 0 is redundant as we check totPopSize > 0
      // FIXME do we want tolerance around that value?
      double max_popSize_fixation =
	*std::max_element(popSize_fixation.begin(), popSize_fixation.end());
      if( (max_popSize_fixation >= fixation_min_size ) &&
	  (max_popSize_fixation >= (totPopSize * (1 - fixation_tolerance) )) ) {
	++num_successive_fixation;
	// DP1("increased num_successive_fixation");
	if( num_successive_fixation >= min_successive_fixation) fixated = true;
      } else {
	// DP1("zeroed num_successive_fixation");
	num_successive_fixation = 0;
      }
    }
  }

  // // DEBUG
  // if(fixated) {
  //   // print fixation_l
  //   // print popSize_fixation
  //   // print totPopSize
  //   DP1("popSize_fixation");
  //   for(size_t fc = 0; fc != popSize_fixation.size(); ++fc) {
  //     DP2(fc);
  //     DP2(popSize_fixation[fc]);
  //   }
  //   DP2(totPopSize);

  // }
  
  if (keepEvery < 0) {
    storeThis = false;
  } else if( currentTime >= (lastStoredSample + keepEvery) ) {
    storeThis = true;
  }

  if( (totPopSize <= 0.0) || (currentTime >= finalTime)  ) {
    simulsDone = true;
  }

  
  // FIXME
  // this is the usual exit condition
  // (totPopSize >= detectionSize) ||
  // 	  ( (lastMaxDr >= detectionDrivers) &&
  // 	    (popSizeOverDDr >= minDetectDrvCloneSz)

  // Now add the prob. of exiting.

  // Doing this is cheaper than drawing unnecessary runifs.
  // Equality, below, leads to suprises with floating point arith.

  // Operates the same as we do with keepEvery, but here we
  // compute the jump in each accepted sample. And here we use >, not
  // >=
  if(currentTime > nextCheckSizeP) {
    checkSizePNow = true;
    nextCheckSizeP = currentTime + checkSizePEvery;
    // Nope; minimal jump can be smaller than checkSizePEvery
    // nextCheckSizeP += checkSizePEvery;
  } else {
    checkSizePNow = false;
  }

  // We do not verify that conditions for exiting are also satisfied
  // at the exit time when extraTime > 0. We could do that,
  // checking again for the conditions (or the reasonable conditions, so
  // probably not detectSizeP). For instance, with fixated.

  // For fixated in particular, note that we evaluate fixation always, but
  // we might be exiting when there is no longer fixation.  But the logic
  // with fixation is probably to use as large a min_successive_fixation
  // as desired and no extraTime.

  // Probably would not need to check lastMaxDr and popSizeOverDDr
  // as those should never decrease. Really?? FIXME

  
  if(AND_DrvProbExit) {
    // The AND of detectionProb and drivers
    // fixated plays no role here, and cannot be passed from R
    if(extraTime > 0) {
      if(done_at <  0) {
	if( (lastMaxDr >= detectionDrivers) &&
	    (popSizeOverDDr >= minDetectDrvCloneSz) &&
	    checkSizePNow  &&
	    detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen) ) {
	  done_at = currentTime + extraTime;
	}
      } else if (currentTime >= done_at) {
  	simulsDone = true;
  	reachDetection = true;
      }
    } else if( (lastMaxDr >= detectionDrivers) &&
	       (popSizeOverDDr >= minDetectDrvCloneSz) &&
	       checkSizePNow  &&
	       detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen) ) {
      simulsDone = true;
      reachDetection = true; 
    }
  } else {
    // The usual OR mechanism of each option
    if(extraTime > 0) {
      if(done_at <  0) {
	if( (fixated) ||
	    (totPopSize >= detectionSize) ||
	    ( (lastMaxDr >= detectionDrivers) &&
	      (popSizeOverDDr >= minDetectDrvCloneSz) ) ||
	    ( checkSizePNow  &&
	      detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen))) {
	  done_at = currentTime + extraTime;
	}
      } else if (currentTime >= done_at) {
	// if(fixated) {
	  simulsDone = true;
	  reachDetection = true;
	// } else {
	//   done_at = -1;
	// }
      }
    } else if( (fixated) ||
	       (totPopSize >= detectionSize) ||
	       ( (lastMaxDr >= detectionDrivers) &&
		 (popSizeOverDDr >= minDetectDrvCloneSz) ) ||
	       ( checkSizePNow  &&
		 detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen)) ) {
      simulsDone = true;
      reachDetection = true;
    }
  }


  // if( checkSizePNow && (lastMaxDr >= detectionDrivers) &&
  // 	       detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen) )  {
  //   simulsDone = true;
  //   reachDetection = true;
  // }


  if(totPopSize >= fatalPopSize) {
    Rcpp::Rcout << "\n\totPopSize > " << fatalPopSize
		<<". You are likely to loose precision and run into numerical issues\n";
       }

  if(simulsDone)
    storeThis = true;




  // Reuse some info for POM

  if( storeThis ) {
    lastStoredSample = currentTime;
    outNS_i++;
    int ndr_lp = 0;
    double l_pop_s = 0.0;
    int largest_clone = -99;

    time_out.push_back(currentTime);

    for(size_t i = 0; i < popParams.size(); ++i) {
      genot_out.push_back(Genotypes[i]);
      popSizes_out.push_back(popParams[i].popSize);
      index_out.push_back(outNS_i);

      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	ndr_lp = getGenotypeDrivers(Genotypes[i], drv).size();
	largest_clone = i;
      }
    }
    sampleTotPopSize.push_back(totPopSize);
    sampleLargestPopSize.push_back(l_pop_s);
    sampleMaxNDr.push_back(max_ndr);
    sampleNDrLargestPop.push_back(ndr_lp);

    if(l_pop_s > 0) {
      if (largest_clone < 0)
    	throw std::logic_error("largest_clone < 0");
      addToPOM(pom, Genotypes[largest_clone], intName, genesInFitness);
    } else {
      addToPOM(pom, "_EXTINCTION_");
    }
  }

  if( !std::isfinite(totPopSize) ) {
    throw std::range_error("totPopSize not finite");
  }
  if( std::isnan(totPopSize) ) {
    throw std::range_error("totPopSize is NaN");
  }
  // For POM
  if( !storeThis ) {
    double l_pop_s = 0.0;
    int largest_clone = -99;
    for(size_t i = 0; i < popParams.size(); ++i) {
      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	largest_clone = i;
      }
    }
    if(l_pop_s > 0) {
      if (largest_clone < 0)
	throw std::logic_error("largest_clone < 0");
      addToPOM(pom, Genotypes[largest_clone], intName, genesInFitness);
    } else {
      addToPOM(pom, "_EXTINCTION_");
    }
  }
}

// FIXME: I might want to return the actual drivers in each period
// and the actual drivers in the population with largest popsize
// Something like what we do now with whichDrivers
// and count_NumDrivers



inline void nr_reshape_to_outNS(Rcpp::NumericMatrix& outNS,
				const vector<vector<int> >& uniqueGenotV,
				const vector<vector<int> >& genot_out_v,
				const vector<double>& popSizes_out,
				const vector<int>& index_out,
				const vector<double>& time_out){

  vector<vector<int> >::const_iterator fbeg = uniqueGenotV.begin();
  vector<vector<int> >::const_iterator fend = uniqueGenotV.end();

  int column;

  for(size_t i = 0; i < genot_out_v.size(); ++i) {
    column = std::distance(fbeg, lower_bound(fbeg, fend, genot_out_v[i]) );
    outNS(index_out[i], column + 1) =  popSizes_out[i];
  }

  for(size_t j = 0; j < time_out.size(); ++j)
    outNS(j, 0) = time_out[j];
}


Rcpp::NumericMatrix create_outNS(const vector<vector<int> >& uniqueGenotypes,
				 const vector<vector<int> >& genot_out_v,
				 const vector<double>& popSizes_out,
				 const vector<int>& index_out,
				 const vector<double>& time_out,
				 const int outNS_i, const int maxram) {
  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)
  // here("after uniqueGenotypes_to_vector");

  int outNS_r, outNS_c, create_outNS;
  if( ( (uniqueGenotypes.size() + 1) *  (outNS_i + 1) ) > ( pow(2, 31) - 1 ) ) {
    Rcpp::Rcout << "\nWARNING: Return outNS object > 2^31 - 1. Not created.\n";
    outNS_r = 1;
    outNS_c = 1;
    create_outNS = 0;
  } else if (
	     static_cast<long>((uniqueGenotypes.size()+1) * (outNS_i+1)) * 8 >
	     (maxram * (1024*1024) ) ) {
    Rcpp::Rcout << "\nWARNING: Return outNS object > maxram. Not created.\n";
    outNS_r = 1;
    outNS_c = 1;
    create_outNS = 0;
  } else {
    outNS_r = outNS_i + 1;
    outNS_c = uniqueGenotypes.size() + 1;
    create_outNS = 1;
  }
  Rcpp::NumericMatrix outNS(outNS_r, outNS_c);
  if(create_outNS) {
    nr_reshape_to_outNS(outNS, uniqueGenotypes,
			genot_out_v,
			popSizes_out,
			index_out, time_out);

  } else {
    outNS(0, 0) = -99;
  }
  return outNS;
}



// FIXME: when creating the 0/1, collapse those that are the same


vector< vector<int> > uniqueGenot_vector(vector<vector<int> >& genot_out) {
  // From genot_out we want the unique genotypes, but each as a single
  // vector. Convert to the vector, then use a set to give unique sorted
  // vector.
  std::set<std::vector<int> > uniqueGenotypes_nr(genot_out.begin(),
						 genot_out.end());
  std::vector<std::vector<int> > uniqueGenotypes_vector_nr (uniqueGenotypes_nr.begin(),
							    uniqueGenotypes_nr.end());
  return uniqueGenotypes_vector_nr;
}



std::vector<std::vector<int> > genot_to_vectorg(const std::vector<Genotype>& go) {
  std::vector<std::vector<int> > go_l;
  std::transform(go.begin(), go.end(), back_inserter(go_l), genotypeSingleVector);
  return go_l;
}



std::string driversToNameString(const std::vector<int>& presentDrivers,
			    const std::map<int, std::string>& intName) {
  std::string strDrivers;
  std::string comma = "";
  for(auto const &g : presentDrivers) {
    strDrivers += (comma + intName.at(g));
    comma = ", ";
  }
  return strDrivers;
}

// No longer used.
// std::string genotypeToIntString(const std::vector<int>& genotypeV,
// 				   const fitness_as_genes& fg) {

//   // The genotype vectors are returned as a string of ints.

//   std::string strGenotype;

//   std::vector<int> order_int;
//   std::vector<int> rest_int;

//   for(auto const &g : genotypeV) {
//     if( binary_search(fg.orderG.begin(), fg.orderG.end(), g)) {
//       order_int.push_back(g);
//     } else {
//       rest_int.push_back(g);
//     }
//   }

//   std::string order_sep = "_";
//   std::string order_part;
//   std::string rest;
//   std::string comma = "";


//   for(auto const &g : order_int) {
// #ifdef _WIN32
//      order_part += (comma + SSTR(g));
// #endif
// #ifndef _WIN32
//     order_part += (comma + std::to_string(g));
// #endif
//     comma = ", ";
//   }
//   comma = "";
//   for(auto const &g : rest_int) {
// #ifdef _WIN32
//      rest += (comma + SSTR(g));
// #endif
// #ifndef _WIN32
//     rest += (comma + std::to_string(g));
// #endif
//     comma = ", ";
//   }
//   if(fg.orderG.size()) {
//     strGenotype = order_part + order_sep + rest;
//   } else {
//     strGenotype = rest;
//   }
//   return strGenotype;
// }


std::string genotypeToNameString(const std::vector<int>& genotypeV,
				 const fitness_as_genes& fg,
				 const std::map<int, std::string>& intName) {

  // The genotype vectors are returned as a string of names. Similar to
  // the Int version, but we map here to names.

  // As the fitness is stored in terms of modules, not genes, we need to
  // check if a _gene_ is in the order part or not by mapping back to
  // modules. That is the fitness_as_genes argument.

  std::string strGenotype;

  std::vector<int> order_int;
  std::vector<int> rest_int;

  for(auto const &g : genotypeV) {
    if( binary_search(fg.orderG.begin(), fg.orderG.end(), g)) {
      order_int.push_back(g);
    } else {
      rest_int.push_back(g);
    }
  }

  std::string order_sep = " _ ";
  std::string order_part;
  std::string rest;
  std::string comma = "";

  // FIXME: when sure no problems, remove at if needed for speed.
  for(auto const &g : order_int) {
    order_part += (comma + intName.at(g));
    comma = " > "; // comma = ", ";
  }
  comma = "";
  for(auto const &g : rest_int) {
    rest += (comma + intName.at(g));
    comma = ", ";
  }
  if(fg.orderG.size()) {
    strGenotype = order_part + order_sep + rest;
  } else {
    strGenotype = rest;
  }
  return strGenotype;
}


std::vector<std::string> genotypesToNameString(const std::vector< vector<int> >& uniqueGenotypesV,
					       const fitness_as_genes fg,
					       // const fitnessEffectsAll& F,
					       const std::map<int, std::string>& intName) {
  //fitness_as_genes fg = fitnessAsGenes(F); // I use this before;
  std::vector<std::string> gs;
  for(auto const &v: uniqueGenotypesV )
      gs.push_back(genotypeToNameString(v, fg, intName));
  return gs;
}


// std::vector<std::string> genotypesToString(const std::vector< vector<int> >& uniqueGenotypesV,
// 					   const fitnessEffectsAll& F,
// 					   bool names = true) {
//   fitness_as_genes fg = fitnessAsGenes(F);
//   std::vector<std::string> gs;

//   if(names) {
//     std::map<int, std::string> intName = mapGenesIntToNames(F);
//     for(auto const &v: uniqueGenotypesV )
//       gs.push_back(genotypeToNameString(v, fg, intName));
//   } else {
//       for(auto const &v: uniqueGenotypesV )
// 	gs.push_back(genotypeToIntString(v, fg));
//   }

//   // exercise: do it with lambdas
//   // std::transform(uniqueGenotypesV.begin(), uniqueGenotypesV.end(),
//   // 		 back_inserter(gs), vectorGenotypeToString);
//   return gs;
// }

Rcpp::IntegerMatrix nr_create_returnGenotypes(const int& numGenes,
					      const std::vector< vector<int> >& uniqueGenotypesV){
  // We loose order here. Thus, there might be several identical columns.
  Rcpp::IntegerMatrix returnGenotypes(numGenes, uniqueGenotypesV.size());
  for(size_t i = 0; i < uniqueGenotypesV.size(); ++i) {
    for(int j : uniqueGenotypesV[i]) {
      returnGenotypes(j - 1, i) = 1;
    }
  }
  return returnGenotypes;
}






static void nr_sample_all_pop_P(std::vector<int>& sp_to_remove,
				std::vector<spParamsP>& popParams,
				// next only used with DEBUGV
				const std::vector<Genotype>& Genotypes,
				const double& tSample,
				const int& mutationPropGrowth){
  sp_to_remove.clear();

  for(size_t i = 0; i < popParams.size(); i++) {
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n     ********* 5.9 ******\n "
	      << "     Species  = " << i
		<< "\n      Genotype = ";
    print_Genotype(Genotypes[i]); //genotypeSingleVector(Genotypes[i])
    //	      << "\n      sp_id = " << genotypeSingleVector(Genotypes[i]) // sp_id[i]
    Rcpp::Rcout << "\n      pre-update popSize = "
	      << popParams[i].popSize
	      << "\n      time of sample = " << tSample
	      << "\n      popParams[i].timeLastUpdate = "
	      << popParams[i].timeLastUpdate
	      << ";\n     t for Algo2 = "
	      << tSample - popParams[i].timeLastUpdate
	      << " \n     species R " << popParams[i].R
	      << " \n     species W " << popParams[i].W
	      << " \n     species death " << popParams[i].death
	      << " \n     species birth " << popParams[i].birth;
#endif

    // Account for forceSampling. When
    // forceSampling, popSize for at least one species
    // was updated in previous loop, so we skip that one
    if(tSample > popParams[i].timeLastUpdate) {
      popParams[i].popSize =
	Algo2_st(popParams[i], tSample, mutationPropGrowth);
    }
    if( popParams[i].popSize <=  0.0 ) {
      // this i has never been non-zero in any sampling time
      // eh??
      // If it is 0 here, remove from _current_ population. Anything that
      // has had a non-zero size at sampling time is preserved (if it
      // needs to be preserved, because it is keepEvery time).
      sp_to_remove.push_back(i);

#ifdef DEBUGV
      Rcpp::Rcout << "\n\n     Removing species i = " << i
		  << " with genotype = ";
      print_Genotype(Genotypes[i]); //genotypeSingleVector(Genotypes[i]);
#endif
    }
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n   post-update popSize = "
	      << popParams[i].popSize << "\n";
#endif
  }
}

// zz: add population size of parent, to get the true LOD
// as in Szendro
void addToPhylog(PhylogName& phylog,
		 const Genotype& parent,
		 const Genotype& child,
		 const double time,
		 const std::map<int, std::string>& intName,
		 const fitness_as_genes& fg,
		 const double pop_size_child) {
  phylog.time.push_back(time);
  phylog.parent.push_back(genotypeToNameString(genotypeSingleVector(parent),
					       fg, intName));
  phylog.child.push_back(genotypeToNameString(genotypeSingleVector(child),
					      fg, intName));
  phylog.pop_size_child.push_back(pop_size_child);
}

// // Only called when the child has pop size of 0
// // so true LOD
// void addToLOD(LOD& lod,
// 	      const Genotype& parent,
// 	      const Genotype& child,
// 	      // const double time,
// 	      const std::map<int, std::string>& intName,
// 	      const fitness_as_genes& fg) {
//   // lod.time.push_back(time);
//   lod.parent.push_back(genotypeToNameString(genotypeSingleVector(parent),
// 					       fg, intName));
//   lod.child.push_back(genotypeToNameString(genotypeSingleVector(child),
// 					      fg, intName));
// }


// Only called when the child has pop size of 0
// so true LOD
// Use a map for LOD, and overwrite the parent:
// we only add when the size of the child is 0
// The key of the map is the child.

// FIXME: we might want to store the time? Not really clear even if that
// makes sense. We would be storing the last time the child (which had 0
// size at that time) arose from the parent.
// A simple kludge is to have two maps, the second with child and time.
// Or do it properly as map<int, genot_time_struct>
// genot_time_struct {string parent; double time}

void addToLOD(std::map<std::string, std::string>& lod,
	      const Genotype& parent,
	      const Genotype& child,
	      // const double time,
	      const std::map<int, std::string>& intName,
	      const fitness_as_genes& fg) {
  std::string parent_str = genotypeToNameString(genotypeSingleVector(parent),
					 fg, intName);
  std::string child_str = genotypeToNameString(genotypeSingleVector(child),
					       fg, intName);
  lod[child_str] = parent_str;
  // // lod.time.push_back(time);
  // lod.parent.push_back(genotypeToNameString(genotypeSingleVector(parent),
  // 					       fg, intName));
  // lod.child.push_back(genotypeToNameString(genotypeSingleVector(child),
  // 					      fg, intName));
}


void addToPOM(POM& pom,
	      const Genotype& genotype,
	      const std::map<int, std::string>& intName,
	      const fitness_as_genes& fg) {

  if (pom.genotypes.empty()) {
    std::string g = genotypeToNameString(genotypeSingleVector(genotype),
				       fg, intName);
    pom.genotypesString.push_back(g);
    pom.genotypes.push_back(genotype);
  } else if ( !(pom.genotypes.back() == genotype) ) {
    // Insert only if different from previous
    std::string g = genotypeToNameString(genotypeSingleVector(genotype),
				       fg, intName);
    pom.genotypesString.push_back(g);
    pom.genotypes.push_back(genotype);
  }
}

// to explicitly signal extinction
void addToPOM(POM& pom,
	      const std::string string) {
  pom.genotypesString.push_back(string);
}



static void nr_innerBNB (const fitnessEffectsAll& fitnessEffects,
			const double& initSize,
			const double& K,
			// const double& alpha,
			// const double& genTime,
			const TypeModel typeModel,
			const int& mutationPropGrowth,
			const std::vector<double>& mu,
			// const double& mu,
			const double& death,
			const double& keepEvery,
			const double& sampleEvery,
			const std::vector<int>& initMutant,
			const time_t& start_time,
			const double& maxWallTime,
			const double& finalTime,
			const double& detectionSize,
			const int& detectionDrivers,
			const double& minDetectDrvCloneSz,
			const double& extraTime,
			const int& verbosity,
			double& totPopSize,
			double& em1,
			double& em1sc,
			// double& n_1,
			// double& en1,
			double& ratioForce,
			double& currentTime,
			int& speciesFS,
			int& outNS_i,
			int& iter,
			std::vector<Genotype>& genot_out,
			std::vector<double>& popSizes_out,
			std::vector<int>& index_out,
			std::vector<double>& time_out,
			std::vector<double>& sampleTotPopSize,
			std::vector<double>& sampleLargestPopSize,
			std::vector<int>& sampleMaxNDr,
			std::vector<int>& sampleNDrLargestPop,
			bool& reachDetection,
			std::mt19937& ran_gen,
			// randutils::mt19937_rng& ran_gen,
			double& runningWallTime,
			bool& hittedWallTime,
			const std::map<int, std::string>& intName,
			const fitness_as_genes& genesInFitness,
			PhylogName& phylog,
			bool keepPhylog,
			const fitnessEffectsAll& muEF,
			const std::vector<int>& full2mutator,
			const double& cPDetect,
			const double& PDBaseline,
			const double& checkSizePEvery,
			const bool& AND_DrvProbExit,
			const std::vector< std::vector<int> >& fixation_l,
			 const double& fixation_tolerance,
			 const int& min_successive_fixation,
			 const double& fixation_min_size,
			 int& ti_dbl_min,
			 int& ti_e3,
			 std::map<std::string, std::string>& lod,
			 // LOD& lod,
			 POM& pom) {

  double nextCheckSizeP = checkSizePEvery;
  const int numGenes = fitnessEffects.genomeSize;

  double mymindummy = 1.0e-11; //1e-10
  double targetmindummy = 1.0e-10; //1e-9
  double minmu = *std::min_element(mu.begin(), mu.end());
  // Very small, but no less than mymindummy, for numerical issues.
  // We can probably go down to 1e-13. 1e-16 is not good as we get lots
  // of pE.f not finite. 1e-15 is probably too close, and even if no pE.f
  // we can get strange behaviors.
  double dummyMutationRate = std::max(std::min(minmu/1.0e4, targetmindummy),
				      mymindummy);
  // This should very rarely happen:
  if(minmu <= mymindummy) { // 1e-9
    double newdd = minmu/100.0;
    Rcpp::Rcout << "WARNING: the smallest mutation rate is "
		<< "<= " << mymindummy << ". That is a really small value"
		<< "(per-base mutation rate in the human genome is"
		<< " ~ 1e-11 to 1e-9). "
		<< "Setting dummyMutationRate to your min/100 = "
		<< newdd
		<< ". There can be numerical problems later.\n";
    dummyMutationRate = newdd;
  }
  // double dummyMutationRate = 1e-10;
  // ALWAYS initialize this here, or reinit or rezero
  genot_out.clear();

  phylog = PhylogName();
  // lod = LOD();
  lod.clear();
  pom = POM();

  popSizes_out.clear();
  index_out.clear();
  time_out.clear();
  totPopSize = 0.0;
  sampleTotPopSize.clear();
  currentTime = 0.0;
  iter = 0;

  outNS_i = -1;

  sampleTotPopSize.clear();
  sampleLargestPopSize.clear();
  sampleMaxNDr.clear();
  sampleNDrLargestPop.clear();
  // end of rezeroing.


  // }
  // anyForceRerunIssues = false;

  bool forceSample = false;
  bool simulsDone = false;
  double lastStoredSample = 0.0;


  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double timeNextPopSample;
  double tSample;

  std::vector<int> newMutations;
  int nextMutant;
  unsigned int numSpecies = 0;
  int numMutablePosParent = 0;


  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 100;
  //int timeL = 1000;

  int iterInterrupt = 50000; //how large should we make this?

  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;

  std::vector<int>sp_to_remove(1);
  sp_to_remove.reserve(10000);

  // those to update
  int to_update = 1; //1: one species; 2: 2 species; 3: all.
  int u_1 = -99;
  int u_2 = -99;

  Genotype newGenotype;
  std::vector<Genotype> Genotypes(1);
  Genotypes[0] = wtGenotype(); //Not needed, but be explicit.

  std::vector<spParamsP> popParams(1);

      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 1 ";
      // print_spP(popParams[0]);
      // // end debug


  const int sp_per_period = 5000;
  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 01 ";
      // print_spP(popParams[0]);
      // // end debug


  spParamsP tmpParam;
  init_tmpP(tmpParam);
  init_tmpP(popParams[0]);
  popParams[0].popSize = initSize;
  totPopSize = initSize;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10000 ";
      // print_spP(popParams[0]);
      // // end debug



  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

    // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;
  //std::multimap<double, int>::iterator m1pos;

  // int ti_dbl_min = 0;
  // int ti_e3 = 0;
  ti_dbl_min = 0;
  ti_e3 = 0;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10002 ";
      // print_spP(popParams[0]);
      // // end debug



  // // Beerenwinkel
  // double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  //McFarland
  double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  // for McFarland error
  em1 = 0.0;
  em1sc = 0.0;
  // n_0 = 0.0;
  // n_1 = 0.0;
  // double tps_0; //, tps_1;
  // tps_0 = totPopSize;
  // tps_1 = totPopSize;

  // en1 = 0;
  double totPopSize_previous = totPopSize;
  double DA_previous = log1p(totPopSize_previous/K);

      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10004 ";
      // print_spP(popParams[0]);
      // // end debug



  int lastMaxDr = 0;
  double done_at = -9;


  int num_successive_fixation = 0; // none so far


#ifdef MIN_RATIO_MUTS_NR
  g_min_birth_mut_ratio_nr = DBL_MAX;
  g_min_death_mut_ratio_nr = DBL_MAX;
  g_tmp1_nr = DBL_MAX;
#endif

      // // FIXME debug
      // Rcpp::Rcout << " popSize[0]  at 1b ";
      // print_spP(popParams[0]);
      // // end debug

    // This long block, from here to X1, is ugly and a mess!
  // This is what takes longer to figure out whenever I change
  // anything. FIXME!!
  if(initMutant.size() > 0) {
    Genotypes[0] = createNewGenotype(wtGenotype(),
				     initMutant,
				     fitnessEffects,
				     ran_gen,
				     false);
    int numGenesInitMut = Genotypes[0].orderEff.size() +
      Genotypes[0].epistRtEff.size() + Genotypes[0].rest.size();
    int numGenesGenotype = fitnessEffects.allGenes.size();
    popParams[0].numMutablePos = numGenesGenotype - numGenesInitMut;
    // Next two are unreachable since caught in R.
    // But just in case, since it would lead to seg fault.
    if(popParams[0].numMutablePos < 0)
      throw std::invalid_argument("initMutant's genotype has more genes than are possible.");
    if(popParams[0].numMutablePos == 0)
      throw std::invalid_argument("initMutant has no mutable positions: genotype with all genes mutated.");
    // popParams[0].numMutablePos = numGenes - 1;
    // From obtainMutations, but initMutant an int vector. But cumbersome.
    // std::vector<int> sortedg = convertGenotypeFromInts(initMutant);
    // sort(sortedg.begin(), sortedg.end());
    // std::vector<int> nonmutated;
    // set_difference(fitnessEffects.allGenes.begin(), fitnessEffects.allGenes.end(),
    // 		   sortedg.begin(), sortedg.end(),
    // 		   back_inserter(nonmutated));
    // popParams[0].numMutablePos = nonmutated.size();



    // Commenting out the unused models!
    // if(typeModel == TypeModel::beerenwinkel) {

    //   popParams[0].death = 1.0; //note same is in McFarland.
    //   // But makes sense here; adjustment in beerenwinkel is via fitness

    //   // initialize to prevent birth/mutation warning with Beerenwinkel
    //   // when no mutator. O.w., the defaults
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = prodFitness(evalGenotypeFitness(Genotypes[0],
    // 								fitnessEffects));
    //   updateRatesBeeren(popParams, adjust_fitness_B, initSize,
    // 			currentTime, alpha, initSize,
    // 			mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland0) {
    //   // death equal to birth of a non-mutant.
    //   popParams[0].death = log1p(totPopSize/K); // log(2.0), except rare cases
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = prodFitness(evalGenotypeFitness(Genotypes[0],
    // 								fitnessEffects));
    //   updateRatesMcFarland0(popParams, adjust_fitness_MF, K,
    // 			    totPopSize,
    // 			    mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland) {
    //   popParams[0].death = totPopSize/K;
    //   popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
    // 								fitnessEffects));
    // } else       if(typeModel == TypeModel::mcfarlandlog) {

    if(typeModel == TypeModel::mcfarlandlog) {
      popParams[0].death = log1p(totPopSize/K);
      popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects, Genotypes, popParams));
    } else if(typeModel == TypeModel::mcfarlandlog_d) {
      popParams[0].death = std::max(1.0, log1p(totPopSize/K));
      popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects, Genotypes, popParams));
    } else if(typeModel == TypeModel::bozic1) {
      tmpParam.birth =  1.0;
      tmpParam.death = -99.9;
    // } else if (typeModel == TypeModel::bozic2) {
    //   tmpParam.birth =  -99;
    //   tmpParam.death = -99;
    } else if (typeModel == TypeModel::exp) {
      tmpParam.birth =  -99;
      tmpParam.death = death;
    } else {
      // caught in R, so unreachable here
      throw std::invalid_argument("this ain't a valid typeModel");
    }
    // if( (typeModel != TypeModel::beerenwinkel) && (typeModel != TypeModel::mcfarland0)
    // 	&& (typeModel != TypeModel::mcfarland) && (typeModel != TypeModel::mcfarlandlog)) // wouldn't matter
    //   nr_fitness(popParams[0], tmpParam,
    // 		 Genotypes[0],
    // 		 fitnessEffects,
    // 		 typeModel, genTime,
    // 		 adjust_fitness_B, adjust_fitness_MF);
    if( (typeModel != TypeModel::mcfarlandlog) &&
	(typeModel != TypeModel::mcfarlandlog_d) ) // wouldn't matter
      nr_fitness(popParams[0], tmpParam,
		 Genotypes[0],
		 fitnessEffects,
		 typeModel, Genotypes, popParams);
    // , genTime);
    //		 adjust_fitness_B, adjust_fitness_MF);
    // we pass as the parent the tmpParam; it better initialize
    // everything right, or that will blow. Reset to init
    init_tmpP(tmpParam);
    addToPOM(pom, Genotypes[0], intName, genesInFitness);
	} else { //no initMutant
    popParams[0].numMutablePos = numGenes;

    // if(typeModel == TypeModel::beerenwinkel) {
    //   popParams[0].death = 1.0;
    //   // initialize to prevent birth/mutation warning with Beerenwinkel
    //   // when no mutator. O.w., the defaults
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = 1.0;
    //   updateRatesBeeren(popParams, adjust_fitness_B, initSize,
    // 			currentTime, alpha, initSize,
    // 			mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland0) {
    //   popParams[0].death = log1p(totPopSize/K);
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = 1.0;
    //   updateRatesMcFarland0(popParams, adjust_fitness_MF, K,
    // 			    totPopSize,
    // 			    mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland) {
    //   popParams[0].birth = 1.0;
    //   popParams[0].death = totPopSize/K;
    //   // no need to call updateRates
    // } else if(typeModel == TypeModel::mcfarlandlog) {
    if(typeModel == TypeModel::mcfarlandlog) {
      if(fitnessEffects.frequencyDependentFitness){
	popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
							     fitnessEffects, Genotypes, popParams));
      } else {
	popParams[0].birth = 1.0;
      }
      popParams[0].death = log1p(totPopSize/K);
      // no need to call updateRates
    } else if(typeModel == TypeModel::mcfarlandlog_d) {
      if(fitnessEffects.frequencyDependentFitness){
	popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
							     fitnessEffects, Genotypes, popParams));
      }else{
	popParams[0].birth = 1.0;
      }
      popParams[0].death = std::max(1.0, log1p(totPopSize/K));
      
    } else if(typeModel == TypeModel::bozic1) {

			if(fitnessEffects.frequencyDependentFitness){
 				popParams[0].birth = prodDeathFitness(evalGenotypeFitness(Genotypes[0],
 					fitnessEffects, Genotypes, popParams));
 			}else{
				popParams[0].birth = 1.0;

			}
      popParams[0].death = 1.0;
    // } else if (typeModel == TypeModel::bozic2) {
    //   popParams[0].birth = 0.5/genTime;
    //   popParams[0].death = 0.5/genTime;
    } else if (typeModel == TypeModel::exp) {

			if(fitnessEffects.frequencyDependentFitness){
				popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
					fitnessEffects, Genotypes, popParams));
			}else{
				popParams[0].birth = 1.0;
			}
      popParams[0].death = death;
    } else {
      throw std::invalid_argument("this ain't a valid typeModel");
    }
  }

  // // these lines (up to, and including, R_F_st)
  // // not needed with mcfarland0 or beerenwinkel
  // if(mutationPropGrowth)
  //   popParams[0].mutation = mu * popParams[0].birth * popParams[0].numMutablePos;
  // else
  //   popParams[0].mutation = mu * popParams[0].numMutablePos;

  popParams[0].mutation = mutationFromScratch(mu, popParams[0], Genotypes[0],
					      fitnessEffects, mutationPropGrowth,
					      full2mutator, muEF,
								Genotypes, popParams);
  W_f_st(popParams[0]);
  R_f_st(popParams[0]);


  // X1: end of mess of initialization block

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  if( keepEvery > 0 ) {
    // We keep the first ONLY if we are storing more than one.
    outNS_i++;
    time_out.push_back(currentTime);

    genot_out.push_back(Genotypes[0]);
    popSizes_out.push_back(popParams[0].popSize);
    index_out.push_back(outNS_i);

    sampleTotPopSize.push_back(popParams[0].popSize);
    sampleLargestPopSize.push_back(popParams[0].popSize);
    sampleMaxNDr.push_back(getGenotypeDrivers(Genotypes[0],
					      fitnessEffects.drv).size());
    sampleNDrLargestPop.push_back(sampleMaxNDr[0]);
  }
  // FIXME: why next line and not just genot_out.push_back(Genotypes[i]);
  // if keepEvery > 0? We do that already.
  // It is just ugly to get a 0 in that first genotype when keepEvery < 0
  // uniqueGenotypes.insert(Genotypes[0].to_ullong());
  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;


#ifdef DEBUGV
  Rcpp::Rcout << "\n the initial species\n";
  print_spP(popParams[0]);
#endif





  while(!simulsDone) {
    // Check how we are doing with time as first thing.
    runningWallTime = difftime(time(NULL), start_time);
    if( runningWallTime > maxWallTime ) {
      hittedWallTime = true;
      forceSample = true;
      simulsDone = true;
    }

    iter++;

    if( !(iter % iterInterrupt))
      Rcpp::checkUserInterrupt();

    if(verbosity > 1) {
      if(! (iter % iterL) ) {
	Rcpp::Rcout << "\n\n    ... iteration " << iter;
	Rcpp::Rcout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	Rcpp::Rcout << "\n\n    ... iteration " << iter;
	Rcpp::Rcout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      Rcpp::Rcout <<"\n\n\n*** Looping through 5.2. Iter = " << iter
		  << ".  Current time " << currentTime <<	" \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    Rcpp::Rcout << " DEBUGV\n";
    Rcpp::Rcout << "\n ForceSample? " << forceSample
	      << "  tSample " << tSample
	      << "  currentTime " << currentTime;
#endif

    if(iter == 1) { // handle special case of first iter
    tmpdouble1 = ti_nextTime_tmax_2_st(popParams[0],
					 currentTime,
					 tSample,
					 ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, 0, tmpdouble1);
      //popParams[0].Flag = false;
      popParams[0].timeLastUpdate = currentTime;
    } else { // any other iter
      if(to_update == 1) {
	// we did not sample or mutate to a different species in previous period
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample,
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	popParams[u_1].timeLastUpdate = currentTime;

#ifdef DEBUGV
	detect_ti_duplicates(mapTimes, tmpdouble1, u_1);
#endif

#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, update one ******\n For to_update = \n "
		  << "     tSample  = " << tSample

		  << "\n\n**   Species  = " << u_1
		    << "\n       genotype =  ";
	print_Genotype(Genotypes[u_1]);
	Rcpp::Rcout << "\n       popSize = " << popParams[u_1].popSize
		  << "\n       currentTime = " << currentTime
		  << "\n       popParams[i].nextMutationTime = "
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth;
#endif

      } else if(to_update == 2) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	tmpdouble2 = ti_nextTime_tmax_2_st(popParams[u_2],
					   currentTime,
					   tSample, ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_2, tmpdouble2);
	popParams[u_1].timeLastUpdate = currentTime;
	popParams[u_2].timeLastUpdate = currentTime;

#ifdef DEBUGV
	detect_ti_duplicates(mapTimes, tmpdouble1, u_1);
	detect_ti_duplicates(mapTimes, tmpdouble2, u_2);
#endif



#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, update two ******\n "
		  << "     tSample  = " << tSample

		  << "\n\n**   Species  = " << u_1
		    << "\n       genotype =  ";
	print_Genotype(Genotypes[u_1]);
	Rcpp::Rcout << "\n       popSize = " << popParams[u_1].popSize
		  << "\n       currentTime = " << currentTime
		  << "\n       popParams[i].nextMutationTime = "
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth


		  << "\n\n**     Species  = " << u_2
		    << "\n       genotype =  ";
	print_Genotype(Genotypes[u_2]);
	Rcpp::Rcout << "\n       popSize = " << popParams[u_2].popSize
		  << "\n       currentTime = " << currentTime
		  << "\n       popParams[i].nextMutationTime = "
		  << tmpdouble2
		  << " \n     species R " << popParams[u_2].R
		  << " \n     species W " << popParams[u_2].W
		  << " \n     species death " << popParams[u_2].death
		  << " \n     species birth " << popParams[u_2].birth;
#endif

      } else { // we sampled, so update all: i.e. to_update == 3
	for(size_t i = 0; i < popParams.size(); i++) {
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample, ti_dbl_min, ti_e3);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].timeLastUpdate = currentTime;
#ifdef DEBUGV
	  detect_ti_duplicates(mapTimes, tmpdouble1, i);
#endif

#ifdef DEBUGV
	  Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, update all ******\n "
		    << "     Species  = " << i
		      << "\n       genotype =  ";
	  print_Genotype(Genotypes[i]);
	  Rcpp::Rcout << "\n       popSize = " << popParams[i].popSize
		    << "\n       currentTime = " << currentTime
		      << "\n       popParams[i].nextMutationTime = "
		      << tmpdouble1
		    << " \n     species R " << popParams[i].R
		    << " \n     species W " << popParams[i].W
		    << " \n     species death " << popParams[i].death
		    << " \n     species birth " << popParams[i].birth;

#endif
	}
      }
    }
    if(forceSample) {
      // A VERY ugly hack. Resetting tSample to jump to sampling.
      tSample = currentTime;
      // Need this, o.w. would skip a sampling.
      timeNextPopSample = currentTime;
    }


    // ******************** 5.3 and do we sample? ***********
    // Find minimum to know if we need to sample the whole pop
    // We also obtain the nextMutant
    getMinNextMutationTime4(nextMutant, minNextMutationTime,
			    mapTimes);

    if(verbosity >= 2) {
      Rcpp::Rcout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime
		<< "; timeNextPopSample = " << timeNextPopSample
		<< "; popParams.size() = " << popParams.size() << "\n";
    }

    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;
      // ************   5.4   ***************
      mutantTimeSinceLastUpdate = currentTime -
	popParams[nextMutant].timeLastUpdate;

      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
					       mutantTimeSinceLastUpdate);

      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
#ifdef DEBUGV
	//if(verbosity > -2) {
	// We always warn about this, since interaction with ti==0
	Rcpp::Rcout << "\n Forced sampling triggered for next loop: \n    " <<
	  " popParams[nextMutant].popSize = " <<
	  popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	Rcpp::Rcout << " when nextMutant = " << nextMutant <<
	  " at iteration " << iter << "\n";
	//}
#endif
      }
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet
      // know total number of different species
      // This is a protection against things going wild. Should
      // not happen in regular usage.
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2;
#ifdef DEBUGV
      	//if(verbosity > -2) // we always warn about this

	Rcpp::Rcout << "\n Forced sampling triggered for next loop "
		  << " when numSpecies = " <<
	  numSpecies << " at iteration " << iter << "\n";
#endif
      }

      if(popParams[nextMutant].numMutablePos != 0) {
	// this is the usual case. The alternative is the dummy or null mutation


	// ************   5.5   ***************

	newMutations.clear();
	// FIXME: nonmutated also returned here
	obtainMutations(Genotypes[nextMutant],
			fitnessEffects,
			numMutablePosParent,
			newMutations,
			ran_gen,
			mu);
	//DP2(newMutations);
	// nr_change
	// getMutatedPos_bitset(mutatedPos, numMutablePosParent, // r,
	// 		     ran_gen,
	// 		     mutablePos,
	// 		     Genotypes[nextMutant],
	// 		     numGenes);

	// ************   5.6   ***************
	newGenotype = createNewGenotype(Genotypes[nextMutant],
					newMutations,
					fitnessEffects,
					ran_gen,
					true);
	// nr_change
	// newGenotype = Genotypes[nextMutant];
	// newGenotype.set(mutatedPos);
	// newGenotype[mutatedPos] = 1;

	// FIXME
	// any speed diff between a) and b)?
	// a)
	new_sp_v(sp, newGenotype, Genotypes);
	// b)
	// sp = 0;
	// sp = new_sp(newGenotype, Genotypes);

	// nr_change
	// new_sp_bitset(sp, newGenotype, Genotypes);

	if(sp == numSpecies) {// New species
	  ++numSpecies;
	  init_tmpP(tmpParam);

	  if(verbosity >= 2) {
	    Rcpp::Rcout <<"\n     Creating new species   " << (numSpecies - 1)
			<< "         from species "  <<   nextMutant;
	  }

#ifdef DEBUGW
	  if( (currentTime - popParams[nextMutant].timeLastUpdate) < 0.0) {
	    DP2(currentTime); //this is set to minNextMutationTime above
	    DP2(minNextMutationTime);
	    DP2(tSample);
	    DP2(popParams[nextMutant].timeLastUpdate);
	    DP2( (currentTime -  popParams[nextMutant].timeLastUpdate) );
	    DP2( (currentTime <  popParams[nextMutant].timeLastUpdate) );
	    DP2( (currentTime ==  popParams[nextMutant].timeLastUpdate) );
	    DP2(nextMutant);
	    DP2(u_1);
	    DP2(u_2);
	    DP2(tmpdouble1);
	    DP2(tmpdouble2);
	    DP2(popParams[nextMutant].timeLastUpdate);
	    DP2(popParams[u_1].timeLastUpdate);
	    DP2(popParams[u_2].timeLastUpdate);
	    DP2( (popParams[u_1].timeLastUpdate - popParams[u_2].timeLastUpdate) );
	    DP2( (popParams[u_1].timeLastUpdate - popParams[nextMutant].timeLastUpdate) );
	    DP2( (popParams[u_1].timeLastUpdate - popParams[0].timeLastUpdate) );
	    print_spP(popParams[nextMutant]);
	    throw std::out_of_range("new species: currentTime - timeLastUpdate[sp] out of range. ***###!!!Serious bug!!!###***");
	  }
#endif
	  tmpParam.popSize = 1;

	  nr_fitness(tmpParam, popParams[nextMutant],
		     newGenotype,
		     fitnessEffects,
		     typeModel, Genotypes, popParams);// , genTime,
		     // adjust_fitness_B, adjust_fitness_MF);

	  if(tmpParam.birth > 0.0) {
	    // if(keepMutationTimes)
	    //   update_mutation_freqs(newMutation, currentTime, mutation_freq_at);
	    //FIXME: phylog
	    // if(keepPhylog)
	    //   addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
	    // 		  intName, genesInFitness);

	    tmpParam.numMutablePos = numMutablePosParent - 1;
	    tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
					       fitnessEffects,
					       mutationPropGrowth, full2mutator,
						    muEF, Genotypes, popParams);
	    // tmpParam.mutation = mutationFromParent(mu, tmpParam, popParams[nextMutant],
	    // 					   newMutations, mutationPropGrowth,
	    // 					   newGenotype, full2mutator,
	    // 					   muEF);


	    //tmpParam.mutation = mu * (numMutablePosParent - 1);
	    if (tmpParam.mutation > 1 )
	      Rcpp::Rcout << "WARNING: mutation > 1\n";
	    if (numMutablePosParent == 1) {
	      if(verbosity >= 1)
		Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
	      // FIXME:varmutrate: give the value of dummy here.
	      tmpParam.mutation = dummyMutationRate; // dummy mutation here. Set some mu.
	    }
	    W_f_st(tmpParam);
	    R_f_st(tmpParam);
	    tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP does what it should.
	    // as this is a new species
	    popParams.push_back(tmpParam);
	    Genotypes.push_back(newGenotype);
	    to_update = 2;
#ifdef MIN_RATIO_MUTS_NR
	    g_tmp1_nr = tmpParam.birth/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;

	    g_tmp1_nr = tmpParam.death/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;
#endif

	    // LOD:
	    // here first call to addToPhylog, with popSize popParams[sp].popSize
	    // and it is 0
	    if(keepPhylog)
	      addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
			  intName, genesInFitness, 0);
	    // LOD, as LOD sensu stricto, always done now
	    addToLOD(lod, Genotypes[nextMutant], newGenotype, // currentTime,
			intName, genesInFitness);

	  } else {// fitness is 0, so we do not add it
	    --sp;
	    --numSpecies;
	    to_update = 1;
	  }
	  // #ifdef DEBUGV
	  if(verbosity >= 3) {
	    Rcpp::Rcout << " \n\n\n Looking at NEW species " << sp << " at creation";
	    Rcpp::Rcout << "\n New Genotype :";
	    print_Genotype(newGenotype);
	    Rcpp::Rcout << "\n Parent Genotype :";
	    print_Genotype(Genotypes[nextMutant]);
	    // Rcpp::Rcout << "\n Genotype = " << genotypeSingleVector(newGenotype); //Genotypes[sp];
	    //Genotypes[sp].to_ullong();
	    Rcpp::Rcout << "\n birth of sp = " << tmpParam.birth;
	    Rcpp::Rcout << "\n death of sp = " << tmpParam.death;
	    // Rcpp::Rcout << "\n s = " << s;
	    Rcpp::Rcout << "\n parent birth = " << popParams[nextMutant].birth;
	    Rcpp::Rcout << "\n parent death = " << popParams[nextMutant].death;
	    // Rcpp::Rcout << "\n parent Genotype = " << genotypeSingleVector(Genotypes[nextMutant]);
	    Rcpp::Rcout << "\n\n popParams parent: \n";
	    print_spP(popParams[nextMutant]);
	    Rcpp::Rcout << "\n\npopParams child: \n";
	    print_spP(tmpParam);
	    }
	  // #endif
	} else {	// A mutation to pre-existing species

	  // What we do here is step 6 of Algorithm 5, in the "Otherwise",
	  // in p. 5 of suppl mat. We will update both, and only these
	  // two.
	  to_update = 2;

#ifdef DEBUGW
	  if( (currentTime - popParams[sp].timeLastUpdate) < 0.0) {
	    // Yes, the difference could be 0 if two next mutation times are identical.
	    // You enable detect_ti_duplicates and use trigger-duplicated-ti.R
	    // to see it.
	    // Often the involved culprits (nextMutant and the other, say sp)
	    // were lastUpdated with tiny difference and they were, when updated
	    // given an identical ti, each in its own run.
	    // Key is not timeLastUpdate. This is a possible sequence of events:
	    //    - at time t0, species that will become nextMutant is updated and gets ti = tinm
	    //    - t1: species u1 gets ti = tinm
	    //    - t2: species u2 gets some ti > tinm
	    //    - tinm becomes minimal, so we mutate u1, and it mutates to u2
	    //    - (so now the timeLastUpdate of u1 = u2 = tinm)
	    //    - nextMutant is now mutated, and it mutates to u2, which becomes sp
	    //    - tinm = timeLastUpdate of u1 and u2.
	    //    - You will also see that number of mutations, or genotypes are such
	    //      that, in this case, u2 is the most mutated, etc.
	    //    - If you enable the detect_ti_duplicates, you would have seen duplicated ti
	    //      for nextMutant and u1

	    //   Even simpler is if above, nextMutant will mutate to u1 (not u2) so u1 becomes sp.
	    DP2(currentTime); //this is set to minNextMutationTime above
	    DP2(minNextMutationTime);
	    DP2(tSample);
	    DP2(popParams[sp].timeLastUpdate);
	    DP2( (currentTime -  popParams[sp].timeLastUpdate) );
	    DP2( (currentTime <  popParams[sp].timeLastUpdate) );
	    DP2( (currentTime ==  popParams[sp].timeLastUpdate) );
	    DP2(sp);
	    DP2(nextMutant);
	    DP2(u_1);
	    DP2(u_2);
	    DP2(tmpdouble1);
	    DP2(tmpdouble2);
	    DP2(popParams[sp].timeLastUpdate);
	    DP2(popParams[nextMutant].timeLastUpdate);
	    DP2(popParams[u_1].timeLastUpdate);
	    DP2(popParams[u_2].timeLastUpdate);
	    DP2( (popParams[u_1].timeLastUpdate - popParams[u_2].timeLastUpdate) );
	    DP2( (popParams[u_1].timeLastUpdate - popParams[nextMutant].timeLastUpdate) );
	    DP2( (popParams[u_1].timeLastUpdate - popParams[0].timeLastUpdate) );
	    print_spP(popParams[sp]);
	    print_spP(popParams[nextMutant]);
	    throw std::out_of_range("currentTime - timeLastUpdate[sp] out of range.  ***###!!!Serious bug!!!###***");
	  }
	  if( (currentTime - popParams[nextMutant].timeLastUpdate) < 0.0) {
	    DP2(currentTime); //this is set to minNextMutationTime above
	    DP2(minNextMutationTime);
	    DP2(tSample);
	    DP2(popParams[nextMutant].timeLastUpdate);
	    DP2( (currentTime -  popParams[nextMutant].timeLastUpdate) );
	    DP2( (currentTime <  popParams[nextMutant].timeLastUpdate) );
	    DP2( (currentTime ==  popParams[nextMutant].timeLastUpdate) );
	    DP2(sp);
	    DP2(nextMutant);
	    DP2(u_1);
	    DP2(u_2);
	    DP2(tmpdouble1);
	    DP2(tmpdouble2);
	    DP2(popParams[sp].timeLastUpdate);
	    DP2(popParams[nextMutant].timeLastUpdate);
	    DP2(popParams[u_1].timeLastUpdate);
	    DP2(popParams[u_2].timeLastUpdate);
	    DP2( (popParams[u_1].timeLastUpdate - popParams[u_2].timeLastUpdate) );
	    DP2( (popParams[u_1].timeLastUpdate - popParams[nextMutant].timeLastUpdate) );
	    DP2( (popParams[u_1].timeLastUpdate - popParams[0].timeLastUpdate) );
	    print_spP(popParams[sp]);
	    print_spP(popParams[nextMutant]);
	    throw std::out_of_range("currentTime - timeLastUpdate[nextMutant] out of range. ***###!!!Serious bug!!!###***");
	  }
#endif
	  // if(verbosity >= 2) {
#ifdef DEBUGV
	    Rcpp::Rcout <<"\n     Mutated to existing species " << sp
			<< " (Genotype = ";
	    print_Genotype(Genotypes[sp]);
	      // << "; sp_id = " << Genotypes[sp].to_ullong()
	    Rcpp::Rcout << ")"
			<< "\n from species "  <<   nextMutant
			<< " (Genotypes = ";
	    print_Genotype(Genotypes[nextMutant]);
	      // << "; sp_id = " << Genotypes[sp].to_ullong()
	    Rcpp::Rcout	<< ")";
	    // }
#endif
	  // FIXME00: the if can be removed??
	    // Possibly. But note that the popParams[sp].popSize can be >
	    // 0, but when updated via Algo2 and added to 1.0 we can end
	    // in 1. Why? Because Algo2 can return a 0. The species
	    // "exist" in the sense that it had non-zero pop size when we
	    // last sampled/updated it.

	    // What we do here is step 6 of Algorithm 5, in the
	    // "Otherwise", in p. 5 of suppl mat.

	    if(popParams[sp].popSize > 0.0) {
	      popParams[sp].popSize = 1.0 +
		Algo2_st(popParams[sp], currentTime, mutationPropGrowth);
	      if(verbosity >= 2) {
		Rcpp::Rcout << "\n New popSize = " << popParams[sp].popSize << "\n";
	      }
	    } else {
	      throw std::range_error("\n popSize == 0 but existing? \n");
	    }
#ifdef DEBUGW
	  // This is wrong!!! if we set it to -999999, then the time to
  	  // next mutation will not be properly updated.  In fact, the
  	  // mapTimes map becomes a mess because the former pv in the
  	  // popParams is not removed so we end up inserting another pair
  	  // for the same species.
	  // popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	  //popParams[sp].Flag = true;

	    //zz: LOD:
	    // here one of the calls to addToPhylog, with popSize popParams[sp].popSize
	    if(keepPhylog)
	      addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
			  intName, genesInFitness, popParams[sp].popSize);


	}
	//   ***************  5.7 ***************
	// u_2 irrelevant if to_update = 1;
	u_1 = nextMutant;
	u_2 = static_cast<int>(sp);
      } else { // the null or dummy mutation case
	// We increase size by 1, as we already called Algo3. And then
	// update the ti.
	++popParams[nextMutant].popSize;
	to_update = 1;
	u_1 = nextMutant;
	u_2 = -99;
	if(verbosity >= 1)
	  Rcpp::Rcout << "Note: updating in null mutation\n";
      }
    } else { //       *********** We are sampling **********
      to_update = 3; //short_update = false;
      if(verbosity >= 2) {
	Rcpp::Rcout <<"\n We are SAMPLING";
	if(tSample < finalTime) {
	  Rcpp::Rcout << " at time " << tSample << "\n";
	} else
	  Rcpp::Rcout <<". We reached finalTime " << finalTime << "\n";
      }

      currentTime = tSample;
      if(verbosity >= 3)
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() << "\n";

      nr_sample_all_pop_P(sp_to_remove,
			  popParams, Genotypes, tSample,
			  mutationPropGrowth);
      timeNextPopSample += sampleEvery;
      // When we call nr_totPopSize ... species that existed between
      // their creation and sampling time are never reflected. That is OK.
      // This is on purpose, but if you track the phylogeny, you might see
      // in the phylogeny things that never get reflected in the pops.by.time
      // object.
      if(sp_to_remove.size())
	remove_zero_sp_nr(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();

      nr_totPopSize_and_fill_out_crude_P(outNS_i, totPopSize,
					 lastStoredSample,
					 genot_out,
					 //sp_id_out,
					 popSizes_out, index_out,
					 time_out,
					 sampleTotPopSize,sampleLargestPopSize,
					 sampleMaxNDr, sampleNDrLargestPop,
					 simulsDone,
					 reachDetection,
					 lastMaxDr,
					 done_at,
					 Genotypes, popParams,
					 currentTime,
					 keepEvery,
					 detectionSize,
					 finalTime,
					 //endTimeEvery,
					 detectionDrivers,
					 verbosity,
					 minDetectDrvCloneSz,
					 extraTime,
					 fitnessEffects.drv,
					 cPDetect,
					 PDBaseline,
					 checkSizePEvery,
					 nextCheckSizeP,
					 ran_gen,
					 AND_DrvProbExit,
					 fixation_l,
					 fixation_tolerance,
					 min_successive_fixation,
					 fixation_min_size,
					 num_successive_fixation,
					 pom, intName,
					 genesInFitness); //keepEvery is for thinning
      if(verbosity >= 3) {
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size()
		  << "\n totPopSize after sampling " << totPopSize << "\n";
      }

      // computeMcFarlandError(e1, n_0, tps_0,
      // 			    typeModel, totPopSize, K); //, initSize);
      computeMcFarlandError_new(em1, em1sc, totPopSize_previous, DA_previous,
				typeModel, totPopSize, K);

      if(simulsDone)
	break; //skip last updateRates

      // if( (typeModel == TypeModel::beerenwinkel) ) {
      // 	updateRatesBeeren(popParams, adjust_fitness_B,
      // 			  initSize, currentTime, alpha, totPopSize,
      // 			  mutationPropGrowth, mu);
      // } else if( (typeModel == TypeModel::mcfarland0) ) {
      // 	updateRatesMcFarland0(popParams, adjust_fitness_MF,
      // 			     K, totPopSize,
      // 			     mutationPropGrowth, mu);
      // } else if( (typeModel == TypeModel::mcfarland) ) {
      // 	updateRatesMcFarland(popParams, adjust_fitness_MF,
      // 			     K, totPopSize);
      // } else if( (typeModel == TypeModel::mcfarlandlog) ) {
      if (typeModel == TypeModel::mcfarlandlog && !fitnessEffects.frequencyDependentFitness){
	
	updateRatesMcFarlandLog(popParams, adjust_fitness_MF, K, totPopSize);

      } else if(typeModel == TypeModel::mcfarlandlog_d && !fitnessEffects.frequencyDependentFitness ) {

	updateRatesMcFarlandLog_D(popParams, adjust_fitness_MF, K, totPopSize);

      } else if (fitnessEffects.frequencyDependentFitness){
	
	if( (typeModel == TypeModel::mcfarlandlog) ) {
	  
	  updateRatesFDFMcFarlandLog(popParams, Genotypes, fitnessEffects,
				     adjust_fitness_MF, K, totPopSize);
	  
	} else if( (typeModel == TypeModel::mcfarlandlog_d) ) {
	  
	  updateRatesFDFMcFarlandLog_D(popParams, Genotypes, fitnessEffects,
				     adjust_fitness_MF, K, totPopSize);
	  
	} else if(typeModel == TypeModel::exp){
	  
	  updateRatesFDFExp(popParams, Genotypes, fitnessEffects);
	  
	}else if(typeModel == TypeModel::bozic1){
	  
	  updateRatesFDFBozic(popParams, Genotypes, fitnessEffects);
	  
	} else {
	  throw std::invalid_argument("this ain't a valid typeModel");
	}
      }

#ifdef MIN_RATIO_MUTS_NR
      // could go inside sample_all_pop but here we are sure death, etc, current
      // But I catch them when they are created. Is this really needed?
      for(size_t i = 0; i < popParams.size(); i++) {
	g_tmp1_nr = popParams[i].birth/popParams[i].mutation;
	if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;

	g_tmp1_nr = popParams[i].death/popParams[i].mutation;
	if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;
      }
#endif

      forceSample = false;
    }
  }
}

/*____ Funciones Modelo espacial _____ */

void printMigrant(const Migrant &deme, const fitnessEffectsAll &fe, const int verbosity)
{
  std::map<int, std::string> intName = mapGenesIntToNames(fe);
  fitness_as_genes genesInFitness = fitnessAsGenes(fe);

  std::vector<Genotype> genot_out;

  Rcpp::Rcout << "\n\tMigrant Population of {" << deme.posOrigen[0] << "," << deme.posOrigen[1] << "," << deme.posOrigen[2] << "} a {" << deme.posDestino[0] << "," << deme.posDestino[1] << "," << deme.posDestino[2] << "}"
              << " tamaño: " << deme.popParams.size();
  for (int i = 0; i < deme.popParams.size(); i++)
  {
    genot_out.clear();
    genot_out.push_back(deme.Genotypes[i]);
    std::vector<std::vector<int>> genot_out_v = genot_to_vectorg(genot_out);
    std::vector<std::string> genotypesAsStrings = genotypesToNameString(genot_out_v, genesInFitness, intName);

    Rcpp::Rcout << "\n\t\t - Clon " << i << " ---- size: " << deme.popParams[i].popSize;
    Rcpp::Rcout << " ---- genotype: " << genotypesAsStrings[0];
    if (verbosity >= 6)
    {
      Rcpp::Rcout << "\n\t\t   birth: " << deme.popParams[i].birth << " ---- death: " << deme.popParams[i].death;
      Rcpp::Rcout << " ---- numMutablePos: " << deme.popParams[i].numMutablePos;
      Rcpp::Rcout << "\n\t\t   W: " << deme.popParams[i].W << " ---- R: " << deme.popParams[i].R;
      Rcpp::Rcout << " ---- mutation: " << deme.popParams[i].mutation;
      //Rcpp::Rcout << " ---- fitness: " << dm.listClones[i].absfitness;
    }
  }
}

void printDeme(const Deme &deme, const fitnessEffectsAll &fe, const int verbosity)
{
  std::map<int, std::string> intName = mapGenesIntToNames(fe);
  fitness_as_genes genesInFitness = fitnessAsGenes(fe);

  std::vector<Genotype> genot_out;

  Rcpp::Rcout << "\n\tDeme population {" << deme.pos[0] << "," << deme.pos[1] << "," << deme.pos[2] << "} : " << deme.totPopSize << " ---- especies: " << deme.numSpecies << " ----- tamaño: " << deme.popParams.size();
  for (int i = 0; i < deme.popParams.size(); i++)
  {
    genot_out.clear();
    genot_out.push_back(deme.Genotypes[i]);
    std::vector<std::vector<int>> genot_out_v = genot_to_vectorg(genot_out);
    std::vector<std::string> genotypesAsStrings = genotypesToNameString(genot_out_v, genesInFitness, intName);

    Rcpp::Rcout << "\n\t\t - Clon " << i << " ---- size: " << deme.popParams[i].popSize;
    Rcpp::Rcout << " ---- genotype: " << genotypesAsStrings[0];
    if (verbosity >= 6)
    {
      Rcpp::Rcout << "\n\t\t   birth: " << deme.popParams[i].birth << " ---- death: " << deme.popParams[i].death;
      Rcpp::Rcout << " ---- numMutablePos: " << deme.popParams[i].numMutablePos;
      Rcpp::Rcout << "\n\t\t   W: " << deme.popParams[i].W << " ---- R: " << deme.popParams[i].R;
      Rcpp::Rcout << " ---- mutation: " << deme.popParams[i].mutation;
      //Rcpp::Rcout << " ---- fitness: " << dm.listClones[i].absfitness;
    }
  }
}

void printMapaDemes(MapaDemes mapaDemes, const fitnessEffectsAll &fe, const int verbosity)
{

  MapaDemes::iterator demes_i;

  Rcpp::Rcout << "\nTotal population: " << mapaDemes.size() << " demes:";
  for (demes_i = mapaDemes.begin(); demes_i != mapaDemes.end(); ++demes_i)
  {
    printDeme(demes_i->second, fe, verbosity);
  }
}

static void spatial_intrademe(const fitnessEffectsAll &fitnessEffects,
                              const double &initSize,
                              const double &K,
                              const TypeModel typeModel,
                              const int &mutationPropGrowth,
                              const std::vector<double> &mu,
                              const double &death,
                              const double &keepEvery,
                              const double &sampleEvery,
                              const std::vector<int> &initMutant,
                              const time_t &start_time,
                              const double &maxWallTime,
                              const double &finalTime,
                              const double &detectionSize,
                              const int &detectionDrivers,
                              const double &minDetectDrvCloneSz,
                              const double &extraTime,
                              const int &verbosity,
                              double &em1,
                              double &em1sc,
                              double &ratioForce,
                              double &currentTime,
                              int &speciesFS,
                              int &outNS_i,
                              int &iter,
                              std::vector<Genotype> &genot_out,
                              std::vector<double> &popSizes_out,
                              std::vector<int> &index_out,
                              std::vector<double> &time_out,
                              std::vector<double> &sampleTotPopSize,
                              std::vector<double> &sampleLargestPopSize,
                              std::vector<int> &sampleMaxNDr,
                              std::vector<int> &sampleNDrLargestPop,
                              bool &reachDetection,
                              std::mt19937 &ran_gen,
                              double &runningWallTime,
                              bool &hittedWallTime,
                              const std::map<int, std::string> &intName,
                              const fitness_as_genes &genesInFitness,
                              PhylogName &phylog,
                              bool keepPhylog,
                              const fitnessEffectsAll &muEF,
                              const std::vector<int> &full2mutator,
                              const double &cPDetect,
                              const double &PDBaseline,
                              const double &checkSizePEvery,
                              const bool &AND_DrvProbExit,
                              const std::vector<std::vector<int>> &fixation_l,
                              const double &fixation_tolerance,
                              const int &min_successive_fixation,
                              const double &fixation_min_size,
                              int &ti_dbl_min,
                              int &ti_e3,
                              std::map<std::string, std::string> &lod,
                              POM &pom,
                              const double &maxFitness,
                              std::vector<spParamsP> &popParams,
                              std::vector<Genotype> &Genotypes,
                              std::multimap<double, int> &mapTimes,
                              unsigned int &numSpecies,
                              double &totPopSize,
                              const double cteSize,
                              Deme &deme)
{

  // Variables de la funcion nr_totPopSize_and_fill_out_crude_P
  double nextCheckSizeP = checkSizePEvery;
  genot_out.clear();
  popSizes_out.clear();
  index_out.clear();
  time_out.clear();
  sampleTotPopSize.clear();
  sampleLargestPopSize.clear();
  sampleMaxNDr.clear();
  sampleNDrLargestPop.clear();
  outNS_i = -1;
  double lastStoredSample = 0.0;
  int lastMaxDr = 0;
  double done_at = -9;
  int num_successive_fixation = 0; // none so far

  //Variables que yo creo que cumplen funciones que no necesito de momento
  phylog = PhylogName();

  // Para Mcfarland:
  double adjust_fitness_MF = -std::numeric_limits<double>::infinity();
  em1 = 0.0;
  em1sc = 0.0;
  double totPopSize_previous = totPopSize;
  double DA_previous = log1p(totPopSize_previous / K);

  /*Variables globales necesarias */
  Genotype newGenotype;
  const int numGenes = fitnessEffects.genomeSize;
  currentTime = 0.0;
  iter = 0;
  bool forceSample = false;
  bool simulsDone = false;
  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double timeNextPopSample;
  double tSample;
  std::vector<int> newMutations;
  int nextMutant;
  int numMutablePosParent = 0;
  unsigned int sp = 0;
  int iterInterrupt = 50000; //how large should we make this?
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  std::vector<int> sp_to_remove(1);
  sp_to_remove.reserve(500);
  spParamsP tmpParam;
  ti_dbl_min = 0;
  ti_e3 = 0;

  int to_update = 3; //1: one species; 2: 2 species; 3: all.
  int u_1 = -99;
  int u_2 = -99;

  double mymindummy = 1.0e-11;     //1e-10
  double targetmindummy = 1.0e-10; //1e-9
  double minmu = *std::min_element(mu.begin(), mu.end());
  // Very small, but no less than mymindummy, for numerical issues.
  // We can probably go down to 1e-13. 1e-16 is not good as we get lots
  // of pE.f not finite. 1e-15 is probably too close, and even if no pE.f
  // we can get strange behaviors.
  double dummyMutationRate = std::max(std::min(minmu / 1.0e4, targetmindummy),
                                      mymindummy);
  // This should very rarely happen:
  if (minmu <= mymindummy)
  { // 1e-9
    double newdd = minmu / 100.0;
    Rcpp::Rcout << "WARNING: the smallest mutation rate is "
                << "<= " << mymindummy << ". That is a really small value"
                << "(per-base mutation rate in the human genome is"
                << " ~ 1e-11 to 1e-9). "
                << "Setting dummyMutationRate to your min/100 = "
                << newdd
                << ". There can be numerical problems later.\n";
    dummyMutationRate = newdd;
  }

  // EMPIEZA SIMULACION:
  timeNextPopSample = currentTime + sampleEvery; // 0 + algo
  numSpecies = popParams.size();                 //solo hemos inicializado un tipo de genotipo en toda la poblacion

  //Se deben resetear todos los parametros death, mutation de cada célula, se ponen 0. Esto se hace sobre todo
  // porque si se ha creado un deme nuevo con parametros death de otros demes, se descompenan y rompen
  // la selectividad genética produvida por el fitness. Recordamos que updateRatesConstante actualiza
  // este valor para mantener la población estable y mantiene todos los deat de todas las células iguales
  for (size_t i = 0; i < popParams.size(); ++i)
  {
    popParams[i].death = 1;
    if (popParams[i].numMutablePos == 0)
    {
      tmpParam.mutation = dummyMutationRate; // dummy mutation here. Set some mu.
    }
    else
    {
      popParams[i].mutation = mutationFromScratch(mu, popParams[i], Genotypes[i],
                                                  fitnessEffects, mutationPropGrowth,
                                                  full2mutator, muEF, Genotypes, popParams);
    }

    //Se calculan parametros W y R de popParams
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
    popParams[i].timeLastUpdate = -99999.99999; //mapTimes_updateP hara su funcion
  }

  //Mientras no se decida parar la simulacion:
  while (!simulsDone)
  {
    //if (iter%1 ==0){
    //  Rcpp::Rcout << "\nPARRA: Iteracion " << iter << " --- tiempo: " << currentTime;
    //  printDeme(deme, fitnessEffects);
    //  print_spP(deme.popParams[0]);
    //  //print_spP(demes[0].popParams[1]);
    //}

    // Miramos como vamos con el tiempo. Start_time se inicializa
    // a 0 en la funcion que llama a  inner_spatial (BNB_Algo5).
    // Este codigo no me interesa
    //runningWallTime = difftime(time(NULL), start_time);
    //// Se mira si es mayor que maxWallTime que marca el final de la simulacion
    //if( runningWallTime > maxWallTime ) {
    //  hittedWallTime = true; // ni idea de para que es este
    //  forceSample = true;
    //  simulsDone = true;
    //}

    // No estoy seguro para que sirven estas interrupciones
    // Creo que para salir de la ejecucion facilmente desde R
    iter++;
    if (!(iter % iterInterrupt))
      Rcpp::checkUserInterrupt();

    //Se calcula cunado se va a mostrar resultados la siguiente vez
    tSample = std::min(timeNextPopSample, finalTime);

    //En la primera iteracion se pasa popParams[0] en vez uno random
    if (to_update == 1)
    { // we did not sample or mutate to a different species in previous period
      // Se llama a la funcion de antes pero con una poblacion calculada la iter anterior
      tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
                                         currentTime,
                                         tSample,
                                         ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
      popParams[u_1].timeLastUpdate = currentTime;
    }
    else if (to_update == 2)
    { // we did not sample in previous period.
      // Lo mismo pero se calculan tiemps para todos los clones
      tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
                                         currentTime,
                                         tSample, ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
      tmpdouble2 = ti_nextTime_tmax_2_st(popParams[u_2],
                                         currentTime,
                                         tSample, ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, u_2, tmpdouble2);
      popParams[u_1].timeLastUpdate = currentTime;
      popParams[u_2].timeLastUpdate = currentTime;
    }
    else
    { // we sampled, so update all: i.e. to_update == 3
      // Lo mismo pero se calculan tiempos para todas los clones
      for (size_t i = 0; i < popParams.size(); i++)
      {
        tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
                                           currentTime,
                                           tSample, ti_dbl_min, ti_e3);
        mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
        popParams[i].timeLastUpdate = currentTime;
      }
    }

    // No comprendo mucho este if que papel juega, de momento ahi se queda.
    // Parece que fuerza que se sample la simulacion
    if (forceSample)
    {
      // A VERY ugly hack. Resetting tSample to jump to sampling.
      tSample = currentTime;
      // Need this, o.w. would skip a sampling.
      timeNextPopSample = currentTime;
    }

    // Obtiene el tiempo de la siguiente mutacion y en el clon que se producira
    getMinNextMutationTime4(nextMutant, minNextMutationTime,
                            mapTimes);

    // Ahora miramos a ver si queremos muestrear la poblacion:
    // Si la siguiente mutacion es antes del siguiente tiempo de muestreo...
    // ... seguimos con el algoritmo
    if (minNextMutationTime <= tSample)
    {                                    // We are not sampling
      currentTime = minNextMutationTime; // Se actualiza el tiempo actual
      // Se calcula el tiempo que hacia que ese clon no se mutaba
      mutantTimeSinceLastUpdate = currentTime - popParams[nextMutant].timeLastUpdate;
      // Se modifica el tamaño de la poblacion de ese clon
      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
                                               mutantTimeSinceLastUpdate);

      // Condicion de parada  o de muestreo que no comprendo muy bien
      if (popParams[nextMutant].popSize > (ratioForce * detectionSize))
      {
        forceSample = true;
        ratioForce = std::min(1.0, 2 * ratioForce);
      }
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet
      // know total number of different species
      // This is a protection against things going wild. Should
      // not happen in regular usage.
      // Ni idea de que es esto pero ahi se queda. numspecies empieza siendo 1,
      // speciesFS se pasa por parametro y nos da igual
      if (!(numSpecies % speciesFS))
      {
        forceSample = true;
        speciesFS *= 2;
      }

      // Entonces, si el clon que se ha decidido mutar tiene alguna celula que mutar continuamos:
      if (popParams[nextMutant].numMutablePos != 0)
      {
        // A continuacion se obtienen las mutaciones y se genera el nuevo genotipo
        newMutations.clear();
        obtainMutations(Genotypes[nextMutant],
                        fitnessEffects,
                        numMutablePosParent, // se sobreescribe con el valor que tenga
                        newMutations,
                        ran_gen,
                        mu);
        newGenotype = createNewGenotype(Genotypes[nextMutant],
                                        newMutations,
                                        fitnessEffects,
                                        ran_gen,
                                        true);

        // Se comprueba si ese genotipo ya existia
        new_sp_v(sp, newGenotype, Genotypes);

        // Si no existia se añade una especie mas:
        if (sp == numSpecies)
        { // New species
          ++numSpecies;

          // Y se inicia el proceso de crear una nueva poblacion con su popParams y tal
          init_tmpP(tmpParam);
          tmpParam.popSize = 1;
          nr_fitness(tmpParam, popParams[nextMutant],
                     newGenotype,
                     fitnessEffects,
                     typeModel, Genotypes, popParams); // , genTime,

          // Si la nueva especie tiene tasa de nacimiento positiva:
          if (tmpParam.birth > 0.0)
          {
            // Se reduce en uno el numero de mutaciones posibles
            tmpParam.numMutablePos = numMutablePosParent - 1;
            tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
                                                    fitnessEffects,
                                                    mutationPropGrowth, full2mutator,
                                                    muEF, Genotypes, popParams);

            //Control sobre el parametro mutation que no termino de saber su uso
            if (tmpParam.mutation > 1)
              Rcpp::Rcout << "WARNING: mutation > 1\n";
            // Si ya no hay mutaciones que hacer se pone ese dummyMutationRate
            if (numMutablePosParent == 1)
            {
              tmpParam.mutation = dummyMutationRate; // dummy mutation here. Set some mu.
            }
            //Se inicializan el resto de campos de la nueva poblacion
            W_f_st(tmpParam);
            R_f_st(tmpParam);
            tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP does what it should.

            // Y se añaden a la lista los nuevos genotipos y la pblacion creada.
            popParams.push_back(tmpParam);
            Genotypes.push_back(newGenotype);
            to_update = 2;
          }
          else
          { // Si el fitness es 0 o menos este clon muere instantaneamente
            --sp;
            --numSpecies;
            to_update = 1;
          }
          // Si la especie ya existia  se hace otra cosa
        }
        else
        {

          to_update = 2;
          // Si esa poblacion tiene mas de un individuo se actualiza la poblacion
          if (popParams[sp].popSize > 0.0)
          {
            popParams[sp].popSize = 1.0 +
                                    Algo2_st(popParams[sp], currentTime, mutationPropGrowth);
          }
          else
          {
            throw std::range_error("\n popSize == 0 but existing? \n");
          }
          // Esto es para mantener unos datos que de momento yo no hago, asi que se puede quitar
          if (keepPhylog)
            addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
                        intName, genesInFitness, popParams[sp].popSize);
        }

        // Se buscaran nuevos tiempos de siguiente mutacion en el clon que acaba de ser mutado y del nuevo generado
        u_1 = nextMutant;
        u_2 = static_cast<int>(sp);

        // En el caso de que el clon elegido no tenga celulas que mutar entonces
        // estamos en el caso de null o dummyMutation
      }
      else
      {
        // Actualizamos la poblacion de ese clon
        ++popParams[nextMutant].popSize;
        to_update = 1;
        u_1 = nextMutant; //calcularemos tiempos de ese clon en la siguiente iteracion
        u_2 = -99;
      }

      // Si vamos por aqui es que es momento de muestrear
    }
    else
    {                        //       *********** We are sampling **********
      to_update = 3;         //Se actualizaran los tiempos minimos de tooodos los individuos
      currentTime = tSample; //El tiempo actual lo paramos en el tiempo de muestreo

      // Esta funcion coge cada poblacion y actualiza su poblacion total.
      // Guarda en sp_to_remove aquellas poblaciones que mueren.
      nr_sample_all_pop_P(sp_to_remove,
                          popParams, Genotypes, tSample,
                          mutationPropGrowth);

      timeNextPopSample += sampleEvery;

      //Eliminamos las especies que han muerto en esta ronda.
      if (sp_to_remove.size())
      {
        remove_zero_sp_nr(sp_to_remove, Genotypes, popParams, mapTimes);
      }
      numSpecies = popParams.size();

      /* COMENTAMOS ESTE CODIGO PORQUE QUEREMOS PRESCINDIR DE ÉL
      // Esta funcion hace mierdas como:
      // detectar condiciones de parada en base a ver si poblaciones se han fijado.
      // sacar datos totales de poblacion y cosas así
      // Ver si se ha superado el tiempo maximo de ejecucion..
      // pero NO TOCA NADA DE LA SIMLACION, solo decide condiciones de parada y salida
      nr_totPopSize_and_fill_out_crude_P(outNS_i, totPopSize, 
           lastStoredSample,
           genot_out, 
           //sp_id_out,
           popSizes_out, index_out,
           time_out, 
           sampleTotPopSize,sampleLargestPopSize,
           sampleMaxNDr, sampleNDrLargestPop,
           simulsDone,
           reachDetection,
           lastMaxDr,
           done_at,
           Genotypes, popParams, 
           currentTime,
           keepEvery,
           detectionSize,
           finalTime,
           //endTimeEvery,
           detectionDrivers,
           verbosity,
           minDetectDrvCloneSz,
           extraTime,
           fitnessEffects.drv,
           cPDetect,
           PDBaseline,
           checkSizePEvery,
           nextCheckSizeP,
           ran_gen,
           AND_DrvProbExit,
           fixation_l,
           fixation_tolerance,
           min_successive_fixation,
           fixation_min_size,
           num_successive_fixation,
           pom, intName,
           genesInFitness); //keepEvery is for thinning

      // Saca un dato de error que me da igual de momento
      computeMcFarlandError_new(em1, em1sc, totPopSize_previous, DA_previous, 
        typeModel, totPopSize, K); 
*/

      //Calculo de la poblacion total:
      totPopSize = 0;
      for (size_t i = 0; i < popParams.size(); ++i)
      {
        totPopSize += popParams[i].popSize;
      }

      // Condiciones de parada:
      if ((totPopSize <= 0.0) || (currentTime >= finalTime))
      {
        simulsDone = true;
      }

      // Si hemos terminado la simulacion salimos sin updatear
      if (simulsDone)
        break; //skip last updateRates

      // Updateamos ratios de nacimiento y muerte para la siguiente iteracion
      if ((typeModel == TypeModel::mcfarlandlog))
      {
        updateRatesMcFarlandLog(popParams, adjust_fitness_MF,
                                K, totPopSize);
      }
      if (typeModel == TypeModel::constante)
      {
        updateRatesConstante(popParams, Genotypes, totPopSize, initSize, cteSize, fitnessEffects, maxFitness, sampleEvery);
      }

      forceSample = false;
    }
  }
}

static void spatial_addDemeFromImmigrant(Deme &deme,
                                         Migrant &migrant,
                                         const double initSize,
                                         const int numGenes,
                                         const fitnessEffectsAll &fitnessEffects,
                                         const int &mutationPropGrowth,
                                         const fitnessEffectsAll &muEF,
                                         const std::vector<int> &full2mutator,
                                         const TypeModel typeModel,
                                         const double &K,
                                         const double &death,
                                         const std::vector<double> &mu)
{

  size_t i, numClonWithoutMuts;
  size_t sp_per_period = 10;
  spParamsP tmpParam;
  Genotype newGenotype;
  bool clonWithoutMuts = false;

  deme.popParams.clear();
  deme.Genotypes.clear();
  deme.numSpecies = 0;
  deme.totPopSize = 0;
  deme.popParams.reserve(sp_per_period);
  deme.Genotypes.reserve(sp_per_period);
  deme.pos[0] = migrant.posDestino[0];
  deme.pos[1] = migrant.posDestino[1];
  deme.pos[2] = migrant.posDestino[2];
  // Creamos un popParam por cada especie que lleve consigo el migrante
  for (i = 0; i < migrant.popParams.size(); i++)
  {
    deme.popParams.push_back(migrant.popParams[i]);
    deme.Genotypes.push_back(migrant.Genotypes[i]);
    deme.popParams[i].timeLastUpdate = -99999.99999;
    deme.numSpecies++;
    deme.totPopSize += deme.popParams[i].popSize;
    if (migrant.popParams[i].numMutablePos == numGenes)
    {
      numClonWithoutMuts = i;
      clonWithoutMuts = true;
    }
  }

  //Ahora añadimos initSize celulas sin mutar al deme (las initSize originales de antes de llegar los migrantes)
  if (clonWithoutMuts)
  {
    deme.popParams[numClonWithoutMuts].popSize += initSize;
    deme.totPopSize += initSize;
  }
  else
  {
    init_tmpP(tmpParam);
    tmpParam.popSize = initSize;
    tmpParam.numMutablePos = numGenes;
    newGenotype = wtGenotype();

    // Calcula los parametros birth and death
    if (typeModel == TypeModel::mcfarlandlog)
    {
      tmpParam.birth = 1.0;
      tmpParam.death = log1p(initSize / K);
    }
    else if (typeModel == TypeModel::bozic1)
    {
      tmpParam.birth = 1.0;
      tmpParam.death = 1.0;
    }
    else if (typeModel == TypeModel::exp)
    {
      tmpParam.birth = 1.0;
      tmpParam.death = death;
    }
    else if (typeModel == TypeModel::constante)
    {
      tmpParam.birth = 1.0;
      tmpParam.death = 1.0;
    }
    else
    {
      throw std::invalid_argument("this ain't a valid typeModel");
    }

    // Inicializa el parametro mutation, pero no se para que sirve este parametro
    tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
                                            fitnessEffects, mutationPropGrowth,
                                            full2mutator, muEF, deme.Genotypes, deme.popParams);


    //Se calculan parametros W y R de popParams
    W_f_st(tmpParam);
    R_f_st(tmpParam);
    tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP hara su funcion

    // Se añaden al deme los genotipos y spParams
    deme.Genotypes.push_back(newGenotype);
    deme.popParams.push_back(tmpParam);

    //Se actualiza poblacion total y numero de especies
    deme.numSpecies++;
    deme.totPopSize += initSize;
  }
}

static void spatial_chooseMigrants(MapaDemes &mapaDemes,
                                   std::vector<Migrant> &migrants,
                                   std::mt19937 &ran_gen,
                                   const double migrationProb,
                                   const double maxMigrationPercentage,
                                   const double cteNumMigrations,
                                   const TypeModel typeModel)
{

  MapaDemes::iterator deme_i;
  Migrant tmpMig;
  std::uniform_real_distribution<double> runif;
  bool cellAded;
  size_t i, j, k, nonRepeated, numInmi;
  int numMigrationCells, aux;
  double currentSize;

  //Elegimos celulas que migran:
  // La cantidad sera un numero aleatorio entre 1 y maxMigrationPorcentage multiplicado por la cteNumMigrations
  numMigrationCells = std::max((double)1, ceil(1 * maxMigrationPercentage * cteNumMigrations));
  //numMigrationCells = std::max((double)1,ceil(runif(ran_gen)*maxMigrationPercentage*cteNumMigrations));
  // Se eligen aleatoriamente un vector de numeros entre 1 y cte de Moran.
  std::vector<int> migrationCells(numMigrationCells);
  for (i = 0; i < numMigrationCells; i++)
  {
    cellAded = false;
    aux = ceil(runif(ran_gen) * (cteNumMigrations - 1)); //Se escoge una al azar entre 1 y cteNumMigrations-1
    while (!cellAded)
    { //Se comprueba si esta repetida
      if (std::find(migrationCells.begin(), migrationCells.end(), aux) != migrationCells.end())
      { //Si lo esta se le suma 1
        aux++;
        continue;
      }
      else
      { //Cuando deje de estar repetida se añade
        migrationCells[i] = aux;
        cellAded = true;
      }
    }
  }
  std::sort(migrationCells.begin(), migrationCells.end()); //Se ordenan en orden creciente

  //DEBUG
  //Rcpp::Rcout << "\n NumMigrationCells " << migrationCells.size() << ": ";
  //for (i=0; i<numMigrationCells; i++){
  //  Rcpp::Rcout << migrationCells[i] << ", " ;
  //}

  //Escogemos las celulas que migraran de cada deme:
  // Por cada deme se recorre el vector de migrationCells, mientras que el numero de celula
  // sea menor que el popSize del primer popParam de ese deme, se van añadiendo las migradas a un popParam
  // nuevo. Cuando el numero es mayor se suma el popSize del siguiente popParam y se van añadiendo las migradas
  // al nuevo popParam. Y así sucesivamente.
  for (deme_i = mapaDemes.begin(); deme_i != mapaDemes.end(); ++deme_i)
  {
    // Primero vemos si vamos a realizar migracion de ese deme o no:
    if (runif(ran_gen) > migrationProb)
    {
      continue;
    }
    currentSize = deme_i->second.popParams[0].popSize;
    k = 0;
    j = 0;
    numInmi = 0;
    nonRepeated = true;
    while (j < numMigrationCells)
    { //Si quedan mas celulas del mismo popParam se van añadiendo
      if (migrationCells[j] <= currentSize)
      {
        if (nonRepeated)
        {
          tmpMig.popParams.push_back(deme_i->second.popParams[k]);
          tmpMig.Genotypes.push_back(deme_i->second.Genotypes[k]);
          tmpMig.posOrigen[0] = deme_i->second.pos[0];
          tmpMig.posOrigen[1] = deme_i->second.pos[1];
          tmpMig.posOrigen[2] = deme_i->second.pos[2];
          tmpMig.popParams[numInmi].popSize = 1;
          nonRepeated = false;
          numInmi++;
        }
        else
        {
          tmpMig.popParams[numInmi - 1].popSize++;
        }
        j++;
      }
      else
      { // si no quedan mas celulas de ese popParam, se cierra y se añade al migrant
        nonRepeated = true;
        k++;
        if (k < deme_i->second.popParams.size())
        { //Si  quedan mas popParams en ese deme se va sumando su popSize
          currentSize += deme_i->second.popParams[k].popSize;
        }
        else
        { // si no ya no entran mas celulas y terminamos la eleccion de migrantes de ese deme
          break;
        }
      }
    }
    // Por ultimo añadimos el migrante si este tiene alguna celula
    if (tmpMig.popParams.size() > 0)
    {
      migrants.push_back(tmpMig);
      tmpMig.popParams.clear();
      tmpMig.Genotypes.clear();
    }
  }
}

static void spatial_chooseDestination(std::vector<Migrant> &migrants,
                                      const TypeMigration typeMigration,
                                      std::mt19937 &ran_gen)
{

  std::uniform_real_distribution<double> runif;
  size_t i, neighbour, x, y, z;

  // Elegir destino
  if (typeMigration == TypeMigration::shortDistanceMoore)
  {
    int mooreNeighborhood[26][3] = {{-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1}, {0, -1, -1}, {0, -1, 0}, {0, -1, 1}, {0, 0, -1}, {0, 0, 1}, {0, 1, -1}, {0, 1, 0}, {0, 1, 1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0}, {1, 0, 1}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};

    for (i = 0; i < migrants.size(); i++)
    {
      neighbour = floor(runif(ran_gen) * 26);
      migrants[i].posDestino[0] = migrants[i].posOrigen[0] + mooreNeighborhood[neighbour][0];
      migrants[i].posDestino[1] = migrants[i].posOrigen[1] + mooreNeighborhood[neighbour][1];
      migrants[i].posDestino[2] = migrants[i].posOrigen[2] + mooreNeighborhood[neighbour][2];
    }
  }
  else if (typeMigration == TypeMigration::largeDistanceMoore)
  {
    int min = 10;
    int interval = 25;
    for (i = 0; i < migrants.size(); i++)
    {
      x = (floor(runif(ran_gen) * interval) + min) * ((runif(ran_gen) < 0.5) ? 1 : -1);
      y = (floor(runif(ran_gen) * interval) + min) * ((runif(ran_gen) < 0.5) ? 1 : -1);
      z = (floor(runif(ran_gen) * interval) + min) * ((runif(ran_gen) < 0.5) ? 1 : -1);
      migrants[i].posDestino[0] = migrants[i].posOrigen[0] + x;
      migrants[i].posDestino[1] = migrants[i].posOrigen[1] + y;
      migrants[i].posDestino[2] = migrants[i].posOrigen[2] + z;
    }
  }
  else
  {
    throw std::logic_error("This TypeMigration don't exist");
  }
}

static void spatial_DeleteEmptyDemes(MapaDemes &mapaDemes, const bool spatialDeleteDemesWithoutMutations, const int numGenes, const size_t verbosity)
{
  size_t i;
  MapaDemes::iterator deme_i;
  std::vector<std::vector<int>> demesToRemove;

  // Apuntamos la posicion de los demes que debemos eliminar, todos aquellos con poblacion = 0.
  for (deme_i = mapaDemes.begin(); deme_i != mapaDemes.end(); ++deme_i)
  {
    if ((deme_i->second.totPopSize == 0) && (deme_i->second.popParams.size() == 0))
    {
      demesToRemove.push_back(deme_i->first);
    } // Ahora miramos si hay demes sin celulas mutadas para eliminarlos tambien
    else if ((deme_i->second.popParams.size() == 1) &&
             (deme_i->second.popParams[0].numMutablePos == numGenes) &&
             (spatialDeleteDemesWithoutMutations))
    {
      demesToRemove.push_back(deme_i->first);
    }
  }

  // Los eliminamos
  for (i = 0; i < demesToRemove.size(); i++)
  {
    if (verbosity >= 4)
      Rcpp::Rcout << "\n Se ha eliminado el deme (" << demesToRemove[i][0] << "," << demesToRemove[i][2] << "," << demesToRemove[i][2] << ")";
    mapaDemes.erase(demesToRemove[i]);
  }
}

static void spatial_Migration(MapaDemes &mapaDemes,
                              std::mt19937 &ran_gen,
                              const double migrationProb,
                              const TypeMigration typeMigration,
                              const double maxMigrationPercentage,
                              const double cteNumMigrations,
                              const bool spatialDeleteDemesWithoutMutations,
                              const double initSize,
                              const int numGenes,
                              const fitnessEffectsAll &fitnessEffects,
                              const int &mutationPropGrowth,
                              const fitnessEffectsAll &muEF,
                              const std::vector<int> &full2mutator,
                              const TypeModel typeModel,
                              const double &K,
                              const double &death,
                              const std::vector<double> &mu,
                              const int verbosity)
{

  MapaDemes::iterator deme_i;
  Deme tmpDeme;
  Migrant tmpMig;
  std::vector<Migrant> migrants;
  size_t i, j, k;
  unsigned int sp = 0;
  int x, y, z;

  spatial_chooseMigrants(mapaDemes,
                         migrants,
                         ran_gen,
                         migrationProb,
                         maxMigrationPercentage,
                         cteNumMigrations,
                         typeModel);

  spatial_chooseDestination(migrants,
                            typeMigration,
                            ran_gen);

  //DEBUG
  if (verbosity >= 5)
  {
    for (i = 0; i < migrants.size(); i++)
    {
      printMigrant(migrants[i], fitnessEffects, 1);
    }
  }

  // Introducir celulas en destino
  for (i = 0; i < migrants.size(); i++)
  {
    x = migrants[i].posDestino[0];
    y = migrants[i].posDestino[1];
    z = migrants[i].posDestino[2];
    deme_i = mapaDemes.find({x, y, z});
    if (deme_i == mapaDemes.end())
    { //No hay deme en esa posicion
      spatial_addDemeFromImmigrant(tmpDeme,
                                   migrants[i],
                                   initSize,
                                   numGenes,
                                   fitnessEffects,
                                   mutationPropGrowth,
                                   muEF,
                                   full2mutator,
                                   typeModel,
                                   K,
                                   death,
                                   mu);
      mapaDemes.insert(MapaDemes::value_type({migrants[i].posDestino[0], migrants[i].posDestino[1], migrants[i].posDestino[2]}, tmpDeme));
    }
    else
    { // Hay deme en esa posicion
      for (j = 0; j < migrants[i].popParams.size(); j++)
      {
        new_sp_v(sp, migrants[i].Genotypes[j], deme_i->second.Genotypes);
        if (sp == deme_i->second.numSpecies)
        { // Si genotipo no esta en deme crearlo
          deme_i->second.numSpecies++;
          deme_i->second.popParams.push_back(migrants[i].popParams[j]);
          deme_i->second.Genotypes.push_back(migrants[i].Genotypes[j]);
          deme_i->second.popParams.back().timeLastUpdate = -99999.99999;
          deme_i->second.totPopSize += deme_i->second.popParams.back().popSize;
        }
        else
        { // Si genotipo esta en deme añadir poblacion
          deme_i->second.popParams[j].popSize += migrants[i].popParams[j].popSize;
          deme_i->second.totPopSize += migrants[i].popParams[j].popSize;
        }
      }
    }
  }

  std::vector<int> sp_to_remove;
  // Restar celulas de origen
  for (i = 0; i < migrants.size(); i++)
  {
    x = migrants[i].posOrigen[0];
    y = migrants[i].posOrigen[1];
    z = migrants[i].posOrigen[2];
    deme_i = mapaDemes.find({x, y, z});
    if (deme_i != mapaDemes.end())
    {
      sp_to_remove.clear();
      for (j = 0; j < migrants[i].popParams.size(); j++)
      {
        new_sp_v(sp, migrants[i].Genotypes[j], deme_i->second.Genotypes);
        if (sp != deme_i->second.numSpecies)
        { // Si genotipo esta en deme restamos poblacion
          deme_i->second.popParams[sp].popSize -= migrants[i].popParams[j].popSize;
          deme_i->second.totPopSize -= migrants[i].popParams[j].popSize;
          if (deme_i->second.popParams[sp].popSize == 0)
          { //Si la poblacion de ese popParam se queda a 0 la eliminamos
            sp_to_remove.push_back(sp);
          }
          else if (deme_i->second.popParams[sp].popSize < 0)
          {
            throw std::logic_error("There is a migrant with more cells than the popParams deme popSize?");
          }
        }
        else
        { // Si genotipo no esta en deme excepcion
          throw std::logic_error("There is a migrant from inexistent popParams of deme?");
        }
      }
      remove_zero_sp_nr(sp_to_remove, deme_i->second.Genotypes, deme_i->second.popParams, deme_i->second.mapTimes);
      deme_i->second.numSpecies = deme_i->second.popParams.size();
    }
    else
    {
      throw std::logic_error("There is a migrant from inexistent deme?");
    }
  }

  //Eliminar demes vacíos.
  spatial_DeleteEmptyDemes(mapaDemes, spatialDeleteDemesWithoutMutations, numGenes, verbosity);
}

static void spatial_interdeme(MapaDemes &mapaDemes,
                              std::mt19937 &ran_gen,
                              const double migrationProb,
                              const double largeDistanceMigrationProb,
                              const double maxMigrationPercentage,
                              const bool spatialDeleteDemesWithoutMutations,
                              const double cteSize,
                              const double initSize,
                              const int numGenes,
                              const fitnessEffectsAll &fitnessEffects,
                              const int &mutationPropGrowth,
                              const fitnessEffectsAll &muEF,
                              const std::vector<int> &full2mutator,
                              const TypeModel typeModel,
                              const double &K,
                              const double &death,
                              const std::vector<double> &mu,
                              const int verbosity)
{

  MapaDemes::iterator deme_i;
  double cteNumMigrations;
  int modeChooseDestination;
  TypeMigration typeMigration;
  double migrationProbability;

  // Calculamos una constante utilizada para calcular la cantidad de celulas que migran.
  // Esta cte es la poblacion media de todos los demes
  if (typeModel == TypeModel::constante)
  {
    cteNumMigrations = cteSize;
  }
  else
  {
    cteNumMigrations = 0;
    for (deme_i = mapaDemes.begin(); deme_i != mapaDemes.end(); ++deme_i)
    {
      cteNumMigrations += deme_i->second.totPopSize;
    }
    cteNumMigrations = ceil(cteNumMigrations / mapaDemes.size());
  }

  // Migracion a corta distancia
  typeMigration = TypeMigration::shortDistanceMoore;
  migrationProbability = migrationProb;
  spatial_Migration(mapaDemes,
                    ran_gen,
                    migrationProbability,
                    typeMigration,
                    maxMigrationPercentage,
                    cteNumMigrations,
                    spatialDeleteDemesWithoutMutations,
                    initSize,
                    numGenes,
                    fitnessEffects,
                    mutationPropGrowth,
                    muEF,
                    full2mutator,
                    typeModel,
                    K,
                    death,
                    mu,
                    verbosity);

  // Migracion a larga distancia
  typeMigration = TypeMigration::largeDistanceMoore;
  migrationProbability = largeDistanceMigrationProb;
  spatial_Migration(mapaDemes,
                    ran_gen,
                    migrationProbability,
                    typeMigration,
                    maxMigrationPercentage,
                    cteNumMigrations,
                    spatialDeleteDemesWithoutMutations,
                    initSize,
                    numGenes,
                    fitnessEffects,
                    mutationPropGrowth,
                    muEF,
                    full2mutator,
                    typeModel,
                    K,
                    death,
                    mu,
                    verbosity);
}

void spatial_initialization(MapaDemes &mapaDemes,
                            const double initSize,
                            const double K,
                            const fitnessEffectsAll &fitnessEffects,
                            const std::vector<int> &initMutant,
                            const double initMutantPercentage,
                            std::mt19937 &ran_gen,
                            const TypeModel typeModel,
                            const int &mutationPropGrowth,
                            const std::vector<double> &mu,
                            const double &death,
                            const fitnessEffectsAll &muEF,
                            const std::vector<int> &full2mutator)
{

  /* Parametros obligatorios globales*/
  const int sp_per_period = 10; //Reservamos memoria para vector de hasta sp_per_period popParams y genotipos
  const int numGenes = fitnessEffects.genomeSize;
  spParamsP tmpParam;
  Genotype newGenotype;
  std::vector<spParamsP>::iterator popParams_begin;
  std::vector<Genotype>::iterator Genotypes_begin;

  /* Inicializamos un deme con todas las celulas sin mutar menos una mutada */
  Deme tmpDeme;
  tmpDeme.popParams.reserve(sp_per_period);
  tmpDeme.Genotypes.reserve(sp_per_period);
  tmpDeme.totPopSize = initSize;
  tmpDeme.pos[0] = 0;
  tmpDeme.pos[1] = 0;
  tmpDeme.pos[2] = 0;
  tmpDeme.numSpecies = 2; //Se inicializan a continuacion

  /* ____ Primero las celulas sin mutar */
  init_tmpP(tmpParam);
  tmpParam.popSize = initSize;
  tmpParam.numMutablePos = numGenes;
  newGenotype = wtGenotype();

  // Calcula los parametros birth and death
  if (typeModel == TypeModel::mcfarlandlog)
  {
    tmpParam.birth = 1.0;
    tmpParam.death = log1p(tmpDeme.totPopSize / K);
  }
  else if (typeModel == TypeModel::bozic1)
  {
    tmpParam.birth = 1.0;
    tmpParam.death = 1.0;
  }
  else if (typeModel == TypeModel::exp)
  {
    tmpParam.birth = 1.0;
    tmpParam.death = death;
  }
  else if (typeModel == TypeModel::constante)
  {
    tmpParam.birth = 1.0;
    tmpParam.death = 1.0;
  }
  else
  {
    throw std::invalid_argument("this ain't a valid typeModel");
  }

  // Inicializa el parametro mutation, pero no se para que sirve este parametro
  tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
                                          fitnessEffects, mutationPropGrowth,
                                          full2mutator, muEF, tmpDeme.Genotypes, tmpDeme.popParams);


  //Se calculan parametros W y R de popParams
  W_f_st(tmpParam);
  R_f_st(tmpParam);
  tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP hara su funcion

  // Se añaden al deme los genotipos y spParams
  tmpDeme.Genotypes.push_back(newGenotype);
  tmpDeme.popParams.push_back(tmpParam);

  /* _____Ahora inicializamos la celula mutada */
  init_tmpP(tmpParam);
  tmpParam.popSize = std::max((double)1, ceil(initMutantPercentage * initSize));
  tmpDeme.popParams[0].popSize = initSize - tmpParam.popSize; //Actualizamos la poblacion de las celulas sin mutar
  tmpParam.numMutablePos = numGenes - initMutant.size();      //Numero de posibles mutaciones en ese genotipo

  /*Se crea un genotipo con un conjunto de celulas mutadas (initMutant) */
  newGenotype = createNewGenotype(wtGenotype(),
                                  initMutant,
                                  fitnessEffects,
                                  ran_gen,
                                  false);

  // Esta funcion actualiza birth y death del primer argumento con lo calculado en el segundo, teniendo en cuenta alguna cosa mas
  nr_fitness(tmpParam, tmpDeme.popParams[0], //El padre del nuevo mutante es la celula sin mutar
             newGenotype,
             fitnessEffects,
             typeModel, tmpDeme.Genotypes, tmpDeme.popParams);

  // Inicializa el parametro mutation, pero no se para que sirve este parametro
  tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
                                          fitnessEffects, mutationPropGrowth,
                                          full2mutator, muEF, tmpDeme.Genotypes, tmpDeme.popParams);

  //Se calculan parametros W y R de popParamsP
  W_f_st(tmpParam);
  R_f_st(tmpParam);
  tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP hara su funcion

  // Se añaden al deme los genotipos y spParamsP
  tmpDeme.Genotypes.push_back(newGenotype);
  tmpDeme.popParams.push_back(tmpParam);

  // Miramos si la primera poblacion tiene poblacion <= 0
  if (tmpDeme.popParams[0].popSize <= 0)
  {
    popParams_begin = tmpDeme.popParams.begin();
    Genotypes_begin = tmpDeme.Genotypes.begin();
    tmpDeme.popParams.erase(popParams_begin);
    tmpDeme.Genotypes.erase(Genotypes_begin);
    tmpDeme.numSpecies = tmpDeme.popParams.size();
    Rcpp::Rcout << "\n Solo creamos una poblacion mutada de tamaño " << tmpDeme.popParams[0].popSize;
  }

  /*Queda inicializado el primer deme*/
  /* _______________________________________________*/

  // Se inserta el deme creado en el mapa de demes
  mapaDemes.insert(MapaDemes::value_type({tmpDeme.pos[0], tmpDeme.pos[1], tmpDeme.pos[2]}, tmpDeme));
}

static void spatial_Output(SpatialOutput &so, MapaDemes &mapaDemes, SpatialSimulation &ss,
                           const fitnessEffectsAll &fe, const int verbosity)
{
  std::map<int, std::string> intName = mapGenesIntToNames(fe);
  fitness_as_genes genesInFitness = fitnessAsGenes(fe);
  std::vector<std::vector<int>> genot_out_v;
  std::vector<std::string> genotypesAsStrings;
  ofstream myfile;

  //SpatialOutput_aux1 so_aux1;
  //std::vector<SpatialOutput_aux1> vector_so_aux1;
  MapaSimulacion::iterator ss_i;

  std::vector<Genotype> vector_Genotypes;
  std::vector<int> vector_popSize;
  std::vector<int> vector_numDemes;
  std::vector<int> vector_numDemesGanador;

  MapaDemes::iterator deme_i;
  unsigned int sp = 0;
  size_t i, j, numSpecies, maxPop, maxGenot = 0;

  //________Preparar salida y simulacion
  // Si la poblacion se extingue en la primera iteracion ss.ms.size() = 0
  i = 0;
  if (ss.ms.size() != 0)
  {
    std::string file = "OncoSimulR.spatial";
    myfile.open(file);
    myfile << "'X'  'Y'  'Z'  ";
    for (j = 0; j < ss.keepEveryIters; j++)
    {
      myfile << "'g_" << j << "'  'e_" << j << "'  ";
    }
    myfile << "\n";

    //Ahora escribimos la informacion de cada deme:
    for (ss_i = ss.ms.begin(); ss_i != ss.ms.end(); ++ss_i)
    {
      i++;
      myfile << "'" << i << "'  ";
      myfile << ss_i->first[0] << "  " << ss_i->first[1] << "  " << ss_i->first[2] << "  ";
      for (j = 0; j < ss_i->second.comienzo; j++)
      {
        myfile << "-1  0  ";
      }
      for (j = 0; j < ss_i->second.genotipos.size(); j++)
      {
        if (ss_i->second.genotipos[j] == -1)
        {
          myfile << "-1  0  ";
        }
        else
        {
          myfile << ss_i->second.genotipos[j] + 1 << "  1  ";
        }
      }
      myfile << "\n";
    }
    myfile.close();

    // Genotipos en orden para escribir la leyenda de la simulacion
    genot_out_v = genot_to_vectorg(ss.vector_Genotypes);
    genotypesAsStrings = genotypesToNameString(genot_out_v, genesInFitness, intName);

    //Infomacion de demes, demes ganadores y popSize por iteracion
    numSpecies = ss.vector_Genotypes.size();
    for (i = 0; i < ss.vector_iter_numDemes.size(); i++)
    {
      for (j = ss.vector_iter_numDemes[i].size(); j < numSpecies; j++)
      {
        ss.vector_iter_numDemes[i].push_back(0);
      }
    }
    so.vector_iter_numDemes = ss.vector_iter_numDemes;

    for (i = 0; i < ss.vector_iter_numDemesGanador.size(); i++)
    {
      for (j = ss.vector_iter_numDemesGanador[i].size(); j < numSpecies; j++)
      {
        ss.vector_iter_numDemesGanador[i].push_back(0);
      }
    }
    so.vector_iter_numDemesGanador = ss.vector_iter_numDemesGanador;

    for (i = 0; i < ss.vector_iter_popSize.size(); i++)
    {
      for (j = ss.vector_iter_popSize[i].size(); j < numSpecies; j++)
      {
        ss.vector_iter_popSize[i].push_back(0);
      }
    }
    so.vector_iter_popSize = ss.vector_iter_popSize;

    // Resto de salida
    so.vector_Genotypes = genotypesAsStrings;
    so.vector_popSize = so.vector_iter_popSize[ss.keepEveryIters - 1];
    so.vector_numDemes = so.vector_iter_numDemes[ss.keepEveryIters - 1];
    so.vector_numDemesGanador = so.vector_iter_numDemesGanador[ss.keepEveryIters - 1];

    // Poblacion total de demes y de celulas
    so.numDemes = mapaDemes.size();
    so.totPopSize = 0;
    for (deme_i = mapaDemes.begin(); deme_i != mapaDemes.end(); ++deme_i)
    {
      for (j = 0; j < deme_i->second.Genotypes.size(); j++)
      {
        so.totPopSize += deme_i->second.popParams[j].popSize;
      }
    }
  }
  else
  {
    so.vector_iter_numDemesGanador = {};
    so.vector_iter_numDemes = {};
    so.vector_iter_popSize = {};
    so.numDemes = 0;
    so.totPopSize = 0;
    so.vector_Genotypes = {};
    so.vector_popSize = {};
    so.vector_numDemes = {};
    so.vector_numDemesGanador = {};
  }

  //______________________________________________________________________
}

static void spatial_KeepSimulation(MapaDemes &mapaDemes, SpatialSimulation &ss)
{

  MapaDemes::iterator deme_i;
  MapaSimulacion::iterator ss_i;
  SpatialSimulation_aux1 ss_aux1;

  std::vector<int> vector_numDemesGanador = {};
  std::vector<int> vector_numDemes = {};
  std::vector<int> vector_popSize = {};
  size_t i, maxPop = 0, maxGenot = 0, numSpecies;
  unsigned int sp = 0;

  numSpecies = ss.vector_Genotypes.size();
  //Inicializamos estos vectores temporales para guardar popSize y demes ganadores.
  // Esta inicializacion permite tambien que funcione si llega una poblacion extinta
  for (i = 0; i < numSpecies; i++)
  {
    vector_numDemesGanador.push_back(0);
    vector_numDemes.push_back(0);
    vector_popSize.push_back(0);
  }

  // En este bucle calculamos genotipo mas grande, controamos la cantidad de genotipos
  // encontrados y añadimos la poblacion a cada genotipo.
  for (deme_i = mapaDemes.begin(); deme_i != mapaDemes.end(); ++deme_i)
  {

    // Escogemos el genotipo mas grande del deme
    maxPop = 0;
    for (i = 0; i < deme_i->second.Genotypes.size(); i++)
    {

      //Miramos si ese genotipo esta en la lista de genotipos actual y si no lo añadimos
      new_sp_v(sp, deme_i->second.Genotypes[i], ss.vector_Genotypes);
      if (sp == numSpecies)
      {
        numSpecies++;
        ss.vector_Genotypes.push_back(deme_i->second.Genotypes[i]);
        vector_numDemesGanador.push_back(0);
        vector_numDemes.push_back(0);
        vector_popSize.push_back(0);
      }

      //Añadimos poblacion de ese genotipo:
      vector_popSize[sp] += deme_i->second.popParams[i].popSize;
      vector_numDemes[sp] += 1;

      //Para calcular el genotipo mas grande
      if (deme_i->second.popParams[i].popSize >= maxPop)
      {
        maxPop = deme_i->second.popParams[i].popSize;
        maxGenot = i;
      }
    }

    //Miramos si ese genotipo esta en la lista de genotipos actual y si no lo añadimos
    new_sp_v(sp, deme_i->second.Genotypes[maxGenot], ss.vector_Genotypes);
    if (sp == numSpecies)
    {
      throw std::logic_error("No puede haber genotipos nuevos ya que han sido añadidos antes");
    }

    //Añadimos deme ganador
    vector_numDemesGanador[sp] += 1;

    // Ahora actualizamos la estrucutura donde guardamos por cada deme e iteracion su genotipo maximo
    //Miramos si la posicion de ese deme ya existia
    ss_i = ss.ms.find({deme_i->second.pos[0], deme_i->second.pos[1], deme_i->second.pos[2]});
    if (ss_i == ss.ms.end())
    { // Si no existia lo añadimos
      ss_aux1.comienzo = ss.keepEveryIters;
      ss_aux1.currentIter = ss.keepEveryIters;
      ss_aux1.genotipos.push_back(sp);
      ss.ms.insert(MapaSimulacion::value_type({deme_i->second.pos[0], deme_i->second.pos[1], deme_i->second.pos[2]}, ss_aux1));
      ss_aux1.genotipos.clear();
    }
    else
    { // Si ya existia añadimos el genotipo nuevo
      ss_i->second.genotipos.push_back(sp);
      ss_i->second.currentIter++;
    }
  }

  //Ahora actualizamos el resto de demes:
  for (ss_i = ss.ms.begin(); ss_i != ss.ms.end(); ++ss_i)
  {
    if (ss_i->second.currentIter != ss.keepEveryIters)
    {
      ss_i->second.currentIter++;
      ss_i->second.genotipos.push_back(-1);
    }
  }

  // Actualizamos las poblaciones y demes con genotipo maximo
  ss.vector_iter_numDemesGanador.push_back(vector_numDemesGanador);
  ss.vector_iter_numDemes.push_back(vector_numDemes);
  ss.vector_iter_popSize.push_back(vector_popSize);
}

static SpatialOutput nr_inner_spatial(const fitnessEffectsAll &fitnessEffects,
                                      const double &initSize,
                                      const double &K,
                                      const TypeModel typeModel,
                                      const int &mutationPropGrowth,
                                      const std::vector<double> &mu,
                                      const double &death,
                                      const double &keepEvery,
                                      const double &sampleEvery,
                                      const std::vector<int> &initMutant,
                                      const time_t &start_time,
                                      const double &maxWallTime,
                                      const double &finalTime,
                                      const double &detectionSize,
                                      const int &detectionDrivers,
                                      const double &minDetectDrvCloneSz,
                                      const double &extraTime,
                                      const int &verbosity,
                                      double &totPopSize,
                                      double &em1,
                                      double &em1sc,
                                      double &ratioForce,
                                      double &currentTime,
                                      int &speciesFS,
                                      int &outNS_i,
                                      int &iter,
                                      std::vector<Genotype> &genot_out,
                                      std::vector<double> &popSizes_out,
                                      std::vector<int> &index_out,
                                      std::vector<double> &time_out,
                                      std::vector<double> &sampleTotPopSize,
                                      std::vector<double> &sampleLargestPopSize,
                                      std::vector<int> &sampleMaxNDr,
                                      std::vector<int> &sampleNDrLargestPop,
                                      bool &reachDetection,
                                      std::mt19937 &ran_gen,
                                      double &runningWallTime,
                                      bool &hittedWallTime,
                                      const std::map<int, std::string> &intName,
                                      const fitness_as_genes &genesInFitness,
                                      PhylogName &phylog,
                                      bool keepPhylog,
                                      const fitnessEffectsAll &muEF,
                                      const std::vector<int> &full2mutator,
                                      const double &cPDetect,
                                      const double &PDBaseline,
                                      const double &checkSizePEvery,
                                      const bool &AND_DrvProbExit,
                                      const std::vector<std::vector<int>> &fixation_l,
                                      const double &fixation_tolerance,
                                      const int &min_successive_fixation,
                                      const double &fixation_min_size,
                                      int &ti_dbl_min,
                                      int &ti_e3,
                                      std::map<std::string, std::string> &lod,
                                      POM &pom,
                                      const double maxFitness,
                                      const double initMutantPercentage,
                                      double cteSize,
                                      const double migrationProb,
                                      const double largeDistanceMigrationProb,
                                      const double maxMigrationPercentage,
                                      const int spatialDemesMax,
                                      const int spatialIterMax,
                                      const int spatialKeepEvery,
                                      const int spatialVerbosity,
                                      const bool spatialDeleteDemesWithoutMutations)
{

  /* Parametros obligatorios globales*/
  const int sp_per_period = 10; //Reservamos memoria para vector de hasta sp_per_period popParams y genotipos
  size_t i, j, k, numIter;
  const int numGenes = fitnessEffects.genomeSize;
  clock_t intra1, intra2, inter1, inter2;

  // Variables agregadas para el modelo espacial
  MapaDemes mapaDemes; //Mapa de los demes, identifica una posicion en el espacio con el deme que la ocupa
  MapaDemes::iterator deme_i;
  bool doneall = false;
  SpatialOutput so;
  SpatialSimulation ss;
  ss.keepEveryIters = 0;

  
  // En el modelo Moran si la cte no se indica se considera que es initSize.
  // De momento este if es redundante pues debe ser el typemodel = Constante para entrar al modelo espacial
  if (typeModel == TypeModel::constante)
  { 
    if (cteSize == -1)
    { // Por defecto se establece initSize como parametro de Moran
      cteSize = initSize;
    }
  }

  // Ininicalizamos el modelo espacial
  spatial_initialization(mapaDemes,
                         initSize,
                         K,
                         fitnessEffects,
                         initMutant,
                         initMutantPercentage,
                         ran_gen,
                         typeModel,
                         mutationPropGrowth,
                         mu,
                         death,
                         muEF,
                         full2mutator);

  //Mensajes para controlar la evolucion de la poblacion
  if (spatialVerbosity >= 3)
  {
    Rcpp::Rcout << "\n Comienzo: ";
    printMapaDemes(mapaDemes, fitnessEffects, spatialVerbosity);
  }
  numIter = 0;

  while (!doneall)
  {
    //__________ INTRADEME  _________
    intra1 = clock();
    for (deme_i = mapaDemes.begin(); deme_i != mapaDemes.end(); ++deme_i)
    {
      spatial_intrademe(
          fitnessEffects,
          initSize,
          K,
          typeModel,
          mutationPropGrowth,
          mu,
          death,
          keepEvery,
          sampleEvery,
          initMutant,
          start_time,
          maxWallTime,
          finalTime,
          detectionSize,
          detectionDrivers,
          minDetectDrvCloneSz,
          extraTime,
          verbosity,
          em1,
          em1sc,
          ratioForce,
          currentTime,
          speciesFS,
          outNS_i,
          iter,
          genot_out,
          popSizes_out,
          index_out,
          time_out,
          sampleTotPopSize,
          sampleLargestPopSize,
          sampleMaxNDr,
          sampleNDrLargestPop,
          reachDetection,
          ran_gen,
          runningWallTime,
          hittedWallTime,
          intName,
          genesInFitness,
          phylog,
          keepPhylog,
          muEF,
          full2mutator,
          cPDetect,
          PDBaseline,
          checkSizePEvery,
          AND_DrvProbExit,
          fixation_l,
          fixation_tolerance,
          min_successive_fixation,
          fixation_min_size,
          ti_dbl_min,
          ti_e3,
          lod,
          pom,
          maxFitness,
          deme_i->second.popParams,
          deme_i->second.Genotypes,
          deme_i->second.mapTimes,
          deme_i->second.numSpecies,
          deme_i->second.totPopSize,
          cteSize,
          deme_i->second);
    }
    intra2 = clock();

    //Eliminar demes vacíos y controlar que quedan demes. Los demes sin mutaciones solo son borrados despues
    // de las migraciones (por eso se pone false)
    spatial_DeleteEmptyDemes(mapaDemes, false, numGenes, spatialVerbosity);
    if (mapaDemes.size() == 0)
    {
      if (spatialVerbosity >= 1)
        Rcpp::Rcout << "\n FINAL: La población se ha extinguido en la iteracion: " << numIter;
      doneall = true;
      continue;
    }

    //Mensajes para controlar la evolucion de la poblacion
    if (spatialVerbosity >= 4)
    {
      Rcpp::Rcout << "\n Intrademe " << numIter;
      printMapaDemes(mapaDemes, fitnessEffects, spatialVerbosity);
    }
    //______________________________________________________

    // Guardado del estado de la simulacion
    if (spatialKeepEvery > 0)
    {
      if (!(numIter % spatialKeepEvery))
      {
        spatial_KeepSimulation(mapaDemes, ss);
        ss.keepEveryIters++;
      }
    }

    // Control de parada
    if (numIter >= spatialIterMax)
    {
      if (spatialVerbosity >= 1)
        Rcpp::Rcout << "\n FINAL: Numero máximo de iteraciones alcanzado: " << numIter;
      doneall = true;
      continue;
    }

    // Control de parada
    if (mapaDemes.size() >= spatialDemesMax)
    {
      //Con esta funcion eliminamos demes que tengan el genotipo vacio, por si así evitamos parar aquí
      // Es redundante el if, porque se hace dentro de la funcion pero es para recordar que esta funcion solo va a eliminar demes si el if es TRUE
      if (spatialDeleteDemesWithoutMutations)
      {
        spatial_DeleteEmptyDemes(mapaDemes, spatialDeleteDemesWithoutMutations, numGenes, spatialVerbosity);
      }
      // Si despues de eliminarlos sigue habiendo más demes de la cuenta, ya paramos
      if (mapaDemes.size() >= spatialDemesMax)
      {
        if (spatialVerbosity >= 1)
          Rcpp::Rcout << "\n FINAL: Numero máximo de demes alcanzado: " << mapaDemes.size();
        doneall = true;
        continue;
      }
    }

    //________________ INTERDEME ________________________
    // La funcion de interdeme hace todo internamente: migracion, migracion, comprueba demes vacios.
    // luego hace migracion a larga distancia y vuelve a comprobar demes vacios
    inter1 = clock();
    spatial_interdeme(mapaDemes,
                      ran_gen,
                      migrationProb,
                      largeDistanceMigrationProb,
                      maxMigrationPercentage,
                      spatialDeleteDemesWithoutMutations,
                      cteSize,
                      initSize,
                      numGenes,
                      fitnessEffects,
                      mutationPropGrowth,
                      muEF,
                      full2mutator,
                      typeModel,
                      K,
                      death,
                      mu,
                      spatialVerbosity);
    inter2 = clock();

    //Mensajes para controlar la evolucion de la poblacion
    if (spatialVerbosity >= 4)
    {
      Rcpp::Rcout << "\n Interdeme " << numIter;
      printMapaDemes(mapaDemes, fitnessEffects, spatialVerbosity);
    }
    //______________________________________________________

    //Mensajes para controlar la evolucion de la simulacion
    if (spatialVerbosity >= 2)
    {
      Rcpp::Rcout << "\n Iteracion " << numIter + 1 << " (" << mapaDemes.size() << " demes) --- intrademe: " << (intra2 - intra1) / (double)CLOCKS_PER_SEC << " ---- interdeme:" << (inter2 - inter1) / (double)CLOCKS_PER_SEC;
    }

    numIter++;
  } //Bucle de iteraciones

  // Se eliminan demes una vez más. Es redundante el if, porque se hace dentro de la funcion pero
  // es para recordar que esta funcion solo va a eliminar demes si el if es TRUE
  if (spatialDeleteDemesWithoutMutations)
  {
    spatial_DeleteEmptyDemes(mapaDemes, spatialDeleteDemesWithoutMutations, numGenes, spatialVerbosity);
  }

  // Guardamos por ultima vez el estado de la simulacion
  spatial_KeepSimulation(mapaDemes, ss);
  ss.keepEveryIters++;

  // Calculamos output del programa
  so.numIter = numIter;
  so.timeForIter = finalTime;
  so.keepEveryIters = ss.keepEveryIters;
  spatial_Output(so, mapaDemes, ss, fitnessEffects, spatialVerbosity);

  // Información final de la poblacion:
  if (spatialVerbosity >= 1)
  {
    Rcpp::Rcout << "\n\n Final Population: " << so.totPopSize << " (" << so.numDemes << " demes " << so.numIter << " iterations)";
    for (j = 0; j < so.vector_Genotypes.size(); j++)
    {
      Rcpp::Rcout << "\n Genotype: " << so.vector_Genotypes[j] << " --- Population: " << so.vector_popSize[j]
                  << " --- Demes: " << so.vector_numDemes[j] << " --- WinDemes: " << so.vector_numDemesGanador[j];
    }
  }

  if (spatialVerbosity >= 3)
  {
    printMapaDemes(mapaDemes, fitnessEffects, spatialVerbosity);
  }

  return so;
}

// [[Rcpp::export]]
Rcpp::List nr_BNB_Algo5(Rcpp::List rFE,
                        Rcpp::NumericVector mu_,
                        double death,
                        double initSize,
                        double sampleEvery,
                        double detectionSize,
                        double finalTime,
                        int initSp,
                        int initIt,
                        double seed,
                        int verbosity,
                        int speciesFS,
                        double ratioForce,
                        Rcpp::CharacterVector typeFitness_,
                        int maxram,
                        int mutationPropGrowth,
                        Rcpp::IntegerVector initMutant_,
                        double maxWallTime,
                        double keepEvery,
                        double K,
                        int detectionDrivers,
                        bool onlyCancer,
                        bool errorHitWallTime,
                        int maxNumTries,
                        bool errorHitMaxTries,
                        double minDetectDrvCloneSz,
                        double extraTime,
                        bool keepPhylog,
                        Rcpp::List MMUEF,
                        Rcpp::IntegerVector full2mutator_,
                        double n2,
                        double p2,
                        double PDBaseline,
                        double cPDetect_i,
                        double checkSizePEvery,
                        bool AND_DrvProbExit,
                        Rcpp::List fixation_i,
                        double maxFitness,
                        double initMutantPercentage,
                        double M,
                        double migrationProb,
                        double largeDistanceMigrationProb,
                        double maxMigrationPercentage,
                        bool spatialModel,
                        int spatialDemesMax,
                        int spatialIterMax,
                        int spatialKeepEvery,
                        int spatialVerbosity,
                        bool spatialDeleteDemesWithoutMutations){
  // double cPDetect){
  // double n2,
  // double p2,
  // double PDBaseline) {


  precissionLoss();
  const std::vector<double> mu = Rcpp::as<std::vector<double> >(mu_);
  const std::vector<int> initMutant = Rcpp::as<std::vector<int> >(initMutant_);
  const TypeModel typeModel = stringToModel(Rcpp::as<std::string>(typeFitness_));

  // A simple, vector-indexed way to map from numeric ids in full to
  // numeric ids in mutator. Recall all genes start with 1. So full2mutator[i-1];
  const std::vector<int> full2mutator = Rcpp::as<std::vector<int> >(full2mutator_);
  // A consistency check

  // const double genTime = 4.0; // should be a parameter. For Bozic only.

  //If seed is -9, then use automatic seed.


  // Code when using randutils
  // randutils::mt19937_rng ran_gen;
  // if(seed == 0)
  //   ran_gen.seed();
  // else {
  //   ran_gen.seed(static_cast<unsigned int>(seed));
  //   // The next does not solve the differences between clang and gcc. So
  //   // keep it simple.
  //   // std::seed_seq s1{static_cast<unsigned int>(seed)};
  //   // ran_gen.seed(s1);
  // }

  unsigned int rseed = static_cast<unsigned int>(seed);
  if(seed == 0) {
    rseed = std::random_device{}();
  }
  std::mt19937 ran_gen(rseed);

  double cPDetect = cPDetect_i;
  if( (n2 > 0) && (p2 > 0) ) {
    if (PDBaseline <= 0) throw std::range_error("PDBaseline <= 0");
    cPDetect = set_cPDetect(n2, p2, PDBaseline);
    if(verbosity >= 1)
      Rcpp::Rcout << "  cPDetect set at " << cPDetect << "\n";
  }

  if( (K < 1 ) && ( (typeModel ==   TypeModel::mcfarlandlog) ||
		    (typeModel ==   TypeModel::mcfarlandlog) ))
    throw std::range_error("K < 1.");
  
  fitnessEffectsAll fitnessEffects =  convertFitnessEffects(rFE);
  //Used at least twice
  std::map<int, std::string> intName = mapGenesIntToNames(fitnessEffects);
  fitness_as_genes genesInFitness = fitnessAsGenes(fitnessEffects);
  PhylogName phylog;
  // LOD lod;
  std::map<std::string, std::string> lod;
  POM pom;

  // Mutator effects
	fitnessEffectsAll muEF;
  if( (full2mutator.size() != 0) ){
		muEF = convertFitnessEffects(MMUEF);
	} else {
		muEF = nullFitnessEffects();
		muEF.frequencyDependentFitness = fitnessEffects.frequencyDependentFitness;
			}
  // Paranoia. We should never end up here.
  if( (full2mutator.size() != 0) && (muEF.genomeSize == 0))
    throw std::logic_error("full2mutator > 0 with mutatorEffects.genomesize 0");
  if( (full2mutator.size() == 0) && (muEF.genomeSize != 0)) {
    throw std::logic_error("full2mutator 0 with mutatorEffects.genomesize != 0");
  }

  // fixation: run until some genotype combinations fixed

  double fixation_tolerance = -9;
  int min_successive_fixation = 100;
  double fixation_min_size = 0.0;
  std::vector < std::vector<int> > fixation_l;

  if( fixation_i.size() != 0 ) {
    Rcpp::List fggl = fixation_i["fixation_list"] ;
    fixation_l = list_to_vector_of_int_vectors(fggl); // FIXME
    fixation_tolerance = Rcpp::as<double>(fixation_i["fixation_tolerance"]);
    min_successive_fixation = Rcpp::as<int>(fixation_i["min_successive_fixation"]);
    fixation_min_size = Rcpp::as<double>(fixation_i["fixation_min_size"]);
  } else {
    fixation_l.resize(0); // explicit
  }


  bool runAgain = true;
  bool reachDetection = false;
  //Output
  std::vector<Genotype> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out;
  std::vector<double> time_out; //only one entry per period!
  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  double totPopSize = 0;
  std::vector<double> sampleTotPopSize;
  std::vector<double> sampleLargestPopSize;
  std::vector<int> sampleMaxNDr; //The largest number of drivers in any
				 //genotype or clone at each  time sample
  std::vector<int> sampleNDrLargestPop; //Number of drivers in the clone
					// with largest size (at each time
					// sample)
  sampleTotPopSize.reserve(initIt);
  sampleLargestPopSize.reserve(initIt);
  sampleMaxNDr.reserve(initIt);
  sampleNDrLargestPop.reserve(initIt);

  int outNS_i = -1; // the column in the outNS
  time_t start_time = time(NULL);
  double runningWallTime = 0;
  bool  hittedWallTime = false;
  bool hittedMaxTries = false;

  // spParamsP tmpParam;
  // std::vector<spParamsP> popParams(1);
  // const int sp_per_period = 5000;

  // popParams.reserve(sp_per_period);
  // Genotypes.reserve(sp_per_period);

  // std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset


  // // multimap to hold nextMutationTime
  // std::multimap<double, int> mapTimes;
  // //std::multimap<double, int>::iterator m1pos;


  // // count troublesome tis
  // int ti_dbl_min = 0;
  // int ti_e3 = 0;



  // // Beerenwinkel
  // double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  // //McFarland
  // double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  // double e1, n_0; //n_1; // for McFarland error
  // double tps_0, tps_1; // for McFarland error
  // tps_0 = 0.0;
  // tps_1 = 0.0;
  // e1 = 0.0;
  // n_0 = 0.0;
  // n_1 = 0.0;

  double em1, em1sc; // new computation of McFarland error
  em1 = 0.0;
  em1sc = 0.0;


  // // For totPopSize_and_fill and bailing out
  // // should be static vars inside funct,
  // // but they keep value over calls in same R session.
  // int lastMaxDr = 0;
  // double done_at = -9;
  // // totalPopSize at time t, at t-1 and the max error.

  // 5.1 Initialize

  int numRuns = 0;
  int numRecoverExcept = 0;
  bool forceRerun = false;

  double currentTime = 0;
  int iter = 0;

  int ti_dbl_min = 0;
  int ti_e3 = 0;

  int accum_ti_dbl_min = 0;
  int accum_ti_e3 = 0;
  // bool AND_DrvProbExit = ( (cpDetect >= 0) &&
  // 			     (detectionDrivers < 1e9) &&
  // 			     (detectionSize < std::numeric_limits<double>::infinity()));

  //Variables modelo espacial
  SpatialOutput so;
  if (!spatialModel)
  {

    while (runAgain)
    {

      if (numRuns >= maxNumTries)
      {
        //  hittedMaxTries This we want here to avoid an extra run and
        //  confusing output
        hittedMaxTries = true;
        Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
        runAgain = false;
        if (errorHitMaxTries)
        {
          Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
          return List::create(Named("HittedWallTime") = false,
                              Named("HittedMaxTries") = true,
                              Named("other") =
                                  List::create(Named("UnrecoverExcept") = false));
        }
        break;
      }

      try
      {
        Rcpp::checkUserInterrupt();

        // it is CRUCIAL that several entries are zeroed (or -1) at the
        // start of innerBNB now that we do multiple runs if onlyCancer = true.

        nr_innerBNB(
            fitnessEffects,
            initSize,
            K,
            // alpha,
            // genTime,
            typeModel,
            mutationPropGrowth,
            mu,
            death,
            keepEvery,
            sampleEvery,
            initMutant,
            start_time,
            maxWallTime,
            finalTime,
            detectionSize,
            detectionDrivers,
            minDetectDrvCloneSz,
            extraTime,
            verbosity,
            totPopSize,
            em1,
            em1sc,
            // n_0,
            // 	  // n_1,
            // 	  en1,
            ratioForce,
            currentTime,
            speciesFS,
            outNS_i,
            iter,
            genot_out,
            popSizes_out,
            index_out,
            time_out,
            sampleTotPopSize,
            sampleLargestPopSize,
            sampleMaxNDr,
            sampleNDrLargestPop,
            reachDetection,
            ran_gen,
            runningWallTime,
            hittedWallTime,
            intName,
            genesInFitness,
            phylog,
            keepPhylog,
            muEF,
            full2mutator,
            cPDetect,
            PDBaseline,
            checkSizePEvery,
            AND_DrvProbExit,
            fixation_l,
            fixation_tolerance,
            min_successive_fixation,
            fixation_min_size,
            ti_dbl_min,
            ti_e3,
            lod,
            pom);
        ++numRuns;
        forceRerun = false;
        accum_ti_dbl_min += ti_dbl_min;
        accum_ti_e3 += ti_e3;
      } catch (rerunExcept &e) {
      Rcpp::Rcout << "\n Recoverable exception " << e.what()
		  << ". Rerunning.";
      forceRerun = true;
      ++numRecoverExcept;
      ++numRuns; // exception should count here!
      accum_ti_dbl_min += ti_dbl_min;
      accum_ti_e3 += ti_e3;
    } catch (const std::exception &e) {
      Rcpp::Rcout << "\n Unrecoverable exception: " << e.what()
		  << ". Aborting. \n";
      return
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = e.what()));
    } catch (...) {
      Rcpp::Rcout << "\n Unknown unrecoverable exception. Aborting."
		  << "(User interrupts also generate this).\n";
      return
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = "Unknown exception"));
    }
    if(hittedWallTime) {
      Rcpp::Rcout << "\n Hitted wall time. Exiting.";
      runAgain = false;
      if(errorHitWallTime) {
	Rcpp::Rcout << "\n Hitting wall time is regarded as an error. \n";
	return
	  List::create(Named("HittedWallTime") = true,
		       Named("HittedMaxTries") = false, // yes, for
							// coherent return
							// objects
		       Named("other") =
		       List::create(Named("UnrecoverExcept") = false));
      }
    // } else if(numRuns > maxNumTries) {
    //   //  hittedMaxTries FIXME this is very, very confusing in limit
    //   // cases.  suppose maxNumTries = 1. We will run two times, and the
    //   // second might have reached cancer, but we will bail out here, as
    //   // numRuns is actually 2. However, we report the value. And we run
    //   // once more than needed.
    //   hittedMaxTries = true;
    //   Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
    //   runAgain = false;
    //   if(errorHitMaxTries) {
    // 	Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
    // 	return
    // 	  List::create(Named("HittedWallTime") = false,
    // 		       Named("HittedMaxTries") = true,
    // 		       Named("other") =
    // 		       List::create(Named("UnrecoverExcept") = false));
    //   }
    } else if(forceRerun) {
      runAgain = true;
      forceRerun = false;
    } else {
      if(onlyCancer) {
	runAgain = !reachDetection;
      } else {
	runAgain = false;
      }
    }
#ifdef DEBUGV
      Rcpp::Rcout << "\n reachDetection = " << reachDetection;
      Rcpp::Rcout << "\n forceRerun =  " << forceRerun  << "\n";

#endif

  } // runAgain loop


  std::vector<std::vector<int> > genot_out_v = genot_to_vectorg(genot_out);
  std::vector<std::vector<int> > uniqueGenotypes_vector_nr  =
    uniqueGenot_vector(genot_out_v);
  IntegerMatrix returnGenotypes =
    nr_create_returnGenotypes(fitnessEffects.genomeSize,
  			      uniqueGenotypes_vector_nr);
  Rcpp::NumericMatrix outNS = create_outNS(uniqueGenotypes_vector_nr,
  					   genot_out_v,
  					   popSizes_out,
  					   index_out, time_out,
  					   outNS_i, maxram);

  int maxNumDrivers = 0;
  int totalPresentDrivers = 0;
  std::vector<int> countByDriver(fitnessEffects.drv.size(), 0);
  std::vector<int> presentDrivers;
  driverCounts(maxNumDrivers, totalPresentDrivers,
	       countByDriver, presentDrivers,
	       returnGenotypes, fitnessEffects.drv);


  std::vector<std::string> genotypesAsStrings =
    genotypesToNameString(uniqueGenotypes_vector_nr, genesInFitness, intName);
  std::string driversAsString =
    driversToNameString(presentDrivers, intName);

  // // // zz: debugging
  // // // Correct too
  // DP1("intName");
  // for(auto mmm: intName) {
  //   Rcpp::Rcout << mmm.first << " :" ;
  //   Rcpp::Rcout << mmm.second << std::endl;
  // }


  // // wrong
  // DP1("genotypesAsStrings");
  // for(auto gas: genotypesAsStrings) {
  //   Rcpp::Rcout << gas;
  //   Rcpp::Rcout << std::endl;
  // }


  std::vector<double> sampleLargestPopProp(outNS_i + 1);
  if((outNS_i + 1) != static_cast<int>(sampleLargestPopSize.size()))
    throw std::length_error("outNS_i + 1 != sampleLargestPopSize.size");
  std::transform(sampleLargestPopSize.begin(), sampleLargestPopSize.end(),
  		 sampleTotPopSize.begin(),
  		 sampleLargestPopProp.begin(),
  		 std::divides<double>());
  NumericMatrix perSampleStats(outNS_i + 1, 5);
  fill_SStats(perSampleStats, sampleTotPopSize, sampleLargestPopSize,
  	      sampleLargestPopProp, sampleMaxNDr, sampleNDrLargestPop);

  // create the lod return pieces. Move to a function later
  std::vector<std::string> lod_parent;
  std::vector<std::string> lod_child;
  for (const auto &l : lod) {
    lod_child.push_back(l.first);
    lod_parent.push_back(l.second);
  }

  return
    List::create(Named("pops.by.time") = outNS,
		 Named("NumClones") = uniqueGenotypes_vector_nr.size(),
		 Named("TotalPopSize") = totPopSize,
		 Named("Genotypes") = returnGenotypes,
		 Named("GenotypesWDistinctOrderEff") = Rcpp::wrap(uniqueGenotypes_vector_nr),
		 Named("GenotypesLabels") = Rcpp::wrap(genotypesAsStrings),
		 Named("MaxNumDrivers") = maxNumDrivers,
		 Named("MaxDriversLast") = sampleMaxNDr[outNS_i],
		 Named("NumDriversLargestPop") =  sampleNDrLargestPop[outNS_i],
		 Named("LargestClone") = sampleLargestPopSize[outNS_i],
		 Named("PropLargestPopLast") = sampleLargestPopProp[outNS_i],
		 Named("FinalTime") = currentTime,
		 Named("NumIter") = iter,
		 Named("HittedWallTime") = hittedWallTime,
		 Named("HittedMaxTries") = hittedMaxTries,
		 Named("TotalPresentDrivers") = totalPresentDrivers,
		 Named("CountByDriver") = countByDriver,
		 Named("OccurringDrivers") = driversAsString,
		 Named("PerSampleStats") = perSampleStats,
		 Named("other") = List::create(Named("attemptsUsed") = numRuns,
					       Named("errorMF") =
					       returnMFE_new(em1sc, typeModel),
					       Named("errorMF_size") =
					       returnMFE_new(em1, typeModel), // Used to be e1, not log
					       // Named("errorMF_n_0") = n_0,
#ifdef MIN_RATIO_MUTS_NR
					       Named("minDMratio") =
					       g_min_death_mut_ratio_nr,
					       Named("minBMratio") =
					       g_min_birth_mut_ratio_nr,
#else
					       Named("minDMratio") = -99,
					       Named("minBMratio") = -99,
#endif
					       //    Named("errorMF_n_1") = n_1,
					       Named("PhylogDF") =  DataFrame::create(
										      Named("parent") = phylog.parent,
										      Named("child") = phylog.child,
										      Named("time") = phylog.time,
										      Named("pop_size_child") = phylog.pop_size_child
										      ),
					       Named("UnrecoverExcept") = false,
					       Named("numRecoverExcept") = numRecoverExcept,
					       Named("accum_ti_dbl_min") = accum_ti_dbl_min,
					       Named("accum_ti_e3") = accum_ti_e3,
					       Named("LOD_DF") = DataFrame::create(
										   Named("parent") = lod_parent, // lod.parent,
										   Named("child") = lod_child //lod.child
										   ),
					       Named("POM") = Rcpp::wrap(pom.genotypesString)
					       )
		 );
  }
  else { //Spatial model
    try {
      Rcpp::checkUserInterrupt();

      so = nr_inner_spatial(
          fitnessEffects,
          initSize,
          K,
          typeModel,
          mutationPropGrowth,
          mu,
          death,
          keepEvery,
          sampleEvery,
          initMutant,
          start_time,
          maxWallTime,
          finalTime,
          detectionSize,
          detectionDrivers,
          minDetectDrvCloneSz,
          extraTime,
          verbosity,
          totPopSize,
          em1,
          em1sc,
          ratioForce,
          currentTime,
          speciesFS,
          outNS_i,
          iter,
          genot_out,
          popSizes_out,
          index_out,
          time_out,
          sampleTotPopSize,
          sampleLargestPopSize,
          sampleMaxNDr,
          sampleNDrLargestPop,
          reachDetection,
          ran_gen,
          runningWallTime,
          hittedWallTime,
          intName,
          genesInFitness,
          phylog,
          keepPhylog,
          muEF,
          full2mutator,
          cPDetect,
          PDBaseline,
          checkSizePEvery,
          AND_DrvProbExit,
          fixation_l,
          fixation_tolerance,
          min_successive_fixation,
          fixation_min_size,
          ti_dbl_min,
          ti_e3,
          lod,
          pom,
          maxFitness,
          initMutantPercentage,
          M,
          migrationProb,
          largeDistanceMigrationProb,
          maxMigrationPercentage,
          spatialDemesMax,
          spatialIterMax,
          spatialKeepEvery,
          spatialVerbosity,
          spatialDeleteDemesWithoutMutations);
    }
    catch (rerunExcept &e)
    {
      Rcpp::Rcout << "\n Recoverable exception " << e.what()
                  << ". Aborting. \n";
      forceRerun = false;
    }
    catch (const std::exception &e)
    {
      Rcpp::Rcout << "\n Unrecoverable exception: " << e.what()
                  << ". Aborting. \n";
      return List::create(Named("other") =
                              List::create(Named("UnrecoverExcept") = true,
                                           Named("ExceptionMessage") = e.what()));
    }
    catch (...)
    {
      Rcpp::Rcout << "\n Unknown unrecoverable exception. Aborting."
                  << "(User interrupts also generate this).\n";
      return List::create(Named("other") =
                              List::create(Named("UnrecoverExcept") = true,
                                           Named("ExceptionMessage") = "Unknown exception"));
    }

    return List::create(Named("SpatialSummary") = DataFrame::create(
                            Named("Genotypes") = so.vector_Genotypes,
                            Named("Population") = so.vector_popSize,
                            Named("Demes") = so.vector_numDemes,
                            Named("WinDemes") = so.vector_numDemesGanador),
                        Named("TotalPopulation") = so.totPopSize,
                        Named("NumIter") = so.numIter,
                        Named("NumDemes") = so.numDemes,
                        Named("TimeForIter") = so.timeForIter,
                        Named("Simulation") = List::create(
                            Named("Genotypes") = so.vector_Genotypes,
                            Named("KeepEveryIters") = so.keepEveryIters,
                            Named("WinDemesPerIter") = so.vector_iter_numDemesGanador,
                            Named("PopSizePerIter") = so.vector_iter_popSize,
                            Named("NumDemesPerIter") = so.vector_iter_numDemes));
  }
}

  // Creating return object:

  // The 0, 1 representation is how most of the work is done in R: do I want
  // to change that?

  // Order: beware of two things: order is important for the "true"
  // genotypes, but is not immediately observable. So for 0,1
  // representation, not needed or used. Thus, maybe I want two
  // representations.

  // Yes, the full Genotye structure is only used when assigning fitness. So
  // could we use a collapsed one with: order + rest? Nope, as whenever I'd
  // create a child from a genotype, I'd need like deconvolve, and go back
  // to the three piece structure. This seems much more expensive than the
  // overloaded == and the usage of the overloaded < (this is only used at
  // the end, when producing the output objects)

  // When you need to alter things that affect a genotype periodically, like
  // the FDF, or the McFarland rates, etc, you will want to do it here:
  // when we call updateRatesMcFarlandLog and similar, e.g.,
  // In BNB_nr.cpp, in nr_innerBNB function
  // when we are sampling.
  // Probably also where nr_fitness is called
  // And generally, in general, were we pass currentTime
