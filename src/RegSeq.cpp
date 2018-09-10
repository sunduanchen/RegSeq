// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadilloExtensions/sample.h>
#include<Rcpp.h>
#include<vector>
#include<string>
#include<cstring>
#include<math.h>
#include<Rmath.h>
#include<Rcpp/stats/random/random.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
RcppExport SEXP TwotailedGSEA(SEXP regulon, SEXP tvalue, SEXP min, SEXP max, SEXP nPerm){

    // Original R-input preprocessing.
    Rcpp::List REGULON(regulon);
    vector<string> TF_name = REGULON.names();
    Rcpp::List TVALUE(tvalue);
    unsigned int MIN = as<int>(min);
    unsigned int MAX = as<int>(max);
    int NPERM = as<int>(nPerm) + 1;

    // Count the number of gene sets satisfying user options: MIN, MAX
    int nGeneset = 0;
    int i, k, PERM;
    for(i=0; i<REGULON.size(); i++){
        Rcpp::List REGULON_tmp(REGULON[i]);
        NumericVector target = as<NumericVector>(REGULON_tmp["tfmode"]);
        if(target.size() >= MIN && target.size() <= MAX){
            nGeneset++;
        }
    }
    if(nGeneset == 0){
        Rcpp::warning("No TF passed the condition for TF-target size. Please check the option min and max for target gene size.");
    }

    NumericMatrix enrich(nGeneset, NPERM);
    vector<int> OUTPUT_SIZE;
    int gs = 0;
    vector<string> GenesetName;
    map<string, int> gene2index;
    vector<int> rankGene;
    double N_R, P_hit, P_miss;
    int left, right, mean;
    bool replace = false;
    NumericVector prob = NumericVector::create();
    vector<int> rankGene2;
    for(i=0; i<REGULON.size(); i++){
        Rcpp::List REGULON_tmp(REGULON[i]);
        NumericVector target = as<NumericVector>(REGULON_tmp["tfmode"]);        // Positive and negative score for each target gene in i-th TF.
        vector<string> target_name = target.names();                            // Positive and negative targets.

        if(target.size() < MIN || target.size() > MAX)
            continue;
        GenesetName.push_back(TF_name[i]);                                      // TFs satisfy the min/max size limits.
        OUTPUT_SIZE.push_back(target.size());

        // Map the gene symbol to integer.
        gene2index.clear();
        NumericVector signature = as<NumericVector>(TVALUE[i]);                 // Sorted gene signature by reversing target genes.
        vector<string> gene_symbol = signature.names();
        int gene_size = signature.size();
        for(k=0; k<gene_size; k++){
            gene2index[gene_symbol[k]] = k;                                     // Store the mapping information of gene symbol and rank.
        }

        // Set rank: The rank of target genes in the sorted signature vector.
        rankGene.clear();
        N_R = 0.0;
        for (k=0; k<target.size(); k++){
            rankGene.push_back(gene2index[target_name[k]]);
            N_R += fabs(signature[rankGene[k]]);                                // Get the sums of all target gene t-value score.
        }
        sort(rankGene.begin(), rankGene.end());                                 // sort the rankGene in ascending order.

        // Enrichment Score calculation
        enrich(gs,0) = P_hit = 0. ;
        for(k=0; k<target.size(); k++){
            P_miss = 1./(gene_size - target.size()) * (rankGene[k]-k);
            if(fabs(enrich(gs,0)) < fabs(P_hit-P_miss)){
                enrich(gs,0) = P_hit-P_miss;
            }
            P_hit += fabs(signature[rankGene[k]])/N_R;
            if(fabs(enrich(gs,0)) < fabs(P_hit - P_miss)){
                enrich(gs,0) = P_hit-P_miss;
            }
        }

        // Do permutation
        for(PERM = 1; PERM<NPERM; PERM++){
            NumericVector permutation_score = clone(signature);
            permutation_score = RcppArmadillo::sample(permutation_score, gene_size, replace, prob);
            NumericVector rank(gene_size);
            for(k=0; k<gene_size; k++){
                left = 0; right = gene_size-1;
                while(left <= right){
                    mean = (left + right) / 2 ;
                    if(permutation_score[k] < signature[mean]){
                        left = mean + 1;
                    }else if(permutation_score[k] > signature[mean]){
                        right = mean - 1;
                    }else{
                        break;
                    }
                }
                rank[k] = mean;                                                 // Obtain the updated rank information of the permutated value.
            }                                                                   // The rank of each permutated value in the original sorted signature vector.
            N_R = 0.0;
            rankGene2.clear();
            for(k=0; k<target.size(); k++){
                rankGene2.push_back(rank(rankGene[k]));                         // New ranks of all target genes.
                N_R += fabs(permutation_score[rankGene[k]]);                    // Same target gene index, but the values are different by permutation.
            }
            sort(rankGene2.begin(), rankGene2.end());
            enrich(gs,PERM) = P_hit = 0.0;
            for(k=0; k<target.size(); k++){
                P_miss = 1./(gene_size - target.size()) * (rankGene2[k] - k);   // Using the updated rank information. Not the previous one.
                if(fabs(enrich(gs,PERM)) < fabs(P_hit - P_miss)){
                    enrich(gs,PERM) = P_hit - P_miss;
                }
                P_hit += fabs(signature[rankGene2[k]])/N_R;
                if(fabs(enrich(gs,PERM)) < fabs(P_hit - P_miss)){
                    enrich(gs,PERM) = P_hit - P_miss;
                }
            }
        }
        gs++;
    }

    // Normalization
    int nPos, nNeg;
    double sumX, sumY;
    for(gs=0; gs<nGeneset; gs++){
        nPos = 0; sumX = 0.;
        nNeg = 0; sumY = 0.;
        for(PERM=1; PERM<NPERM; PERM++){
            if(enrich(gs,PERM) >= 0){
                ++nPos;
                sumX += enrich(gs,PERM);
            }else{
                ++nNeg ;
                sumY += enrich(gs,PERM);
            }
        }
        if(nPos > 0){
            sumX /= nPos;
            for(PERM=0; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) >= 0){
                    enrich(gs,PERM) /= sumX;
                }
            }
        }
        if(nNeg > 0){
            sumY = fabs(sumY/nNeg);
            for(PERM=0; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) < 0){
                    enrich(gs,PERM) /= sumY;
                }
            }
        }
    }

    // Estimating significance: normal p value.
    double nominalP;
    string updown;
    vector<string> OUTPUT_GENESET, DIRECTION;
    vector<double> OUTPUT_NES, OUTPUT_NominalP;
    for(gs=0; gs<nGeneset; gs++){
        nominalP = 0.;
        nPos = nNeg = 0;
        if(enrich(gs,0) >= 0){
            for(PERM=1; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) >= 0){
                    ++nPos;
                    if(enrich(gs,PERM) >= enrich(gs,0)){
                        ++nominalP;
                    }
                }
            }
        }else{
            for(PERM=1; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) < 0){
                    ++nNeg;
                    if(enrich(gs,PERM) <= enrich(gs,0)){
                        ++nominalP;
                    }
                }
            }
        }
        if(enrich(gs,0) >= 0){
            nominalP /= nPos;
            nominalP = Rf_fprec(nominalP, 4);
            updown = "ACTIVATED";
        }else{
            nominalP /= nNeg;
            nominalP = Rf_fprec(nominalP, 4);
            updown = "REPRESSED";
        }

        // Result vectors
        OUTPUT_GENESET.push_back(GenesetName[gs]);
        DIRECTION.push_back(updown);
        OUTPUT_NES.push_back(enrich[gs]);
        OUTPUT_NominalP.push_back(nominalP);
    }
    DataFrame df = DataFrame::create(Named("GenesetName", OUTPUT_GENESET),
                                     Named("Size", OUTPUT_SIZE),
                                     Named("NES", OUTPUT_NES),
                                     Named("Nominal P-value", OUTPUT_NominalP),
                                     Named("Direction", DIRECTION)
    );
    return wrap(df);
}


// [[Rcpp::export]]
RcppExport SEXP OnetailedGSEA(SEXP regulon, SEXP tvalue, SEXP min, SEXP max, SEXP nPerm){

    // Original R-input preprocessing.
    Rcpp::List REGULON(regulon);
    vector<string> TF_name = REGULON.names();
    NumericVector TVALUE = as<NumericVector>(tvalue);
    unsigned int MIN = as<int>(min);
    unsigned int MAX = as<int>(max);
    int NPERM = as<int>(nPerm) + 1;

    // Count the number of gene sets satisfying user options: MIN, MAX
    int nGeneset = 0;
    int i, k, PERM;
    for(i=0; i<REGULON.size(); i++){
        Rcpp::List REGULON_tmp(REGULON[i]);
        NumericVector target = as<NumericVector>(REGULON_tmp["tfmode"]);
        if(target.size() >= MIN && target.size() <= MAX){
            nGeneset++;
        }
    }
    if(nGeneset == 0){
        Rcpp::warning("No TF passed the condition for TF-target size. Please check the option min and max for target gene size.");
    }

    // Store the rank of transcription factor target gene member.
    map<string, int> gene2index;
    std::vector<vector<int> > rankGene;
    std::vector<vector<int> > rankGene2;
    vector<string> GenesetName;
    vector<int> OUTPUT_SIZE;
    for(i=0; i<nGeneset; i++){
        rankGene.push_back(vector<int>());
        rankGene2.push_back(vector<int>());
    }
    int gene_size = TVALUE.size();
    vector<string> gene_symbol = TVALUE.names();
    for(k=0; k<gene_size; k++){
        gene2index[gene_symbol[k]] = k;                                         // Map the gene symbol to integer.
    }
    int gs = 0;
    for(i=0; i<REGULON.size(); i++){
        Rcpp::List REGULON_tmp(REGULON[i]);
        NumericVector target = as<NumericVector>(REGULON_tmp["tfmode"]);
        vector<string> target_name = target.names();
        if(target.size() < MIN || target.size() > MAX){
            continue;
        }
        GenesetName.push_back(TF_name[i]);
        OUTPUT_SIZE.push_back(target.size());
        for (k=0; k<target.size(); k++){
            rankGene[gs].push_back(gene2index[target_name[k]]);                 // Set rank: The rank of target genes in the sorted signature vector.
        }
        gs++;
    }

    // Enrichment Score calculation.
    double N_R, P_hit, P_miss;
    NumericMatrix enrich(nGeneset, NPERM);
    for(gs=0; gs<nGeneset; gs++){
        N_R = 0.0;
        for(k=0; k<rankGene[gs].size(); k++){
            N_R += fabs(TVALUE[rankGene[gs][k]]);
        }
        sort(rankGene[gs].begin(), rankGene[gs].end());
        enrich(gs,0) = P_hit = 0.;
        for(k=0; k<rankGene[gs].size(); k++){
            P_miss = 1./(gene_size - rankGene[gs].size()) * (rankGene[gs][k]-k);
            if(enrich(gs,0) < (P_hit-P_miss)){
                enrich(gs,0) = P_hit-P_miss;
            }
            P_hit += fabs(TVALUE[rankGene[gs][k]])/N_R;
            if(enrich(gs,0) < (P_hit - P_miss)){
                enrich(gs,0) = P_hit-P_miss;
            }
        }
    }

    // Do permutation
    int left, right, mean;
    bool replace = false;
    NumericVector prob = NumericVector::create();
    for(PERM=1; PERM<NPERM; PERM++){
        NumericVector TVALUE_temp = clone(TVALUE);
        TVALUE_temp = RcppArmadillo::sample(TVALUE_temp, gene_size, replace, prob);
        NumericVector rank(gene_size);
        for(k=0; k<gene_size; k++){
            left = 0; right = gene_size-1;
            while(left <= right){
                mean = (left + right) / 2 ;
                if(TVALUE_temp[k] < TVALUE[mean]){
                    left = mean + 1;
                }else if(TVALUE_temp[k] > TVALUE[mean]){
                    right = mean - 1;
                }else{
                    break;
                }
            }
            rank[k] = mean;
        }
        for(gs=0; gs<nGeneset; gs++){
            N_R = 0.0;
            rankGene2[gs].clear();
            for(k=0; k<rankGene[gs].size(); k++){
                rankGene2[gs].push_back(rank(rankGene[gs][k]));
                N_R += fabs(TVALUE_temp[rankGene[gs][k]]);
            }
            sort(rankGene2[gs].begin(), rankGene2[gs].end());
            enrich(gs,PERM) = P_hit = 0.0;
            for(k=0; k<rankGene2[gs].size(); k++){
                P_miss = 1./(gene_size - rankGene2[gs].size()) * (rankGene2[gs][k] - k);
                if(enrich(gs,PERM) < (P_hit - P_miss)){
                    enrich(gs,PERM) = P_hit-P_miss;
                }
                P_hit += fabs(TVALUE[rankGene2[gs][k]])/N_R;
                if(enrich(gs,PERM) < (P_hit - P_miss)){
                    enrich(gs,PERM) = P_hit - P_miss;
                }
            }
        }
    }

    // Normalization
    int orgPos = 0, orgNeg = 0, allPos = 0, allNeg = 0;
    int nPos, nNeg;
    double sumX, sumY;
    for(gs=0; gs<nGeneset; gs++){
        if(enrich(gs,0) >= 0){
            ++orgPos;
        }
        for(PERM=1; PERM<NPERM; PERM++){
            if (enrich(gs,PERM) >= 0){
                ++allPos;
            }
        }
    }
    orgNeg = nGeneset - orgPos;
    allNeg = nGeneset*(NPERM-1) - allPos;

    for(gs=0; gs<nGeneset; gs++){
        nPos = 0; sumX = 0.;
        nNeg = 0; sumY = 0.;
        for(PERM = 1; PERM<NPERM; PERM++){
            if(enrich(gs,PERM) >=0){
                ++nPos;
                sumX += enrich(gs,PERM);
            }else{
                ++nNeg ;
                sumY += enrich(gs,PERM);
            }
        }
        if(nPos > 0){
            sumX /= nPos;
            for(PERM = 0; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) >= 0){
                    enrich(gs,PERM) /= sumX;
                }
            }
        }
        if(nNeg > 0){
            sumY = fabs(sumY/nNeg);
            for(PERM=0; PERM<NPERM; PERM++){
                if(enrich(gs,PERM)<0){
                    enrich(gs,PERM) /= sumY;
                }
            }
        }
    }

    // Estimating significance: normal p value, FDR
    double FDRmo, FDRja, FDR;
    double nominalP;
    vector<string> OUTPUT_GENESET;
    vector<double> OUTPUT_NES, OUTPUT_NominalP, OUTPUT_FDR;
    for(gs=0; gs<nGeneset; gs++){
        nominalP = 0.;
        nPos = nNeg = 0;
        FDRmo = FDRja = 0.;
        if(enrich(gs,0) >= 0){
            for(PERM=1; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) >= 0){
                    ++nPos;
                    if(enrich(gs,PERM) >= enrich(gs,0)){
                        ++nominalP;
                    }
                }
            }
            for(k=0; k<nGeneset; k++){
                if(enrich(gs,0) <= enrich(k,0)){
                    ++FDRmo;
                }
            }
            FDRmo /= orgPos;
            for(k=0; k<nGeneset; k++){
                for(PERM=1; PERM<NPERM; PERM++){
                    if(enrich(gs,0) <= enrich(k,PERM)){
                        ++FDRja;
                    }
                }
            }
            FDRja /= allPos;
        }else{
            for(PERM=1; PERM<NPERM; PERM++){
                if(enrich(gs,PERM) < 0){
                    ++nNeg;
                    if(enrich(gs,PERM) <= enrich(gs,0)){
                        ++nominalP;
                    }
                }
            }
            for(k=0; k<nGeneset; k++){
                if(enrich(gs,0) >= enrich(k,0)){
                    ++FDRmo;
                }
            }
            FDRmo /= orgNeg;
            for(k=0; k<nGeneset; k++){
                for(PERM=1; PERM<NPERM; PERM++){
                    if(enrich(gs,0) >= enrich(k,PERM)){
                        ++FDRja;
                    }
                }
            }
            FDRja /= allNeg;
        }
        if(enrich(gs,0) >= 0){
            nominalP /= nPos;
            nominalP = Rf_fprec(nominalP, 4);
        }else{
            nominalP /= nNeg;
            nominalP = Rf_fprec(nominalP, 4);
        }
        FDR = FDRja/FDRmo;
        FDR = Rf_fprec(FDR, 4);

        OUTPUT_GENESET.push_back(GenesetName[gs]);
        OUTPUT_NES.push_back(enrich[gs]);
        OUTPUT_NominalP.push_back(nominalP);
        OUTPUT_FDR.push_back(FDR);
    }
    DataFrame df = DataFrame::create(Named("GenesetName", OUTPUT_GENESET),
                                     Named("Size", OUTPUT_SIZE),
                                     Named("NES", OUTPUT_NES),
                                     Named("Nominal P-value", OUTPUT_NominalP),
                                     Named("FDR Q-value", OUTPUT_FDR)
    );
    return wrap(df);
}
