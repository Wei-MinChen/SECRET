////////////////////////////////////////////////////////////// ////////
// GDT.cpp
// (c) 2023-2026 Wei-Min Chen
//
// This file is distributed as part of the SECRET source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile SECRET.
//
// All computer programs have bugs. Use this file at your own risk.
//
// May 22, 2026

#include <math.h>
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

//extern int plotRSQ(const char *input, const char *traitName, const char *rpath);

void GDT_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *tuneTrait, bool rplotFlag, int coreCount)
{
    printf("\nOptions in effect:\n");
    printf("\t--gdt\n");
    printf("\t--score %s\n", rsfile);
    printf("\t--pheno %s\n", phefile);
    if (gTraitList[0]) printf("\t--trait %s\n", gTraitList);
    if(dictfile[0])
        printf("\t--diction %s\n", dictfile);
    if(tuneTrait[0])
        printf("\t--strat %s\n", tuneTrait);
    if (rplotFlag)
        printf("\t--rplot\n");
    if (coreCount)
        printf("\t--cpus %d\n", coreCount);
    if (prefix)
        printf("\t--prefix %s\n", prefix);
    printf("\n");
    StringArray familyID(0), sampleID, alltraitNames;
    Matrix alltraits;   // Step 1: Read in traits
    ReadTraits(phefile, alltraitNames, familyID, sampleID, alltraits, coreCount);
    int alltraitCount = alltraitNames.Length();
    int idCount = sampleID.Length();
    IntArray covariates(0), affections(0), quantraits(0);//Step 2: Make traits
    MakeTraits(alltraits, alltraitNames, covariates, affections, quantraits);
    int affectionCount = affections.Length();
    if (affectionCount == 0) { printf("Dichotomous traits not found.\n"); return; }
    Matrix risks;
    StringArray rsNames;    // Step 3: Read in risk scores

    IntArray rsKey;
    ReadRisks(rsfile, rsNames, familyID, sampleID, risks, rsKey, coreCount);//Step 4
    int keyCount = rsKey.Length();
    if (keyCount < idCount) {
        for (int i = 0; i < keyCount; i++) {
            int i0 = rsKey[i];
            familyID[i] = familyID[i0];
            sampleID[i] = familyID[i0];
            for (int t = 0; t < alltraitCount; t++)
                alltraits[t][i] = alltraits[t][i0];
        }   // End of i loop
        idCount = keyCount;
        familyID.Dimension(idCount);
        sampleID.Dimension(idCount);
        alltraits.Dimension(alltraitCount, idCount);
    }

//    ReadRisks(rsfile, rsNames, familyID, sampleID, risks, coreCount);
    int rsCount = rsNames.Length();
    printf("%d risk scores for %d samples loaded.\n", rsCount, idCount);
    StringArray annotation(0);
    if (dictfile[0]) {  // Step 4: Read in risk score annotations
        ReadKeys(dictfile, rsNames, annotation);
        printf("Data dictionary loaded.\n");
    }
    int annotationCount = annotation.Length();
    int famCount = 0;
    IntArray famIndex(0);
    StringIntHash fidLookup;
    fidLookup.Clear();
    for (int i = 0; i < idCount; i++) {
        int k = fidLookup.Integer(familyID[i]);
        if (k == -1) {  // starting a new family
            famIndex.Push(i);
            fidLookup.SetInteger((const char*)familyID[i], famCount);
            famCount++;
        }
    }   // End of i loop 
    IntArray *famMember = new IntArray[famCount];
    for (int f = 0; f < famCount; f++)
        famMember[f].Dimension(0);
    for (int i = 0; i < idCount; i++) {
        int f = fidLookup.Integer(familyID[i]);
        if (f == -1) error("Hash errors for family IDs");
        famMember[f].Push(i);
    }   // End of i loop
    IntArray famSize(famCount);
    for (int f = 0; f < famCount; f++)
        famSize[f] = famMember[f].Length();
    bool popFlag = true;
    for(int f = 0; f < famCount; f++)
        if (famSize[f] > 1) {
            popFlag = false;
            break;
        }
    if(popFlag){
        printf("Families are needed for GDT scans.\n");
        return;
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int r = 0; r < rsCount; r++) {
        double sd0 = sqrt(risks[r].Var());
        if (sd0 < 1E-10)
            for (int i = 0; i < idCount; i++)
                risks[r][i] = _NAN_;
        else
            for (int i = 0; i < idCount; i++)
                risks[r][i] /= sd0;
    }   // End of r loop
    Matrix GR(famCount, rsCount);
    Vector famTrait;
    Matrix sigma(famCount, rsCount);
    Vector gamma(famCount);
    IntArray *sampleSize = new IntArray[famCount];
    for (int f = 0; f < famCount; f++) sampleSize[f].Dimension(rsCount);
    Vector Z(rsCount), Rsq(rsCount);
    String outfile0, outfile;
    Vector oneChisq;
    QuickIndex idx;
    char buffer[0x20000];
    for (int a = 0; a < affectionCount; a++) {
        GR.Zero();
        sigma.Zero();
        int trait = affections[a];
        printf("PGS-wide scan for trait %s starts...\n", (const char*)alltraitNames[trait]);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int f = 0; f < famCount; f++) {
            sampleSize[f].Zero();
            for (int r = 0; r < rsCount; r++) {
                double var = 0.0;
                double sumG = 0.0;
                double sumY = 0.0;
                double gr = 0.0;
                int memberCount = 0;
                for (int i = 0; i < famSize[f]; i++) {
                    int id = famMember[f][i];
                    double y = alltraits[trait][id];
                    double g = risks[r][id];
                    if (y != _NAN_ && y != 0.0 && g != _NAN_) {
                        gr += g * y;
                        var += g * g;
                        sumG += g;
                        sumY += y;
                        memberCount++;
                    }
                }   // End of i loop
                if (memberCount < 2) continue;
                GR[f][r] = gr - sumG * sumY / memberCount;
                sigma[f][r] = (var - sumG * sumG / memberCount) / memberCount;
                sampleSize[f][r] = memberCount;
            }   // End of r loop
        }   // End of f loop

        IntArray ss(rsCount);
        ss.Zero();
        for (int r = 0; r < rsCount; r++)
            for (int f = 0; f < famCount; f++)
                ss[r] += sampleSize[f][r];
        gamma.Zero();
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int f = 0; f < famCount; f++) {
            double sumSigma = sigma[f].Sum();
            if(sumSigma != 0.0)
                gamma[f] = GR[f].InnerProduct(GR[f]) / sumSigma;
        }   // End of f loop
        printf("  The average inflation factor gamma is %.3lf\n", gamma.Average());
        Z.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int r = 0; r < rsCount; r++) {
            double sum = 0.0;
            double z = 0.0;
            for (int f = 0; f < famCount; f++) {
                z += GR[f][r];
                sum += gamma[f] * sigma[f][r];
            }   // End of f loop
            if (sum != 0.0)
                Z[r] = z / sqrt(sum);
        }   // End of r loop
        const double ZCUTOFF = 1.959964;
        IntArray validrs(0);
        for (int r = 0; r < rsCount; r++) {
            double z = Z[r];
            if (z == _NAN_ || z > ZCUTOFF || z < -ZCUTOFF) continue;
            validrs.Push(r);
        }   // End of r loop
        int gammarsCount = validrs.Length();
        for (int f = 0; f < famCount; f++) {
            double sumSigma = 0.0, gr2 = 0.0;
            for (int r = 0; r < gammarsCount; r++) {
                int rs = validrs[r];
                sumSigma += sigma[f][rs];
                gr2 += GR[f][rs] * GR[f][rs];
            }   // End of r loop
            if(sumSigma != 0.0)
                gamma[f] = gr2 / sumSigma;
        }   // End of f loop
        printf("  The average inflation factor gamma is %.3lf (with %d PGS)\n", gamma.Average(), gammarsCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int r = 0; r < rsCount; r++) {
            double sum = 0.0;
            double z = 0.0;
            double rsq = 0.0;
            double sum2 = 0.0;
            for (int f = 0; f < famCount; f++) {
                z += GR[f][r];
                sum += gamma[f] * sigma[f][r];
                rsq += GR[f][r] * GR[f][r];
                sum2 += gamma[f] * sigma[f][r] * sampleSize[f][r];
            }   // End of f loop
            if (sum != 0.0) {
                Z[r] = z / sqrt(sum);
                Rsq[r] = rsq / sum2;
            }
        }   // End of r loop
        outfile0.Copy(prefix);
        outfile0.Add("_");
        outfile0.Add(alltraitNames[trait]);
        outfile = outfile0;
        outfile.Add(".gdt");
        FILE *fp = fopen(outfile, "wb");
        if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
        int annotCount = annotation.Length();
        int pbuffer = 0;
        pbuffer += sprintf(&buffer[pbuffer], "PGS\tN\tZ\tPvalue\tRsq");
        if (annotCount)
            pbuffer += sprintf(&buffer[pbuffer], "\tAnnotation");
        validrs.Dimension(0);
        oneChisq.Dimension(0);
        for (int r = 0; r < rsCount; r++) {
            if (Z[r] == _NAN_) continue;
            oneChisq.Push(Z[r] * Z[r]);
            validrs.Push(r);
        }   // End of r loop
        int validrsCount = validrs.Length();
        idx.Index(oneChisq);
        for (int r = 0; r < validrsCount; r++) {
            int rs = validrs[idx[validrsCount - 1 - r]];
            pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d\t%.4lf\t%.2G\t%.5lf", 
                (const char*)rsNames[rs], ss[rs], Z[rs], ndist(fabs(Z[rs]))*2, Rsq[rs]);
            if (annotCount)
                pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)annotation[rs]);
            if (pbuffer > 0xFFFF) {  // buffer big enough for writing
                fwrite(buffer, 1, pbuffer, fp);
                pbuffer = 0;
            }
        }   // End of r loop
        pbuffer += sprintf(&buffer[pbuffer], "\n");
        fwrite(buffer, 1, pbuffer, fp);
        fclose(fp);
        printf("  Genomic control lambda for PGS-wide scan is %.5lf\n", oneChisq[idx[validrsCount / 2]] / 0.4549364);
        printf("  GDT scan of %d risk scores for trait %s writes to %s\n",
            validrsCount, (const char*)alltraitNames[trait], (const char*)outfile);

        outfile = outfile0;
        outfile.Add(".rs");
        fp = fopen(outfile, "wb");
        if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
        pbuffer = 0;
        pbuffer += sprintf(&buffer[pbuffer], "ID\t%s", (const char*)alltraitNames[trait]);
        IntArray topRS(0);
        for (int r = 0; r < 10; r++) {
            int rs = validrs[idx[validrsCount - 1 - r]];
            topRS.Push(rs);
            pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)rsNames[rs]);
        }   // End of r loop
        Matrix riskAdj(10, idCount);
        riskAdj.Set(_NAN_);
        for (int f = 0; f < famCount; f++) {
            for (int r = 0; r < 10; r++) {
                int rs = topRS[r];
                double sumG = 0.0;
                int memberCount = 0;
                for (int i = 0; i < famSize[f]; i++) {
                    int id = famMember[f][i];
                    double y = alltraits[trait][id];
                    double g = risks[rs][id];
                    if (y != _NAN_ && y != 0.0 && g != _NAN_) {
                        sumG += g;
                        memberCount++;
                    }
                }   // End of i loop
                if (memberCount < 2) continue;
                for (int i = 0; i < famSize[f]; i++) {
                    int id = famMember[f][i];
                    riskAdj[r][id] = risks[rs][id] - sumG / memberCount;
                }   // End of i loop
            }   // End of r loop
        }   // End of f loop
        for (int i = 0; i < idCount; i++) {
            if (alltraits[trait][i] == _NAN_ || alltraits[trait][i] == 0 || risks[topRS[0]][i] == _NAN_ || riskAdj[0][i]==_NAN_) continue;
            pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d", (const char*)sampleID[i], int(alltraits[trait][i]+0.5));
            for (int r = 0; r < 10; r++)
                pbuffer += sprintf(&buffer[pbuffer], "\t%.3lf", riskAdj[r][i]);
            if (pbuffer > 0xFFFF) {  // buffer big enough for writing
                fwrite(buffer, 1, pbuffer, fp);
                pbuffer = 0;
            }
        }   // End of i loop
        pbuffer += sprintf(&buffer[pbuffer], "\n");
        fwrite(buffer, 1, pbuffer, fp);
        fclose(fp);
        printf("  TOP 10 PGS scores for trait %s writes to %s\n", (const char*)alltraitNames[trait], (const char*)outfile);
        //if(rplotFlag) plotRSQ(outfile0, alltraitNames[trait], "");
    }   // End of a loop
}
