////////////////////////////////////////////////////////////// ////////
// TDT.cpp
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
// Oct 2, 2026

#include <math.h>
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

extern double ComputeAUC(Vector & risks, IntArray & diseaseStatus, Vector & stat);

void TDT_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *tuneTrait, int nTop, int coreCount)
{
  printf("\nOptions in effect:\n");
  printf("\t--tdt\n");
  printf("\t--score %s\n", rsfile);
  printf("\t--pheno %s\n", phefile);
  if (gTraitList[0]) printf("\t--trait %s\n", gTraitList);
  if (dictfile[0])
    printf("\t--diction %s\n", dictfile);
  if (tuneTrait[0])
    printf("\t--strat %s\n", tuneTrait);
  if (nTop)
    printf("\t--ntop %d\n", nTop);
  if (coreCount)
    printf("\t--cpus %d\n", coreCount);
  if (prefix)
    printf("\t--prefix %s\n", prefix);
  printf("\n");
  StringArray familyID(0), sampleID, familyID0, sampleID0, FAID, MOID, FAID0,MOID0, alltraitNames;
  Matrix alltraits;   // Step 1: Read in traits
  ReadTraits4(phefile, alltraitNames, familyID, sampleID, FAID, MOID, alltraits, coreCount);
  int alltraitCount = alltraitNames.Length();
  int idCount = sampleID.Length();
  IntArray covariates(0), affections(0), quantraits(0);//Step 2: Make traits
  MakeTraits(alltraits, alltraitNames, covariates, affections, quantraits);
  int affectionCount = affections.Length();
  if (affectionCount == 0) { printf("Dichotomous traits not found.\n"); return; }
  Matrix risks;
  StringArray rsNames;    // Step 3: Read in risk scores
  IntArray rsKey;
  familyID0 = familyID;
  sampleID0 = sampleID;
  FAID0 = FAID;
  MOID0 = MOID;
  ReadRisks(rsfile, rsNames, familyID, sampleID, risks, rsKey, coreCount);//Step 4
  int keyCount = rsKey.Length();
  if (keyCount < idCount) {
    for (int i = 0; i < keyCount; i++) {
      int i0 = rsKey[i];
      familyID[i] = familyID0[i0];
      sampleID[i] = sampleID0[i0];
      FAID[i] = FAID0[i0];
      MOID[i] = MOID0[i0];
      for (int t = 0; t < alltraitCount; t++)
        alltraits[i][t] = alltraits[i0][t];
    }   // End of i loop
    idCount = keyCount;
    familyID.Dimension(idCount);
    sampleID.Dimension(idCount);
    alltraits.Dimension(idCount, alltraitCount);
  }else{
    familyID = familyID0;
    sampleID = sampleID0;
    FAID = FAID0;
    MOID = MOID0;
  }
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
  for (int f = 0; f < famCount; f++)
    if (famSize[f] > 1) {
      popFlag = false;
      break;
      }
  if (popFlag) {
    printf("Families are needed for TDT scans.\n");
    return;
  }
  IntArray trios(0);
  StringArray iid;
  for (int f = 0; f < famCount; f++) {
    iid.Clear();
    for (int i = 0; i < famSize[f]; i++)
      iid.Push(sampleID[famMember[f][i]]);
    for (int i = 0; i < famSize[f]; i++){
      int id = famMember[f][i];
      if (FAID[id] != "0" && MOID[id] != "0") {
        int fa = iid.Find(FAID[id]);
        int mo = iid.Find(MOID[id]);
        if (fa > -1 && mo > -1) {
          trios.Push(id);
          trios.Push(famMember[f][fa]);
          trios.Push(famMember[f][mo]);
        }
      }
    }   // End of i loop
  }   // End of f loop
  int trioCount = trios.Length() / 3;
  if (trioCount == 0) {
    printf("No trios found in the dataset.\n"); return;
  }else
    printf("%d trios in the dataset.\n", trioCount);
  int rsCount = rsNames.Length();
  printf("%d polygenic scores for %d samples loaded.\n", rsCount, idCount);
  StringArray annotation(0);
  if (dictfile[0]) {  // Step 4: Read in annotations
    ReadKeys(dictfile, rsNames, annotation);
    printf("Data dictionary loaded.\n");
  }
  int annotCount = annotation.Length();
  Matrix trioRisk(rsCount, trioCount); trioRisk.Set(_NAN_);
  Matrix adjRisk[3];
  for(int p = 0; p < 3; p++)
    adjRisk[p].Dimension(rsCount, trioCount);
  for (int t = 0; t < trioCount; t++) {
    int id = trios[t * 3];
    int fa = trios[t * 3 + 1];
    int mo = trios[t * 3 + 2];
    for (int r = 0; r < rsCount; r++) {
      if (risks[r][fa] == _NAN_ || risks[r][mo] == _NAN_ || risks[r][id] == _NAN_) continue;
      trioRisk[r][t] = risks[r][id];
      adjRisk[1][r][t] = risks[r][id] - risks[r][fa];
      adjRisk[2][r][t] = risks[r][id] - risks[r][mo];
      adjRisk[0][r][t] = (adjRisk[1][r][t] + adjRisk[2][r][t]) * 0.5;
    }   // End of r loop
  }   // End of t loop
  Vector Z(rsCount), Z_cc(rsCount), AUC[5], pAUC(rsCount), means[3], SDs[3];
  for (int k = 0; k < 3; k++) {
    means[k].Dimension(rsCount);
    SDs[k].Dimension(rsCount);
    SDs[k].Zero();
  } // End of k loop
  for (int p = 0; p < 5; p++)
    AUC[p].Dimension(rsCount);
  String outfile0, outfile;
  Vector oneChisq;
  QuickIndex idx;
  char buffer[0x20000];
  IntArray validrs, affTrio, unaffTrio, A[5];
  Vector Y[3][5], Y0;
  int sextrait = alltraitNames.Find("SEX");
  if (sextrait == -1) { printf("More coding on SEX is needed.\n"); return;}
  for (int a = 0; a < affectionCount; a++) {
    int trait = affections[a];
    if (trait == sextrait) continue;
    affTrio.Dimension(0);
    unaffTrio.Dimension(0);
    for (int i = 0; i < trioCount; i++)
      if (alltraits[trios[i * 3]][trait] == 2.0) affTrio.Push(i);
      else if (alltraits[trios[i * 3]][trait] == 1.0) unaffTrio.Push(i);
    int afftrioCount = affTrio.Length();
    int unafftrioCount = unaffTrio.Length();
    if(unafftrioCount==0) printf("No controls. TDT only.\n");
    printf("PGS-wide scan for trait %s starts...\n", (const char*)alltraitNames[trait]);
    Z.Set(_NAN_);
    Vector stat;
    for (int r = 0; r < rsCount; r++) {
      for (int k = 0; k < 3; k++)
        for(int p = 0; p < 5; p++)
          Y[k][p].Dimension(0);
      for(int p = 0; p < 5; p++)
        A[p].Dimension(0);
      for (int i = 0; i < afftrioCount; i++) {
        for (int p = 0; p < 3; p++) {
          Y[2][p].Push(adjRisk[p][r][affTrio[i]]);
          A[p].Push(2);
        } // End of p loop
        if (int(alltraits[trios[affTrio[i] * 3]][sextrait] + 0.5) == 1) {
          Y[2][3].Push(adjRisk[0][r][affTrio[i]]);
          A[3].Push(2);
        }
        else if (int(alltraits[trios[affTrio[i] * 3]][sextrait] + 0.5) == 2) {
          Y[2][4].Push(adjRisk[0][r][affTrio[i]]);
          A[4].Push(2);
        }
      } // End of i loop
      for (int i = 0; i < unafftrioCount; i++) {
        for (int p = 0; p < 3; p++) {
          Y[1][p].Push(adjRisk[p][r][unaffTrio[i]]);
          A[p].Push(1);
        } // End of p loop
        if (int(alltraits[trios[unaffTrio[i] * 3]][sextrait] + 0.5) == 1) {
          Y[1][3].Push(adjRisk[0][r][unaffTrio[i]]);
          A[3].Push(1);
        }
        else if (int(alltraits[trios[unaffTrio[i] * 3]][sextrait] + 0.5) == 2) {
          Y[1][4].Push(adjRisk[0][r][unaffTrio[i]]);
          A[4].Push(1);
        }
      }   // End of i loop
      for (int p = 0; p < 5; p++) {
        Y[0][p] = Y[2][p];
        if(unafftrioCount) {
          Y[0][p].Stack(Y[1][p]);
          AUC[p][r] = ComputeAUC(Y[0][p], A[p], stat);
        }
      }   // End of p loop
      for (int k = 0; k < 3; k++) {
        if(Y[k][0].Length()==0) continue;
        means[k][r] = Y[k][0].Average();
        SDs[k][r] = sqrt(Y[k][0].Var());
      }   // End of k loop
      if (SDs[2][r] < 1E-100) continue;
      Z[r] = sqrt(afftrioCount)*means[2][r] / SDs[2][r];
      if (unafftrioCount==0 || SDs[1][r] < 1E-100) continue;
      pAUC[r] = ComputeAUC(Y0, A[0], stat);
      Z_cc[r] = (means[2][r] - means[1][r]) / sqrt(SDs[2][r] * SDs[2][r] / afftrioCount + SDs[1][r] * SDs[1][r] / unafftrioCount);
    }   // End of r loop
    outfile0.Copy(prefix);
    outfile0.Add("_");
    outfile0.Add(alltraitNames[trait]);
    outfile = outfile0;
    outfile.Add(".tdt");
    FILE *fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    int pbuffer = 0;
    if(unafftrioCount)
      pbuffer += sprintf(&buffer[pbuffer], "PGS\tN_Case\tPGSCase\tSD_Case\tZ_TDT\tP_TDT\tN_Cont\tPGSCont\tSD_Cont\tZ_CC\tP_CC\tAUC_MP\tAUC_ETM\tAUC_ETF\tAUC_ET");
    else
      pbuffer += sprintf(&buffer[pbuffer], "PGS\tN_Case\tPGSCase\tSD_Case\tZ_TDT\tP_TDT");
    if (annotCount)
      pbuffer += sprintf(&buffer[pbuffer], "\tAnnotation");
    validrs.Dimension(0);
    oneChisq.Dimension(0);
    for (int r = 0; r < rsCount; r++) {
      if (Z[r] == _NAN_) continue;
      oneChisq.Push(Z[r]);
      validrs.Push(r);
    }   // End of r loop
    int validrsCount = validrs.Length();
    idx.Index(oneChisq);
    for (int r = 0; r < validrsCount; r++) {
      int rs = validrs[idx[validrsCount - 1 - r]];
      if(unafftrioCount)
        pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d\t%.4lf\t%.5lf\t%.3lf\t%.2G\t%d\t%.4lf\t%.5lf\t%.3lf\t%.2G\t%.5lf\t%.5lf\t%.5lf\t%.5lf",
        (const char*)rsNames[rs], afftrioCount, means[2][rs], SDs[2][rs], Z[rs], ndist(fabs(Z[rs]))*2,
        unafftrioCount, means[1][rs], SDs[1][rs], Z_cc[rs], ndist(fabs(Z_cc[rs]))*2, pAUC[rs], AUC[3][rs], AUC[4][rs], AUC[0][rs]);
      else
        pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d\t%.4lf\t%.5lf\t%.3lf\t%.2G",
        (const char*)rsNames[rs], afftrioCount, means[2][rs], SDs[2][rs], Z[rs], ndist(fabs(Z[rs]))*2);
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
    printf("  TDT scan of %d polygenic scores for trait %s writes to %s\n",
      validrsCount, (const char*)alltraitNames[trait], (const char*)outfile);
    const int topCOUNT = nTop==0? (rsCount<100?rsCount:100): nTop;
    IntArray topRS(0);
    for (int r = 0; r < topCOUNT; r++) {
      int rs = validrs[idx[validrsCount - 1 - r]];
      topRS.Push(rs);
      pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)rsNames[rs]);
    }   // End of r loop
    outfile = outfile0;
    outfile.Add("_offspringPGS.txt");
    fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "ID\t%s", (const char*)alltraitNames[trait]);
    for (int r = 0; r < topCOUNT; r++)
      pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)rsNames[topRS[r]]);
    for (int i = 0; i < trioCount; i++) {
      int id = trios[i * 3];
      pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d", (const char*)sampleID[id], int(alltraits[id][trait]+0.5));
      for (int r = 0; r < topCOUNT; r++)
        pbuffer += sprintf(&buffer[pbuffer], "\t%.3lf", trioRisk[topRS[r]][i]);
      if (pbuffer > 0xFFFF) {  // buffer big enough for writing
        fwrite(buffer, 1, pbuffer, fp);
        pbuffer = 0;
      }
    }   // End of i loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
    fclose(fp);
    printf("  Top %d offspring PGS scores for trait %s write to %s\n", topCOUNT, (const char*)alltraitNames[trait], (const char*)outfile);
    outfile = outfile0;
    outfile.Add("_ETPGS.txt");
    fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "ID\t%s", (const char*)alltraitNames[trait]);
    for (int r = 0; r < topCOUNT; r++)
      pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)rsNames[topRS[r]]);
    for (int i = 0; i < trioCount; i++) {
      int id = trios[i * 3];
      pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d", (const char*)sampleID[id], int(alltraits[id][trait] + 0.5));
      for (int r = 0; r < topCOUNT; r++)
        pbuffer += sprintf(&buffer[pbuffer], "\t%.3lf", adjRisk[0][topRS[r]][i]);
      if (pbuffer > 0xFFFF) {  // buffer big enough for writing
        fwrite(buffer, 1, pbuffer, fp);
        pbuffer = 0;
      }
    }   // End of i loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
    fclose(fp);
    printf("  Top %d Excess Transmission scores for trait %s write to %s\n", topCOUNT, (const char*)alltraitNames[trait], (const char*)outfile);
    outfile = outfile0;
    outfile.Add("_POOPGS.txt");
    fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "ID\t%s", (const char*)alltraitNames[trait]);
    for (int r = 0; r < topCOUNT; r++)
      pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)rsNames[topRS[r]]);
    for (int i = 0; i < trioCount; i++) {
      int id = trios[i * 3];
      pbuffer += sprintf(&buffer[pbuffer], "\n%s\t%d", (const char*)sampleID[id], int(alltraits[id][trait]+0.5));
      for (int r = 0; r < topCOUNT; r++)
        pbuffer += sprintf(&buffer[pbuffer], "\t%.3lf", (adjRisk[1][topRS[r]][i] - adjRisk[2][topRS[r]][i])*0.5);
      if (pbuffer > 0xFFFF) {  // buffer big enough for writing
        fwrite(buffer, 1, pbuffer, fp);
        pbuffer = 0;
      }
    } // End of i loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
    fclose(fp);
    printf("  Top %d POO PGS scores for trait %s write to %s\n", topCOUNT, (const char*)alltraitNames[trait], (const char*)outfile);
  } // End of a loop
}

