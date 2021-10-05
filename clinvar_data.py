#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import sys
import requests
import pandas as pd
import xml.etree.ElementTree as ET
Type_of_variant = []
nucleotide_change = []
OMIM = []
MIM = []
dbSNP = []
MedGen = []
HGNC = []
Genecard = []
Chr_pos = []
ACC = []
Gene = []
GRCh37_Strand = []
GRCh37_Start = []
GRCh37_Stop = []
GRCh37_Chr = []
GRCh38_Strand = []
GRCh38_Start = []
GRCh38_Stop = []
GRCh38_Chr = []
Pathogenicity = []
Disease = []
Review_status = []
description = []
df = open(r'C:\RNA-seq\Clinvar_database\Clinvar.csv','w')
Clinvar_original = ET.parse(r'C:\RNA-seq\Clinvar_database\test.xml')
root = Clinvar_original.getroot()
index_length = len(root)
j = 0
for nucleotide in root:
    print(nucleotide.attrib)
    j = j + 1
    print(j)
    for Ref in nucleotide.iter('ReferenceClinVarAssertion'):
        for MESSETTYPE in Ref.iter('MeasureSet'):
            for Measure in MESSETTYPE.iter('Measure'):
                print(Measure.findall(".//Attribute"))
                if len(Measure.findall(".//Attribute")) == True:
                    print("yes")
                for AttribSet in Measure.iter('AttributeSet'):
                    for Attribute in AttribSet.iter('Attribute'):
                        nucleotide_change.append(Attribute.text)
                        print(nucleotide_change)
                if len(Measure.findall(".//Attribute")) != True:
                    print("no")
                    for Name in Measure.iter('Name'):
                        for ElementValue in Name.iter('ElementValue'):
                            if ElementValue.attrib['Type'] == "Preferred":
                                variant = ElementValue.text.split(',')[1:]
                                if len(variant) == True:
                                    nucleotide_change.append(variant[0].replace(" ","",1))
                                    print(nucleotide_change)
                                else:
                                    continue
                            else:
                                continue
                else:
                    continue
        for Met in Ref.iter('Measure'):
            Type_of_variant.append(Met.attrib['Type'])
            print(Type_of_variant)
            for XRef in Met.iter('XRef'):
                if XRef.attrib['DB'] == "OMIM" and XRef.attrib['Type'] == "MIM":
                    MIM.append(XRef.attrib['ID'])
                    print(MIM)
                elif  XRef.attrib['DB'] == "dbSNP" and XRef.attrib['Type'] == "rs":
                    dbSNP.append("rs" + XRef.attrib['ID'])
                    print(dbSNP)
                elif XRef.attrib['DB'] == "MedGen":
                    MedGen.append(XRef.attrib['ID'])
                    print(MedGen)
                elif XRef.attrib['DB'] == "HGNC":
                    HGNC.append(XRef.attrib['ID'].split(":")[1])
                    print(HGNC)
                elif XRef.attrib['DB'] == "Gene":
                    Genecard.append(XRef.attrib['ID'])
                    print(Genecard)
                else:
                    continue
        for Met in Ref.iter('Measure'):
            for XRef in Met.iter('XRef'):
                if len(MIM)<j:
                    MIM.append("None")
                elif len(dbSNP)<j:
                    dbSNP.append("None")
                    print(dbSNP)
                elif len(MedGen)<j:
                    MedGen.append("None")
                    print(MedGen)
                elif len(HGNC)<j:
                    HGNC.append("None")
                    print(HGNC)
                elif len(Genecard)<j:
                    Genecard.append("None")
                    print(Genecard)
                else:
                    continue
            for chromosome_pos in Met.iter('CytogeneticLocation'):
                Chr_pos.append(chromosome_pos.text)
                print(Chr_pos)
        for MesSet in Ref.iter('MeasureSet'):
            ACC.append(MesSet.attrib['Acc'])
            print(ACC)
            for Mes in MesSet.iter('Measure'):
                for Mes_Rel in Mes.iter('MeasureRelationship'):
                    for Symbol in Mes_Rel.iter('Symbol'):
                        for ElementValue in Symbol.iter('ElementValue'):
                            Gene.append(ElementValue.text)
                            print(Gene)
                    for SequenceLocation in Mes_Rel.iter('SequenceLocation'):
                        if SequenceLocation.attrib['Assembly'] == "GRCh37":
                            GRCh37_Strand.append(SequenceLocation.attrib['Strand'])
                            print(GRCh37_Strand)
                            GRCh37_Start.append(SequenceLocation.attrib['display_start'])
                            print(GRCh37_Start)
                            GRCh37_Stop.append(SequenceLocation.attrib['display_stop'])
                            print(GRCh37_Stop)
                            GRCh37_Chr.append(SequenceLocation.attrib['Chr'])
                            print(GRCh37_Chr)
                        elif SequenceLocation.attrib['Assembly'] == "GRCh38":
                            GRCh38_Strand.append(SequenceLocation.attrib['Strand'])
                            print(GRCh38_Strand)
                            GRCh38_Start.append(SequenceLocation.attrib['display_start'])
                            print(GRCh38_Start)
                            GRCh38_Stop.append(SequenceLocation.attrib['display_stop'])
                            print(GRCh38_Stop)
                            GRCh38_Chr.append(SequenceLocation.attrib['Chr'])
                            print(GRCh38_Chr)
        for Cli in Ref.iter('ClinicalSignificance'):
            for Des in Cli.iter('Description'):
                Pathogenicity.append(Des.text)
                print(Pathogenicity)
    for ClinVarAssertion in nucleotide.iter('ClinVarAssertion'):
        for ClinVarSubmission in ClinVarAssertion.iter('ClinVarSubmissionID'):
            Disease.append(ClinVarSubmission.attrib['localKey'].split('_')[1])
            print(Disease)
    for ClinVarAssertion in nucleotide.iter('ClinVarAssertion'):
        for Omim_ID in ClinVarAssertion.iter('ExternalID'):
            OMIM.append(Omim_ID.attrib['ID'])
            print(OMIM)
        for ClinicalSig in ClinVarAssertion.iter('ClinicalSignificance'):
            for Review_Status in ClinicalSig.iter('ReviewStatus'):
                Review_status.append(Review_Status.text)
                print(Review_status)
            for Description in ClinicalSig.iter('Description'):
                description.append(Description.text)
                print(description)
for i in range(0,(index_length)-1):
    print((Type_of_variant)[0])
    df = pd.DataFrame({'Type_of_variant':Type_of_variant,
                       'nucleotide_change':nucleotide_change, 
                       'OMIM_ID':OMIM,
                       'MIM_ID':MIM,
                       'dbSNP_ID':dbSNP,
                       'Med_Gen':MedGen,
                       'HGNC':HGNC, 
                       'Gene_card':Genecard,
                       'Chr_pos':Chr_pos,
                       'ACC_ID':ACC,
                       'Gene':Gene,
                       'GRCh37_Strand':GRCh37_Strand,
                       'GRCh37_Start':GRCh37_Start,
                       'GRCh37_Stop':GRCh37_Stop,
                       'GRCh37_Chr':GRCh37_Chr,
                       'GRCh38_Strand':GRCh38_Strand,
                       'GRCh38_Start':GRCh38_Start,
                       'GRCh38_Stop':GRCh38_Stop,
                       'GRCh38_Chr':GRCh38_Chr,
                       'Pathogenicity':Pathogenicity,
                       'Disease':Disease,
                       'Description':description,
                       'Preview_Status':Review_status})
df.to_csv(r'C:\RNA-seq\Clinvar_database\Clinvar.csv')
