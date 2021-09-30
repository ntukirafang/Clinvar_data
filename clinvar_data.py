#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import sys
import requests
import pandas as pd
import xml.etree.ElementTree as ET
Clinvar_original = ET.parse(r'C:\RNA-seq\Clinvar_database\test.xml')
root = Clinvar_original.getroot()
for nucleotide in root:
    for Ref in nucleotide.iter('ReferenceClinVarAssertion'):
        for Met in Ref.iter('Measure'):
            Type_of_variant = Met.attrib['Type']
            print(Type_of_variant)
            for Name in Met.iter('Name'):
                for ElementValue in Name.iter('ElementValue'):
                    if len((ElementValue.text).split(",")) > 1:
                        emty_list = []
                        for element in ElementValue.text.split(',')[1:]:
                            emty_list.append(element.strip())
                        nucleotide_change = emty_list
                        print(nucleotide_change)
                    else:
                        continue
            for XRef in Met.iter('XRef'):
                if XRef.attrib['DB'] == "OMIM" and XRef.attrib['Type'] != "MIM":
                    OMIM = XRef.attrib['ID']
                    print(OMIM)
                elif XRef.attrib['DB'] == "OMIM" and XRef.attrib['Type'] == "MIM":
                    MIM = XRef.attrib['ID']
                    print(MIM)
                elif  XRef.attrib['DB'] == "dbSNP" and XRef.attrib['Type'] == "rs":
                    dbSNP = "rs" + XRef.attrib['ID']
                    print(dbSNP)
                elif XRef.attrib['DB'] == "MedGen":
                    MedGen = XRef.attrib['ID']
                    print(MedGen)
                elif XRef.attrib['DB'] == "HGNC":
                    HGNC = XRef.attrib['ID'].split(":")[1]
                    print(HGNC)
                elif XRef.attrib['DB'] == "Gene":
                    Genecard = XRef.attrib['ID']
                    print(Genecard)
            for chromosome_pos in Met.iter('CytogeneticLocation'):
                Chr_pos = chromosome_pos.text
                print(Chr_pos)
        for MesSet in Ref.iter('MeasureSet'):
            ACC = MesSet.attrib['Acc']
            print(ACC)
            for Mes in MesSet.iter('Measure'):
                for Mes_Rel in Mes.iter('MeasureRelationship'):
                    for Symbol in Mes_Rel.iter('Symbol'):
                        for ElementValue in Symbol.iter('ElementValue'):
                            Gene = ElementValue.text
                            print(Gene)
                    for SequenceLocation in Mes_Rel.iter('SequenceLocation'):
                        if SequenceLocation.attrib['Assembly'] == "GRCh37":
                            GRCh37_Strand = SequenceLocation.attrib['Strand']
                            print(GRCh37_Strand)
                            GRCh37_Start = SequenceLocation.attrib['display_start']
                            print(GRCh37_Start)
                            GRCh37_Stop = SequenceLocation.attrib['display_stop']
                            print(GRCh37_Stop)
                            GRCh37_Chr = SequenceLocation.attrib['Chr']
                            print(GRCh37_Chr)
                        elif SequenceLocation.attrib['Assembly'] == "GRCh38":
                            GRCh38_Strand = SequenceLocation.attrib['Strand']
                            print(GRCh38_Strand)
                            GRCh38_Start = SequenceLocation.attrib['display_start']
                            print(GRCh38_Start)
                            GRCh38_Stop = SequenceLocation.attrib['display_stop']
                            print(GRCh38_Stop)
                            GRCh38_Chr = SequenceLocation.attrib['Chr']
                            print(GRCh38_Chr)
        for Cli in Ref.iter('ClinicalSignificance'):
            for Des in Cli.iter('Description'):
                Pathogenicity = print(Des.text)
                
        for TraitSet in Ref.iter('TraitSet'):
            for Trait in TraitSet.iter('Trait'):
                for Name in Trait.iter('Name'):
                    for ElementValue in Name.iter('ElementValue'):
                        if ElementValue.attrib['Type'] == "Alternate":
                            Alternate_disease = print("Alternate" + ":" + ElementValue.text)
                        elif ElementValue.attrib['Type'] == "Preferred":
                            Preferred_disease = print("Preferred" + ":" + ElementValue.text)
    for ClinVarAssertion in nucleotide.iter('ClinVarAssertion'):
        for Omim_ID in ClinVarAssertion.iter('ExternalID'):
            OMIM = print(Omim_ID.attrib['ID'])
        for ClinicalSig in ClinVarAssertion.iter('ClinicalSignificance'):
            for Review_Status in ClinicalSig.iter('ReviewStatus'):
                Review_status = print(Review_Status.text)
            for Description in ClinicalSig.iter('Description'):
                description = print(Description.text)        
