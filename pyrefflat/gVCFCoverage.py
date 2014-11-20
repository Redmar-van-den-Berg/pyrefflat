__author__ = 'ahbbollen'

import argparse
import math
import vcf
import json

import parser as refparser

def asBedgraph(inputf, outputf, reff, mode):
    vcfReader = vcf.Reader(filename=inputf)
    bedWriter = open(outputf, 'wb')
    refReader = open(reff, 'rb')

    for sample in vcfReader.samples:
        trackline = "track type=bedGraph name='{0}' description='gVCFCoverage' graphType='bar'".format(sample)
        bedWriter.write(bytes(trackline + "\n"))
        for line in refReader:
            record = refparser.Record(line, reff)
            for exon in record.exons:
                try:
                    vcf_records = vcfReader.fetch(str(exon.chr), int(exon.start), int(exon.stop))
                except ValueError:
                    vcf_records = []
                DPS = []
                for v_record in vcf_records:
                    DPS.append(v_record.genotype(sample)['DP'])
                bedline = [exon.chr, exon.start, exon.stop]
                if len(DPS) > 0:
                    bedline.append(sum(DPS)/len(DPS))
                    bedline = map(str, bedline)
                    bedWriter.write(bytes("\t".join(bedline) + "\n"))
                else:
                    bedline.append(0)
                    bedline = map(str, bedline)
                    bedWriter.write(bytes("\t".join(bedline) + "\n"))

    bedWriter.close()
    refReader.close()

def asCSV(inputf, outputf, reff, sep, mode):
    vcfReader = vcf.Reader(filename=inputf)
    bedWriter = open(outputf, 'wb')
    refReader = open(reff, 'rb')

    for sample in vcfReader.samples:
        headerline = "Sample {0}; Gene \t Transcript \t Chr \t Exon number \t Exon start \t Exon stop \t Coverage \n".format(sample)
        bedWriter.write(bytes(headerline))

        for line in refReader:
            record = refparser.Record(line, reff)
            gene = record.gene
            transcript = record.transcript
            for exon in record.exons:
                chr = exon.chr
                start = exon.start
                stop = exon.stop
                n = exon.number
                try:
                    vcf_records = vcfReader.fetch(str(chr), int(start), int(stop))
                except ValueError:
                    vcf_records = []
                #print vcf_records
                DPS = []
                for v_record in vcf_records:
                    DPS.append(v_record.genotype(sample)['DP'])
                line = [gene, transcript, chr, n, start, stop]
                if len(DPS) > 0:
                    line.append(sum(DPS)/len(DPS))
                    #print DPS
                else:
                    line.append(0)
                line = map(str, line)
                bedWriter.write(bytes(sep.join(line) + "\n"))

    bedWriter.close()
    refReader.close()


def asJSON(inputf, outputf, reff, mode):
    vcfReader = vcf.Reader(filename=inputf)
    jsonWriter = open(outputf, 'wb')
    refReader = open(reff, 'rb')

    jdict = {}

    for sample in vcfReader.samples:
        jdict[sample] = {}

        for line in refReader:
            record = refparser.Record(line, reff)
            gene = record.gene
            transcript = record.transcript
            jdict[sample][gene] = {"transcript": transcript, "exons" : []}
            for exon in record.exons:
                chr = exon.chr
                start = exon.start
                stop = exon.stop
                n = exon.number
                try:
                    vcf_records = vcfReader.fetch(str(chr), int(start), int(stop))
                except:
                    vcf_records = []
                DPS = []
                for v_record in vcf_records:
                    DPS.append(v_record.genotype(sample)['DP'])
                if len(DPS) > 0:
                    coverage = sum(DPS)/len(DPS)
                else:
                    coverage = 0
                exondict = {"chr": chr, "start": start, "stop" : stop, "number": n, "coverage" : coverage}
                jdict[sample][gene]["exons"].append(exondict)

    jsonWriter.write(json.dumps(jdict))
    jsonWriter.close()
    refReader.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #generic arguments
    parser.add_argument('-I', '--input', help="Your input gzipped gVCF")
    parser.add_argument('-O', '--output', help="Path to output file")
    parser.add_argument('-R', '--refflat', help="Path to refflat file")
    parser.add_argument('-m', '--mode', choices=['BP_RESOLUTION', 'GVCF'], help="gVCF mode")
    parser.add_argument('-om', '--output-mode', choices=['bedgraph', 'csv', 'json'], help="Output file type")

    # optional arguments
    parser.add_argument('--tab-delimited', action='store_true', help='Output CSV as TSV')

    args = parser.parse_args()

    if args.output_mode == 'bedgraph':
        asBedgraph(args.input, args.output, args.refflat, args.mode)
    elif args.output_mode == 'csv':
        if args.tab_delimited:
            sep = "\t"
        else:
            sep = ","
        asCSV(args.input, args.output, args.refflat, sep, args.mode)
    elif args.output_mode == 'json':
        asJSON(args.input, args.output, args.refflat, args.mode)

