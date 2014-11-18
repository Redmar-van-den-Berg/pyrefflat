__author__ = 'ahbbollen'

import argparse
import pysam
# should be
# from pyrefflat import parser as...
import parser as refparser


def asBedgraph(args):
    inputf, outputf, reff, trackline = args.input, args.output, args.refflat, args.track_line
    bReader = pysam.Samfile(inputf, 'rb')
    rReader = open(reff, 'rb')
    bedWriter = open(outputf, 'wb')
    if trackline:
        bedWriter.write(bytes('track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n', 'UTF-8'))
    for line in rReader:
        record = refparser.Record(line, reff)
        for exon in record.exons:
            n_reads = bReader.count(exon.chr, int(exon.start), int(exon.stop))
            coverage = int(n_reads/(int(exon.stop) - int(exon.start)))
            bedWriter.write(bytes("{0}\t{1}\t{2}\t{3}\n".format(exon.chr, exon.start, exon.stop, coverage), 'UTF-8'))
    return True

def asCSV(args):
    inputf, outputf, reff, sep = args.input, args.output, args.refflat, args.separator
    bReader = pysam.Samfile(inputf, 'rb')
    rReader = open(reff, 'rb')
    csvWriter = open(outputf, 'wb')
    header = ["Gene", "Transcript", "Chr", "Exon", "Exon start", "Exon stop", "Mean coverage"]
    csvWriter.write(bytes(sep.join(header) + "\n", 'UTF-8'))
    for line in rReader:
        record = refparser.Record(line, reff)
        gene = record.gene
        transcript = record.transcript
        for exon in record.exons:
            chr, start, stop, n = exon.chr, exon.start, exon.stop, exon.number
            n_reads = bReader.count(chr, int(start), int(stop))
            coverage = int(n_reads/(int(stop) - int(start)))
            line = [gene, transcript, chr, n, start, stop, coverage]
            line = map(str, line)
            csvWriter.write(bytes(sep.join(line) + "\n", 'UTF-8'))

    return True

def asJson(inputf, outputf, reff):
    return True

def info(inputf, outputf, reff):
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)

    # generic arguments
    parser.add_argument('-I', '--input', help="Input bam file")
    parser.add_argument('-O', '--output', help="Output file")
    parser.add_argument('-R', '--refflat', help="RefFlat file")

    # subcommands
    subparsers = parser.add_subparsers()

    parser_bedgraph = subparsers.add_parser('bedgraph', help="Output as bedgraph", parents=[parser])
    parser_bedgraph.add_argument('--track-line', action='store_true', help="Show track line")
    parser_bedgraph.set_defaults(func=asBedgraph)

    parser_csv = subparsers.add_parser('csv', help="Output as CSV", epilog="Useful for loading in excel", parents=[parser])
    parser_csv.add_argument('-s', '--separator', help="Which field separator to use")
    parser_csv.set_defaults(func=asCSV)

    parser_json = subparsers.add_parser('json', help="Output as JSON", parents=[parser])
    parser_json.set_defaults(func=asJson)

    parser_info = subparsers.add_parser('info', help="Show some info", parents=[parser])
    parser_info.set_defaults(func=info)

    # parse
    args = parser.parse_args()
    args.func(args)
