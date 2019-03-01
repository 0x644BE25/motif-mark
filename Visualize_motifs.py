#!/usr/bin/env python3

#########################################
# Visualize_Motifs.py                   #
# Carolyn Brewster cbrewste@uoregon.edu #
# 2019.02.28                            #
#                                       #
# Creates a visualization of motifs     #
# in provided FASTA sequences.          #
#########################################


#---------------- IMPORTS --------------#

import cairo
import re
import argparse

#---------------- GLOBAL VARS ----------#

global ntDict, colors
ntDict = {'T':'[tT]','A':'[aA]','C':'[cC]','U':'[tT]','G':'[gG]','Y':'[ctCT]', 
          'R':'[agAG]','W':'[atAT]','S':'[gcGC]','K':'[tgTG]','M':'[caCA]','N':'[acugACUG]'}

#--------------- METHODS ---------------#

def getArgs():
    """ get command line input """
    parser=argparse.ArgumentParser(description="Visualize motifs in FASTA sequences")
    parser.add_argument("-f", "--fasta", help="path to FASTA file", required=True, type=str)
    parser.add_argument("-m", "--motif", help="path to motif text file", required=True, type=str)
    parser.add_argument("-o", "--output", help="(optional) output file, default is motifs.svg", default='motifs.svg', type=str)
    parser.add_argument("-c", "--colors", help="(optional) text file of colors in (r,g,b) format, 1 per line", type=str)
    parser.add_argument("-a", "--allowOverlaps", help="allow overlaps, defaults to True", default=True, type=bool)
    
    return(parser.parse_args())

def makePattern(motif,allowOverlaps):
    """ From nucleotide sequence, create appropriate regex query, with any T --> U. """
    regex = '('
    if allowOverlaps: 
        regex = regex+'?='
    for char in motif:
        regex = regex+ntDict[char.upper()]
    regex = regex+')'
    return(regex)

def getMotifs(filename,colors,allowOverlaps):
    """ Parse file of motifs, return list of tuples of (motif, regex, color)."""
    file = open(filename, 'r')
    motifs = []
    i = 0
    for line in file.readlines():
        motifs.append((line.strip(), makePattern(line.strip(),allowOverlaps), colors[i]))
        i += 1
        
    file.close()
    
    return(motifs)

def parseFile(filename):
    """ Parse FASTA or FASTQ file and return list of geneID, sequence tuples. """
    file = open(filename, 'r')
    res = []
    maxlen = 0
    name,seq = "",""
    
    # iterate over lines
    for line in file.readlines():
        if line.startswith('>'):
            if len(name) >= 1:
                res.append((name,seq))
                maxlen = max(maxlen,len(seq))
                seq = ""
            name = line[1:].strip()
        else:
            seq = seq+line.strip()
    
    # take care of the last one
    res.append((name,seq))
    maxlen = max(maxlen,len(seq))
            
    file.close()
    
    return(res, maxlen)

def drawTranscript(context,seq,i):
    """ Draw exon and flanking intronson existing SVG. """
    name, bases = seq[0],seq[1]
    
    context.set_source_rgb(.5,.5,.5)
    x,y = 50, i*150
    context.move_to(x,y)
    context.show_text(name)
    y += 50
    context.move_to(x,y)
    
    for char in bases:
        x += 1
        if char.islower():
            context.set_line_width(10)
        else:
            context.set_line_width(30)
        context.line_to(x,y)
        context.stroke()
        context.move_to(x,y)
    
    return(True)

def hexToColor(hex):
    """ Convert hexadecimal color codes to PyCairo's rgb format. """
    hex = hex.strip("#")
    r = int(hex[0:2], 16)/255
    g = int(hex[2:4], 16)/255
    b = int(hex[4:6], 16)/255

    return(r,g,b)

    
def rgbToColor(rgb):
    """ Convert standard RGB to PyCairo's format. """
    r,g,b = rgb
    
    return(int(r)/255,int(g)/255,int(b)/255)


def getColors(filename):
    """ Read a file of RGB and/or hex colors and convert to PyCairo's RGB format. """
    res = []
    file = open(filename, 'r')
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        if line == '':
            pass
        else:
            line = line.strip("()")
            line = line.split(',')
            if len(line) > 1:
                res.append(rgbToColor(line))
            else:
                res.append(hexToColor(line[0]))
            
    return(res)


def drawMotifs(context,motif,seq,i):
    """ Draw motifs on existing transcript. """
    sites = [s.start() for s in re.finditer(motif[1],seq)]
    x,y = 50, (i*150)+50
    context.set_line_width(30)
    color = motif[2]
    context.set_source_rgb(color[0],color[1],color[2])
    for site in sites:
        context.move_to(site+50,y)
        context.line_to(site+50+len(motif[0]),y)
        context.stroke()
    return(True)


def doSeq(context,seq,motifs,i):
    """ Find and mark all instances of each motif in a sequence. """
    drawTranscript(context,seq,i)
    for motif in motifs:
        pass
        drawMotifs(context,motif,seq[1],i)
    return(True)


def main():
    
    args = getArgs()
    fasta,motiffile,allowOverlaps,output = args.fasta, args.motif,args.allowOverlaps,args.output
    colors = [(0.933, 0.376, 0.333), (0.376, 0.827, 0.580), 
          (1.0, 0.851, 0.490), (1.0, 0.608, 0.522), (0.659, 0.855, 0.863), 
          (0.271, 0.482, 0.596), (0.114, 0.208, 0.341), (0.667, 0.965, 0.514), (1.0, 0.624, 0.110)]
    if args.colors != None:
        print(hmm)
        colors = getColors(args.colors)

    seqs,maxlen = parseFile(fasta)
    motifs = getMotifs(motiffile,colors,allowOverlaps)
    
    surface = cairo.SVGSurface(output, maxlen+100, (len(seqs)+1)*150)
    context = cairo.Context(surface)
    context.set_font_size(25)
    context.move_to(50,50)
    
    for i in range(len(motifs)):
        name,color = motifs[i][0],motifs[i][2]
        context.move_to(100*(i+1)-50, 50)
        context.set_source_rgb(color[0],color[1],color[2])
        context.show_text(name)
    context.move_to(50,50)
    
    i = 1
    for seq in seqs:
        doSeq(context,seq,motifs,i)
        i += 1
    surface.finish()
    
    return(True)

if __name__ == "__main__":
    main()
