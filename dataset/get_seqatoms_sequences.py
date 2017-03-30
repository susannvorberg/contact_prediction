import argparse
import wget
import os
import sys
from Bio import pairwise2
import Bio.PDB as BP
from Bio.PDB.Polypeptide import PPBuilder
import math
import utils.io_utils as io

def getPDBInfo(pdb, pdbpath):
    # extract information from the PDB file using Biopython
    # Biopython parsers
    parser = BP.PDBParser()
    ppb = PPBuilder()

    pdbseq = ""
    # PDB descriptors
    # name = pdb.split("/")[-1].split(".")[0].split("_")[0]
    if len(pdb.split("/")[-1].split(".")[0].split("_")) > 1:
        chain = pdb.split("/")[-1].split(".")[0].split("_")[1]
    else:
        # not "_" delimiters in pdb name
        chain = pdb.split("/")[-1].split(".")[0][4]

    if os.path.exists(os.path.join(pdbpath, pdb + "_ren.pdb")):
        structure = parser.get_structure(pdb, os.path.join(pdbpath, pdb + "_ren.pdb"))
    residue_list = BP.Selection.unfold_entities(structure[0][chain], 'R')

    # the build_peptides method has an option for treating non-standard amino acids
    for pp in ppb.build_peptides(structure[0][chain], aa_only=False):
        pdbseq += (pp.get_sequence())
    return structure, residue_list, pdbseq


def eucl(v1, v2):
    xx = (v1[0] - v2[0]) ** 2
    yy = (v1[1] - v2[1]) ** 2
    zz = (v1[2] - v2[2]) ** 2
    return math.sqrt(xx + yy + zz)


def getNAtom(residue):
    if not residue.has_id('N'):
        return None
    else:
        return residue['N']


def getNDistances(res1, res2):
    a1 = getNAtom(res1)
    a2 = getNAtom(res2)
    if a1 != None and a2 != None:
        d1 = eucl(a1.get_vector(), a2.get_vector())
    else:
        d1 = 1000
    return d1



def map_sequences(key, fasta_sequence, mmcif_sequence, pdb_dir):
    """
    Input: fasta_sequence (string), mmCif_sequence (string)
    Output: mapping (dict)
    Example: {1: 1, 2: 3, 3: 4, 4: 5, 5: 7} (position 1 of fasta sequence
    is aligned to 1 of the mmCif sequence and so on)

    Assumes that fasta_sequence and mmCif_sequence are strings and returns a
    alignment that maps position i from the fasta_sequence to position j from
    mmCif_sequence.
    """
    # print fasta_sequence
    # print mmcif_sequence
    match = 2
    mismatch = -2
    template_gap_open = -3
    template_gap_extend = -0.1
    query_gap_open = -3
    query_gap_extend = -0.1
    alns = pairwise2.align.globalmd(fasta_sequence, mmcif_sequence, match, mismatch,
                                    query_gap_open, query_gap_extend, template_gap_open, template_gap_extend)
    if (len(alns) == 0):
        raise Exception("Error: Sequences could not be aligned!\n\t" + fasta_sequence + "\n\t" + mmcif_sequence + "\n")
        pass

    # take = len(alns)-1
    top_aln = alns[0]
    aln_from, aln_to, score, begin, end = top_aln
    # print fasta_sequence
    # print mmcif_sequence
    # print aln_from
    # print aln_to
    if len(alns) > 1:
        aln_to = correct_alignment(key, pdb_dir, aln_from, aln_to)
    # print aln_from
    # print aln_to
    return get_mapping(aln_from, aln_to, 1, 1)


def correct_alignment(pdb, pdb_dir, aln_from, aln_to):
    try:
        structure, residue_list, pdbseq = getPDBInfo(pdb, pdb_dir)
    except UnboundLocalError:
        print "no PDB"
        return aln_to
    # print structure
    # print residue_list
    # print len(residue_list)
    # print len(aln_to.replace("-", ""))

    # find where gaps are
    gaps = find_gaps(pdb, aln_to)
    # print gaps

    # now see if the aminoacids at the edges of the gap are
    # problematic, if yes shift the gap
    for gap in gaps:
        # print gap, aln_from, aln_to, residue_list
        bad_gap, move = check_edge(gap, aln_from, aln_to, residue_list)
        if bad_gap:
            aln_to = correct_gap(move, aln_to)

    return aln_to


def correct_gap(move, aln_to):
    """
    change positions between the characters in positions move[0] and move[1]
    """
    a = aln_to[move[0]]
    b = aln_to[move[1]]

    # print move
    if move[0] == -1:
        for i in range(move[1], len(aln_to)):
            aln_to = aln_to[0:i] + "-" + aln_to[i + 1:]
        return aln_to

    aln_to = aln_to[0:move[1]] + a + aln_to[move[1] + 1:]
    aln_to = aln_to[0:move[0]] + b + aln_to[move[0] + 1:]
    # print aln_to
    return aln_to


def check_edge(gap, aln_from, aln_to, residue_list):
    # print gap
    """
    this is where the fun happens:
    1. check if there is a same aminoacid conflict at the edges of the
    gap
    2. if yes, then check with the structure; where does the ambiguous
    aminoacid belong to in the ATOM record?
    3. if the offending aa is in the correct place everything is rosy.
    If not we gotta put it back in its place --> save the move in "move"
    """
    # two pairs to check in the COMBS sequence:
    # gap-1 --- gap_end
    # gap --- gap_end+1
    # print len(residue_list)
    if aln_from[gap[0] - 1] == aln_from[gap[1]]:
        # go back in the sequence as long as there are other aminoacids:
        # I have to screen for the case where there is a gap in
        # the beginning of the sequence.
        # AAASY...
        # AA-SY...
        # should obviously be -AASY...

        pos_next = 1
        pos_prev = -1

        i = gap[0] - 1
        diff = 0
        while i > 0 and aln_from[i] == aln_from[gap[1]]:
            aln_to = aln_to[0:i] + "-" + aln_to[i + 1:]
            i = i - 1
            diff = diff + 1
            pos_next = pos_next + 1
        # print i
        gap[0] = i + 1
        # print gap
        # print "gap[0]-1 --- gap[1]"
        # print aln_to[0:gap[1]+2]
        # print aln_from[0:gap[1]+2]
        # print aln_to[0:gap[1]+2].replace("-", "")
        # we need to see if this is aligned correctly
        # find which ATOM position gap[0]-1 corresponds to
        # and see if it is closer to gap[0]-2 or to gap[1]+1
        truncated = aln_to[0:gap[0] + 1].replace("-", "")
        # print truncated
        pos = len(truncated) - 1  # <-- this is the position where we need to look in the structure
        # print pos
        # print truncated[pos]
        # print aln_from[pos]
        # print len(residue_list)

        # if pos+1 == len(residue_list):
        # 	return False, ""

        # is res[pos] closer to res[pos+1] or res[pos-1]?
        # if the former, then the alignment has been done correctly and there is nothing to
        # change here --> return False, ""
        # if the latter, then we need to move the gap one aminoacid.

        if pos == 0:
            prev = []
        else:
            prev = residue_list[pos + pos_prev]

        curr = residue_list[pos]

        try:
            next = residue_list[pos + pos_next]
        except IndexError:
            next = []

        # print prev
        # print curr
        # print next

        # prev_curr = getNDistances(prev, curr)
        try:
            prev_curr = getNDistances(prev, curr)
        except AttributeError:
            prev_curr = 10

        try:
            curr_next = getNDistances(curr, next) - diff * 3.5
        except AttributeError:
            curr_next = 10

        # print "prev_curr", prev_curr
        # print "curr_next", curr_next



        if prev_curr < 4.0 and prev_curr < curr_next:
            # ok folks nothing to see here
            return False, ""
        elif curr_next < 4.0 and curr_next < prev_curr:
            # we need to go from
            # ...AVsvsndapv...  to
            # ...AvsvsndapV...
            # or from
            # ...AV--------...  to
            # ...A--------V...
            # so we move position V, or gap[0]
            return True, [gap[0] - 1, gap[1]]

    if aln_from[gap[0]] == aln_from[gap[1] + 1]:
        # print "gap[0] --- gap[1]+1"
        # print aln_to[0:gap[1]+2]
        # print aln_from[0:gap[1]+2]
        # print aln_to[0:gap[1]+2].replace("-", "")
        # we need to see if this is aligned correctly
        # find which ATOM position gap[1]+1 corresponds to
        # and see if it is closer to gap[1]+2 or to gap[0]-1
        truncated = aln_to[0:gap[1] + 2].replace("-", "")
        # print truncated
        pos = len(truncated) - 1  # <-- this is the position where we need to look in the structure
        # print pos
        # print truncated[pos]
        # print aln_from[pos]
        # print len(residue_list)

        if pos >= len(residue_list):
            print "reslist too short"
            return True, [-1, len(residue_list)]

        # is res[pos] closer to res[pos+1] or res[pos-1]?
        # if the former, then the alignment has been done correctly and there is nothing to
        # change here --> return False, ""
        # if the latter, then we need to move the gap one aminoacid.
        prev = residue_list[pos - 1]
        curr = residue_list[pos]
        try:
            next = residue_list[pos + 1]
        except IndexError:
            next = []

        # print prev
        # print curr
        # print next

        prev_curr = getNDistances(prev, curr)
        try:
            curr_next = getNDistances(curr, next)
        except AttributeError:
            curr_next = 10

        # print "prev_curr", prev_curr
        # print "curr_next", curr_next
        if curr_next < 4.0 and curr_next < prev_curr:
            # ok folks nothing to see here
            return False, ""
        elif prev_curr < 4.0 and prev_curr < curr_next:
            # we need to go from
            # ...AVsvsndapv...  to
            # ...AvsvsndapV...
            # or from
            # ...AV--------...  to
            # ...A--------V...
            # so we move position V, or gap[0]
            return True, [gap[0], gap[1] + 1]
    return False, ""


def find_gaps(pdb, aln):
    """
    go through a sequence with gaps and return the positions at which the gaps start
    and end. Also check gaps for sanity - gaps of the form ----G-g---g-g--g-- cannot
    be dealt with (see 1bf5_A_04)
    """
    gaps = []
    start = -1
    end = -1
    for pos in range(len(aln)):
        # print aln[pos], start, end
        if aln[pos] != "-" and start >= 0:
            end = pos - 1
            gaps.append([start, end])
            start = -1
            end = -1
        elif aln[pos] == "-" and start < 0:
            start = pos
            # print start, end
            # print gaps

    for gap in gaps:
        print pdb, str(gap[0]) + "-" + str(gap[1]), gap[1] - gap[0] + 1

    check_sanity(pdb, gaps)
    return gaps


def check_sanity(pdb, gaps):
    """
    check that the regions between gaps are not of length 1, as this is
    something that we should take a look at
    """
    for index in range(len(gaps) - 1):
        gap1 = gaps[index]
        gap2 = gaps[index + 1]
        if gap2[0] - gap1[1] == 2:
            print "check", pdb
    return


def get_mapping(query_aln, template_aln, query_start, template_start):
    mapping = dict()

    query_pos = query_start - 1
    template_pos = template_start - 1

    for i in range(len(query_aln)):
        if query_aln[i] != "-":
            query_pos += 1
        if template_aln[i] != "-":
            template_pos += 1

        if query_aln[i] != "-" and template_aln[i] != "-":
            mapping[query_pos] = template_pos

    return mapping


def interpret_map(alignment, combs, atom):
    """
    takes the combs sequence, the atom sequence and the mapping between them
    and returns a sequence like those returned by SEQATOMS, where the gaps in
    the ATOM sequence are replaced by the corresponding COMBS aminoacids in
    lowercase.
    """

    new_atom = ""
    for key in range(1, len(combs) + 1):
        try:
            new_atom = new_atom + atom[alignment[key] - 1]
        except KeyError:
            new_atom = new_atom + combs[key - 1].lower()

    return new_atom

def retrieve_seqatom_from_database(id, seqatoms_filename):

    url="http://www.bioinformatics.nl/tools/seqatoms/cgi-bin/getseqs?db=cath&id="+id
    seqatoms_filename = wget.download(url, out=seqatoms_filename)

    #if protein is not in seqatoms database
    #seqatoms will be empty
    seqatoms = io.read_fasta(seqatoms_filename)

    return seqatoms


def main():

    parser = argparse.ArgumentParser(description="Generate SEQATOM sequences from deprecated database or recompute")
    parser.add_argument("-c", "--combs", dest="combs", help="COMBS sequences")
    parser.add_argument("-a", "--atom", dest="atom", help="ATOM sequences of the same chains")
    parser.add_argument("-o", "--output", dest="output", help="output dir for new seq files")
    parser.add_argument("-p", "--pdb", dest="pdb", help="dompdb directory")

    args = parser.parse_args()
    combs_fasta_file    = args.combs
    atom_fasta_file     = args.atom
    out_dir             = args.output
    pdb_dir             = args.pdb

    #read fasta files
    combs_id, combs_fasta   = io.read_fasta(combs_fasta_file)
    atom_id, atom_fasta     = io.read_fasta(atom_fasta_file)

    seqatoms_filename   = out_dir + "/" + combs_id + ".seqatom"

    #try ot get seqatom fasta from seqatom server
    seqatoms = retrieve_seqatom_from_database(combs_id, seqatoms_filename)

    #if protein is not in seqatoms database, generate seqatoms
    if seqatoms is None:
        with open(seqatoms_filename, 'w') as out:
            out.write(">" + combs_id + "\n")
            ali = map_sequences(combs_id, combs_fasta, atom_fasta, pdb_dir)
            getseq = interpret_map(ali, combs_fasta, atom_fasta)
            out.write(getseq + "\n")
        print("Generated SEQATOMS file for " + combs_id)
    else:
        print("Retrieved " + combs_id + " from SEQATOMS Cath database")
    sys.stdout.flush()




if __name__ == "__main__":
    main()