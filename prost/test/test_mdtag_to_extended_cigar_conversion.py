# Copyright (C) 2015-2019 University of Oregon
#
# This file is part of Prost!.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the License with this program.
#
# Written by Jason Sydes and Peter Batzel.
# Advised by Thomas Desvignes.

from pytest import fixture
import re

# For iterating!
import itertools

def get_mdtag_tokens(mdtag_str):
    """Parse the mdtag string into tokens.

    Returns:
        (int|String): A tuple of with elements of type int or String.  The
            As does the MD tag does in the SAM file, the mdtag_tokens tuple
            starts and ends with an int, and alternates between ints and
            Strings.  For example, an mdtag string of '14^TCCAC7A0A10' is
            parsed into (14, '^TCCAC', 7, 'A', 0, 'A', 10).

    """
    # mdtag_str = '14^TCCAC7AA10'
    # l = ['14', '^TCCAC', '7', 'AA', '10']
    l = re.split('([0-9]+)', mdtag_str)[1:-1]
    for i in xrange(0, len(l), 2):
        l[i] = int(l[i])
    # l = [14, '^TCCAC', 7, 'AA', 10]
    return tuple(l)


def get_cigar_tokens(cigar_str):
    """Parse the cigar string into tokens.

    Returns:
        ((int,String),): A tuple of (int, String) tuples, where each
            (int, String) tuple represents one cigar token. For example, a
            cigar string of "17=1X5=" will be parsed into
            ((17, '='), (1, 'X'), (5, '=')).
    """
    # cigar_str = "10=1X9="
    # l = ['10', '=', '1', 'X', '9', '=', '']
    l = re.split('([MIDNSHP=X])', cigar_str)
    # Yup, this works.
    return tuple((int(n),s) for n,s in itertools.izip(l[::2], l[1::2]))

def hit(cigar_str, mdtag_str, ext_cigar_str):
    """Returns a simulated hit, where mdtag and cigar string are expected to
    be converted to ext_cigar.  Returns the tuplized version of each."""
    return (
        get_cigar_tokens(cigar_str),
        get_mdtag_tokens(mdtag_str),
        get_cigar_tokens(ext_cigar_str)
    )

def convert(cigar_tokens, mdtag_tokens):

    # The new extended cigar tokens
    ext = []

    # make a copy of cigar_tokens and mdtag_tokens, as we'll be modifying them.
    cigar = [list(e) for e in cigar_tokens]
    mdtag = list(mdtag_tokens)

    print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    print "in cigar: {}".format(cigar)
    print "in mdtag: {}".format(mdtag)
    print "in ext  : {}".format(ext)

    while(len(cigar)):

        # Treat mdtag leading '0's as a noop
        if mdtag[0] == 0:
            mdtag.pop(0)

        # Convert INSH - Convert directly without "decrementing" mdtag
        if cigar[0][1] in ('I', 'N', 'S', 'H'):

            # Convert cigar to ext directly
            ext.append(tuple(cigar[0]))

            # Shift off the first token that is now converted
            cigar.pop(0)

        # Convert D - 'decrement' mdtag
        #       Just convert directly to ext_cigar, and "decrement" the
        #       MDtag by removing the '^' and the letters following the '^'.
        elif cigar[0][1] == 'D':
            # The Integer part of the cigar should equal the number of
            # letters following the '^' in the MDTag.
            assert(mdtag[0][0]) == '^'
            assert(len(mdtag[0]) - 1 == cigar[0][0])

            # Convert cigar to ext directly
            ext.append(tuple(cigar[0]))

            # Shift off the '^NNNN' from the mdtag, and (len, 'D') from cigar
            mdtag.pop(0)
            cigar.pop(0)

        # Convert 'M' (depends on MDtag)
        elif cigar[0][1] == 'M':

            # mdtag has a leading '0' or an NT.
            #   The second case (mdtag has a leading NT) is due to some aligners
            #   not following the SAM spec (e.g "3C4AT4" instead of "3C4A0T4").
            # Example cases:
            #       3C4AT4
            #       3C4A0T4
            #       0A0T0A55M
            if (mdtag[0] == 0 or (type(mdtag[0]) == str and mdtag[0].isalpha())):
                mismatches = 0
                while (mdtag[0] == 0 or (type(mdtag[0]) == str and mdtag[0].isalpha())):
                    print "DEBUG BBB, mdtag={}".format(mdtag)
                    if mdtag[0] == 0:
                        # shift off the leading 0, don't increment
                        mdtag.pop(0)

                    if (len(mdtag)):
                        # We assume that now the next token is just NTs. Assert!
                        assert(mdtag[0].isalpha())
                        assert(mdtag[0][0] != '^')

                        # Shift off leading NTs, increment count
                        mismatches += len(mdtag.pop(0))
                    else:
                        # That was the last '0'
                        break

                # Update old cigar, finalize extended cigar
                ext.append((mismatches, 'X'))
                cigar[0][0] -= mismatches
                if cigar[0][0] == 0:
                    cigar.pop(0)

            # int_of(cigar) <= int_of(mdtag)
            elif cigar[0][0] <= mdtag[0]:
                # Get the number of matching NTs and update mdtag, cigar, and ext
                count = cigar[0][0]
                mdtag[0] -= count
                ext.append((count, '='))
                cigar.pop(0)
                # End of the road?
                # NEEDS MORE EXAMINATION....
                #       Basically, please examine what happens at the termination of conversion...
                if (len(cigar) == 0):
                    assert(len(mdtag) == 1 and mdtag[0] == 0)
                elif (len(cigar) > 0):
                    # Hrm...
                    assert(len(mdtag) > 1)
                else:
                    raise Exception, 'Not possible...'

            # int_of(mdtag) < int_of(cigar)
            elif mdtag[0] < cigar[0][0]:
                # Get the number of matching NTs and update mdtag, cigar, and ext
                count = mdtag[0]
                mdtag[0] -= count
                cigar[0][0] -= count
                ext.append((count, '='))

            else:
                raise Exception, "Impossible to be here"

        print "during cigar: {}".format(cigar)
        print "during mdtag: {}".format(mdtag)
        print "during ext  : {}".format(ext)
        print

    return tuple(tuple(e) for e in ext)



class TestConversion:

    @fixture
    def hit_A(self):
        return hit('14M', '14', '14=')

    @fixture
    def hit_B(self):
        return hit('14M', '3C10', '3=1X10=')

    @fixture
    def hit_C1(self):
        return hit('14M', '3C4A0T4', '3=1X4=2X4=')

    @fixture
    def hit_C2(self):
        return hit('14M', '3C4AT4', '3=1X4=2X4=')

    @fixture
    def hit_D(self):
        return hit('14M', '8AT4', '8=2X4=')

    @fixture
    def hit_E(self):
        return hit('14M', '0C12T0', '1X12=1X')

    #@fixture
    #def hit_ZZZZZZZZZZZ(self):
    #    return hit('56M1D45M', '56^A45', '')





    @fixture
    def hit_F(self):
        """Adapted from http://kimplove.blogspot.com/2011/03/md-tag.html.
                     123456789012345678901234
        REF:         ATCGTAGCTAATTTGGACATCGGT
        READ:        ATCGTAGCTATTTTGG--ATCGGT
        MD TAG:      10        A5   ^AC6            Full: 10A5^AC6
        CIGAR:       16M             2D6M           Full: 16M2D6M
        EXTCIGAR:    10=       1X                   Full: 10=1X5=2D6=
                                5=   2D6=
        """
        return hit('16M2D6M', '10A5^AC6', '10=1X5=2D6=')

    @fixture
    def hit_G(self):
        """Adapted from http://kimplove.blogspot.com/2011/03/md-tag.html.
                     123456789012345678901234
        REF:         ATCGTAGCTAATTTGGACATCGGT       **
        READ:        ATCGTAGCTAATTTGGACATCGGT (ATCGTGGAGCTAATTTGGACATCGGT)
        MD TAG:      24                             Full: 24
        CIGAR:       5M   2I19M                     Full: 5M2I19M
        EXTCIGAR:    5=   2I19=                     Full: 5=2I19=

        REF:         AGCTAATTTGGACATCGGT       **
        READ:        AGCTAATTTGGACATCGGT (ATCGTGGAGCTAATTTGGACATCGGT)
        MD TAG:      19
        CIGAR:       2I19M
        EXTCIGAR:    2I19=
        """
        return hit('16M2D6M', '10A5^AC6', '10=1X5=2D6=')

    # Assertions noticed
    #       1. Len of reference = In MDTAG: Count of letters + sum(Integers, though definitely not those preceded by 'S' or 'H', and careful about 'P' too!)
    #       2. Len of reference =

    # So...
    #   1. I think that 'I's can be directly converted.
    #   2. I think that 'N's can maybe be directly converted (without examining the )

    # Example:
    #       See https://www.biostars.org/p/109333/#115106
    #
    # As an example, here's part of a bam file with a read pair containing a chimeric read. The top hit is soft clipped and the second top hit is hard clipped and marked as secondary by BWA (-M option).
    #
    # 20692128    97    viral_genome    21417    60    69M32S    chr7    101141242    0    TACATCTTCTCCCTCTCTCACGACACAAGAATTAGTCACATAGGGATGTTCTCGTAAATCTACATTATCTTACAAAAACATTTTTTAAAAATTTGCTAGGT    GGGGGGGGGGGGGGEGGEGGGGGGGGGFGGGGGGGGGGGGGEGFFGGGGGGGFGGFGGGGEGGGGGGGGGGGEGEFFGGGFEGGGGGFGCGGGFBGGGBG@    NM:i:4    MD:Z:6G34G6C5C14    AS:i:49    XS:i:0    SA:Z:chr7,101141091,+,66S35M,60,0;
    # 20692128    353    chr7    101141091    60    66H35M    =    101141242    252    ATCTTACAAAAACATTTTTTAAAAATTTGCTAGGT    GGGGGGEGEFFGGGFEGGGGGFGCGGGFBGGGBG@    NM:i:0    MD:Z:35    AS:i:35    XS:i:23    SA:Z:gi|224020395|ref|NC_001664.2|,21417,+,69M32S,60,4;
    # 20692128    145    chr7    101141242    60    101M    gi|224020395|ref|NC_001664.2|    21417    0    GCAACAGAGCGAGACCCTATATTCATGAGTGTTGCAATGAGCCAAGTAGTGGAGGTTGGCTTTTGAAGGCAGAAAAGGACTGAGAAAAGCTAACACAGAGA    FEGCGGGGGCGEFCDEEEEGGGGGGGGGGGGGGGEGGGGGGFGGGEGGG

    # So...
    #   * I: Convert directly from cigar to ext_cigar without decrementing
    #        MDtag int.
    #   * N: Same as 'I' (?)
    #   * S: Same as 'I' (?)
    #   * H: Same as 'I' (?)
    #   * P: ???
    #   * D: Works with MDTag:
    #           The Integer part of the cigar should equal the number of
    #           letters following the '^' in the MDTag (assert this!).
    #           Just convert directly to ext_cigar, and "decrement" the
    #           MDtag by removing the '^' and the letters following the '^'.
    #   * M: Depends on MDTag.
    #           if str_of(cigar) == 'M':
    #               if int_of(mdtag) < int_of(cigar):
    #
    #               elif int_of(cigar) < int_of(mdtag):
    #
    #               else: # int_of(cigar) == int_of(mdtag)



    @fixture
    def hit_G(self):
        return hit('', '', '')


    def test_A(self, hit_A):
        cigar, mdtag, ext = hit_A
        assert(convert(cigar, mdtag) == ext)

    def test_B(self, hit_B):
        _hit = hit_B
        cigar, mdtag, ext = _hit
        print "exp ext  : {}".format(_hit[2])
        ans = convert(cigar, mdtag)
        print "out ext  : {}".format(ans)
        assert(ans == ext)
        #assert(convert(cigar, mdtag) == ext)

    def test_C1(self, hit_C1):
        _hit = hit_C1
        cigar, mdtag, ext = _hit
        print "exp ext  : {}".format(_hit[2])
        ans = convert(cigar, mdtag)
        print "out ext  : {}".format(ans)
        assert(ans == ext)
        #assert(convert(cigar, mdtag) == ext)

    def test_C2(self, hit_C2):
        _hit = hit_C2
        cigar, mdtag, ext = _hit
        print "exp ext  : {}".format(_hit[2])
        ans = convert(cigar, mdtag)
        print "out ext  : {}".format(ans)
        assert(ans == ext)
        #assert(convert(cigar, mdtag) == ext)

    def test_D(self, hit_D):
        _hit = hit_D
        cigar, mdtag, ext = _hit
        print "exp ext  : {}".format(_hit[2])
        ans = convert(cigar, mdtag)
        print "out ext  : {}".format(ans)
        assert(ans == ext)
        #assert(convert(cigar, mdtag) == ext)

    def test_E(self, hit_E):
        _hit = hit_E
        cigar, mdtag, ext = _hit
        print "exp ext  : {}".format(_hit[2])
        ans = convert(cigar, mdtag)
        print "out ext  : {}".format(ans)
        assert(ans == ext)
        #assert(convert(cigar, mdtag) == ext)

    def test_F(self, hit_F):
        _hit = hit_F
        cigar, mdtag, ext = _hit
        print "exp ext  : {}".format(_hit[2])
        ans = convert(cigar, mdtag)
        print "out ext  : {}".format(ans)
        assert(ans == ext)
        #assert(convert(cigar, mdtag) == ext)


# vim: softtabstop=4:shiftwidth=4:expandtab
