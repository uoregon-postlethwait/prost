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
from prost.alignment import BBMapAlignmentExecutionHit
from prost.prost import ModificationThing

def hit(line):
    """Helper function. Converts a space separated representation of a SAM
    alignment line to a SamAlignmentExecutionHit (current BBMap version)."""
    return BBMapAlignmentExecutionHit('\t'.join(line.split()))

class TestModificationThing:
    """Test the functionality that helps to calculate the percent columns in
    the by_genomic_location tab (e.g. %_other_edited).  Specifically, we are
    testing ModificationThing.

    Examples:
        All examples assume that the BinStarter is a '1'.

        Example1:
          +-strand:
            1234567890123456789012345
            CCCCCTTTTTGGGGGAAAAACCCCC    main_hit : 25=
                      X    XX            mismatches wrt main_hit
         AAACCCCCTTTTTCGGGGTGAAACCCCCAA
         12345678901234567890123456789  member   : 1=1X10=1X4=2X10=
          X           X    XX           mismatches wrt genome

            mismatches_coords_ref_main   = (0,11,16,17)
            mismatches_coords_ref_member = (2,13,18,19)  <-- Not currently produced or
            offset = -2                                      used in the code, for educational
                                                             purposes only.
        Example2:
          +-strand:
            1234567890123456789012345
            CCCCCTTTTTGGGGGAAAAACCCCC    main_hit : 25=
                      X    XX            mismatches wrt main_hit
              CCCTTTTTCGGGGTGAAACCC
              123456789012345678901      member   : 8=1X4=2X6=
                      X    XX            mismatches wrt genome

            mismatches_coords_ref_main   = (11,16,17)
            mismatches_coords_ref_member = ( 9,14,15)
            offset = 2
    """

    ########################################
    ### MainHit (BinStarterHit) Fixtures ###
    ########################################

    @fixture
    def main_hit_A(self):
        # main hit A, plus strand
        return hit("CCCCCTTTTTGGGGGAAAAACCCCC 0 lg5 1000011 * 25= * * * * *")

    @fixture
    def main_hit_A_(self):
        # main hit A, minus strand
        return hit("CCCCCTTTTTGGGGGAAAAACCCCC 16 lg5 1000011 * 25= * * * * *")

    @fixture
    def main_hit_B(self):
        # main hit B, plus strand
        return hit("AAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 18= * * * * *")

    @fixture
    def main_hit_B_(self):
        # main hit B, minus strand
        return hit("AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *")

    ##################################
    ### Mismatch Coordinates tests ###
    ##################################

    # See Example1 and Example2 in the class' doc string above.

    @fixture
    def ex_1(self):
        # Example 1
        return hit("AACCCCCTTTTTCGGGGTGAAACCCCCAA 0 lg5 1000009 * 1=1X10=1X4=2X10= * * * * *")
    @fixture
    def ex_1_(self):
        # Example 1, minus strand
        return hit("AACCCCCTTTTTCGGGGTGAAACCCCCAA 16 lg5 1000009 * 10=2X4=1X10=1X1= * * * * *")

    @fixture
    def ex_2(self):
        # Example 2
        return hit("CCCTTTTTCGGGGTGAAACCC 0 lg5 1000013 * 8=1X4=2X8= * * * * *")
    @fixture
    def ex_2_(self):
        # Example 2, minus strand
        return hit("CCCTTTTTCGGGGTGAAACCC 16 lg5 1000013 * 8=2X4=1X8= * * * * *")

    def test_mismatches_coords_ex_1(self, main_hit_A, ex_1):
        m = ModificationThing(main_hit_A, ex_1)
        assert (0, 11, 16, 17) == m._mismatches_coords(main_hit_A, ex_1)

    def test_mismatches_coords_ex_1_minus(self, main_hit_A_, ex_1_):
        m = ModificationThing(main_hit_A_, ex_1_)
        assert (0, 11, 16, 17) == m._mismatches_coords(main_hit_A_, ex_1_)

    def test_mismatches_coords_ex_2(self, main_hit_A, ex_2):
        m = ModificationThing(main_hit_A, ex_2)
        assert (11, 16, 17) == m._mismatches_coords(main_hit_A, ex_2)

    def test_mismatches_coords_ex_2_minus(self, main_hit_A_, ex_2_):
        m = ModificationThing(main_hit_A_, ex_2_)
        assert (9, 14, 15) == m._mismatches_coords(main_hit_A_, ex_2_)

    #################################################
    ### Property tests: 3p_alt_cut, No seed_shift ###
    #################################################

    # Example3a: 3p_alt_cut (1nt 3p main overhang) (no seed shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 17=

    # Example3b: 3p_alt_cut (1nt 3p member overhang) (no seed shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 19=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_3a(self):
        return hit("AAAAAAAAAAAAAAAAA 0 lg5 1000011 * 17= * * * * *")
    @fixture
    def ex_3a_(self):
        return hit("AAAAAAAAAAAAAAAAA 16 lg5 1000012 * 17= * * * * *")

    @fixture
    def ex_3b(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 19= * * * * *")
    @fixture
    def ex_3b_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 19= * * * * *")

    def test_3p_alt_cut_ex_3a(self, main_hit_B, ex_3a):
        m = ModificationThing(main_hit_B, ex_3a)
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_ex_3a_minus(self, main_hit_B_, ex_3a_):
        m = ModificationThing(main_hit_B_, ex_3a_)
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_ex_3b(self, main_hit_B, ex_3b):
        m = ModificationThing(main_hit_B, ex_3b)
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_ex_3b_minus(self, main_hit_B_, ex_3b_):
        m = ModificationThing(main_hit_B_, ex_3b_)
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted
        assert not m.is_3p_mismatch

    ##################################################
    ### Property tests: 3p_mismatch, No seed_shift ###
    ##################################################

    # Example4a: 1nt 3p_mismatch (no 3p_alt_cut, no seed_shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 18=1X

    # Example4b: 3nt 3p_mismatch (no 3p_alt_cut, no seed_shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 18=3X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_4a(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 18=1X * * * * *")
    @fixture
    def ex_4a_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 1X18= * * * * *")

    @fixture
    def ex_4b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 18=3X * * * * *")
    @fixture
    def ex_4b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 16 lg5 1000008 * 3X18= * * * * *")

    def test_3p_mismatch_ex_4a(self, main_hit_B, ex_4a):
        m = ModificationThing(main_hit_B, ex_4a)
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_ex_4a_minus(self, main_hit_B_, ex_4a_):
        m = ModificationThing(main_hit_B_, ex_4a_)
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_ex_4b(self, main_hit_B, ex_4b):
        m = ModificationThing(main_hit_B, ex_4b)
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_ex_4b_minus(self, main_hit_B_, ex_4b_):
        m = ModificationThing(main_hit_B_, ex_4b_)
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted
        assert not m.is_3p_alt_cut

    #################################################################
    ### Property tests: 3p_alt_cut and 3p_mismatch, No seed_shift ###
    #################################################################

    # Example5a: 3nt 3p_mismatch and 3p_alt_cut (1nt 3p main overhang) (no seed shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 17=3X

    # Example5b: 3nt 3p_mismatch and 3p_alt_cut (1nt 3p member overhang) (no seed shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 19=3X

    # Example5c: 1nt 3p_mismatch and 3p_alt_cut (1nt 3p main overhang) (no seed shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 17=1X

    # Example5d: 5nt 3p_mismatch and 3p_alt_cut (5nt 3p main overhang) (no seed shift)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 13=5X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_5a(self):
        return hit("AAAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 17=3X * * * * *")
    @fixture
    def ex_5a_(self):
        return hit("AAAAAAAAAAAAAAAAAAAA 16 lg5 1000009 * 3X17= * * * * *")

    @fixture
    def ex_5b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 19=3X * * * * *")
    @fixture
    def ex_5b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAA 16 lg5 1000007 * 3X19= * * * * *")

    @fixture
    def ex_5c(self):
        return hit("AAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 17=1X * * * * *")
    @fixture
    def ex_5c_(self):
        return hit("AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 1X17= * * * * *")

    @fixture
    def ex_5d(self):
        return hit("AAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 13=5X * * * * *")
    @fixture
    def ex_5d_(self):
        return hit("AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 5X13= * * * * *")

    def test_3p_mismatch_and_3p_alt_cut_ex_5a(self, main_hit_B, ex_5a):
        m = ModificationThing(main_hit_B, ex_5a)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5a_minus(self, main_hit_B_, ex_5a_):
        m = ModificationThing(main_hit_B_, ex_5a_)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5b(self, main_hit_B, ex_5b):
        m = ModificationThing(main_hit_B, ex_5b)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5b_minus(self, main_hit_B_, ex_5b_):
        m = ModificationThing(main_hit_B_, ex_5b_)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5c(self, main_hit_B, ex_5c):
        m = ModificationThing(main_hit_B, ex_5c)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5c_minus(self, main_hit_B_, ex_5c_):
        m = ModificationThing(main_hit_B_, ex_5c_)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5d(self, main_hit_B, ex_5d):
        m = ModificationThing(main_hit_B, ex_5d)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    def test_3p_mismatch_and_3p_alt_cut_ex_5d_minus(self, main_hit_B_, ex_5d_):
        m = ModificationThing(main_hit_B_, ex_5d_)
        assert m.is_3p_mismatch
        assert m.is_3p_alt_cut
        assert not m.is_seed_shifted

    ###################################################
    ### Property tests: 3p_alt_cut, With seed_shift ###
    ###################################################

    # Example6a: 3p_alt_cut (1nt 3p main overhang) with seed_shift (2nt 5p main overhang)
    #     Offset = 2
    #     MainSeq: 18=
    #     MembSeq: 15=

    # Example6b: 3p_alt_cut (1nt 3p member overhang) with seed_shift (2nt 5p main overhang)
    #     Offset = 2
    #     MainSeq: 18=
    #     MembSeq: 17=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_6a(self):
        return hit("AAAAAAAAAAAAAAA 0 lg5 1000013 * 15= * * * * *")
    @fixture
    def ex_6a_(self):
        return hit("AAAAAAAAAAAAAAA 16 lg5 1000012 * 15= * * * * *")

    @fixture
    def ex_6b(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000013 * 17= * * * * *")
    @fixture
    def ex_6b_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 17= * * * * *")

    def test_3p_alt_cut_with_seed_shift_ex_6a(self, main_hit_B, ex_6a):
        m = ModificationThing(main_hit_B, ex_6a)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_with_seed_shift_ex_6a_minus(self, main_hit_B_, ex_6a_):
        # minus strand
        m = ModificationThing(main_hit_B_, ex_6a_)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_with_seed_shift_ex_6b(self, main_hit_B, ex_6b):
        m = ModificationThing(main_hit_B, ex_6b)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_with_seed_shift_ex_6b_minus(self, main_hit_B_, ex_6b_):
        # minus strand
        m = ModificationThing(main_hit_B_, ex_6b_)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    # Example7a: 3p_alt_cut (1nt 3p main overhang) with seed_shift (2nt 5p member overhang)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 19=

    # Example7b: 3p_alt_cut (1nt 3p member overhang) with seed_shift (2nt 5p member overhang)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 21=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_7a(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 19= * * * * *")
    @fixture
    def ex_7a_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000012 * 19= * * * * *")

    @fixture
    def ex_7b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 21= * * * * *")
    @fixture
    def ex_7b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 21= * * * * *")

    def test_3p_alt_cut_with_seed_shift_ex_7a(self, main_hit_B, ex_7a):
        m = ModificationThing(main_hit_B, ex_7a)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_with_seed_shift_ex_7a_minus(self, main_hit_B_, ex_7a_):
        # minus strand
        m = ModificationThing(main_hit_B_, ex_7a_)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_with_seed_shift_ex_7b(self, main_hit_B, ex_7b):
        m = ModificationThing(main_hit_B, ex_7b)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    def test_3p_alt_cut_with_seed_shift_ex_7b_minus(self, main_hit_B_, ex_7b_):
        # minus strand
        m = ModificationThing(main_hit_B_, ex_7b_)
        assert m.is_3p_alt_cut
        assert m.is_seed_shifted
        assert not m.is_3p_mismatch

    ####################################################
    ### Property tests: 3p_mismatch, With seed_shift ###
    ####################################################

    # Example8a: 1nt 3p_mismatch (no 3p_alt_cut) with seed_shift (2nt 5p main overhang)
    #     Offset = 2
    #     MainSeq: 18=
    #     MembSeq: 16=1X

    # Example8b: 3nt 3p_mismatch (no 3p_alt_cut) with seed_shift (2nt 5p main overhang)
    #     Offset = 2
    #     MainSeq: 18=
    #     MembSeq: 16=3X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_8a(self):
        return hit("AAAAAAAAAAAAAAAAA 0 lg5 1000013 * 16=1X * * * * *")
    @fixture
    def ex_8a_(self):
        return hit("AAAAAAAAAAAAAAAAA 16 lg5 1000010 * 1X16= * * * * *")

    @fixture
    def ex_8b(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000013 * 16=3X * * * * *")

    @fixture
    def ex_8b_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000008 * 3X16= * * * * *")

    def test_3p_mismatch_with_seed_shift_ex_8a(self, main_hit_B, ex_8a):
        m = ModificationThing(main_hit_B, ex_8a)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_with_seed_shift_ex_8a_minus(self, main_hit_B_, ex_8a_):
        m = ModificationThing(main_hit_B_, ex_8a_)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_with_seed_shift_ex_8b(self, main_hit_B, ex_8b):
        m = ModificationThing(main_hit_B, ex_8b)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_with_seed_shift_ex_8b_minus(self, main_hit_B_, ex_8b_):
        m = ModificationThing(main_hit_B_, ex_8b_)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    # Example9a: 1nt 3p_mismatch (no 3p_alt_cut) with seed_shift (2nt 5p member overhang)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 20=1X

    # Example9b: 3nt 3p_mismatch (no 3p_alt_cut) with seed_shift (2nt 5p member overhang)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 20=3X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_9a(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 20=1X * * * * *")
    @fixture
    def ex_9a_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 1X20= * * * * *")

    @fixture
    def ex_9b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 20=3X * * * * *")

    @fixture
    def ex_9b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAAA 16 lg5 1000008 * 3X20= * * * * *")

    def test_3p_mismatch_with_seed_shift_ex_9a(self, main_hit_B, ex_9a):
        m = ModificationThing(main_hit_B, ex_9a)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_with_seed_shift_ex_9a_minus(self, main_hit_B_, ex_9a_):
        m = ModificationThing(main_hit_B_, ex_9a_)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_with_seed_shift_ex_9b(self, main_hit_B, ex_9b):
        m = ModificationThing(main_hit_B, ex_9b)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    def test_3p_mismatch_with_seed_shift_ex_9b_minus(self, main_hit_B_, ex_9b_):
        m = ModificationThing(main_hit_B_, ex_9b_)
        assert m.is_3p_mismatch
        assert m.is_seed_shifted
        assert not m.is_3p_alt_cut

    #################################################################################
    ### Property tests: Other Edited, With and Without seed_shift, No 3p_mismatch ###
    #################################################################################

    # Example20a: other edited at nt17 (no seed shift, no 3p_mismatch)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 16=1X1=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_20a(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 16=1X1= * * * * *")
    @fixture
    def ex_20a_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 1=1X16= * * * * *")

    def test_other_edited_n17_no_seed_shift_ex_20a(self, main_hit_B, ex_20a):
        m = ModificationThing(main_hit_B, ex_20a)
        assert m.is_other_edited
        assert not m.is_3p_mismatch
        assert not m.is_seed_shifted

    def test_other_edited_n17_no_seed_shift_ex_20a_minus(self, main_hit_B_, ex_20a_):
        m = ModificationThing(main_hit_B_, ex_20a_)
        assert m.is_other_edited
        assert not m.is_3p_mismatch
        assert not m.is_seed_shifted

    # Example20b: other edited at nt17 with seed_shift (2nt 5p main overhang) (no 3p_mismatch)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 18=1X1=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_20b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 18=1X1= * * * * *")
    @fixture
    def ex_20b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 1=1X18= * * * * *")

    def test_other_edited_n17_with_seed_shift_ex_20b(self, main_hit_B, ex_20b):
        m = ModificationThing(main_hit_B, ex_20b)
        assert m.is_other_edited
        assert not m.is_3p_mismatch
        assert m.is_seed_shifted

    def test_other_edited_n17_with_seed_shift_ex_20b_minus(self, main_hit_B_, ex_20b_):
        m = ModificationThing(main_hit_B_, ex_20b_)
        assert m.is_other_edited
        assert not m.is_3p_mismatch
        assert m.is_seed_shifted

    ##############################################################################
    ### Property tests: Other Edited, 3p_mismatch, With and Without seed_shift ###
    ##############################################################################

    # Example21a: other edited at nt17 (no seed shift) and 3p_mismatch
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 16=1X1=1X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_21a(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 16=1X1=1X * * * * *")
    @fixture
    def ex_21a_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 1X1=1X16= * * * * *")

    def test_other_edited_n17_and_3p_mismatch_no_seed_shift_ex_21a(self, main_hit_B, ex_21a):
        m = ModificationThing(main_hit_B, ex_21a)
        assert m.is_other_edited
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted

    def test_other_edited_n17_and_3p_mismatch_no_seed_shift_ex_21a_minus(self, main_hit_B_, ex_21a_):
        m = ModificationThing(main_hit_B_, ex_21a_)
        assert m.is_other_edited
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted

    # Example21b: other edited at nt17 with seed_shift (2nt 5p main overhang) and 3p_mismatch
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 18=1X1=1X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_21b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 18=1X1=1X * * * * *")
    @fixture
    def ex_21b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 16 lg5 1000010 * 1X1=1X18= * * * * *")

    def test_other_edited_n17_and_3p_mismatch_with_seed_shift_ex_21b(self, main_hit_B, ex_21b):
        m = ModificationThing(main_hit_B, ex_21b)
        assert m.is_other_edited
        assert m.is_3p_mismatch
        assert m.is_seed_shifted

    def test_other_edited_n17_and_3p_mismatch_with_seed_shift_ex_21b_minus(self, main_hit_B_, ex_21b_):
        m = ModificationThing(main_hit_B_, ex_21b_)
        assert m.is_other_edited
        assert m.is_3p_mismatch
        assert m.is_seed_shifted


    #####################################################################################
    ### Property tests: Not Other Edited, With and Without seed_shift, No 3p_mismatch ###
    #####################################################################################

    # Example22a: NOT other edited (3p-suppl-edited at nt16) (no seed shift) (no 3p_mismatch)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 15=1X1=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_22a(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 15=1X1= * * * * *")
    @fixture
    def ex_22a_(self):
        return hit("AAAAAAAAAAAAAAAAAAA 16 lg5 1000012 * 1=1X15= * * * * *")

    def test_not_other_edited_n16_no_seed_shift_ex_22a(self, main_hit_B, ex_22a):
        m = ModificationThing(main_hit_B, ex_22a)
        assert not m.is_other_edited
        assert not m.is_3p_mismatch
        assert not m.is_seed_shifted

    def test_not_other_edited_n16_no_seed_shift_ex_22a_minus(self, main_hit_B_, ex_22a_):
        m = ModificationThing(main_hit_B_, ex_22a_)
        assert not m.is_other_edited
        assert not m.is_3p_mismatch
        assert not m.is_seed_shifted

    # Example22b: NOT other edited (3p-suppl-edited at nt16) with seed_shift (2nt 5p main overhang) (no 3p_mismatch)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 17=1X1=

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_22b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 17=1X1= * * * * *")
    @fixture
    def ex_22b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAA 16 lg5 1000012 * 1=1X17= * * * * *")

    def test_not_other_edited_n16_seed_shift_ex_22b(self, main_hit_B, ex_22b):
        m = ModificationThing(main_hit_B, ex_22b)
        assert not m.is_other_edited
        assert not m.is_3p_mismatch
        assert m.is_seed_shifted

    def test_not_other_edited_n16_seed_shift_ex_22b_minus(self, main_hit_B_, ex_22b_):
        m = ModificationThing(main_hit_B_, ex_22b_)
        assert not m.is_other_edited
        assert not m.is_3p_mismatch
        assert m.is_seed_shifted

    ######################################################################################
    ### Property tests: Not Other Edited, With and Without seed_shift, yes 3p_mismatch ###
    ######################################################################################

    # Example23a: NOT other edited (3p-suppl-edited at nt16) (no seed shift) (yes 3p_mismatch)
    #     Offset = 0
    #     MainSeq: 18=
    #     MembSeq: 15=1X1=3X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_23a(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 15=1X1=3X * * * * *")
    @fixture
    def ex_23a_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAA 16 lg5 1000009 * 3X1=1X15= * * * * *")

    def test_not_other_edited_n16_no_seed_shift_3p_mismatch_ex_23a(self, main_hit_B, ex_23a):
        m = ModificationThing(main_hit_B, ex_23a)
        assert not m.is_other_edited
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted

    def test_not_other_edited_n16_no_seed_shift_3p_mismatch_ex_23a_minus(self, main_hit_B_, ex_23a_):
        m = ModificationThing(main_hit_B_, ex_23a_)
        assert not m.is_other_edited
        assert m.is_3p_mismatch
        assert not m.is_seed_shifted

    # Example24b: NOT other edited (3p-suppl-edited at nt16) with seed_shift (2nt 5p main overhang) (yes 3p_mismatch)
    #     Offset = -2
    #     MainSeq: 18=
    #     MembSeq: 17=1X1=3X

    # recall: main_hit_B:  "AAAAAAAAAAAAAAAAAA  0 lg5 1000011 * 18= * * * * *"
    # recall: main_hit_B_: "AAAAAAAAAAAAAAAAAA 16 lg5 1000011 * 18= * * * * *"

    @fixture
    def ex_24b(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAAAA 0 lg5 1000009 * 17=1X1=3X * * * * *")
    @fixture
    def ex_24b_(self):
        return hit("AAAAAAAAAAAAAAAAAAAAAAAA 16 lg5 1000012 * 3X1=1X17= * * * * *")

    def test_not_other_edited_n16_seed_shift_3p_mismatch_ex_24b(self, main_hit_B, ex_24b):
        m = ModificationThing(main_hit_B, ex_24b)
        assert not m.is_other_edited
        assert m.is_3p_mismatch
        assert m.is_seed_shifted

    def test_not_other_edited_n16_seed_shift_3p_mismatch_ex_24b_minus(self, main_hit_B_, ex_24b_):
        m = ModificationThing(main_hit_B_, ex_24b_)
        assert not m.is_other_edited
        assert m.is_3p_mismatch
        assert m.is_seed_shifted


    ##########################################################
    ### Test 'disagreements() method of ModificationThing. ###
    ##########################################################

    @fixture
    def hit_yes_seed_shifted_yes_3p_alt_cut(self):
        return hit("AAAAAAAAAAAAAAAAAA 0 lg5 1000010 * 18= * * * * *")

    #def test_verify_hit_above_1(self, main_hit_B, hit_yes_seed_shifted_yes_3p_alt_cut):
    #    m = ModificationThing(main_hit_B, hit_yes_seed_shifted_yes_3p_alt_cut)
    #    assert m.is_seed_shifted
    #    assert m.is_3p_alt_cut
    #    assert not m.is_other_edited
    #    assert not m.is_3p_mismatch

    @fixture
    def hit_no_seed_shifted_yes_3p_alt_cut(self):
        return hit("AAAAAAAAAAAAAAAAA 0 lg5 1000011 * 17= * * * * *")

    #def test_verify_hit_above_2(self, main_hit_B, hit_no_seed_shifted_yes_3p_alt_cut):
    #    m = ModificationThing(main_hit_B, hit_no_seed_shifted_yes_3p_alt_cut)
    #    assert not m.is_seed_shifted
    #    assert m.is_3p_alt_cut
    #    assert not m.is_other_edited
    #    assert not m.is_3p_mismatch

    @fixture
    def hit_yes_seed_shifted_no_3p_alt_cut(self):
        return hit("AAAAAAAAAAAAAAAAAAA 0 lg5 1000010 * 19= * * * * *")

    #def test_verify_hit_above_3(self, main_hit_B, hit_yes_seed_shifted_no_3p_alt_cut):
    #    m = ModificationThing(main_hit_B, hit_yes_seed_shifted_no_3p_alt_cut)
    #    assert m.is_seed_shifted
    #    assert not m.is_3p_alt_cut
    #    assert not m.is_other_edited
    #    assert not m.is_3p_mismatch

    @fixture
    def hit_no_seed_shifted_no_3p_alt_cut(self):
        return hit("AAAAAAAAAAAAAAAAAA 0 lg5 1000011 * 18= * * * * *")

    #def test_verify_hit_above_4(self, main_hit_B, hit_no_seed_shifted_no_3p_alt_cut):
    #    m = ModificationThing(main_hit_B, hit_no_seed_shifted_no_3p_alt_cut)
    #    assert not m.is_seed_shifted
    #    assert not m.is_3p_alt_cut
    #    assert not m.is_other_edited
    #    assert not m.is_3p_mismatch

    @fixture
    def mt1_yes_seed_shifted_yes_3p_alt_cut(self,
            main_hit_B, hit_yes_seed_shifted_yes_3p_alt_cut):
        return ModificationThing(
            main_hit_B, hit_yes_seed_shifted_yes_3p_alt_cut)

    @fixture
    def mt2_no_seed_shifted_yes_3p_alt_cut(self,
            main_hit_B, hit_no_seed_shifted_yes_3p_alt_cut):
        return ModificationThing(
            main_hit_B, hit_no_seed_shifted_yes_3p_alt_cut)

    @fixture
    def mt3_yes_seed_shifted_no_3p_alt_cut(self,
            main_hit_B, hit_yes_seed_shifted_no_3p_alt_cut):
        return ModificationThing(
            main_hit_B, hit_yes_seed_shifted_no_3p_alt_cut)

    @fixture
    def mt4_no_seed_shifted_no_3p_alt_cut(self,
            main_hit_B, hit_no_seed_shifted_no_3p_alt_cut):
        return ModificationThing(
            main_hit_B, hit_no_seed_shifted_no_3p_alt_cut)

    def test_mt1_mt1_no_disagree(self,
            mt1_yes_seed_shifted_yes_3p_alt_cut):
        d = ModificationThing.disagreements([
                mt1_yes_seed_shifted_yes_3p_alt_cut,
                mt1_yes_seed_shifted_yes_3p_alt_cut])
        assert d == ()

    def test_mt2_mt2_no_disagree(self,
            mt2_no_seed_shifted_yes_3p_alt_cut):
        d = ModificationThing.disagreements([
                mt2_no_seed_shifted_yes_3p_alt_cut,
                mt2_no_seed_shifted_yes_3p_alt_cut])
        assert d == ()

    def test_mt1_mt2_disagree_seed_shifted(self,
            mt1_yes_seed_shifted_yes_3p_alt_cut,
            mt2_no_seed_shifted_yes_3p_alt_cut):
        d = ModificationThing.disagreements([
                mt1_yes_seed_shifted_yes_3p_alt_cut,
                mt2_no_seed_shifted_yes_3p_alt_cut])
        assert d == ('is_seed_shifted',)

    def test_mt1_mt3_disagree_3p_alt_cut(self,
            mt1_yes_seed_shifted_yes_3p_alt_cut,
            mt3_yes_seed_shifted_no_3p_alt_cut):
        d = ModificationThing.disagreements([
                mt1_yes_seed_shifted_yes_3p_alt_cut,
                mt3_yes_seed_shifted_no_3p_alt_cut])
        assert d == ('is_3p_alt_cut',)

    def test_mt1_mt2_mt3_disagree_seed_shifted_3p_alt_cut(self,
            mt1_yes_seed_shifted_yes_3p_alt_cut,
            mt2_no_seed_shifted_yes_3p_alt_cut,
            mt3_yes_seed_shifted_no_3p_alt_cut):
        d = ModificationThing.disagreements([
                mt1_yes_seed_shifted_yes_3p_alt_cut,
                mt2_no_seed_shifted_yes_3p_alt_cut,
                mt3_yes_seed_shifted_no_3p_alt_cut])
        assert d == ('is_3p_alt_cut', 'is_seed_shifted')

    def test_mt2_mt3_mt4_disagree_seed_shifted_3p_alt_cut(self,
            mt2_no_seed_shifted_yes_3p_alt_cut,
            mt3_yes_seed_shifted_no_3p_alt_cut,
            mt4_no_seed_shifted_no_3p_alt_cut):
        d = ModificationThing.disagreements([
                mt2_no_seed_shifted_yes_3p_alt_cut,
                mt3_yes_seed_shifted_no_3p_alt_cut,
                mt4_no_seed_shifted_no_3p_alt_cut])
        assert d == ('is_3p_alt_cut', 'is_seed_shifted')

# vim: softtabstop=4:shiftwidth=4:expandtab
