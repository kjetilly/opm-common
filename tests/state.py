import unittest
import sunbeam

class TestState(unittest.TestCase):

    spe3 = None
    FAULTS_DECK = """
RUNSPEC

DIMENS
 10 10 10 /
GRID
DX
1000*0.25 /
DY
1000*0.25 /
DZ
1000*0.25 /
TOPS
100*0.25 /
FAULTS
  'F1'  1  1  1  4   1  4  'X' /
  'F2'  5  5  1  4   1  4  'X-' /
/
MULTFLT
  'F1' 0.50 /
  'F2' 0.50 /
/
EDIT
MULTFLT /
  'F2' 0.25 /
/
OIL

GAS

TITLE
The title

START
8 MAR 1998 /

PROPS
REGIONS
SWAT
1000*1 /
SATNUM
1000*2 /
\
"""

    def setUp(self):
        if self.spe3 is None:
            self.spe3 = sunbeam.parse('spe3/SPE3CASE1.DATA')
        self.state = self.spe3

    def test_repr_title(self):
        self.assertTrue('EclipseState' in repr(self.state))
        self.assertEqual('SPE 3 - CASE 1', self.state.title)

    def test_state_nnc(self):
        self.assertFalse(self.state.has_input_nnc())

    def test_grid(self):
        grid = self.state.grid()
        self.assertTrue('EclipseGrid' in repr(grid))
        self.assertEqual(9, grid.getNX())
        self.assertEqual(9, grid.getNY())
        self.assertEqual(4, grid.getNZ())
        self.assertEqual(9*9*4, grid.nactive())
        self.assertEqual(9*9*4, grid.cartesianSize())
        g,i,j,k = 295,7,5,3
        self.assertEqual(g, grid.globalIndex(i,j,k))
        self.assertEqual((i,j,k), grid.getIJK(g))

    def test_config(self):
        cfg = self.spe3.cfg()
        smry = cfg.summary()
        self.assertTrue('EclipseConfig' in repr(cfg))
        self.assertTrue('SummaryConfig' in repr(smry))
        self.assertTrue('WOPR' in smry) # hasKeyword
        self.assertFalse('NONO' in smry) # hasKeyword

        init = cfg.init()
        self.assertTrue(init.hasEquil())
        self.assertFalse(init.restartRequested())
        self.assertEqual(0, init.getRestartStep())

        rst = cfg.restart()
        self.assertFalse(rst.getWriteRestartFile(0))
        self.assertEqual(8, rst.getFirstRestartStep())

        sim = cfg.simulation()
        self.assertFalse(sim.hasThresholdPressure())
        self.assertFalse(sim.useCPR())
        self.assertTrue(sim.hasDISGAS())
        self.assertTrue(sim.hasVAPOIL())

    def test_faults(self):
        self.assertEquals([], self.spe3.faultNames())
        self.assertEquals({}, self.spe3.faults())
        faultdeck = sunbeam.parse(self.FAULTS_DECK)
        self.assertEqual(['F1', 'F2'], faultdeck.faultNames())
        # 'F2'  5  5  1  4   1  4  'X-' / \n"
        f2 = faultdeck.faultFaces('F2')
        self.assertTrue((4,0,0,'X-') in f2)
        self.assertFalse((3,0,0,'X-') in f2)
