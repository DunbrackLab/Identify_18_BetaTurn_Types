# Identify_18_BetaTurn_Types
Identifies beta turns in proteins according to Shapovalov, Vucetic, and Dunbrack, PLOSCompBio 2019. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006844

Requires mkdssp from DSSP4.5 : https://github.com/PDB-REDO/dssp
It uses the command line version of DSSP4.5 not the python module. The way to install it in /usr/local/bin/mkdssp is:

    git clone https://github.com/PDB-REDO/dssp.git
    cd dssp
    cmake -S . -B build
    cmake --build build
    cmake --install build

# Usage

    python3 Identify_18_BetaTurn_Types.py filename.cif > outfilename
    python3 Identify_18_BetaTurn_Types.py filename.pdb > outfilename


# Output

    python3 Identify_18_BetaTurn_Types.py 3q4z.cif
    turn  num chn  res1 res4    seq  dssp    type  prev_name          Dist DistAng CA1-CA4     omega2    phi2    psi2   omega3    phi3     psi3  omega4   filename
    turn    1 A     129  132    ALED CGGG    AD    I                0.1491   22.26    5.59     177.38  -49.71  -46.55   176.06  -46.67   -13.00  178.89   3e5a
    turn    2 A     130  133    LEDF GGGE    AD    I                0.0383   11.23    5.71     176.06  -46.67  -13.00   178.89  -82.81    -9.11 -177.51   3e5a
    turn    3 A     141  144    KGKF SGGG    pD    II'              0.0670   14.88    5.57     170.63   54.60 -132.45   179.20  -57.92   -15.22  176.87   3e5a
    turn    4 A     142  145    GKFG GGGT    AD    I                0.0440   12.05    5.23     179.20  -57.92  -15.22   176.87 -103.69    10.00 -173.87   3e5a
    turn    5 A     144  147    FGNV GTTE    AD    I                0.0549   13.45    5.92    -173.87  -62.87    0.73  -156.65 -102.92   -12.84 -164.28   3e5a
    turn    6 A     152  155    EKQS ETTT    AD    I                0.0470   12.44    5.82    -166.98  -58.12  -32.18  -171.95  -91.98   -46.88  179.84   3e5a
    turn    7 A     153  156    KQSK TTTC    AD    I                0.0733   15.56    5.38    -171.95  -91.98  -46.88   179.84  -79.56   -23.87 -173.57   3e5a
    turn    8 A     190  193    HPNI CTTB    AD    I                0.1066   18.79    5.26     176.32  -51.59  -36.82  -178.19 -105.12    27.26  176.33   3e5a
    turn    9 A     202  205    DATR CSSE    AD    I                0.0405   11.56    5.34    -175.29  -53.82  -31.56  -170.02 -123.16   -24.81 -171.69   3e5a
    turn   10 A     213  216    APLG CTTC    AD    I                0.0317   10.21    6.16    -168.14  -69.94  -12.27  -174.51  -98.76     3.43  172.87   3e5a
    turn   11 A     248  251    HSKR HTTT    AD    I                0.0302    9.97    4.79     172.32  -61.66  -10.94  -176.61 -109.75    -2.38 -179.14   3e5a
    turn   12 A     254  257    HRDI CCCC    dD    new              0.1841   24.77    6.82     166.89   78.66   11.03   171.83 -151.58    44.81  176.55   3e5a
    turn   13 A     258  261    KPEN SGGG    AD    I                0.0531   13.24    5.58    -168.83  -51.82  -34.16   179.70  -71.36   -29.86  178.93   3e5a
    turn   14 A     259  262    PENL GGGE    AD    I                0.0570   13.71    5.39     179.70  -71.36  -29.86   178.93  -80.98    11.92 -179.12   3e5a
    turn   15 A     265  268    GSAG CTTS    AD    I                0.1287   20.67    5.11     176.88  -32.72  -55.49  -172.09  -90.42     9.97 -177.37   3e5a
    turn   16 A     275  278    FGWS CTTC    AD    I                0.0393   11.38    4.99     179.12  -57.36  -21.14   174.95 -116.79     2.26 -176.09   3e5a
    turn   17 A     281  284    APSS CSSS    AD    I                0.1712   23.88    6.09    -170.22  -78.63  -32.65  -179.43 -132.21   -70.39 -178.60   3e5a
    turn   18 A     292  295    TLDY CGGG    AD    I                0.1458   22.02    5.40    -175.45  -28.17  -59.24  -168.10  -74.97    -7.03 -173.65   3e5a
    turn   19 A     293  296    LDYL GGGC    AD    I                0.0386   11.27    5.61    -168.10  -74.97   -7.03  -173.65 -113.31    -2.33 -177.08   3e5a
    turn   20 A     307  310    DEKV CTTH    AD    I                0.0469   12.44    6.42    -179.02  -61.71  -14.04  -178.29  -69.54   -14.54  170.42   3e5a
    turn   21 A     327  330    PPFE CTTC    AZ    new_prev_VIII    0.0788   16.14    6.03    -171.03  -72.30  -19.14  -174.07 -114.37    23.77 -178.54   3e5a
    turn   22 A     331  334    ANTY CSSH    AB1   new_prev_VIII    0.0414   11.68    6.77     177.27  -69.86  -46.55  -176.60 -112.91   147.79  173.15   3e5a
    turn   23 A     349  352    PDFV CTTS    AD    I                0.0482   12.60    5.78    -175.41  -52.09  -33.03  -177.37  -69.67   -15.52 -175.53   3e5a
    turn   24 A     365  368    KHNP CSSG    AB2   VIII             0.0704   15.25    6.35     177.12  -55.50  -44.16   178.27  -84.08   115.85 -176.81   3e5a
    turn   25 A     367  370    NPSQ SGGG    AD    I                0.0815   16.41    5.35    -176.81  -62.49  -34.82   175.12  -57.76   -26.21 -179.37   3e5a
    turn   26 A     368  371    PSQR GGGS    AD    I                0.0231    8.71    5.33     175.12  -57.76  -26.21  -179.37  -83.35    -8.55 -172.47   3e5a
    turn   27 B      18   21    NFSS CTTC    AG    new_prev_VIII    0.2072   26.31    6.56    -170.89  -73.24  -26.65   179.65  -83.56    32.22  126.36   3e5a
   
The output gives the residues of each beta turn (res1-res4), the sequence of the 4-residue turn, the dssp assignment of the turn ("C" is coil when DSSP does not report a secondary structure letter), the new turn type (e.g. "AD", "Pa", "Pd", etc.), the classical turn type (e.g. "I", "II"; "new_prev_VIII" indicates it is a new turn type but would formerly have been close to a type VIII turn). Then the distance in our metric, which the average of D=2(1-cos(d_theta)), where theta are the angles given on each line: omega2, phi2, psi2, omega3, phi3, psi3, omega4, which connect CA of the first residue to CA of the 4th residue of each turn. d_theta is the difference between the PDB dihedral angle and the medoid for that turn type, determined by the clustering described in Shapovalov et al. Following the distance in our metric ("Dist"), is the distance in degrees ("DistAng"), which is just the average angle distance converted back into an angle in degrees (theta = arccos(1 - D/2)). The CA1-CA4 distance is given next followed by all the dihedral angles.
