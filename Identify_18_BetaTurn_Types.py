# BetaTurnIdentifier.py
# usage:
#     python3 Identify_18_BetaTurn_Types.py filename.cif > outfilename
#     python3 Identify_18_BetaTurn_Types.py filename.pdb > outfilename
# Maxim Shapovalov, Slobadan Vucetic, Roland Dunbrack
# Sept 5, 2025
# Correspondence: roland.dunbrack@fccc.edu

# Please cite:
# Shapovalov M, Vucetic S, Dunbrack RL Jr (2019).
# A new clustering and nomenclature for beta turns derived from high-resolution protein structures.
# PLoS Comput Biol 15(3): e1006844.
# https://doi.org/10.1371/journal.pcbi.1006844

import gzip
import math
import sys, os
from Bio.PDB import *
import subprocess
import shlex
import tempfile


TURN_DIST_CUTOFF = 0.2359   # this it the cutoff for the average D (distance from medoid of each cluster) over 7 dihedral angles where D=2(1-cos(dtheta)).

def angle_distance(angle1, angle2):
    delta = angle1 - angle2
    delta_rad = math.radians(delta)
    delta_cos = math.cos(delta_rad)
    distance = 2.0*(1.0-delta_cos)
    return distance

def find_theta_in_degrees(D):
    if D < 0 or D > 4:
        raise ValueError("D must be in the range [0, 4] for a real solution.")
    
    cos_theta = 1.0 - D / 2.0
    theta_rad = math.acos(cos_theta)
    theta_deg = math.degrees(theta_rad)
    
    return theta_deg


def getres1(res3):
    three21={   'ALA' : 'A', 'CYS' : 'C', 'ASP' : 'D', 'GLU' : 'E', 'PHE' : 'F', 'GLY' : 'G',
	        'HIS' : 'H', 'ILE' : 'I', 'LYS' : 'K', 'LEU' : 'L', 'MET' : 'M', 'ASN' : 'N',
	        'PRO' : 'P', 'GLN' : 'Q', 'ARG' : 'R', 'SER' : 'S', 'THR' : 'T', 'VAL' : 'V', 
	        'TRP' : 'W', 'TYR' : 'Y'}
    if res3 in three21: return three21[res3]
    else:
        return 'X'

def define_turn_library():
    # definitions for 18 beta turn types from Shapovalov, Vucetic, and Dunbrack, 2019.
    result={}
    result['AD']=     {'no_by_size':  1,   'cluster_size':  6413,  'frequency':  0.49217,  'bturn_name':  'AD',     'prev_name':  'I',
                       'median_below_7A_ca1_ca4':  5.48,  'mean_below_7A_ca1_ca4':  5.52,  'median_any_ca1_ca4':  5.51,  'mean_any_ca1_ca4':  5.65,
                       'mode_omega2':  185.12,  'mode_phi2':  -62.25,   'mode_psi2':  -23.48,   'mode_omega3':  181.88,  'mode_phi3':  -96.25,
                       'mode_psi3':  -2.44,   'mode_omega4':  179.02,  'mode_aa1':  'D',  'mode_aa2':  'P',  'mode_aa3':  'D',  'mode_aa4':  'G',
                       'mode_ss1':  'C',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'C',  'mode_pdb_id':  '2EAB',  'mode_chain_id':  'A',
                       'mode_res1_id':  598,  'medoid_omega2':  184.88,  'medoid_phi2':  -65.27,   'medoid_psi2':  -23.52,   'medoid_omega3':  182.83,
                       'medoid_phi3':  -98.73,   'medoid_psi3':  -18.06,  'medoid_omega4':  181.10,  'medoid_aa1':  'P',  'medoid_aa2':  'Y',
                       'medoid_aa3':  'A',  'medoid_aa4':  'R',  'medoid_ss1':  'T',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'G',
                       'medoid_pdb_id':  '5A71',  'medoid_chain_id':  'A',  'medoid_res1_id':  299}
    result['Pd']=     {'no_by_size':  2,   'cluster_size':  1556,  'frequency':  0.11942,  'bturn_name':  'Pd',     'prev_name':  'II',
                       'median_below_7A_ca1_ca4':  5.61,  'mean_below_7A_ca1_ca4':  5.62,  'median_any_ca1_ca4':  5.61,  'mean_any_ca1_ca4':  5.68,
                       'mode_omega2':  179.57,  'mode_phi2':  -55.26,   'mode_psi2':  133.10,   'mode_omega3':  179.60,  'mode_phi3':  91.18,
                       'mode_psi3':  -6.22,   'mode_omega4':  180.91,  'mode_aa1':  'T',  'mode_aa2':  'K',  'mode_aa3':  'G',  'mode_aa4':  'T',
                       'mode_ss1':  'C',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'C',  'mode_pdb_id':  '1W23',  'mode_chain_id':  'A',
                       'mode_res1_id':  89,   'medoid_omega2':  178.69,  'medoid_phi2':  -55.28,   'medoid_psi2':  132.76,   'medoid_omega3':  175.53,
                       'medoid_phi3':  82.38,    'medoid_psi3':  -2.88,   'medoid_omega4':  179.60,  'medoid_aa1':  'P',  'medoid_aa2':  'D',
                       'medoid_aa3':  'G',  'medoid_aa4':  'D',  'medoid_ss1':  'E',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '3QZB',  'medoid_chain_id':  'A',  'medoid_res1_id':  62}
    result['Pa']=     {'no_by_size':  3,   'cluster_size':  802,   'frequency':  0.06155,  'bturn_name':  'Pa',     'prev_name':  'new_prev_II',
                       'median_below_7A_ca1_ca4':  6.00,  'mean_below_7A_ca1_ca4':  5.96,  'median_any_ca1_ca4':  6.12,  'mean_any_ca1_ca4':  6.18,
                       'mode_omega2':  176.22,  'mode_phi2':  -59.93,   'mode_psi2':  135.21,   'mode_omega3':  178.77,  'mode_phi3':  58.69,
                       'mode_psi3':  27.54,   'mode_omega4':  177.53,  'mode_aa1':  'G',  'mode_aa2':  'A',  'mode_aa3':  'L',  'mode_aa4':  'D',
                       'mode_ss1':  'C',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'H',  'mode_pdb_id':  '2FFY',  'mode_chain_id':  'A',
                       'mode_res1_id':  214,  'medoid_omega2':  180.97,  'medoid_phi2':  -79.76,   'medoid_psi2':  149.11,   'medoid_omega3':  170.92,
                       'medoid_phi3':  61.47,    'medoid_psi3':  32.21,   'medoid_omega4':  180.14,  'medoid_aa1':  'G',  'medoid_aa2':  'D',
                       'medoid_aa3':  'N',  'medoid_aa4':  'V',  'medoid_ss1':  'C',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '3KWE',  'medoid_chain_id':  'A',  'medoid_res1_id':  95}
    result['ad']=     {'no_by_size':  4,   'cluster_size':  748,   'frequency':  0.05741,  'bturn_name':  'ad',     'prev_name':  "I'",
                       'median_below_7A_ca1_ca4':  5.44,  'mean_below_7A_ca1_ca4':  5.52,  'median_any_ca1_ca4':  5.45,  'mean_any_ca1_ca4':  5.58,
                       'mode_omega2':  182.81,  'mode_phi2':  47.08,    'mode_psi2':  44.65,    'mode_omega3':  174.36,  'mode_phi3':  82.66,
                       'mode_psi3':  1.09,    'mode_omega4':  182.92,  'mode_aa1':  'Y',  'mode_aa2':  'K',  'mode_aa3':  'G',  'mode_aa4':  'R',
                       'mode_ss1':  'E',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'E',  'mode_pdb_id':  '3AWU',  'mode_chain_id':  'B',
                       'mode_res1_id':  49,   'medoid_omega2':  177.14,  'medoid_phi2':  54.08,    'medoid_psi2':  37.14,    'medoid_omega3':  174.19,
                       'medoid_phi3':  77.07,    'medoid_psi3':  7.66,    'medoid_omega4':  182.45,  'medoid_aa1':  'F',  'medoid_aa2':  'D',
                       'medoid_aa3':  'G',  'medoid_aa4':  'K',  'medoid_ss1':  'E',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'E',
                       'medoid_pdb_id':  '3GNE',  'medoid_chain_id':  'A',  'medoid_res1_id':  34}
    result['AB1']=    {'no_by_size':  5,   'cluster_size':  648,   'frequency':  0.04973,  'bturn_name':  'AB1',    'prev_name':  'new_prev_VIII',
                       'median_below_7A_ca1_ca4':  6.72,  'mean_below_7A_ca1_ca4':  6.62,  'median_any_ca1_ca4':  7.52,  'mean_any_ca1_ca4':  7.50,
                       'mode_omega2':  184.02,  'mode_phi2':  -67.20,   'mode_psi2':  -30.96,   'mode_omega3':  172.54,  'mode_phi3':  -136.03,
                       'mode_psi3':  162.29,  'mode_omega4':  182.93,  'mode_aa1':  'P',  'mode_aa2':  'R',  'mode_aa3':  'V',  'mode_aa4':  'P',
                       'mode_ss1':  'C',  'mode_ss2':  'S',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '4ATE',  'mode_chain_id':  'A',
                       'mode_res1_id':  219,  'medoid_omega2':  178.29,  'medoid_phi2':  -76.75,   'medoid_psi2':  -33.07,   'medoid_omega3':  180.05,
                       'medoid_phi3':  -137.67,  'medoid_psi3':  155.99,  'medoid_omega4':  179.53,  'medoid_aa1':  'D',  'medoid_aa2':  'F',
                       'medoid_aa3':  'Y',  'medoid_aa4':  'G',  'medoid_ss1':  'C',  'medoid_ss2':  'S',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '3ARX',  'medoid_chain_id':  'A',  'medoid_res1_id':  392}
    result['AZ']=     {'no_by_size':  6,   'cluster_size':  625,   'frequency':  0.04797,  'bturn_name':  'AZ',     'prev_name':  'new_prev_VIII',
                       'median_below_7A_ca1_ca4':  5.80,  'mean_below_7A_ca1_ca4':  5.84,  'median_any_ca1_ca4':  5.90,  'mean_any_ca1_ca4':  6.00,
                       'mode_omega2':  181.01,  'mode_phi2':  -74.11,   'mode_psi2':  -27.71,   'mode_omega3':  182.08,  'mode_phi3':  -140.41,
                       'mode_psi3':  75.09,   'mode_omega4':  190.54,  'mode_aa1':  'A',  'mode_aa2':  'G',  'mode_aa3':  'T',  'mode_aa4':  'P',
                       'mode_ss1':  'T',  'mode_ss2':  'T',  'mode_ss3':  'B',  'mode_ss4':  'T',  'mode_pdb_id':  '4RFU',  'mode_chain_id':  'A',
                       'mode_res1_id':  76,   'medoid_omega2':  181.51,  'medoid_phi2':  -83.81,   'medoid_psi2':  -17.76,   'medoid_omega3':  182.70,
                       'medoid_phi3':  -128.88,  'medoid_psi3':  61.72,   'medoid_omega4':  177.81,  'medoid_aa1':  'G',  'medoid_aa2':  'N',
                       'medoid_aa3':  'V',  'medoid_aa4':  'P',  'medoid_ss1':  'T',  'medoid_ss2':  'S',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '3LO8',  'medoid_chain_id':  'A',  'medoid_res1_id':  62}
    result['AB2']=    {'no_by_size':  7,   'cluster_size':  603,   'frequency':  0.04628,  'bturn_name':  'AB2',    'prev_name':  'VIII',
                       'median_below_7A_ca1_ca4':  6.42,  'mean_below_7A_ca1_ca4':  6.32,  'median_any_ca1_ca4':  7.95,  'mean_any_ca1_ca4':  7.72,
                       'mode_omega2':  175.48,  'mode_phi2':  -69.29,   'mode_psi2':  -30.16,   'mode_omega3':  169.41,  'mode_phi3':  -120.25,
                       'mode_psi3':  128.00,  'mode_omega4':  177.81,  'mode_aa1':  'G',  'mode_aa2':  'L',  'mode_aa3':  'I',  'mode_aa4':  'K',
                       'mode_ss1':  'T',  'mode_ss2':  'S',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '4YPO',  'mode_chain_id':  'A',
                       'mode_res1_id':  112,  'medoid_omega2':  179.88,  'medoid_phi2':  -76.37,   'medoid_psi2':  -33.97,   'medoid_omega3':  184.30,
                       'medoid_phi3':  -116.48,  'medoid_psi3':  120.10,  'medoid_omega4':  183.38,  'medoid_aa1':  'V',  'medoid_aa2':  'A',
                       'medoid_aa3':  'E',  'medoid_aa4':  'K',  'medoid_ss1':  'S',  'medoid_ss2':  'S',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '4UAS',  'medoid_chain_id':  'A',  'medoid_res1_id':  147}
    result['pD']=     {'no_by_size':  8,   'cluster_size':  402,   'frequency':  0.03085,  'bturn_name':  'pD',     'prev_name':  "II'",
                       'median_below_7A_ca1_ca4':  5.44,  'mean_below_7A_ca1_ca4':  5.60,  'median_any_ca1_ca4':  5.66,  'mean_any_ca1_ca4':  6.13,
                       'mode_omega2':  179.73,  'mode_phi2':  57.44,    'mode_psi2':  -130.18,  'mode_omega3':  181.81,  'mode_phi3':  -95.43,
                       'mode_psi3':  11.23,   'mode_omega4':  181.96,  'mode_aa1':  'K',  'mode_aa2':  'G',  'mode_aa3':  'S',  'mode_aa4':  'R',
                       'mode_ss1':  'T',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'C',  'mode_pdb_id':  '3BS2',  'mode_chain_id':  'A',
                       'mode_res1_id':  129,  'medoid_omega2':  177.42,  'medoid_phi2':  56.16,    'medoid_psi2':  -135.60,  'medoid_omega3':  183.17,
                       'medoid_phi3':  -90.92,   'medoid_psi3':  4.81,    'medoid_omega4':  180.04,  'medoid_aa1':  'K',  'medoid_aa2':  'G',
                       'medoid_aa3':  'Y',  'medoid_aa4':  'D',  'medoid_ss1':  'E',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '4G9S',  'medoid_chain_id':  'B',  'medoid_res1_id':  45}
    result['AG']=     {'no_by_size':  9,   'cluster_size':  166,   'frequency':  0.01274,  'bturn_name':  'AG',     'prev_name':  'new_prev_VIII',
                       'median_below_7A_ca1_ca4':  6.37,  'mean_below_7A_ca1_ca4':  6.19,  'median_any_ca1_ca4':  7.24,  'mean_any_ca1_ca4':  7.21,
                       'mode_omega2':  179.98,  'mode_phi2':  -66.42,   'mode_psi2':  -19.26,   'mode_omega3':  176.01,  'mode_phi3':  -82.48,
                       'mode_psi3':  63.16,   'mode_omega4':  183.48,  'mode_aa1':  'V',  'mode_aa2':  'N',  'mode_aa3':  'R',  'mode_aa4':  'A',
                       'mode_ss1':  'H',  'mode_ss2':  'T',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '5HB7',  'mode_chain_id':  'A',
                       'mode_res1_id':  199,  'medoid_omega2':  184.32,  'medoid_phi2':  -71.78,   'medoid_psi2':  -15.84,   'medoid_omega3':  187.19,
                       'medoid_phi3':  -87.70,   'medoid_psi3':  74.65,   'medoid_omega4':  181.63,  'medoid_aa1':  'N',  'medoid_aa2':  'A',
                       'medoid_aa3':  'N',  'medoid_aa4':  'A',  'medoid_ss1':  'T',  'medoid_ss2':  'T',  'medoid_ss3':  'C',  'medoid_ss4':  'H',
                       'medoid_pdb_id':  '4L8A',  'medoid_chain_id':  'A',  'medoid_res1_id':  119}
    result['BcisP']=  {'no_by_size':  10,  'cluster_size':  135,   'frequency':  0.01036,  'bturn_name':  'BcisP',  'prev_name':  'VIb',
                       'median_below_7A_ca1_ca4':  5.70,  'mean_below_7A_ca1_ca4':  5.70,  'median_any_ca1_ca4':  6.07,  'mean_any_ca1_ca4':  6.04,
                       'mode_omega2':  178.76,  'mode_phi2':  -137.53,  'mode_psi2':  119.38,   'mode_omega3':  359.11,  'mode_phi3':  -66.49,
                       'mode_psi3':  163.88,  'mode_omega4':  179.76,  'mode_aa1':  'P',  'mode_aa2':  'S',  'mode_aa3':  'P',  'mode_aa4':  'A',
                       'mode_ss1':  'C',  'mode_ss2':  'S',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '2CIW',  'mode_chain_id':  'A',
                       'mode_res1_id':  228,  'medoid_omega2':  176.27,  'medoid_phi2':  -127.05,  'medoid_psi2':  120.42,   'medoid_omega3':  358.00,
                       'medoid_phi3':  -65.63,   'medoid_psi3':  161.27,  'medoid_omega4':  182.14,  'medoid_aa1':  'N',  'medoid_aa2':  'N',
                       'medoid_aa3':  'P',  'medoid_aa4':  'K',  'medoid_ss1':  'G',  'medoid_ss2':  'S',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '2BMO',  'medoid_chain_id':  'B',  'medoid_res1_id':  110}
    result['dD']=     {'no_by_size':  11,  'cluster_size':  131,   'frequency':  0.01005,  'bturn_name':  'dD',     'prev_name':  'new',
                       'median_below_7A_ca1_ca4':  6.63,  'mean_below_7A_ca1_ca4':  6.48,  'median_any_ca1_ca4':  7.79,  'mean_any_ca1_ca4':  7.62,
                       'mode_omega2':  180.94,  'mode_phi2':  94.08,    'mode_psi2':  -1.19,    'mode_omega3':  186.08,  'mode_phi3':  -127.95,
                       'mode_psi3':  15.44,   'mode_omega4':  170.94,  'mode_aa1':  'D',  'mode_aa2':  'G',  'mode_aa3':  'Y',  'mode_aa4':  'H',
                       'mode_ss1':  'S',  'mode_ss2':  'S',  'mode_ss3':  'C',  'mode_ss4':  'S',  'mode_pdb_id':  '4EZI',  'mode_chain_id':  'A',
                       'mode_res1_id':  282,  'medoid_omega2':  182.57,  'medoid_phi2':  99.89,    'medoid_psi2':  -17.34,   'medoid_omega3':  185.02,
                       'medoid_phi3':  -113.78,  'medoid_psi3':  9.87,    'medoid_omega4':  180.65,  'medoid_aa1':  'P',  'medoid_aa2':  'G',
                       'medoid_aa3':  'Y',  'medoid_aa4':  'E',  'medoid_ss1':  'T',  'medoid_ss2':  'T',  'medoid_ss3':  'C',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '3JQ0',  'medoid_chain_id':  'A',  'medoid_res1_id':  521}
    result['PcisD']=  {'no_by_size':  12,  'cluster_size':  125,   'frequency':  0.00959,  'bturn_name':  'PcisD',  'prev_name':  'VIa1',
                       'median_below_7A_ca1_ca4':  5.74,  'mean_below_7A_ca1_ca4':  5.76,  'median_any_ca1_ca4':  5.76,  'mean_any_ca1_ca4':  5.80,
                       'mode_omega2':  175.04,  'mode_phi2':  -59.70,   'mode_psi2':  144.37,   'mode_omega3':  9.21,    'mode_phi3':  -92.77,
                       'mode_psi3':  8.32,    'mode_omega4':  175.56,  'mode_aa1':  'D',  'mode_aa2':  'S',  'mode_aa3':  'P',  'mode_aa4':  'L',
                       'mode_ss1':  'C',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'C',  'mode_pdb_id':  '5AOZ',  'mode_chain_id':  'A',
                       'mode_res1_id':  443,  'medoid_omega2':  175.04,  'medoid_phi2':  -59.70,   'medoid_psi2':  144.37,   'medoid_omega3':  9.21,
                       'medoid_phi3':  -92.77,   'medoid_psi3':  8.32,    'medoid_omega4':  175.56,  'medoid_aa1':  'D',  'medoid_aa2':  'S',
                       'medoid_aa3':  'P',  'medoid_aa4':  'L',  'medoid_ss1':  'C',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '5AOZ',  'medoid_chain_id':  'A',  'medoid_res1_id':  443}
    result['dN']=     {'no_by_size':  13,  'cluster_size':  102,   'frequency':  0.00783,  'bturn_name':  'dN',     'prev_name':  'new',
                       'median_below_7A_ca1_ca4':  6.42,  'mean_below_7A_ca1_ca4':  6.37,  'median_any_ca1_ca4':  7.83,  'mean_any_ca1_ca4':  7.60,
                       'mode_omega2':  183.90,  'mode_phi2':  69.32,    'mode_psi2':  8.93,     'mode_omega3':  177.80,  'mode_phi3':  -131.73,
                       'mode_psi3':  -63.28,  'mode_omega4':  185.11,  'mode_aa1':  'I',  'mode_aa2':  'E',  'mode_aa3':  'H',  'mode_aa4':  'G',
                       'mode_ss1':  'T',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'E',  'mode_pdb_id':  '5AGD',  'mode_chain_id':  'A',
                       'mode_res1_id':  232,  'medoid_omega2':  178.53,  'medoid_phi2':  76.06,    'medoid_psi2':  -3.15,    'medoid_omega3':  179.31,
                       'medoid_phi3':  -122.82,  'medoid_psi3':  -50.46,  'medoid_omega4':  179.84,  'medoid_aa1':  'P',  'medoid_aa2':  'G',
                       'medoid_aa3':  'V',  'medoid_aa4':  'T',  'medoid_ss1':  'B',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'T',
                       'medoid_pdb_id':  '4E9X',  'medoid_chain_id':  'A',  'medoid_res1_id':  1076}
    result['Dd']=     {'no_by_size':  14,  'cluster_size':  87,    'frequency':  0.00668,  'bturn_name':  'Dd',     'prev_name':  'new',
                       'median_below_7A_ca1_ca4':  6.81,  'mean_below_7A_ca1_ca4':  6.70,  'median_any_ca1_ca4':  7.60,  'mean_any_ca1_ca4':  7.65,
                       'mode_omega2':  182.75,  'mode_phi2':  -115.16,  'mode_psi2':  15.69,    'mode_omega3':  186.11,  'mode_phi3':  100.65,
                       'mode_psi3':  -12.30,  'mode_omega4':  176.58,  'mode_aa1':  'L',  'mode_aa2':  'D',  'mode_aa3':  'G',  'mode_aa4':  'S',
                       'mode_ss1':  'S',  'mode_ss2':  'S',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '2RFR',  'mode_chain_id':  'A',
                       'mode_res1_id':  141,  'medoid_omega2':  177.54,  'medoid_phi2':  -99.08,   'medoid_psi2':  19.56,    'medoid_omega3':  173.80,
                       'medoid_phi3':  108.72,   'medoid_psi3':  -14.53,  'medoid_omega4':  182.16,  'medoid_aa1':  'D',  'medoid_aa2':  'E',
                       'medoid_aa3':  'G',  'medoid_aa4':  'G',  'medoid_ss1':  'H',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'G',
                       'medoid_pdb_id':  '4HZ8',  'medoid_chain_id':  'A',  'medoid_res1_id':  167}
    result['PcisP']=  {'no_by_size':  15,  'cluster_size':  59,    'frequency':  0.00453,  'bturn_name':  'PcisP',  'prev_name':  'new_prev_VIb',
                       'median_below_7A_ca1_ca4':  6.33,  'mean_below_7A_ca1_ca4':  6.24,  'median_any_ca1_ca4':  6.74,  'mean_any_ca1_ca4':  6.79,
                       'mode_omega2':  175.26,  'mode_phi2':  -66.01,   'mode_psi2':  147.89,   'mode_omega3':  0.25,    'mode_phi3':  -75.72,
                       'mode_psi3':  142.11,  'mode_omega4':  176.61,  'mode_aa1':  'D',  'mode_aa2':  'A',  'mode_aa3':  'P',  'mode_aa4':  'Y',
                       'mode_ss1':  'T',  'mode_ss2':  'C',  'mode_ss3':  'S',  'mode_ss4':  'E',  'mode_pdb_id':  '5DZE',  'mode_chain_id':  'A',
                       'mode_res1_id':  189,  'medoid_omega2':  183.55,  'medoid_phi2':  -78.86,   'medoid_psi2':  145.45,   'medoid_omega3':  353.82,
                       'medoid_phi3':  -78.06,   'medoid_psi3':  140.34,  'medoid_omega4':  178.16,  'medoid_aa1':  'M',  'medoid_aa2':  'K',
                       'medoid_aa3':  'P',  'medoid_aa4':  'Q',  'medoid_ss1':  'H',  'medoid_ss2':  'C',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '4UDX',  'medoid_chain_id':  'X',  'medoid_res1_id':  35}
    result['cisDA']=  {'no_by_size':  16,  'cluster_size':  52,    'frequency':  0.00399,  'bturn_name':  'cisDA',  'prev_name':  'new',
                       'median_below_7A_ca1_ca4':  5.36,  'mean_below_7A_ca1_ca4':  5.42,  'median_any_ca1_ca4':  5.36,  'mean_any_ca1_ca4':  5.42,
                       'mode_omega2':  3.33,    'mode_phi2':  -94.18,   'mode_psi2':  7.96,     'mode_omega3':  187.90,  'mode_phi3':  -61.40,
                       'mode_psi3':  -38.41,  'mode_omega4':  184.42,  'mode_aa1':  'Y',  'mode_aa2':  'P',  'mode_aa3':  'D',  'mode_aa4':  'D',
                       'mode_ss1':  'E',  'mode_ss2':  'T',  'mode_ss3':  'T',  'mode_ss4':  'T',  'mode_pdb_id':  '3VGI',  'mode_chain_id':  'A',
                       'mode_res1_id':  56,   'medoid_omega2':  0.14,    'medoid_phi2':  -97.53,   'medoid_psi2':  3.93,     'medoid_omega3':  184.13,
                       'medoid_phi3':  -66.86,   'medoid_psi3':  -35.89,  'medoid_omega4':  182.02,  'medoid_aa1':  'A',  'medoid_aa2':  'P',
                       'medoid_aa3':  'W',  'medoid_aa4':  'F',  'medoid_ss1':  'T',  'medoid_ss2':  'T',  'medoid_ss3':  'T',  'medoid_ss4':  'E',
                       'medoid_pdb_id':  '1UAI',  'medoid_chain_id':  'A',  'medoid_res1_id':  42}
    result['pG']=     {'no_by_size':  17,  'cluster_size':  32,    'frequency':  0.00246,  'bturn_name':  'pG',     'prev_name':  'new',
                       'median_below_7A_ca1_ca4':  6.58,  'mean_below_7A_ca1_ca4':  6.43,  'median_any_ca1_ca4':  8.16,  'mean_any_ca1_ca4':  8.06,
                       'mode_omega2':  180.06,  'mode_phi2':  73.58,    'mode_psi2':  -162.38,  'mode_omega3':  180.27,  'mode_phi3':  -78.96,
                       'mode_psi3':  77.17,   'mode_omega4':  188.34,  'mode_aa1':  'L',  'mode_aa2':  'G',  'mode_aa3':  'R',  'mode_aa4':  'Y',
                       'mode_ss1':  'T',  'mode_ss2':  'C',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '4V28',  'mode_chain_id':  'A',
                       'mode_res1_id':  95,   'medoid_omega2':  188.70,  'medoid_phi2':  68.28,    'medoid_psi2':  -139.60,  'medoid_omega3':  178.07,
                       'medoid_phi3':  -79.92,   'medoid_psi3':  121.91,  'medoid_omega4':  175.54,  'medoid_aa1':  'H',  'medoid_aa2':  'G',
                       'medoid_aa3':  'T',  'medoid_aa4':  'Q',  'medoid_ss1':  'S',  'medoid_ss2':  'S',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '1KQP',  'medoid_chain_id':  'A',  'medoid_res1_id':  81}
    result['cisDP']=  {'no_by_size':  18,  'cluster_size':  25,    'frequency':  0.00192,  'bturn_name':  'cisDP',  'prev_name':  'new',
                       'median_below_7A_ca1_ca4':  6.36,  'mean_below_7A_ca1_ca4':  6.29,  'median_any_ca1_ca4':  6.52,  'mean_any_ca1_ca4':  6.60,
                       'mode_omega2':  11.10,   'mode_phi2':  -86.47,   'mode_psi2':  4.44,     'mode_omega3':  169.99,  'mode_phi3':  -71.39,
                       'mode_psi3':  157.80,  'mode_omega4':  175.12,  'mode_aa1':  'A',  'mode_aa2':  'P',  'mode_aa3':  'L',  'mode_aa4':  'V',
                       'mode_ss1':  'T',  'mode_ss2':  'T',  'mode_ss3':  'S',  'mode_ss4':  'C',  'mode_pdb_id':  '3B4U',  'mode_chain_id':  'A',
                       'mode_res1_id':  268,  'medoid_omega2':  11.10,   'medoid_phi2':  -86.47,   'medoid_psi2':  4.44,     'medoid_omega3':  169.99,
                       'medoid_phi3':  -71.39,   'medoid_psi3':  157.80,  'medoid_omega4':  175.12,  'medoid_aa1':  'A',  'medoid_aa2':  'P',
                       'medoid_aa3':  'L',  'medoid_aa4':  'V',  'medoid_ss1':  'T',  'medoid_ss2':  'T',  'medoid_ss3':  'S',  'medoid_ss4':  'C',
                       'medoid_pdb_id':  '3B4U',  'medoid_chain_id':  'A',  'medoid_res1_id':  268}
    return result


def parse_dssp_struct_summary(mmCIF_path):
    """
    Parse the _dssp_struct_summary loop from an mmCIF file produced by mkdssp:

    mkdssp input.cif input_dssp.cif

    Returns
    -------
    dict
        {(label_asym_id, label_seq_id): {col_name: value, ...}, ...}
        Column names are taken from the item names after the dot,
        e.g., '_dssp_struct_summary.secondary_structure' -> 'secondary_structure'.
        The key fields 'label_asym_id' and 'label_seq_id' are excluded from the value dict.

    example dssp mmCIF file:
    
    loop_
    _dssp_struct_summary.entry_id
    _dssp_struct_summary.label_asym_id
    _dssp_struct_summary.label_seq_id
    _dssp_struct_summary.label_comp_id
    _dssp_struct_summary.secondary_structure
    _dssp_struct_summary.ss_bridge
    _dssp_struct_summary.helix_3_10
    _dssp_struct_summary.helix_alpha
    _dssp_struct_summary.helix_pi
    _dssp_struct_summary.helix_pp
    _dssp_struct_summary.bend
    _dssp_struct_summary.chirality
    _dssp_struct_summary.sheet
    _dssp_struct_summary.strand
    _dssp_struct_summary.ladder_1
    _dssp_struct_summary.ladder_2
    _dssp_struct_summary.accessibility
    _dssp_struct_summary.TCO
    _dssp_struct_summary.kappa
    _dssp_struct_summary.alpha
    _dssp_struct_summary.phi
    _dssp_struct_summary.psi
    _dssp_struct_summary.x_ca
    _dssp_struct_summary.y_ca
    _dssp_struct_summary.z_ca
    1B9O A 1   LYS . . . . . . . . . . . . ?      .     .      .      .  141.4   7.2 28.4 41.0
    1B9O A 2   GLN B . . . . . . - A A A . ? -0.965     . -156.0 -111.0  110.9   7.1 24.6 40.2
    1B9O A 3   PHE . . . . . . . - . . . . ? -0.478  10.2 -126.4  -82.6  158.5   3.5 23.4 40.3
    1B9O A 4   THR . . . > . . . - . . . . ? -0.588  30.6 -106.4  -90.9  162.3   2.3 19.8 41.0
    1B9O A 5   LYS H . . > . . S + . . . . ?  0.910 120.9   48.3  -54.6  -45.2   0.0 18.1 38.6
    1B9O A 6   CYS H 1 . > . . S + . . . . ?  0.859 108.4   54.5  -70.7  -33.8  -2.9 18.4 41.0
    1B9O A 7   GLU H . . > . . S + . . . . ?  0.920 111.2   44.0  -60.7  -45.4  -2.3 22.1 41.6
    1B9O A 8   LEU H . . X . . S + . . . . ?  0.871 107.4   60.4  -68.8  -36.5  -2.3 22.9 37.9
    1B9O A 9   SER H . . < . . S + . . . . ?  0.928 107.8   45.3  -55.8  -40.9  -5.4 20.8 37.4
    1B9O A 10  GLN H . > < . . S + . . . . ?  0.883 115.6   46.3  -67.3  -40.6  -7.2 23.0 39.8
    1B9O A 11  LEU H . 3 < . . S + . . . . ?  0.855 119.6   38.9  -70.7  -34.9  -6.0 26.2 38.2
    1B9O A 12  LEU T . > X . . S + . . . . ?  0.290  80.4  110.7  -96.3    4.2  -6.7 25.0 34.7
    1B9O A 13  LYS G . X 4 . . S + . . . . ?  0.856  77.1   53.0  -50.1  -39.3 -10.0 23.3 35.4
    1B9O A 14  ASP G . 3 4 . . S + . . . . ?  0.668  99.4   61.6  -78.5  -12.1 -11.9 25.9 33.5
    1B9O A 15  ILE G . X 4 . . S + . . . . ?  0.572  71.8  126.1  -87.8   -2.8  -9.9 25.7 30.4
    1B9O A 16  ASP T . < < . . S + . . . . ? -0.392  80.5   14.4  -56.8  125.4 -10.9 22.1 29.9
    1B9O A 17  GLY T . > . . . S + . . . . ?  0.172  85.8  146.1   91.5   -7.5 -12.3 21.8 26.3
    1B9O A 18  TYR G . X . . . S + . . . . ? -0.376  83.0    0.1  -59.8  130.6 -11.0 25.2 25.2
    1B9O A 19  GLY G . 3 . . . S - . . . . ?  0.722 130.9  -68.6   62.3   15.2 -10.1 24.8 21.5
    1B9O A 20  GLY G . < . . . S + . . . . ?  0.687 100.4  135.0   77.2   18.2 -11.2 21.2 21.8

    """

    # reading dssp structure summary table only
    category_prefix = "_dssp_struct_summary."

    def norm(val):  # remove enclosing single or double quotes from values
        if len(val) >= 2 and ((val[0] == val[-1] == '"') or (val[0] == val[-1] == "'")):
            return val[1:-1]
        return val

    result = {}

    with open(mmCIF_path, "r", encoding="utf-8") as fh:
        lines = iter(fh)
        for line in lines:
            if line.strip() != "loop_":
                continue

            # Gather item names for this loop
            item_names = []
            for line in lines:
                s = line.strip()
                if s.startswith("_"):
                    item_names.append(s)
                else:
                    first_data_line = line
                    break  # first non-item line enters data section

            if not item_names:
                # Malformed loop; skip to end marker
                for line in lines:
                    if line.strip() == "#":
                        break
                continue

            # Process only loops that include our category
            if not any(n.startswith(category_prefix) for n in item_names):
                # Consume data lines for unrelated loop and move on
                for line in lines:
                    if line.strip() == "#":
                        break
                continue

            ncols = len(item_names)

            # Map column indices to short column names for our category
            idx_to_short = {}
            for i, name in enumerate(item_names):
                if name.startswith(category_prefix):
                    idx_to_short[i] = name.split(".", 1)[1]

            # Accumulate all data tokens in this loop until '#'
            tokens = []
            def feed(s):
                # mmCIF values are whitespace-separated; shlex handles quotes
                return shlex.split(s, posix=True)

            pending_lines = [first_data_line]
            for line in lines:
                pending_lines.append(line)
                if line.strip() == "#":
                    break

            for s in pending_lines:
                if s.strip() in {"", "#"}:
                    continue
                tokens.extend(feed(s))

            if ncols == 0:
                continue

            # Ensure token count is a multiple of ncols (be forgiving if not)
            if len(tokens) % ncols != 0:
                tokens = tokens[: (len(tokens) // ncols) * ncols]

            # Walk rows
            for rstart in range(0, len(tokens), ncols):
                row = tokens[rstart : rstart + ncols]
                if len(row) < ncols:
                    continue

                # Extract only our category's fields
                row_vals = {}
                for idx, short in idx_to_short.items():
                    row_vals[short] = norm(row[idx])

                if not row_vals:
                    continue  # Mixed-category row without our fields

                label_asym_id = row_vals.get("label_asym_id")
                label_seq_id  = row_vals.get("label_seq_id")
                if label_asym_id is None or label_seq_id is None:
                    continue  # Cannot key this row

                # mkdssp does NOT provide auth_asym_id, pdb_strand_id (same as auth_asym_id), or auth_seq_id, so we use label_asym_id and label_seq_id as keys
                # label_seq_id is 1_to_N numbering for each chain; label_asym_id is PDB's chain_id which may differ from author
                key = (label_asym_id, label_seq_id)
                value = {k: v for k, v in row_vals.items()
                         if k not in ("label_asym_id", "label_seq_id")}
                result[key] = value

    return result


def parse_pdbx_poly_seq_scheme(mmCIF_path):
    """
    Parse the _pdbx_poly_seq_scheme loop from an mmCIF file.
    mkdssp mmCIF file may be missing some fields.
    Key should be (pdb_strand_id, pdb_seq_num) but for mkdssp files missing these (usually when converted from PDB-format files,
    asym_id or auth_seq_num may be used instead.

    Returns
    -------
    dict
        {(pdb_strand_id, pdb_seq_num): {col_name: value, ...}, ...}
        All column names are taken from the item names after the dot,
        e.g., '_pdbx_poly_seq_scheme.seq_id' -> 'seq_id'.
        The key fields 'pdb_strand_id' and 'pdb_seq_num' are excluded from the value dict.
    """
    result = {}

    with open(mmCIF_path, "r", encoding="utf-8") as fh:
        lines = iter(fh)
        for line in lines:
            if line.strip() != "loop_":
                continue

            # Collect item names for this loop
            item_names = []
            for line in lines:
                s = line.strip()
                if s.startswith("_"):
                    item_names.append(s)
                else:
                    first_data_line = line
                    break  # first non-item line starts the data section

            if not item_names:
                # malformed loop; skip to the next '#'
                for line in lines:
                    if line.strip() == "#":
                        break
                continue

            # Only process loops that include _pdbx_poly_seq_scheme items
            cat_prefix = "_pdbx_poly_seq_scheme."
            if not any(n.startswith(cat_prefix) for n in item_names):
                # consume data rows of this unrelated loop, then move on
                # Data ends at a line with just '#'
                tokens_needed = len(item_names)
                tokens_accum = []
                # account for the line we already read (first_data_line)
                to_process = [first_data_line]
                for line in lines:
                    to_process.append(line)
                    if line.strip() == "#":
                        break
                # No further action; skip unrelated loop entirely
                continue

            # We will parse all rows in this loop; some loops may (rarely) mix categories,
            # so keep the full item_names but only extract columns with our prefix.
            full_names = item_names
            ncols = len(full_names)

            # Map indices -> short names for our category
            idx_to_short = {}
            for i, name in enumerate(full_names):
                if name.startswith(cat_prefix):
                    idx_to_short[i] = name.split(".", 1)[1]

            # Helper to normalize tokens
            def norm(val):
                v = val
                # strip surrounding quotes if any (shlex already handles most)
                if (len(v) >= 2) and ((v[0] == v[-1] == '"') or (v[0] == v[-1] == "'")):
                    v = v[1:-1]
                return v

            # Accumulate tokens across lines until we hit '#'
            tokens = []
            def feed_line(s):
                # shlex handles quotes and escapes; mmCIF values are whitespace-separated
                return shlex.split(s, posix=True)

            # Start with the line we already read
            pending_lines = [first_data_line]
            for line in lines:
                pending_lines.append(line)
                if line.strip() == "#":
                    break

            for s in pending_lines:
                if s.strip() in {"", "#"}:
                    continue
                tokens.extend(feed_line(s))

            # Split flat token stream into rows of ncols
            if ncols == 0:
                return result  # nothing to do

            if len(tokens) % ncols != 0:
                # Try to be forgiving: truncate extra tokens if any (rare/malformed files)
                tokens = tokens[: (len(tokens) // ncols) * ncols]

            for rstart in range(0, len(tokens), ncols):
                row = tokens[rstart : rstart + ncols]
                if len(row) < ncols:
                    continue

                # Extract our category fields from the row
                row_vals = {}
                for idx, short in idx_to_short.items():
                    row_vals[short] = norm(row[idx])

                if not row_vals:
                    continue  # row didn't contain our category (mixed loop row)

                # Form the key

                pdb_strand_id = row_vals.get("pdb_strand_id") or row_vals.get("asym_id") or " "
                pdb_seq_num   = row_vals.get("pdb_seq_num") or row_vals.get("auth_seq_num")
                if pdb_seq_num is None:
                    continue  # Can't key this row
                
                key = (pdb_strand_id, pdb_seq_num)

                # Build value dict 
                value = {k: v for k, v in row_vals.items()}
                #                         if k not in ("pdb_strand_id", "pdb_seq_num")}

                result[key] = value

    return result

def run_dssp(input_pdb_or_cif, dssp_executable="/usr/local/bin/mkdssp"):
    """
    Run mkdssp on the input PDB/mmCIF file after removing problematic _audit_conform lines
    (which cause mkdssp to fail on AlphaFold models). Returns the path to the DSSP output file.
    """
    if not os.path.isfile(input_pdb_or_cif):
        raise FileNotFoundError(f"Structure file not found: {input_pdb_or_cif}")
    if not os.path.isfile(dssp_executable):
        raise FileNotFoundError(f"DSSP executable not found: {dssp_executable}")

    # Remove _audit_conform lines and write to a temporary cleaned file
    with open(input_pdb_or_cif, 'r') as infile, tempfile.NamedTemporaryFile('w', delete=False, suffix=".cif") as cleaned:
        for line in infile:
            if not line.startswith("_audit_conform."):
                cleaned.write(line)
        cleaned_input_path = cleaned.name

    root = os.path.splitext(os.path.basename(input_pdb_or_cif))[0]
    dssp_output_path = f"{root}_dssp.cif"

    cmd = [dssp_executable, cleaned_input_path, dssp_output_path]

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        os.remove(cleaned_input_path)
        raise RuntimeError(f"Error running DSSP: {e.stderr.decode()}")
    
    # Optionally remove the temporary cleaned input file
    os.remove(cleaned_input_path)

    return dssp_output_path



def compute_phi(structure, model, chain, prev_residue, curr_residue):          
    phi=999.00
    if (prev_residue.has_id('C') and
        curr_residue.has_id('N') and
        curr_residue.has_id('CA') and
        curr_residue.has_id('C') and
        curr_residue.id[1] == prev_residue.id[1]+1):
        prev_c=structure[model.id][chain.id][prev_residue.id]['C'].get_vector()
        curr_n=structure[model.id][chain.id][curr_residue.id]['N'].get_vector() 
        curr_a=structure[model.id][chain.id][curr_residue.id]['CA'].get_vector() 
        curr_c=structure[model.id][chain.id][curr_residue.id]['C'].get_vector()
        phi=round(math.degrees(calc_dihedral(prev_c,curr_n,curr_a,curr_c)),6)
    return phi

def compute_psi(structure, model,chain, curr_residue, next_residue):
    psi=999.00
    if (curr_residue.has_id('N') and
        curr_residue.has_id('CA') and
        curr_residue.has_id('C') and
        next_residue.has_id('N') and
        next_residue.id[1] == curr_residue.id[1]+1):
        curr_n=structure[model.id][chain.id][curr_residue.id]['N'].get_vector() 
        curr_a=structure[model.id][chain.id][curr_residue.id]['CA'].get_vector() 
        curr_c=structure[model.id][chain.id][curr_residue.id]['C'].get_vector()
        next_n=structure[model.id][chain.id][next_residue.id]['N'].get_vector()
        psi=round(math.degrees(calc_dihedral(curr_n,curr_a,curr_c,next_n)),6)
    return psi

def compute_omega(structure,model,chain,prev_residue,curr_residue):
    omega=999.00
    if (prev_residue.has_id('CA') and
        prev_residue.has_id('C') and
        curr_residue.has_id('N') and
        curr_residue.has_id('CA') and
        curr_residue.id[1]==prev_residue.id[1]+1):
        prev_ca=structure[model.id][chain.id][prev_residue.id]['CA'].get_vector() 
        prev_c=structure[model.id][chain.id][prev_residue.id]['C'].get_vector() 
        curr_n=structure[model.id][chain.id][curr_residue.id]['N'].get_vector()
        curr_a=structure[model.id][chain.id][curr_residue.id]['CA'].get_vector()
        omega=round(math.degrees(calc_dihedral(prev_ca,prev_c,curr_n,curr_a)),6)
    return omega

# chi angles are not used but are legacy from another script; might be useful
def compute_chi1(structure,model,chain,curr_residue):
    chi1=999.00
    curr_n=0
    curr_a=0
    curr_b=0
    curr_g=0

    atomname='N'
    if curr_residue.has_id(atomname): curr_n= structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    atomname='CA'
    if curr_residue.has_id(atomname): curr_a=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    atomname='CB'
    if curr_residue.has_id(atomname): curr_b=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    for atomname in ('CG', 'CG1', 'OG', 'OG1', 'SG', 'NG', 'PG'):
        if curr_residue.has_id(atomname): curr_g=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    if curr_n and curr_a and curr_b and curr_g:
        chi1=round(math.degrees(calc_dihedral(curr_n, curr_a, curr_b, curr_g)),6)
    return chi1

def compute_chi2(structure,model,chain,curr_residue):
    chi2=999.00
    curr_a=0
    curr_b=0
    curr_g=0
    curr_d=0

    atomname='CA'
    if curr_residue.has_id(atomname): curr_a=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    atomname='CB'
    if curr_residue.has_id(atomname): curr_b=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    for atomname in ('CG', 'CG1', 'OG', 'OG1', 'SG', 'NG', 'PG'):
        if curr_residue.has_id(atomname): curr_g=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CD', 'CD1', 'OD', 'OD1', 'ND', 'ND1', 'SD', 'PD'):   
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    
    if curr_residue.resname=='MSE':
        atomname='SE'
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    if curr_residue.resname=='SEP' or curr_residue.resname=="TPO":
        atomname='P'
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    if curr_a and curr_b and curr_g and curr_d:
        chi2=round(math.degrees(calc_dihedral(curr_a,curr_b,curr_g,curr_d)),6)
        
    return chi2
    
def compute_chi3(structure,model,chain,curr_residue):
    chi3=999.00
    if curr_residue.resname in ('PHE','TYR','TRP','HIS','PTR'): return chi3
    curr_b=0
    curr_g=0
    curr_d=0
    curr_e=0

    atomname='CB'
    if curr_residue.has_id(atomname): curr_b=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    for atomname in ('CG', 'CG1', 'OG', 'OG1', 'SG', 'NG'):
        if curr_residue.has_id(atomname): curr_g=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CD', 'CD1', 'OD', 'OD1', 'ND', 'ND1', 'SD'):   # SD and SE are for MET and MSE in D position
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CE', 'CE1', 'OE', 'OE1', 'NE', 'NE1'): 
        if curr_residue.has_id(atomname): curr_e=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
           
    if curr_residue.resname=='MSE':
        atomname='SE'
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    if curr_residue.resname=='SEP' or curr_residue.resname=="TPO":
        atomname='P'
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
        atomname='OP1'
        if curr_residue.has_id(atomname): curr_e=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()

    if curr_b and curr_g and curr_d and curr_e:
        chi3=round(math.degrees(calc_dihedral(curr_b,curr_g,curr_d,curr_e)),6)
        
    return chi3
    
def compute_chi4(structure,model,chain,curr_residue):
    chi4=999.00
    if curr_residue.resname in ('PHE','TYR','TRP','PTR','MET','MSE'): return chi4
    curr_g=0
    curr_d=0
    curr_e=0
    curr_z=0

    for atomname in ('CG', 'CG1', 'OG', 'OG1', 'SG', 'NG', 'PG'):
        if curr_residue.has_id(atomname): curr_g=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CD', 'CD1', 'OD', 'OD1', 'ND', 'ND1', 'SD', 'PD'):   
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CE', 'CE1', 'OE', 'OE1', 'NE', 'NE1', 'SE', 'PE', 'SE'):   
        if curr_residue.has_id(atomname): curr_e=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CZ', 'CZ1', 'OZ', 'OZ1', 'NZ', 'NZ1'):   
        if curr_residue.has_id(atomname): curr_z=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
           
    if curr_g and curr_d and curr_e and curr_z:
        chi4=round(math.degrees(calc_dihedral(curr_g,curr_d,curr_e,curr_z)),6)
        
    return chi4

def compute_chi5(structure,model,chain,curr_residue):
    chi5=999.00
    if curr_residue.resname in ('TYR','TRP','PTR','MSE','MET'): return chi5
    curr_d=0
    curr_e=0
    curr_z=0
    curr_h=0

    for atomname in ('CD', 'CD1', 'OD', 'OD1', 'ND', 'ND1', 'SD', 'PD'):   
        if curr_residue.has_id(atomname): curr_d=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CE', 'CE1', 'OE', 'OE1', 'NE', 'NE1', 'SE', 'PE', 'SE'):   
        if curr_residue.has_id(atomname): curr_e=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CZ', 'CZ1', 'OZ', 'OZ1', 'NZ', 'NZ1'):   
        if curr_residue.has_id(atomname): curr_z=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
    for atomname in ('CH', 'CH1', 'OH', 'OH1', 'NH', 'NH1'):   
        if curr_residue.has_id(atomname): curr_h=structure[model.id][chain.id][curr_residue.id][atomname].get_vector()
           
    if curr_d and curr_e and curr_z and curr_h:
        chi5=round(math.degrees(calc_dihedral(curr_d,curr_e,curr_z,curr_h)),6)
        
    return chi5

# Main

def main():
    pdbfilename=sys.argv[1]
    base = os.path.basename(pdbfilename).replace(".gz", "")
    shortfilename = os.path.splitext(base)[0]
    
    
    # run analysis
    dssp_file_path=run_dssp(pdbfilename)  # run mkdssp from DSSP4.5 which outputs an mmCIF file with DSSP data; easier to pass than standard dssp output
    seqscheme_dict = parse_pdbx_poly_seq_scheme(dssp_file_path) # read mmcif poly_seq_scheme to get correspondence of sequence (1 to N) to auth_seq_nums
    dssp_dict = parse_dssp_struct_summary(dssp_file_path)  # 
    turnlibrary=define_turn_library()  # read turn library data of 18 turn types from our 2019 paper
    
    if '.gz' in pdbfilename.lower():      
        handle=gzip.open(pdbfilename, 'rt')
    else:
        handle=open(pdbfilename, 'r')
    
    if '.cif' in pdbfilename.lower():
        parser=MMCIFParser(QUIET=True)
    if '.pdb' in pdbfilename.lower():
        parser=PDBParser(QUIET=True)
    
    structure1=parser.get_structure(pdbfilename, handle)
        
    datadict={}
    
    for model1 in structure1:
        for chain1 in model1:
            aa_residues = [r for r in chain1 if is_aa(r)]
            if len(aa_residues) < 4:
                continue  # Skip short chains — can't form a 4-residue turn
    
            datadict[chain1] = []
    
            for i, res in enumerate(aa_residues):
                phi = psi = omega = chi1 = chi2 = chi3 = chi4 = chi5 = 999.0
    
                if i > 0:
                    phi = compute_phi(structure1, model1, chain1, aa_residues[i - 1], res)
                    omega = compute_omega(structure1, model1, chain1, aa_residues[i - 1], res)
    
                if i + 1 < len(aa_residues):
                    psi = compute_psi(structure1, model1, chain1, res, aa_residues[i + 1])
    
                chi1 = compute_chi1(structure1, model1, chain1, res)
                chi2 = compute_chi2(structure1, model1, chain1, res)
                chi3 = compute_chi3(structure1, model1, chain1, res)
                chi4 = compute_chi4(structure1, model1, chain1, res)
                chi5 = compute_chi5(structure1, model1, chain1, res)
    
                resnum = res.id[1]
                res3letter = res.resname
                res1letter = getres1(res3letter)
    
                # Safe lookup of DSSP label
                auth_key = (chain1.id, str(resnum))
                if auth_key in seqscheme_dict:
                    label_asym_id = seqscheme_dict[auth_key]['asym_id']
                    label_seq_id = str(seqscheme_dict[auth_key]['seq_id'])
                    label_key = (label_asym_id, label_seq_id)
                else:
                    label_key = None
                    # Optional: warn about missing mapping
                    print(f"Warning: Missing seqscheme_dict entry for {auth_key}", file=sys.stderr)
    
                dssp_value = "."
                three10 = "."
    
                if label_key and label_key in dssp_dict:
                    dssp_value = dssp_dict[label_key].get('secondary_structure', '.')
                    three10 = dssp_dict[label_key].get('helix_3_10', '.')
    
                if dssp_value == ".":
                    dssp_value = "C"  # for "coil"
    
                CAcoor = res['CA'].get_vector() if res.has_id('CA') else None
    
                datadict[chain1].append({
                    'model1': model1,
                    'chain_id': chain1.id,
                    'resnum': resnum,
                    'resname': res3letter,
                    'res1': res1letter,
                    'omega': omega,
                    'phi': phi,
                    'psi': psi,
                    'chi1': chi1,
                    'chi2': chi2,
                    'chi3': chi3,
                    'chi4': chi4,
                    'chi5': chi5,
                    'dssp_ss': dssp_value,
                    'three10': three10,
                    'auth_key': auth_key,
                    'label_key': label_key,
                    'CAcoor': CAcoor
                })
    
    
    # Find beta turns
    nturns=0

    # print header
    print(f'{"turn":<4} {"num":>4} {"chn":<4} {"res1":<4} {"res4":<4}    {"seq":<4} {"dssp":<4}    {"type":<5}  {"prev_name":<13}    {"Dist":>6} {"DistAng":>7} {"CA1-CA4":>7}    {"omega2":>7} {"phi2":>7} {"psi2":>7}  {"omega3":>7} {"phi3":>7}  {"psi3":>7} {"omega4":>7}   {"filename"}')
    
    
    for chain in datadict:
        for i in range(len(datadict[chain]) - 3):
    
            vector1=datadict[chain][i]['CAcoor']   
            if vector1 is None: continue  # CA atom is missing
            vector4=datadict[chain][i+3]['CAcoor']  
            if vector4 is None: continue  # CA atom is missing
            CA1_CA4_distance=(vector1 - vector4).norm()
            dssp1=datadict[chain][i]['dssp_ss']
            dssp2=datadict[chain][i+1]['dssp_ss']
            dssp3=datadict[chain][i+2]['dssp_ss']
            dssp4=datadict[chain][i+3]['dssp_ss']
            dssp_string = dssp1 + dssp2 + dssp3 + dssp4

            # get helix_3_10 field from dssp. Three-residue 3-10 helices have "3" in this field for middle residue
            threeten1=datadict[chain][i]['three10']
            threeten2=datadict[chain][i+1]['three10']
            threeten3=datadict[chain][i+2]['three10']
            threeten4=datadict[chain][i+3]['three10']
    
            dssp_string = dssp1 + dssp2 + dssp3 + dssp4
            three10string = threeten1 + threeten2 + threeten3 + threeten4
            if dssp_string == "HGGG": continue  # skip G helices abutting alpha helices
            if dssp_string == "GGGH": continue  # skip G helices abutting alpha helices
            if dssp_string == "GGGG": continue  # skip long G helices
            if dssp2 in ('H','E'): continue  # skip regular secondary structure for residue 2
            if dssp3 in ('H','E'): continue  # skip regular secondary structure for residue 3
            
            # include GGG if they are three residues but not longer; DSSP always has "3" in middle residue of helix_3_10 value
            if dssp_string.startswith("GGG") or dssp_string.endswith("GGG"):
                if "3" not in three10string: continue
    
            omega2=datadict[chain][i+1]['omega']
            if omega2 == 999.0: continue
            phi2=datadict[chain][i+1]['phi']
            if phi2 == 999.0: continue
            psi2=datadict[chain][i+1]['psi']
            if psi2 == 999.0: continue
            omega3=datadict[chain][i+2]['omega']
            if omega3 == 999.0: continue
            phi3=datadict[chain][i+2]['phi']
            if phi3 == 999.0: continue
            psi3=datadict[chain][i+2]['psi']
            if psi3 == 999.0: continue
            omega4=datadict[chain][i+3]['omega']
            if omega4 == 999.0: continue
            if CA1_CA4_distance>7.0: continue
    
            minD=100.0
            for turn in turnlibrary:
                medoid_omega2=turnlibrary[turn]['medoid_omega2']
                medoid_phi2=turnlibrary[turn]['medoid_phi2']
                medoid_psi2=turnlibrary[turn]['medoid_psi2']
                medoid_omega3=turnlibrary[turn]['medoid_omega3']
                medoid_phi3=turnlibrary[turn]['medoid_phi3']
                medoid_psi3=turnlibrary[turn]['medoid_psi3']
                medoid_omega4=turnlibrary[turn]['medoid_omega4']
    
                D = (angle_distance(omega2, medoid_omega2) +
                     angle_distance(phi2, medoid_phi2) +
                     angle_distance(psi2, medoid_psi2) +
                     angle_distance(omega3, medoid_omega3) +
                     angle_distance(phi3, medoid_phi3) +
                     angle_distance(psi3, medoid_psi3) +
                     angle_distance(omega4, medoid_omega4))/7.0
                if D<minD:
                    minD = D
                    minturn=turn
                    
            if minD<= TURN_DIST_CUTOFF:
                nturns += 1
                resnum1 = int(datadict[chain][i]['auth_key'][1])
                resnum4 = resnum1+3
                prevname = turnlibrary[minturn]['prev_name']
                mintheta=find_theta_in_degrees(minD)
                seq = datadict[chain][i]['res1'] + datadict[chain][i+1]['res1'] + datadict[chain][i+2]['res1'] + datadict[chain][i+3]['res1']
                print(f'turn {nturns:4} {chain.id:4} {resnum1:4} {resnum4:4}    {seq:4} {dssp_string:4}    {minturn:5}  {prevname:13}    {minD:6.4f} {mintheta:7.2f} {CA1_CA4_distance:7.2f}    {omega2:7.2f} {phi2:7.2f} {psi2:7.2f}  {omega3:7.2f} {phi3:7.2f}  {psi3:7.2f} {omega4:7.2f}   {shortfilename}')
    
    
if __name__ == "__main__":
    main()
    
