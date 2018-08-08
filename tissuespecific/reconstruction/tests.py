#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 15:27:16 2017

@author: acabbia
"""

import unittest
from cobra.io import read_sbml_model


class Tests(unittest.TestCase):
     
     def test_EX_AA_present(self):
         model=read_sbml_model('/home/acabbia/Documents/Muscle_Model/tissuespecific/tissuespecific/models/SKM_latest_2017-10-26.xml')
         for aa in ['EX_tyr_L(e)',
                    'EX_trp_L(e)',
                    'EX_thr_L(e)',
                    'EX_leu_L(e)',
                    'EX_ile_L(e)',
                    'EX_lys_L(e)',
                    'EX_ala_L(e)',
                    'EX_val_L(e)',
                    'EX_met_L(e)',
                    'EX_ser_L(e)',
                    'EX_asp_L(e)',
                    'EX_gly(e)',
                    'EX_pro_L(e)',
                    'EX_cys_L(e)',
                    'EX_glu_L(e)',
                    'EX_gln_L(e)',
                    'EX_arg_L(e)',
                    'EX_phe_L(e)',
                    'EX_his_L(e)',
                    'EX_asn_L(e)']:
              
              self.assertIn(aa, set(model.reactions))
              

         
         
         





if __name__ == '__main__':
    unittest.main()
                 