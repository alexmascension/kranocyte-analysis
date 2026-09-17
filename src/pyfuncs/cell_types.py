
DICT_MARKERS_MAJOR_POPULATIONS = {
    'Kranocyte': ['Rasgrp2', 'Tenm2', 'Inpp4b', 'Foxd2os', 'Malt1', 'Gfra2', 'Shisa3', 'Malt1', 'Thrsp', 'Gpld1', 'Smim41', 'Plxdc1', 
                  'Dlk1', 'Fetub', 'Saa1', 'Gria1', 'Greb1', 'Col9a2', 'Gli1', 'Cst6'], 
    'FAP': ['Ly6a', 'Nova1', 'Fstl1', 'Ifi205', 'Slfn5', 'Fbln2', 'Dpep1', 'Fn1', 'Pdgfra', 'Lbp', 'Ifi204'],
    'Tenocyte': ['Cilp2', 'Chad', 'Scx', 'Mkx', 'Rflnb', 'Ptx4', 'Gas2', 'Kctd1', 'Edil3', 'Col11a2', 'Sox6', 'Scube2', 'Cdh2', 'Matn4'],
    'Satellite': ['Dmd', 'Pax7', 'Chodl', 'Fgfr4', 'Notch3', 'Tenm4', 'Peg3', 'Fry'], 
    'Endothelial': ['Egfl7', 'Ptprb', 'Cdh5', 'Vwf', 'Esam', 'Cd36', 'Flt1', 'Tspan13', 'Emcn'],
    'Immune': ['Lyz2', 'Ckb', 'Tyrobp', 'Iqgap1', 'Laptm5', 'Ctsc', 'Ehd4', 'Ctss', 'Fcer1g', 'Ccl4', 'C1qb',] 
}


DICT_MARKERS_MAJOR_POPULATIONS = {
 
    # ---------------- las que ya tenias (sin tocar) ----------------
    'Kranocyte': ['Rasgrp2', 'Tenm2', 'Inpp4b', 'Foxd2os', 'Malt1', 'Gfra2', 'Shisa3', 'Thrsp',
                  'Gpld1', 'Smim41', 'Plxdc1', 'Dlk1', 'Fetub', 'Saa1', 'Gria1', 'Greb1',
                  'Col9a2', 'Gli1', 'Cst6'],
 
    'FAP': ['Pdgfra', 'Ly6a', 'Dcn', 'Lum', 'Col3a1', 'Cd34', 'Fbln2', 'Fstl1', 'Gsn',
            'Nova1', 'Ifi205', 'Slfn5', 'Dpep1', 'Fn1', 'Lbp', 'Ifi204', 'Cilp', 'Lox', 'Postn'],
 
    'Tenocyte': ['Scx', 'Tnmd', 'Mkx', 'Thbs4', 'Comp', 'Fmod', 'Cilp2', 'Chad', 'Rflnb',
                 'Ptx4', 'Gas2', 'Kctd1', 'Edil3', 'Col11a2', 'Sox6', 'Scube2', 'Cdh2', 'Matn4'],
 
    'Satellite': ['Pax7', 'Myf5', 'Vcam1', 'Calcr', 'Sdc4', 'Chodl', 'Dmd', 'Fgfr4',
                  'Notch3', 'Tenm4', 'Peg3', 'Fry'],
 
    'Endothelial': ['Pecam1', 'Cdh5', 'Egfl7', 'Ptprb', 'Vwf', 'Esam', 'Flt1', 'Emcn',
                    'Tspan13', 'Cd36'],
 
    'Immune': ['Ptprc', 'Lyz2', 'Tyrobp', 'Fcer1g', 'Laptm5', 'Ctss', 'Ctsc', 'C1qb',
               'Ccl4', 'Ckb', 'Iqgap1', 'Ehd4'],
 
    # ---------------- poblaciones que faltaban ----------------
 
    # Mural / vascular no-endotelial
    'Pericyte': ['Rgs5', 'Kcnj8', 'Abcc9', 'Pdgfrb', 'Notch3', 'Vtn', 'Higd1b',
                 'Ndufa4l2', 'Cspg4', 'Colec11'],                       # [NIC][SOU][OPR]
 
    'SMC': ['Myh11',  'Cnn1', 'Myl9', 'Lmod1', 'Pln', 'Des',
            'Actg2', ],                                            # [SOU][SONG]
 
    'SMMC': ['Myh11', 'Itga7', 'Ntrk3', 'Notch3', 'Pdgfrb', 'Rgs5',
             'Cd248', 'Des'],                                            # [NIC] smooth muscle-mesenchymal cells
 
    'Lymphatic_EC': ['Prox1', 'Lyve1', 'Pdpn', 'Flt4', 'Ccl21a', 'Mmrn1', 'Reln',
                     'Nts', 'Fgl2', 'Thy1'],                             # subtipo EC, ver seccion 5
 
    # Linaje neural / glial
    'Glial_Schwann': ['Plp1', 'Mpz', 'Sox10', 'S100b', 'Mbp', 'Ptn', 'Kcna1', 'Cdh19',
                      'Ngfr', 'Mal'],                                    # [NIC][SOU][SONG]
 
    # # Miofibra / nucleos musculares
    # 'Myonuclei': ['Acta1', 'Ckm', 'Tnnt3', 'Tnnc2', 'Mylpf', 'Actn3', 'Neb', 'Ttn',
    #               'Atp2a1', 'Mb'],                                       # [OPR][NIC]
 
    # 'Myotenocyte': ['Tnnc2', 'Scx', 'Col22a1', 'Ankrd1', 'Tnmd', 'Mkx', 'Acta1',
    #                 'Csrp3', 'Thbs4'],                                   # [NIC] union miotendinosa
 
    # # Mesenquimal terminal
    # 'Adipocyte': ['Adipoq', 'Plin1', 'Cfd', 'Fabp4', 'Lep', 'Car3', 'Cidec', 'Retn',
    #               'Lpl', 'Scd1'],                                        # [SOU]
 
    # Inmunes que no separabas
    'Neutrophil': ['S100a8', 'S100a9', 'Retnlg', 'Csf3r', 'Il1r2', 'Mmp9', 'Lcn2',
                   'Cxcr2', 'Msrb1', 'Wfdc21'],                          # [OPR][SOU]
 
    'Mast_cell': ['Cpa3', 'Cma1', 'Mcpt4', 'Kit', 'Fcer1a', 'Ms4a2', 'Cd200r3',
                  'Il1rl1', 'Hdc', 'Slc18a2'],                           # [SOU][SONG]
 
    'T_NK': ['Cd3e', 'Cd3d', 'Cd3g', 'Lck', 'Il7r', 'Trbc2', 'Nkg7', 'Gzma',
             'Ncr1', 'Klrd1'],                                           # [OPR][SOU][SONG]
 
    'B_cell': ['Cd79a', 'Cd79b', 'Ms4a1', 'Ighm', 'Ighd', 'Cd19', 'Ebf1', 'Bank1',
               'H2-Ob', 'Iglc1'],                                        # [SOU]
 
    'Dendritic': ['Cd74', 'H2-Ab1', 'H2-Aa', 'H2-Eb1', 'Cd209a', 'Napsa', 'Flt3',
                  'Tmem176a', 'Tmem176b', 'Ccr7'],                       # [OPR][SOU]
 
    # # Poblaciones raras / contexto especifico (isquemia, amputacion, osificacion)
    # 'Erythrocyte': ['Hba-a1', 'Hba-a2', 'Hbb-bs', 'Hbb-bt', 'Alas2', 'Slc4a1',
    #                 'Bpgm', 'Snca'],                                     # [SOU]
 
    # 'Osteoclast': ['Acp5', 'Ctsk', 'Mmp9', 'Dcstamp', 'Ocstamp', 'Atp6v0d2', 'Nfatc1',
    #                'Car2', 'Src'],                                       # [SOU]
 
    # 'Epithelial_basal': ['Krt5', 'Krt14', 'Krt10', 'Epcam', 'Sfn', 'Dsp', 'Krt17',
    #                      'Col17a1'],                                     # [SOU][SONG] (piel en amputaciones)
}





DICT_MARKERS_TENO = {
    'Teno_D': ['Sparcl1', 'F2r', 'Col4a2', 'Ptpre', 'Col22a1', 'Fbn2', 'Col18a1', 'Rapgef4', 'Lrrn2', 'Reln', 'Adam23',],
    'Teno_C': ['Ccdc3', 'Bmpr1b', 'Cav1', 'Wif1', 'Kera', 'C3', 'Sema3b', 'Sema3c', 'Chrdl1', 'Grem2', 'Itga2', 'Hpgd', 'Bmp3', 'Clu'],
    'Teno_B': ['Gpx3', 'Col8a1', 'Bicc1', 'Runx2', 'Smoc2', 'Camk4', 'Spon1', 'Wnt16', 'Slit2', 'Ostn', 'H19', 'Igf2'],
    'Teno_A': ['Abca8a', 'Celf2', 'Mfap5', 'Col14a1', 'Dpt', 'Abca8b', 'Lum', 'Cd34']
}

DICT_MARKERS_SATELLITE = {
    'Sat_A1': ['Pde3a', 'Frem1', 'Pde10a', 'Bric2', 'Cdk14', 
               'Rora', 'Esr1', 'Ror1', 'Dmd', 'Ptpng', 
               'Col19a1', 'Tnik', 'Bicc1', 'Pde7b', 'Bnc2', 
               'Kcnma1'],   # No clear markers
    'Sat_A2': ['Rps16', 'Rps3a1', 'Rps10', 'Rps19', 'Rplp0', 'Rps2', 
               'Rpl15', 'Rpl30', 'Rpsa'],   # No clear markers
    'Sat_B': ['Bgn', 'Dcn', 'Col1a2', 'Col1a1', 'Dlc1', 'Pcolce', 'Nid1', 'Rnase4', 'Spry2', 'Serpinf1', 'Clec3b', 'Col6a3', 'Mfap5', 'Lum', 'C1s1', 'Ly6a', 'Pdgfra', 
              'Phldb1', 'Ly6c1', 'Loxl1', 'Meox2', 'Vcan', 'B4galt1', 'Svil', 'Gpc3', 'Tns1', 'Dpysl3', 
              ]
}


DICT_MARKERS_KRANOCYTE = {
    'Krano_A': ['Dlk1', 'Saa1', 'Fndc1', 'Hmcn2', 'Gpx3', 'Aldh1a2', 'Tnfaip2', 'Dkk2', 'Tspan11', 'Cfb', 'Mn1', 'Col14a1', 'Syk', 'Spon1', 'Ly6a'], 
    'Krano_B': ['Wnt6', 'Cldn1', 'Stra6', 'Sbspon', 'Dach1', 'Ndnf', 'Efemp1', 'Fxyd6', 'Prdm1', 'Gfra1', 'Bmp7', 'Rapgef4', 'Negr1'],
    'Krano_C': ['Nipal1', 'Vit', 'Tenm2', 'Apod', 'Thbs4', 'Gpc3', 'Cp', 'Piezo2', 'Spp1', 'Cyp26b1'],
    }

DICT_INITIAL_KRANO_POPS = {
    'A': ['6030408B16Rik', 'Smim41', 'Col9a2', 'Dlk1', 'Shisa3',  'Saa1',  'Nipal1'],
    'B': ['Lypd2', 'Wnt6', 'Cldn1', 'Moxd1', 'Mansc4', 'Dleu7', 'Efnb3', 'Stra6', 'Sbspon', 'Ace2', 'Hcn4', 'Cldn22', 'Wnt10a', 'Ocln'],
}
LIST_INITIAL_KRANO_POPS_COLORS = ['#bcbcbc', '#900C3F', '#286e87']

DICT_MARKERS_FAP = {
    'FAP_A': ['Sbsn', 'Efna5', 'Adamts16', 'Dact2', 'Krtdap', 'Aldh1a3', 'Ildr2', 
              'Car8', 'Rorb', 'Uchl1', 'Bmp3', 'Sema3e', 'Dmkn', 'Ano3', 'Itgb7'],
    'FAP_AB': ['Cmah', 'Stmn4', 'Robo1', 'Fam167a', 'Duox1', 'Gm33624', 
               'Cd55', 'Calm1', 'Ugp2', 'Efoc1', 'Prdm8', 'Kcp', 'Des'],
    'FAP_AF': ['C3', 'Postn', 'Slc4a4', 'Ifi27l2a', 'Il1rl2', 'Smpd3', 
               'Il1r2', 'Sfrp2', 'Vegfd', 'Cxcl13', 'Prkg2', 'C7', 
               'Cpz', 'Cdh9', 'Edn1', 'Sema3b', 'Mmp27'],
    
    'FAP_B1': ['Grm8', 'Rora', 'Plxdc2', 'Ank2', 'Pcdh7', 'Eda', 'Tnik', 'Adgrl3', 'Elmo1'],
    'FAP_B2': ['Cxcl14', 'Hsd11b1', 'Smoc2', 'Col15a1', 'Mme', 'G0s2', 'Col4a1', 'Lamb1', 'Col4a2', 
               'Cldn15', 'Clec14a', 'Nmb', 'Vwa1', 'Crlf1',],
    'FAP_C': ['Gdf10', 'C2', 'Cygb', 'Fmo2', 'Nr2f2', 'Csmd1', 'Prdm6', 'Abcc9', 'Gria4', 'Clec11a', 'Ism1', 
               'Gata6', 'Pltp', 'Hhip', 'Tspan11', 'Tent5c'],
    'FAP_D': ['Steap4', 'Apoe', 'Cyp1b1', 'Aoc3', 'Ccl19-ps4', 'Sfrp1', 'Adam12', 'Emb', 
               'Agt', 'Pparg',  ],
    'FAP_E': ['Mgp', 'Hmcn1', 'Meox1', 'Meox2', 'Clu', 'Etl4', 'Kctd12', 'Boc', 'Daam2', 'Matn2', 
              'Robo2', 'Clec1a', 'Myo10', 'Sgip1', 'Myo1b', 'Ptn', 'Mettl24', 'Rgs7bp'],
    'FAP_F': ['Thbs2', 'Gpx3', 'Cilp', 'Fgl2', 'Lox', 'Aspn', 'Ccdc80', 'Igfbp3', 'Dkk2', 'Ecrg4', 'Tnmd', 
              'Col12a1', 'Pappa2', 'Egfl6', 'Fibin', 'Ccn5', 'Sfrp4', 'Fxyd6', ],
    }


# We apply this 
DICT_RENAMING = {
    'FAP_A': 'FAP.1',
    'FAP_AB': 'FAP.1',
    'FAP_AF': 'FAP.1',
    'FAP_B1': 'FAP.4',    
    'FAP_B2': 'FAP.4',
    'FAP_E': 'FAP.2',
    'FAP_C': 'FAP.3',
    'FAP_D': 'FAP.3',
    'FAP_F': 'FAP.3',
    'Krano_A': 'FAP.5',
    'Krano_B': 'FAP.6',
    'Krano_C': 'FAP.6',
    'Endothelial*': 'ENDO',
    'Sat_A1': 'SAT', 
    'Sat_A2': 'SAT', 
    'Sat_B': 'SAT', 
    'Teno_A': 'TNMD', 
    'Teno_D': 'TNMD',
    'Teno_B': 'TNMD', 
    'Teno_C': 'TNMD', 

    # Other populations
    "B_cell": "Immune", 
    "Dendritic": "Immune",
    "Lymphatic_EC": "Immune",
    "Mast_cell": "Immune",
    "T_NK": "Immune",
    "Neutrophil": "Immune",

    'Sat_U': 'SAT', 
    'Krano_U': 'FAP.5',
    "SMC": "SMC-SMMC",
    "SMMC": "SMC-SMMC",
    "Glial_Schwann": "GLIA",
    'Pericyte': 'ENDO',
                               
}



PALETTE_CELL_TYPE = {
    # ---- Minor (más oscuros) ----
    # FAP (7)
    'FAP.1':"#A1D9F0",
    'FAP.2' :"#4C92AD",
    'FAP.3':"#12485E",

    'FAP.4':"#B2A1F0",
    'FAP.5' :"#694CAD",
    'FAP.6':"#26125E",

    "TNMD":"#6F5214",
    'SAT':"#2E9277",
    'ENDO':   "#C23C64",
    'SMC-SMMC':      "#79DD6F",
    'GLIA':      "#EE2561",
    'Immune':      "#E766E7",
}




# We apply this 
DICT_RENAMING_INTEGRATION = {
    'FAP_A': 'FAP.1',
    'FAP_AB': 'FAP.1',
    'FAP_AF': 'FAP.1',
    'FAP_B1': 'FAP.4',    
    'FAP_B2': 'FAP.4',
    'FAP_E': 'FAP.2',
    'FAP_C': 'FAP.3',
    'FAP_D': 'FAP.3',
    'FAP_F': 'FAP.3',
    'Krano_A': 'FAP.5',
    'Krano_B': 'FAP.6',
    'Krano_C': 'FAP.6',
    'Endothelial*': 'ENDO',
    'Sat_A1': 'SAT', 
    'Sat_A2': 'SAT', 
    'Sat_B': 'SAT', 
    'Teno_A': 'TNMD', 
    'Teno_D': 'TNMD',
    'Teno_B': 'TNMD', 
    'Teno_C': 'TNMD', 

    # Other populations
    "Immune": "IMM.MONO-MAC", 
    "B_cell": "IMM.B", 
    "Dendritic": "IMM.DEN",
    "Mast_cell": "IMM.MAST",
    "T_NK": "IMM.NK",
    "Neutrophil": "IMM.NEU",

    'Sat_U': 'SAT', 
    'Krano_U': 'FAP.5',
    "SMC": "SMC-SMMC",
    "SMMC": "SMC-SMMC",
    "Glial_Schwann": "GLIA",
    "Lymphatic_EC": "ENDO.LYMPH",
    'Pericyte': 'ENDO.PERI',
                               
}




PALETTE_CELL_TYPE_INTEGRATION = {
    # ---- Minor (más oscuros) ----
    # FAP (7)
    'FAP.1':"#A1D9F0",
    'FAP.2' :"#4C92AD",
    'FAP.3':"#12485E",
    'FAP.4':"#B2A1F0",
    'FAP.5' :"#694CAD",
    'FAP.6':"#26125E",

    "TNMD":"#6F5214",
    'SAT':"#2E9277",
    'SMC-SMMC':      "#79DD6F",
    'GLIA':      "#CAC031",

    'ENDO':   "#C23C64",
    'ENDO.LYMPH':   "#F8608D",
    'ENDO.PERI':   "#991F43",

    'IMM.MONO-MAC':      "#331903",
    'IMM.B':      "#D36B17",
    'IMM.NK':      "#F7A765",
    'IMM.DEN':      "#96490A",
    'IMM.MAST':      "#D35217",
    'IMM.NEU':      "#92390F",
}