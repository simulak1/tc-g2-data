# G2-97 Atomization Energies Reference Data
# E_f12: CCSD(T)-F12 atomization energy in kJ/mol
# E_exp: Experimental atomization energy in kJ/mol
#
# Keys ("formula") match ccsd_energies.csv molecule names

G2_97_ATOMIZATION = [
    # Aluminum compounds
    {"formula": "AlCl3", "name": "Aluminum trichloride", "E_f12": 1311.49, "E_exp": 1310.70},
    {"formula": "AlF3", "name": "Aluminum trifluoride", "E_f12": 1808.29, "E_exp": 1807.99},
    # Boron compounds
    {"formula": "BCl3", "name": "Boron trichloride", "E_f12": 1352.15, "E_exp": 1359.87},
    {"formula": "BF3", "name": "Boron trifluoride", "E_f12": 1963.98, "E_exp": 1972.30},
    # Beryllium compounds
    {"formula": "BeH", "name": "Beryllium monohydride", "E_f12": 212.50, "E_exp": 212.50},
    # Carbon tetrachloride/fluoride
    {"formula": "CCl4", "name": "Tetrachloromethane", "E_f12": 1312.24, "E_exp": 1318.96},
    {"formula": "CF4", "name": "Tetrafluoromethane", "E_f12": 2000.29, "E_exp": 2005.71},
    # CH species
    {"formula": "CH", "name": "Methylidyne radical", "E_f12": 350.91, "E_exp": 351.60},
    {"formula": "CH2_singlet", "name": "Singlet carbene", "E_f12": 755.71, "E_exp": 757.45},
    {"formula": "CH2", "name": "Triplet carbene", "E_f12": 793.88, "E_exp": 797.23},
    {"formula": "CH3", "name": "Methyl radical", "E_f12": 1282.77, "E_exp": 1287.21},
    {"formula": "CH4", "name": "Methane", "E_f12": 1752.56, "E_exp": 1757.82},
    # CH2X2 halomethanes
    {"formula": "CH2Cl2", "name": "Dichloromethane", "E_f12": 1548.52, "E_exp": 1554.42},
    {"formula": "CH2F2", "name": "Difluoromethane", "E_f12": 1827.75, "E_exp": 1832.86},
    # CH2O, CH2O2
    {"formula": "CH2O", "name": "Formaldehyde", "E_f12": 1561.51, "E_exp": 1567.43},
    {"formula": "CH2O2", "name": "Formic acid", "E_f12": 2093.28, "E_exp": 2100.94},
    # CH3X halomethanes
    {"formula": "CH3Cl", "name": "Chloromethane", "E_f12": 1650.66, "E_exp": 1656.17},
    # CH3O isomers
    {"formula": "CH3O_methylene_hydroxide", "name": "Hydroxymethyl radical", "E_f12": 1709.38, "E_exp": 1715.47},
    {"formula": "CH3O", "name": "Methoxy radical", "E_f12": 1668.97, "E_exp": 1674.66},
    # CH3S
    {"formula": "CH3S", "name": "Methylthio radical", "E_f12": 1595.13, "E_exp": 1601.03},
    # CH4O, CH4S
    {"formula": "CH4O", "name": "Methanol", "E_f12": 2141.49, "E_exp": 2147.84},
    {"formula": "CH4S", "name": "Thiomethanol", "E_f12": 1981.10, "E_exp": 1987.63},
    # CHX3 trihalomethanes
    {"formula": "CHCl3", "name": "Trichloromethane", "E_f12": 1438.15, "E_exp": 1444.48},
    {"formula": "CHF3", "name": "Trifluoromethane", "E_f12": 1919.52, "E_exp": 1924.71},
    # CHO
    {"formula": "CHO", "name": "Formyl radical", "E_f12": 1164.33, "E_exp": 1169.59},
    # CN, HCN
    {"formula": "CN", "name": "Cyano radical", "E_f12": 753.39, "E_exp": 758.56},
    {"formula": "CNH", "name": "Hydrogen cyanide", "E_f12": 1303.41, "E_exp": 1310.97},
    # CNH3O2 isomers
    {"formula": "CNH3O2_methyl_nitrite", "name": "Methyl nitrite", "E_f12": 2499.90, "E_exp": 2507.94},
    {"formula": "CNH3O2", "name": "Nitromethane", "E_f12": 2508.86, "E_exp": 2518.40},
    # CNH5 (methylamine)
    {"formula": "CNH5", "name": "Methylamine", "E_f12": 2427.06, "E_exp": 2434.89},
    # C2 species
    {"formula": "C2Cl4", "name": "Tetrachloroethylene", "E_f12": 1963.63, "E_exp": 1976.57},
    {"formula": "C2F4", "name": "Tetrafluoroethylene", "E_f12": 2458.80, "E_exp": 2470.91},
    {"formula": "C2H", "name": "Ethynyl radical", "E_f12": 1103.15, "E_exp": 1112.08},
    {"formula": "C2H2", "name": "Acetylene", "E_f12": 1685.16, "E_exp": 1695.79},
    {"formula": "C2H2O", "name": "Ketene", "E_f12": 2220.04, "E_exp": 2232.05},
    {"formula": "C2H2O2", "name": "Glyoxal", "E_f12": 2645.79, "E_exp": 2657.71},
    {"formula": "C2H3", "name": "Vinyl radical", "E_f12": 1855.58, "E_exp": 1864.92},
    {"formula": "C2H3Cl", "name": "Vinyl chloride", "E_f12": 2265.70, "E_exp": 2276.50},
    {"formula": "C2H3F", "name": "Vinyl fluoride", "E_f12": 2389.99, "E_exp": 2400.61},
    {"formula": "C2H3O", "name": "Carbonyl methane", "E_f12": 2426.20, "E_exp": 2436.67},
    {"formula": "C2H3OCl", "name": "Acetyl chloride", "E_f12": 2787.13, "E_exp": 2798.87},
    {"formula": "C2H3OF", "name": "Acetyl fluoride", "E_f12": 2946.26, "E_exp": 2957.91},
    {"formula": "C2H4", "name": "Ethylene", "E_f12": 2348.65, "E_exp": 2358.75},
    # C2H4O isomers
    {"formula": "C2H4O", "name": "Acetaldehyde", "E_f12": 2824.22, "E_exp": 2835.43},
    {"formula": "C2H4O_oxirane", "name": "Oxirane", "E_f12": 2714.69, "E_exp": 2725.51},
    # C2H4O2 isomers
    {"formula": "C2H4O2", "name": "Acetic acid", "E_f12": 3350.83, "E_exp": 3363.76},
    {"formula": "C2H4O2_methyl_formate", "name": "Methyl formate", "E_f12": 3280.39, "E_exp": 3292.61},
    # C2H4S
    {"formula": "C2H4S", "name": "Thiirane", "E_f12": 2609.09, "E_exp": 2620.57},
    # C2H5 species
    {"formula": "C2H5", "name": "Ethyl radical", "E_f12": 2516.37, "E_exp": 2526.01},
    {"formula": "C2H5Cl", "name": "Ethyl chloride", "E_f12": 2887.09, "E_exp": 2897.60},
    {"formula": "C2H5O", "name": "Ethoxy radical", "E_f12": 2908.38, "E_exp": 2919.13},
    # C2H6
    {"formula": "C2H6", "name": "Ethane", "E_f12": 2971.49, "E_exp": 2981.64},
    # C2H6O isomers
    {"formula": "C2H6O_dimethyl_ether", "name": "Dimethyl ether", "E_f12": 3330.63, "E_exp": 3341.56},
    {"formula": "C2H6O", "name": "Ethanol", "E_f12": 3381.59, "E_exp": 3392.98},
    # C2H6OS
    {"formula": "C2H6OS", "name": "Dimethyl sulfoxide", "E_f12": 3574.86, "E_exp": 3587.65},
    # C2H6S isomers
    {"formula": "C2H6S_dimethyl_thioether", "name": "Dimethyl sulfide", "E_f12": 3206.19, "E_exp": 3217.58},
    {"formula": "C2H6S", "name": "Thioethanol", "E_f12": 3210.77, "E_exp": 3222.30},
    # C2N2
    {"formula": "C2N2", "name": "Cyanogen", "E_f12": 2084.80, "E_exp": 2100.20},
    {"formula": "C2NF3", "name": "Trifluoroacetonitrile", "E_f12": 2675.64, "E_exp": 2688.27},
    {"formula": "C2NH3", "name": "Acetonitrile", "E_f12": 2565.94, "E_exp": 2578.87},
    {"formula": "C2NH5", "name": "Aziridine", "E_f12": 3000.47, "E_exp": 3012.94},
    {"formula": "C2NH5O", "name": "Acetamide", "E_f12": 3618.80, "E_exp": 3633.73},
    # C2NH7 isomers
    {"formula": "C2NH7", "name": "Dimethylamine", "E_f12": 3627.84, "E_exp": 3640.55},
    {"formula": "C2NH7_ethyl_amine", "name": "Ethylamine", "E_f12": 3662.09, "E_exp": 3674.87},
    # C3H4 isomers
    {"formula": "C3H4", "name": "Allene", "E_f12": 2928.72, "E_exp": 2944.45},
    {"formula": "C3H4_cyclopropene", "name": "Cyclopropene", "E_f12": 2836.78, "E_exp": 2852.25},
    {"formula": "C3H4_propyne", "name": "Propyne", "E_f12": 2934.48, "E_exp": 2950.50},
    # C3H6 isomers
    {"formula": "C3H6_cyclopropane", "name": "Cyclopropane", "E_f12": 3555.66, "E_exp": 3571.30},
    {"formula": "C3H6", "name": "Propene", "E_f12": 3587.40, "E_exp": 3602.67},
    # C3H6O
    {"formula": "C3H6O", "name": "Acetone", "E_f12": 4080.98, "E_exp": 4097.39},
    # C3H7
    {"formula": "C3H7", "name": "Isopropyl radical", "E_f12": 3755.06, "E_exp": 3769.88},
    {"formula": "C3H7Cl", "name": "1-Chloropropane", "E_f12": 4115.27, "E_exp": 4130.72},
    # C3H8
    {"formula": "C3H8", "name": "Propane", "E_f12": 4199.18, "E_exp": 4214.27},
    # C3H8O isomers
    {"formula": "C3H8O", "name": "Methoxyethane", "E_f12": 4571.15, "E_exp": 4587.12},
    {"formula": "C3H8O_isopropyl_alcohol", "name": "Isopropyl alcohol", "E_f12": 4624.64, "E_exp": 4641.05},
    # C3N
    {"formula": "C3NH3", "name": "Acrylonitrile", "E_f12": 3173.61, "E_exp": 3191.46},
    {"formula": "C3NH9", "name": "Trimethylamine", "E_f12": 4840.35, "E_exp": 4857.95},
    # C4H6 isomers
    {"formula": "C4H6_1_3_butadiene", "name": "1,3-Butadiene", "E_f12": 4217.43, "E_exp": 4237.73},
    {"formula": "C4H6_2_butyne", "name": "2-Butyne", "E_f12": 4179.82, "E_exp": 4201.18},
    {"formula": "C4H6_bicyclo_1_1_0_butane", "name": "Bicyclo[1.1.0]butane", "E_f12": 4105.06, "E_exp": 4125.62},
    {"formula": "C4H6", "name": "Cyclobutene", "E_f12": 4170.00, "E_exp": 4189.71},
    {"formula": "C4H6_methylene_cyclopropane", "name": "Methylenecyclopropane", "E_f12": 4133.88, "E_exp": 4154.81},
    # C4H8 isomers
    {"formula": "C4H8", "name": "Cyclobutane", "E_f12": 4791.64, "E_exp": 4811.19},
    {"formula": "C4H8_isobutene", "name": "Isobutene", "E_f12": 4829.32, "E_exp": 4849.73},
    # C4H9
    {"formula": "C4H9", "name": "tert-Butyl radical", "E_f12": 4996.68, "E_exp": 5016.63},
    # C4H10 isomers
    {"formula": "C4H10_isobutane", "name": "Isobutane", "E_f12": 5433.13, "E_exp": 5453.18},
    {"formula": "C4H10", "name": "n-Butane", "E_f12": 5427.17, "E_exp": 5447.22},
    # C4 heterocycles
    {"formula": "C4H4O", "name": "Furan", "E_f12": 4142.37, "E_exp": 4164.49},
    {"formula": "C4H4S", "name": "Thiophene", "E_f12": 4016.15, "E_exp": 4038.49},
    {"formula": "C4NH5", "name": "Pyrrole", "E_f12": 4464.55, "E_exp": 4488.94},
    # C5
    {"formula": "C5H8", "name": "Spiropentane", "E_f12": 5350.47, "E_exp": 5376.77},
    {"formula": "C5NH5", "name": "Pyridine", "E_f12": 5155.83, "E_exp": 5183.94},
    # C6
    {"formula": "C6H6", "name": "Benzene", "E_f12": 5696.88, "E_exp": 5727.81},
    # Halogens
    {"formula": "Cl2", "name": "Dichlorine", "E_f12": 247.16, "E_exp": 248.22},
    {"formula": "F2", "name": "Difluorine", "E_f12": 162.15, "E_exp": 162.31},
    {"formula": "FCl", "name": "Chlorine monofluoride", "E_f12": 261.80, "E_exp": 262.43},
    {"formula": "F3Cl", "name": "Chlorine trifluoride", "E_f12": 537.00, "E_exp": 537.44},
    # CO species
    {"formula": "CO", "name": "Carbon monoxide", "E_f12": 1083.14, "E_exp": 1087.57},
    {"formula": "CO2", "name": "Carbon dioxide", "E_f12": 1625.87, "E_exp": 1633.95},
    {"formula": "COF2", "name": "Carbonyl fluoride", "E_f12": 1755.19, "E_exp": 1762.10},
    {"formula": "COS", "name": "Carbonyl sulfide", "E_f12": 1398.66, "E_exp": 1406.74},
    # CS species
    {"formula": "CS", "name": "Carbon monosulfide", "E_f12": 715.40, "E_exp": 719.47},
    {"formula": "CS2", "name": "Carbon disulfide", "E_f12": 1165.93, "E_exp": 1174.07},
    # Hydrogen halides
    {"formula": "HCl", "name": "Hydrogen chloride", "E_f12": 448.52, "E_exp": 449.58},
    {"formula": "HF", "name": "Hydrogen fluoride", "E_f12": 592.09, "E_exp": 593.02},
    {"formula": "HOCl", "name": "Hypochlorous acid", "E_f12": 693.67, "E_exp": 695.38},
    # OH, SH
    {"formula": "HO", "name": "Hydroxyl radical", "E_f12": 447.50, "E_exp": 448.30},
    {"formula": "HS", "name": "Mercapto radical", "E_f12": 365.64, "E_exp": 366.55},
    # H2 species
    {"formula": "H2", "name": "Dihydrogen", "E_f12": 457.73, "E_exp": 457.73},
    {"formula": "H2O", "name": "Water", "E_f12": 973.05, "E_exp": 974.94},
    {"formula": "H2O2", "name": "Hydrogen peroxide", "E_f12": 1124.05, "E_exp": 1126.34},
    {"formula": "H2S", "name": "Hydrogen sulfide", "E_f12": 766.92, "E_exp": 768.72},
    # Lithium
    {"formula": "Li2", "name": "Dilithium", "E_f12": 101.24, "E_exp": 101.24},
    {"formula": "LiF", "name": "Lithium fluoride", "E_f12": 583.12, "E_exp": 583.99},
    {"formula": "LiH", "name": "Lithium hydride", "E_f12": 242.27, "E_exp": 242.27},
    # Sodium
    {"formula": "Na2", "name": "Disodium", "E_f12": 71.55, "E_exp": 71.78},
    {"formula": "NaCl", "name": "Sodium chloride", "E_f12": 411.98, "E_exp": 412.96},
    # Nitrogen species
    {"formula": "N2", "name": "Dinitrogen", "E_f12": 951.59, "E_exp": 955.82},
    {"formula": "N2H4", "name": "Hydrazine", "E_f12": 1827.15, "E_exp": 1832.69},
    {"formula": "N2O", "name": "Nitrous oxide", "E_f12": 1127.48, "E_exp": 1133.70},
    {"formula": "NF3", "name": "Trifluoroamine", "E_f12": 862.79, "E_exp": 863.68},
    {"formula": "NH", "name": "Imidogen", "E_f12": 346.38, "E_exp": 347.02},
    {"formula": "NH2", "name": "Amino radical", "E_f12": 761.30, "E_exp": 762.95},
    {"formula": "NH3", "name": "Ammonia", "E_f12": 1242.94, "E_exp": 1245.99},
    {"formula": "NO", "name": "Nitric oxide", "E_f12": 636.75, "E_exp": 639.28},
    {"formula": "NO2", "name": "Nitrogen dioxide", "E_f12": 950.01, "E_exp": 954.10},
    {"formula": "NOCl", "name": "Nitrosyl chloride", "E_f12": 801.10, "E_exp": 803.43},
    # Oxygen species
    {"formula": "O2", "name": "Dioxygen", "E_f12": 504.36, "E_exp": 505.88},
    {"formula": "O2S", "name": "Sulfur dioxide", "E_f12": 1086.86, "E_exp": 1091.61},
    {"formula": "O3", "name": "Ozone", "E_f12": 614.29, "E_exp": 615.78},
    {"formula": "OCl", "name": "Monochlorine monoxide", "E_f12": 270.00, "E_exp": 271.20},
    {"formula": "OF2", "name": "Difluorine monoxide", "E_f12": 392.34, "E_exp": 392.68},
    {"formula": "OS", "name": "Sulfur monoxide", "E_f12": 526.33, "E_exp": 528.72},
    # Phosphorus
    {"formula": "P2", "name": "Diphosphorus", "E_f12": 485.14, "E_exp": 489.29},
    {"formula": "PF3", "name": "Phosphorus trifluoride", "E_f12": 1527.81, "E_exp": 1530.92},
    {"formula": "PH2", "name": "Phosphino radical", "E_f12": 644.21, "E_exp": 645.47},
    {"formula": "PH3", "name": "Phosphane", "E_f12": 1010.33, "E_exp": 1012.24},
    # Sulfur
    {"formula": "S2", "name": "Disulfur", "E_f12": 431.42, "E_exp": 433.94},
    # Silicon species
    {"formula": "Si2", "name": "Disilicon", "E_f12": 307.15, "E_exp": 307.75},
    {"formula": "Si2H6", "name": "Disilane", "E_f12": 2240.41, "E_exp": 2240.54},
    {"formula": "SiCH6", "name": "Methylsilane", "E_f12": 2626.11, "E_exp": 2631.54},
    {"formula": "SiCl4", "name": "Silicon tetrachloride", "E_f12": 1625.41, "E_exp": 1627.56},
    {"formula": "SiF4", "name": "Silicon tetrafluoride", "E_f12": 2416.45, "E_exp": 2419.74},
    {"formula": "SiH2", "name": "Singlet silylene", "E_f12": 642.99, "E_exp": 643.11},
    {"formula": "SiH2_triplet", "name": "Triplet silylene", "E_f12": 557.56, "E_exp": 555.72},
    {"formula": "SiH3", "name": "Silyl radical", "E_f12": 954.30, "E_exp": 953.44},
    {"formula": "SiH4", "name": "Silane", "E_f12": 1358.08, "E_exp": 1357.91},
    {"formula": "SiO", "name": "Silicon monoxide", "E_f12": 804.83, "E_exp": 809.19},
]
