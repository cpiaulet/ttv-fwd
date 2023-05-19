  
import numpy as np

esinw= [0.00217, 0.00001, 0.00267, -0.00461, -0.00051, 0.00128, -0.00002]
esinw_unc = [0.00244, 0.00171, 0.00112, 0.00087, 0.00087, 0.00070, 0.00044]

ecosw = [-0.00215, 0.00055, -0.00496, 0.00433, -0.00840, 0.00380, -0.00365]
ecosw_unc = [0.00332, 0.00232, 0.00186, 0.00149, 0.00130, 0.00112, 0.00077]

T0 = [7257.55044, 7258.58728, 7257.06768, 7257.82771, 7257.07426, 7267.71462, 7249.60676]
P_days = [1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866]
t_ref = 7257.93115525

first_tts = [7322.515310, 7282.805700, 7560.797300, 7312.713000, 7321.525200, 7294.786000, 7662.554360]
Mstar_sun = 0.0898
t_start_integ_days = 7257
t_end_integ_days = 7700


#nsamp = 40000000
nsamp = 1000
# CALCULATION OF THE ECCENTRICITY OF PLANET B
esinw_b = np.random.normal(esinw[0], esinw_unc[0], nsamp)
ecosw_b = np.random.normal(ecosw[0], ecosw_unc[0], nsamp)
ecc_b_s = np.sqrt((esinw_b)**2+(ecosw_b)**2) # Samples
ecc_b_unc = np.std(ecc_b_s)

# OF PLANET C
esinw_c = np.random.normal(esinw[1], esinw_unc[1], nsamp)
ecosw_c = np.random.normal(ecosw[1], ecosw_unc[1], nsamp)
ecc_c_s  = np.sqrt((esinw_c)**2+(ecosw_c)**2)  # Samples
ecc_c_unc = np.std(ecc_c_s)

# OF PLANET D
esinw_d = np.random.normal(esinw[2], esinw_unc[2], nsamp)
ecosw_d = np.random.normal(ecosw[2], ecosw_unc[2], nsamp)
ecc_d_s = np.sqrt((esinw_d)**2+(ecosw_d)**2)  # Samples
ecc_d_unc = np.std(ecc_d_s)

# OF PLANET E
esinw_e = np.random.normal(esinw[3], esinw_unc[3], nsamp)
ecosw_e = np.random.normal(ecosw[3], ecosw_unc[3], nsamp)
ecc_e_s = np.sqrt((esinw_e)**2+(ecosw_e)**2)  # Samples
ecc_e_unc = np.std(ecc_e_s)

# OF PLANET F
esinw_f = np.random.normal(esinw[-3], esinw_unc[-3], nsamp)
ecosw_f = np.random.normal(ecosw[-3], ecosw_unc[-3], nsamp)
ecc_f_s = np.sqrt((esinw_f)**2+(ecosw_f)**2)  # Samples
ecc_f_unc = np.std(ecc_f_s)

# OF PLANET G
esinw_g = np.random.normal(esinw[-2], esinw_unc[-2], nsamp)
ecosw_g = np.random.normal(ecosw[-2], ecosw_unc[-2], nsamp)
ecc_g_s = np.sqrt((esinw_g)**2+(ecosw_g)**2)  # Samples
ecc_g_unc = np.std(ecc_g_s)

# OF PLANET H
esinw_h = np.random.normal(esinw[-1], esinw_unc[-1], nsamp)
ecosw_h = np.random.normal(ecosw[-1], ecosw_unc[-1], nsamp)
ecc_h_s = np.sqrt((esinw_h)**2+(ecosw_h)**2)  # Samples
ecc_h_unc = np.std(ecc_h_s)

#############################################################################

# CALCULATION OF THE ARGUMENT OF PERIASTRON (omega) OF PLANET B
omega_b_s =np.arctan2(esinw_b,ecosw_b) * 180/np.pi
omega_b_unc = np.std(omega_b_s)

# OF PLANET C
omega_c_s =np.arctan2(esinw_c,ecosw_c) * 180/np.pi
omega_c_unc = np.std(omega_c_s)

# OF PLANET D
omega_d_s =np.arctan2(esinw_d,ecosw_d) * 180/np.pi
omega_d_unc = np.std(omega_d_s)

# OF PLANET E
omega_e_s =np.arctan2(esinw_e,ecosw_e) * 180/np.pi
omega_e_unc = np.std(omega_e_s)

# OF PLANET F
omega_f_s =np.arctan2(esinw_f,ecosw_f) * 180/np.pi
omega_f_unc = np.std(omega_f_s)

# OF PLANET G
omega_g_s =np.arctan2(esinw_g,ecosw_g) * 180/np.pi
omega_g_unc = np.std(omega_g_s)

# OF PLANET H
omega_h_s =np.arctan2(esinw_h,ecosw_h) * 180/np.pi
omega_h_unc = np.std(omega_h_s)

#############################################################################

# FUNCTIONS THAT CALCULATE THE MEAN ANOMALY (M) OF EACH PLANET 
def get_E(ecc, omega):
  term_1 = np.sqrt((1.-ecc)/(1.+ecc))
  term_2 = np.tan(0.5 * (np.pi/2. - omega))
  return (2. * np.arctan(term_1 * term_2))

def get_M0(t_ref, P_days, T0, ecc, omega):     # gives value of M0 in radians
  E = get_E(ecc, omega)
  M0 = (2. * np.pi)/P_days * ((t_ref - T0)%P_days) + E - (ecc * np.sin(E))
  return M0

def rad_to_degrees(angle):
  res = angle * 180/np.pi 
  return res

def get_M0_degrees(t_ref, P_days, T0, ecc, omega):
  M0_rad = get_M0(t_ref, P_days, T0, ecc, omega)
  M0_deg = rad_to_degrees(M0_rad)
  return M0_deg

# CALCULATION OF THE MEAN ANOMALY OF PLANET B
get_E(ecc_b_s, omega_b_s)
M0_b_s = get_M0_degrees(t_ref, P_days[0], T0[0], ecc_b_s, omega_b_s)

# OF PLANET C 
get_E(ecc_c_s, omega_c_s)
M0_c_s = get_M0_degrees(t_ref, P_days[1], T0[1], ecc_c_s, omega_c_s)

# OF PLANET D
get_E(ecc_d_s,omega_d_s)
M0_d_s = get_M0_degrees(t_ref, P_days[2], T0[2], ecc_d_s, omega_d_s)

# OF PLANET E
get_E(ecc_e_s,omega_e_s)
M0_e_s = get_M0_degrees(t_ref, P_days[3], T0[3], ecc_e_s, omega_e_s)

# OF PLANET F
get_E(ecc_f_s,omega_f_s)
M0_f_s = get_M0_degrees(t_ref, P_days[4], T0[4], ecc_f_s, omega_f_s)

# OF PLANET G
get_E(ecc_g_s,omega_g_s)
M0_g_s = get_M0_degrees(t_ref, P_days[5], T0[5], ecc_g_s, omega_g_s)

# OF PLANET H
get_E(ecc_h_s,omega_h_s)
M0_h_s = get_M0_degrees(t_ref, P_days[6], T0[6], ecc_h_s, omega_h_s)

#############################################################################

# DICTIONARY 
pla_b = dict()
pla_b["mass_earth"] = 	1.374
pla_b["period_days"] = 1.510826
pla_b["inc_deg"] = 89.728
pla_b["Omega_deg"] = 180
pla_b["ecc"] = np.median(ecc_b_s)
pla_b["omega_deg"] = np.median(omega_b_s)
pla_b["M"] = np.median(M0_b_s)

pla_c = dict()
pla_c["mass_earth"] =   1.308
pla_c["period_days"] = 2.421937
pla_c["inc_deg"] = 89.778
pla_c["Omega_deg"] = 180
pla_c["ecc"] = np.median(ecc_c_s) 
pla_c["omega_deg"] = np.median(omega_c_s)
pla_c["M"] = np.median(M0_c_s)

pla_d = dict()
pla_d["mass_earth"] =   0.388
pla_d["period_days"] = 4.049219
pla_d["inc_deg"] = 89.896
pla_d["Omega_deg"] = 180
pla_d["ecc"] = np.median(ecc_d_s)
pla_d["omega_deg"] = np.median(omega_d_s)
pla_d["M"] = np.median(M0_d_s)

pla_e = dict()
pla_e["mass_earth"] =  0.692
pla_e["period_days"] = 6.101013
pla_e["inc_deg"] = 89.793
pla_e["Omega_deg"] = 180
pla_e["ecc"] = np.median(ecc_e_s)
pla_e["omega_deg"] = np.median(omega_e_s)
pla_e["M"] = np.median(M0_e_s)

pla_f = dict()
pla_f["mass_earth"] =   1.039
pla_f["period_days"] = 9.20750
pla_f["inc_deg"] = 89.793
pla_f["Omega_deg"] = 180
pla_f["ecc"] = np.median(ecc_f_s)
pla_f["omega_deg"] = np.median(omega_f_s)
pla_f["M"] = np.median(M0_f_s)

pla_g = dict()
pla_g["mass_earth"] =   1.321
pla_g["period_days"] = 12.352446
pla_g["inc_deg"] = 89.742
pla_g["Omega_deg"] = 180
pla_g["ecc"] = np.median(ecc_g_s)
pla_g["omega_deg"] = np.median(omega_g_s)
pla_g["M"] = np.median(M0_g_s)

pla_h = dict()
pla_h["mass_earth"] =   0.326
pla_h["period_days"] = 18.772866
pla_h["inc_deg"] = 89.805
pla_h["Omega_deg"] = 180
pla_h["ecc"] = np.median(ecc_h_s)
pla_h["omega_deg"] = np.median(omega_h_s)
pla_h["M"] = np.median(M0_h_s)

print("hehe")