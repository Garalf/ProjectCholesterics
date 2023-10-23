import numpy as np
import pyllama as ll
import cholesteric as ch
from multiprocessing import Pool
import os
import time

# Model Const Parameters :

n_av = 1.481
#biref = 0.1
# n_e = n_av + 0.5 * biref
# n_o = n_av - 0.5 * biref
n_entry = n_av
n_exit = n_av
N_per = 20

pitch = 500

# Ranged Parameters :
# pitch_0 = 250
# pitch_end = 450
# Nb_pitch = 5
# pitches = np.linspace(pitch_0 ,pitch_end,Nb_pitch, dtype=int)

biref = np.linspace(0.001,0.3,10)

Nb_angle = 45
angles =  np.linspace(np.pi/4,np.pi/2,Nb_angle)

# Spectral Range
wl_nm_list = np.arange(400, 800, 1)

# Output directory
out_dir = "SpectrumsBiref"


def cholesteric_spectrums_from_p_and_theta(params):
    t0 = time.time()
    p, th = params
    chole = ch.Cholesteric(pitch360=p, tilt_rad=0, handedness=1)
    spectrum = ll.Spectrum(wl_nm_list,
                           "CholestericModel",
                           dict(chole=chole,
                                n_e=n_e,
                                n_o=n_o,
                                n_entry=n_entry,
                                n_exit=n_exit,
                                N_per=N_per,
                                theta_in_rad=th))
    spectrum.calculate_refl_trans(circ=True, method="SM", talk=False)
    spectrum.export(f"{out_dir}\\spectrum_p{p}_theta{th:.3f}.pck")
    print(f"Spectrum p:{p}, th:{th:.3f} Exported in {time.time()-t0} s\n")



def cholesteric_spectrums_from_biref_and_theta(params):
    t0 = time.time()
    bir, th = params
    n_e = n_av + 0.5 * bir
    n_o = n_av - 0.5 * bir
    chole = ch.Cholesteric(pitch360=pitch, tilt_rad=0, handedness=1)
    spectrum = ll.Spectrum(wl_nm_list,
                           "CholestericModel",
                           dict(chole=chole,
                                n_e=n_e,
                                n_o=n_o,
                                n_entry=n_entry,
                                n_exit=n_exit,
                                N_per=N_per,
                                theta_in_rad=th))
    spectrum.calculate_refl_trans(circ=True, method="SM", talk=False)
    spectrum.export(f"{out_dir}\\spectrum_b{bir:.3f}_theta{th:.3f}.pck")
    print(f"Spectrum b:{bir}, th:{th:.3f} Exported in {time.time()-t0} s\n")


if __name__ == '__main__':

    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass
    #combinaisons = [(p, th) for p in pitches for th in angles]
    combinaisons = [(b, th) for b in biref for th in angles]

    num_processes = 8
    with Pool(num_processes) as pool:
        pool.map(cholesteric_spectrums_from_biref_and_theta,combinaisons)
