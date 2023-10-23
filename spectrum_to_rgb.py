import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import colour as clr
import os


cmfs = clr.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
illuminant = clr.SDS_ILLUMINANTS["D65"]


def piecewise_model(pitch, theta, n, delta_n, wl):
    lambda_0 = n*pitch*np.cos(theta)-(delta_n*pitch/(2*np.cos(theta)))
    lambda_1 = n * pitch*np.cos(theta) + (delta_n * pitch /2)
    theta_max = np.arcsin(n / (n + delta_n/2))  #Trying to fit the saturation to 0.5
    spectrum = {}
    for i in range(len(wl)):

        if theta >= theta_max or lambda_0<=wl[i]<=lambda_1:
            spectrum[str(wl[i])] = 0.5
        else:
            spectrum[str(wl[i])] = 0
    return clr.SpectralDistribution(spectrum, name=f"Model {pitch}, {theta*180/np.pi}:0.2f, biref{delta_n:.3f}")

def load_spectrum_files_and_compute_RGB(spectrum_dir):
    all_spectrum = []
    for filename in os.listdir(spectrum_dir):
        data = np.load(spectrum_dir+"\\"+filename, allow_pickle=True)
        all_spectrum.append(data)
    return all_spectrum

def compute_sds_reflection(data):
    wl = data['wl']
    R_R_to_R = {}
    R_R_to_L = {}
    R_L_to_R = {}
    R_L_to_L = {}
    data_R_unpola = 0.5*(data['R_R_to_R']
               + data['R_R_to_L']
               + data['R_L_to_R']
               + data['R_L_to_L'])
    R_unpola = {}
    for i in range(len(wl)):
        R_R_to_R[str(wl[i])] = data['R_R_to_R'][i]
        R_R_to_L[str(wl[i])] = data['R_R_to_L'][i]
        R_L_to_R[str(wl[i])] = data['R_L_to_R'][i]
        R_L_to_L[str(wl[i])] = data['R_L_to_L'][i]
        R_unpola[str(wl[i])] = data_R_unpola[i]

    theta_deg = data['theta_in_rad'] * 180 / np.pi
    data['sd_R_R_to_R'] = clr.SpectralDistribution(R_R_to_R, name=f"R_R_to_R pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_R_R_to_L'] = clr.SpectralDistribution(R_R_to_L, name=f"R_R_to_L pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_R_L_to_R'] = clr.SpectralDistribution(R_L_to_R, name=f"R_L_to_R pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_R_L_to_L'] = clr.SpectralDistribution(R_L_to_L, name=f"R_L_to_L pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_R_unpola'] = clr.SpectralDistribution(R_unpola, name=f"R_unpola pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")

def compute_sds_transmission(data):
    wl = data['wl']
    T_R_to_R = {}
    T_R_to_L = {}
    T_L_to_R = {}
    T_L_to_L = {}
    data_T_unpola = 0.5*(data['T_R_to_R']
               + data['T_R_to_L']
               + data['T_L_to_R']
               + data['T_L_to_L'])
    T_unpola = {}
    for i in range(len(wl)):
        T_R_to_R[str(wl[i])] = data['T_R_to_R'][i]
        T_R_to_L[str(wl[i])] = data['T_R_to_L'][i]
        T_L_to_R[str(wl[i])] = data['T_L_to_R'][i]
        T_L_to_L[str(wl[i])] = data['T_L_to_L'][i]
        T_unpola[str(wl[i])] = data_T_unpola[i]

    theta_deg = data['theta_in_rad'] * 180 / np.pi
    data['sd_T_R_to_R'] = clr.SpectralDistribution(T_R_to_R, name=f"T_R_to_R pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_T_R_to_L'] = clr.SpectralDistribution(T_R_to_L, name=f"T_R_to_L pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_T_L_to_R'] = clr.SpectralDistribution(T_L_to_R, name=f"T_L_to_R pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_T_L_to_L'] = clr.SpectralDistribution(T_L_to_L, name=f"T_L_to_L pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")
    data['sd_T_unpola'] = clr.SpectralDistribution(T_unpola, name=f"T_unpola pitch = {data['chole'].pitch:3}, theta = {theta_deg:0.2f}")



def sd_to_sRGB(sd):
    xyz = clr.sd_to_XYZ(sd,cmfs,illuminant)
    RGB =  clr.XYZ_to_sRGB(xyz/100)
    return np.clip(RGB,0,1)


def plot_sds_and_model_theta_pitch(spectums):
    def update(val):
        th = slider_th.val
        p = slider_pitch.val
        spectrum = find_spectrum_from_theta_and_p(spectrums, th, p)
        pitch = spectrum['chole'].pitch
        theta = spectrum['theta_in_rad']
        n = (spectrum['n_e'] + spectrum['n_o']) / 2
        delta_n = abs(spectrum['n_e'] - spectrum['n_o'])
        wl = spectrum['wl']
        model = piecewise_model(pitch,theta,n,delta_n,wl)
        ax.cla()
        clr.plotting.plot_multi_sds([model, spectrum['sd_R_unpola']], axes=ax)

        model_swatch = clr.plotting.ColourSwatch(sd_to_sRGB(model))
        spectrum_swatch = clr.plotting.ColourSwatch(sd_to_sRGB(spectrum['sd_R_unpola']))
        clr.plotting.plot_multi_colour_swatches([spectrum_swatch, model_swatch], columns=2, show=False, axes=ax1)
        ax.set_ylim(0, 1)
        fig.subplots_adjust(bottom=0.23)

    fig, [ax,ax1] = plt.subplots(2,1)
    axTheta = fig.add_axes([0.25, 0.10, 0.65, 0.03])
    axPitch = fig.add_axes([0.25, 0.05, 0.65, 0.03])
    slider_th = Slider(axTheta, 'Theta', 0.0, 90)
    slider_pitch = Slider(axPitch, 'Pitch', 250, 500)
    slider_th.on_changed(update)
    slider_pitch.on_changed(update)
    update(0)



def find_spectrum_from_theta_and_p(spectrums, th, p):
    nearest_spectrum = spectrums[0]
    for spectrum in spectrums:
        if abs(nearest_spectrum['chole'].pitch-p)>=abs(spectrum['chole'].pitch-p) and abs(th-nearest_spectrum['theta_in_rad']*180/np.pi)>=abs(th-spectrum['theta_in_rad']*180/np.pi):
            nearest_spectrum = spectrum
    return nearest_spectrum


def plot_sds_and_model_theta_biref(spectums):
    def update(val):
        th = slider_th.val
        b = slider_biref.val
        spectrum = find_spectrum_from_theta_and_b(spectrums, th, b)
        pitch = spectrum['chole'].pitch
        theta = spectrum['theta_in_rad']
        n = (spectrum['n_e'] + spectrum['n_o']) / 2
        delta_n = abs(spectrum['n_e'] - spectrum['n_o'])
        wl = spectrum['wl']
        model = piecewise_model(pitch,theta,n,delta_n,wl)
        ax.cla()
        clr.plotting.plot_multi_sds([model, spectrum['sd_R_unpola']], axes=ax)

        model_swatch = clr.plotting.ColourSwatch(sd_to_sRGB(model))
        spectrum_swatch = clr.plotting.ColourSwatch(sd_to_sRGB(spectrum['sd_R_unpola']))
        clr.plotting.plot_multi_colour_swatches([spectrum_swatch, model_swatch], columns=2, show=False, axes=ax1)
        ax.set_ylim(0, 1)
        fig.subplots_adjust(bottom=0.23)

    fig, [ax,ax1] = plt.subplots(2,1)
    axTheta = fig.add_axes([0.25, 0.10, 0.65, 0.03])
    axBiref = fig.add_axes([0.25, 0.05, 0.65, 0.03])
    slider_th = Slider(axTheta, 'Theta', 45, 90)
    slider_biref= Slider(axBiref, 'Biref', 0, 0.3)
    slider_th.on_changed(update)
    slider_biref.on_changed(update)
    update(0)



def find_spectrum_from_theta_and_b(spectrums, th, b):
    nearest_spectrum = spectrums[0]
    for spectrum in spectrums:
        if abs(abs( nearest_spectrum['n_e'] -  nearest_spectrum['n_o'])-b)>=abs(abs(spectrum['n_e'] - spectrum['n_o'])-b) and abs(th-nearest_spectrum['theta_in_rad']*180/np.pi)>=abs(th-spectrum['theta_in_rad']*180/np.pi):
            nearest_spectrum = spectrum
    return nearest_spectrum

def plot_color_for_theta_and_b(spectrum):
    swatch = []
    model_swatch = []
    for spectrum in spectrums:
        pitch = spectrum['chole'].pitch
        theta = spectrum['theta_in_rad']
        n = (spectrum['n_e'] + spectrum['n_o']) / 2
        delta_n = abs(spectrum['n_e'] - spectrum['n_o'])
        wl = spectrum['wl']
        model = piecewise_model(pitch,theta,n,delta_n,wl)
        swatch.append(clr.plotting.ColourSwatch(sd_to_sRGB(spectrum['sd_R_unpola'])))
        model_swatch.append(clr.plotting.ColourSwatch(sd_to_sRGB(model)))

    fig1, [ax0, ax01] = plt.subplots(2, 1)
    clr.plotting.plot_multi_colour_swatches(swatch, width=1, height=1, columns=45, axes=ax0, show=False)
    ax0.yaxis.set(ticks=[0, 2, 4, 6, 8,10], ticklabels=[0.001, 0.067, 0.134, 0.2, 0.267,0.3])
    ax0.xaxis.set(ticks=[0, 10, 20, 30, 45], ticklabels=[45, 55, 65, 75, 90])
    clr.plotting.plot_multi_colour_swatches(model_swatch,width=1,height=1,columns=45,axes=ax01)
    ax01.yaxis.set(ticks=[0, 2, 4, 6, 8, 10], ticklabels=[0.001, 0.067, 0.134, 0.2, 0.267, 0.3])
    ax01.xaxis.set(ticks=[0, 10, 20, 30, 45], ticklabels=[45, 55, 65, 75, 90])
    plt.xlabel("theta incident starting from 45°")
    plt.ylabel("birefringence")
    plt.legend()


if __name__ == '__main__':
    num_processes = 8
    spectrums = load_spectrum_files_and_compute_RGB("SpectrumsBiref")
    for spectrum in spectrums:
        compute_sds_reflection(spectrum)
    #plot_sds_and_model_theta_biref(spectrums)
    plot_color_for_theta_and_b(spectrum)





